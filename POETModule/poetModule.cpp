#include "YRECIO.h"
#include "Common.h"
#include "OrbitSolver.h"

#include <Python.h>
#include <sstream>

#define STR(x) #x
#define STRING(str) STR(str)
#define DEFAULT_MSTAR				1.0
#define DEFAULT_MPLANET				1.0
#define DEFAULT_RPLANET				1.0
#define DEFAULT_QSTAR				1e7
#define DEFAULT_P0					3.0
#define DEFAULT_WINDK				0.17
#define DEFAULT_WIND_WSAT			2.45
#define DEFAULT_CORE_ENV_COUPLING	0.028
#define DEFAULT_STAR_P0				7.0
#define DEFAULT_DISK_DISSIPATION	0.005
#define DEFAULT_STELLAR_INTERP		"serialized_evolution"

#ifdef __cplusplus
extern "C" {
	static PyObject *serialize_stellar_evolution(PyObject *self,
			PyObject *args);
	static PyObject *evolve_orbit(PyObject *self, PyObject *pos_args,
			PyObject *kw_args);
}
#endif 

static PyObject *serialize_stellar_evolution(PyObject *self,
		PyObject *args)
{
	const char *filename=DEFAULT_STELLAR_INTERP;
	if (!PyArg_ParseTuple(args, "|s", &filename)) {
		PyErr_SetString(PyExc_TypeError, 
				"Expecting exactly one argument: the filename where to "
				"save the serialized interpolation.");
		return NULL;
	}
	std::cout << "Serializing the evolution read from " STRING(YREC_TRACK_PATH)
		" to " << filename << std::endl;
	YRECEvolution toSave(STRING(YREC_TRACK_PATH), 0, 2.0, 2.0);
	toSave.save_state(filename);
	Py_RETURN_NONE;
}

///Fills the given python list with the values in source_list. It is assumed
///that py_list has the correct size already.
int fill_list(PyObject *py_list, const std::list<double> &source_list,
		const std::string &list_name)
{
	Py_ssize_t py_i=0;
	for(std::list<double>::const_iterator src_i=source_list.begin();
			src_i!=source_list.end(); src_i++) {
		if(PyList_SetItem(py_list, py_i, Py_BuildValue("d", *src_i))) {
			std::ostringstream msg;
			msg << "Failed to set " << list_name << "[" << py_i << "]";
			PyErr_SetString(PyExc_IndexError, msg.str().c_str());
			return -1;
		}
		py_i++;
	}
	return 0;
}

///Fills the given python list with strings representing the evolution modes
///in source_list.
int fill_list(PyObject *py_list, const std::list<EvolModeType> &source_list)
{
	Py_ssize_t py_i=0;
	for(std::list<EvolModeType>::const_iterator src_i=source_list.begin();
			src_i!=source_list.end(); src_i++) {
		std::ostringstream mode;
		mode << *src_i;
		if(PyList_SetItem(py_list, py_i,
					Py_BuildValue("s", mode.str().c_str()))) {
			std::ostringstream msg;
			msg << "Failed to set evolution mode [" << py_i << "]";
			PyErr_SetString(PyExc_IndexError, msg.str().c_str());
			return -1;
		}
		py_i++;
	}
	return 0;
}

///Generates the return value of the evolve_orbit function (see docstring
///for details).
PyObject *pythonize_tabulated_orbit(const OrbitSolver &solver)
{
	const std::list<double> *ages=solver.get_tabulated_var(AGE);
	size_t num_steps=ages->size();
	PyObject *t_list=PyList_New(num_steps);
	if(!t_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated ages list.");
		return NULL;
	}
	if(fill_list(t_list, *ages, "tabulated ages")) {
		Py_DECREF(t_list);
		return NULL;
	}

	PyObject *mode_list=PyList_New(num_steps);
	if(!mode_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated evolution mode list.");
		Py_DECREF(t_list); 
		return NULL;
	}
	if(fill_list(mode_list, *(solver.get_tabulated_evolution_mode()))) {
		Py_DECREF(t_list); Py_DECREF(mode_list);
		return NULL;
	}

	PyObject *a_list=PyList_New(num_steps);
	if(!a_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated semimajor axes list.");
		Py_DECREF(t_list); Py_DECREF(mode_list);
		return NULL;
	}
	if(fill_list(a_list, *(solver.get_tabulated_var(SEMIMAJOR)),
				"tabulated semimajor axes")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		return NULL;
	}

	PyObject *a_dot_list=PyList_New(num_steps);
	if(!a_dot_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated semimajor axes derivative "
				"list.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		return NULL;
	}
	if(fill_list(a_dot_list, *(solver.get_tabulated_var_deriv(SEMIMAJOR)),
				"tabulated semimajor axis derivatives")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list);
		return NULL;
	}

	PyObject *Lconv_list=PyList_New(num_steps);
	if(!Lconv_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated convective angular momentum "
				"list.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list);
		return NULL;
	}
	if(fill_list(Lconv_list, *(solver.get_tabulated_var(LCONV)),
				"Tabulated convective angular momenta")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		return NULL;
	}

	PyObject *Lconv_dot_list=PyList_New(num_steps);
	if(!Lconv_dot_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated convective angular momentum "
				"derivative list.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		return NULL;
	}
	if(fill_list(Lconv_dot_list, *(solver.get_tabulated_var_deriv(LCONV)),
				"Tabulated convective angular momentum derivatives")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list);
		return NULL;
	}

	PyObject *Lrad_list=PyList_New(num_steps);
	if(!Lrad_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated radiative angular momentum "
				"list.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list);
		return NULL;
	}
	if(fill_list(Lrad_list, *(solver.get_tabulated_var(LRAD)),
				"Tabulated radiative angulare momenta")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		return NULL;
	}

	PyObject *Lrad_dot_list=PyList_New(num_steps);
	if(!Lrad_dot_list) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated radiative angular momentum "
				"derivative list.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		return NULL;
	}
	if(fill_list(Lrad_dot_list, *(solver.get_tabulated_var_deriv(LRAD)),
				"Tabulated radiative angular momentum derivatives")) {
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list);
		return NULL;
	}
	PyObject *result=PyDict_New();
	if(!result) {
		PyErr_SetString(PyExc_MemoryError,
				"Failed to allocate tabulated orbit dictionary.");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list);
		return NULL;
	}

	if(PyDict_SetItemString(result, "t", t_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['t']");
		Py_DECREF(t_list); Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(t_list);

	if(PyDict_SetItemString(result, "evol_mode", mode_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['evol_mode']");
		Py_DECREF(mode_list); Py_DECREF(a_list);
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(mode_list);

	if(PyDict_SetItemString(result, "a", a_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['a']");
		Py_DECREF(a_list); Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(a_list);

	if(PyDict_SetItemString(result, "a_dot", a_dot_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['a_dot']");
		Py_DECREF(a_dot_list); Py_DECREF(Lconv_list);
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(a_dot_list);

	if(PyDict_SetItemString(result, "Lconv", Lconv_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['Lconv']");
		Py_DECREF(Lconv_list); Py_DECREF(Lconv_dot_list);
		Py_DECREF(Lrad_list); Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(Lconv_list);

	if(PyDict_SetItemString(result, "Lconv_dot", Lconv_dot_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['Lconv_dot']");
		Py_DECREF(Lconv_dot_list); Py_DECREF(Lrad_list);
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(Lconv_dot_list);

	if(PyDict_SetItemString(result, "Lrad", Lrad_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['Lrad']");
		Py_DECREF(Lrad_list); Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(Lrad_list);

	if(PyDict_SetItemString(result, "Lrad_dot", Lrad_dot_list)) {
		PyErr_SetString(PyExc_KeyError,
				"Failed to set tabulated orbital evolution['Lrad_dot']");
		Py_DECREF(Lrad_dot_list); Py_DECREF(result);
		return NULL;
	}
	Py_DECREF(Lrad_dot_list);

	return result;
}

//If argument is a single floating point value store it in value, otherwise
//check that it is a sequence and if num_orbits is zero it gets set to the
//sequence length, otheriwe check that the sequenc is of length num_orbits.
//
//Returns 0 if all check succeed and 1 otherwise.
//
//Sets the correct AssertionError if a check fails.
static int process_evolve_argument(PyObject *argument, double &value,
		Py_ssize_t &num_orbits)
{
	if(argument==NULL) return 0;
	if(PyFloat_Check(argument)) {
		value=PyFloat_AsDouble(argument);
		return 0;
	} else if(!PySequence_Check(argument)) {
		PyErr_SetString(PyExc_AssertionError,
				"star mass is neither a single floating point value nor "
				"a sequence!");
		return 1;
	}
	if(num_orbits==0) {
		num_orbits=PySequence_Size(argument);
		return 0;
	} else if(PySequence_Size(argument)!=num_orbits) {
		PyErr_SetString(PyExc_AssertionError,
				"Not all input lists are of the same length!");
		return 1;
	}
	return 0;
}

//If argument is a sequence sets value to its index-th entry, otherwise
//does not touch value.
//
//Returns non-zero if some python exception occurs, in which case the
//appropriate exception is set.
static int get_from_sequence(PyObject *argument, Py_ssize_t index,
		double &value)
{
	if(argument!=NULL && PySequence_Check(argument)) {
		PyObject *arg_i=PyList_GetItem(argument, index);
		value=PyFloat_AsDouble(arg_i);
		Py_DECREF(arg_i);
		if(PyErr_Occurred()) return 1;
	} 
	return 0;
}

static PyObject *evolve_orbit(PyObject *self, PyObject *pos_args,
		PyObject *kw_args)
{
	double star_mass=DEFAULT_MSTAR,
		   planet_mass=DEFAULT_MPLANET,
		   planet_radius=DEFAULT_RPLANET,
		   Qstar=DEFAULT_QSTAR,
		   P0=NaN, a0=NaN,
		   wind_K=DEFAULT_WINDK,
		   wind_wsat=DEFAULT_WIND_WSAT,
		   core_env_coupling=DEFAULT_CORE_ENV_COUPLING,
		   star_P0=NaN, star_W0=NaN,
		   disk_lifetime=DEFAULT_DISK_DISSIPATION;
	PyObject *star_mass_obj=NULL, 
			 *planet_mass_obj=NULL,
			 *planet_radius_obj=NULL,
			 *Qstar_obj=NULL,
			 *P0_obj=NULL, *a0_obj=NULL,
			 *wind_K_obj=NULL,
			 *wind_wsat_obj=NULL,
			 *core_env_coupling_obj=NULL,
			 *star_P0_obj=NULL, *star_W0_obj=NULL,
			 *disk_lifetime_obj=NULL;
	const char *stellar_evolution_interp=DEFAULT_STELLAR_INTERP;
	char *kwlist[]={const_cast<char*>("star_mass"),
		const_cast<char*>("planet_mass"),
		const_cast<char*>("planet_radius"),
		const_cast<char*>("Qstar"),
		const_cast<char*>("P0"),
		const_cast<char*>("a0"),
		const_cast<char*>("wind_K"),
		const_cast<char*>("wind_wsat"),
		const_cast<char*>("core_env_coupling"),
		const_cast<char*>("star_P0"),
		const_cast<char*>("star_W0"),
		const_cast<char*>("disk_lifetime"), 
		const_cast<char*>("stellar_evolution_interp"), NULL};
	if(!PyArg_ParseTupleAndKeywords(pos_args, kw_args, "|OOOOOOOOOOOs",
				kwlist, &star_mass_obj, &planet_mass_obj, &planet_radius_obj,
				&Qstar_obj, &P0_obj, &a0_obj, &wind_K_obj, &wind_wsat_obj,
				&core_env_coupling_obj, &star_P0_obj, &star_W0_obj,
				&disk_lifetime_obj, &stellar_evolution_interp))
		return NULL;
	Py_ssize_t num_orbits=0;
	if(process_evolve_argument(star_mass_obj, star_mass, num_orbits))
		return NULL;
	if(process_evolve_argument(planet_mass_obj, planet_mass, num_orbits))
		return NULL;
	if(process_evolve_argument(planet_radius_obj, planet_radius, num_orbits))
		return NULL;
	if(process_evolve_argument(Qstar_obj, Qstar, num_orbits))
		return NULL;
	if(P0_obj!=NULL && a0_obj!=NULL) {
		PyErr_SetString(PyExc_AssertionError,
				"Only one of P0 or a0 arguments should be specified.");
		return NULL;
	} else if(P0_obj==NULL) P0=DEFAULT_P0;
	if(process_evolve_argument(P0_obj, P0, num_orbits))
		return NULL;
	if(process_evolve_argument(a0_obj, a0, num_orbits))
		return NULL;
	if(process_evolve_argument(wind_K_obj, wind_K, num_orbits))
		return NULL;
	if(process_evolve_argument(wind_wsat_obj, wind_wsat, num_orbits))
		return NULL;
	if(process_evolve_argument(core_env_coupling_obj, core_env_coupling,
				num_orbits)) return NULL;
	if(star_P0_obj!=NULL && star_W0_obj!=NULL) {
		PyErr_SetString(PyExc_AssertionError,
				"Only one of star_P0 or star_W0 arguments should be "
				"specified.");
		return NULL;
	} else if(star_P0_obj==NULL) star_P0=DEFAULT_STAR_P0;
	if(process_evolve_argument(star_P0_obj, star_P0, num_orbits))
		return NULL;
	if(process_evolve_argument(star_W0_obj, star_W0, num_orbits))
		return NULL;
	if(process_evolve_argument(disk_lifetime_obj, disk_lifetime, num_orbits))
		return NULL;
	if(num_orbits==0) num_orbits=1;
	try {
		StellarEvolution evol;
		evol.load_state(stellar_evolution_interp);
		PyObject *orbit_list=(num_orbits>1 ? PyList_New(num_orbits) : NULL);
		for(Py_ssize_t orb_ind=0; orb_ind<num_orbits; orb_ind++) {
			if(get_from_sequence(star_mass_obj, orb_ind, star_mass))
				return NULL;
			if(get_from_sequence(planet_mass_obj, orb_ind, planet_mass))
				return NULL;
			if(get_from_sequence(planet_radius_obj, orb_ind, planet_radius))
				return NULL;
			if(get_from_sequence(Qstar_obj, orb_ind, Qstar)) return NULL;
			if(get_from_sequence(P0_obj, orb_ind, P0)) return NULL;
			if(get_from_sequence(a0_obj, orb_ind, a0)) return NULL;
			if(a0_obj==NULL)
				a0=std::pow(AstroConst::G*star_mass*AstroConst::solar_mass*
						std::pow(P0*AstroConst::day, 2)/4.0/M_PI/M_PI,
						1.0/3.0)/AstroConst::AU;
			if(get_from_sequence(wind_K_obj, orb_ind, wind_K)) return NULL;
			if(get_from_sequence(wind_wsat_obj, orb_ind, wind_wsat))
				return NULL;
			if(get_from_sequence(core_env_coupling_obj, orb_ind,
						core_env_coupling)) return NULL;
			if(get_from_sequence(star_P0_obj, orb_ind, star_P0)) return NULL;
			if(get_from_sequence(star_W0_obj, orb_ind, star_W0)) return NULL;
			if(star_W0_obj==NULL) star_W0=2.0*M_PI/star_P0;
			if(get_from_sequence(disk_lifetime_obj, orb_ind, disk_lifetime))
				return NULL;
			Star star(star_mass, Qstar, wind_K, wind_wsat, core_env_coupling, 0,
					star_W0, disk_lifetime, evol);
			Planet planet(&star, planet_mass, planet_radius, a0);
			StellarSystem system(&star, &planet);
			OrbitSolver solver(std::min(star.core_formation_age(), disk_lifetime),
					star.get_lifetime(), 1e-5);
			solver(system, Inf, disk_lifetime, a0);
			PyObject *orbit=pythonize_tabulated_orbit(solver);
			if(num_orbits==1) return orbit;
			PyList_SetItem(orbit_list, orb_ind, orbit);
		}
		return orbit_list;
	} catch (Error::General err) {
		std::ostringstream msg;
		msg << err.what() << ": " << err.get_message();
		PyErr_SetString(PyExc_RuntimeError, msg.str().c_str());
		return NULL;
	}
}

static PyMethodDef poetMethods[] = {
	{"serialize_stellar_evolution", serialize_stellar_evolution,
		METH_VARARGS,
		"Prepares and saves in a file the interpolation between YREC "
		"tracks.\n"
		"\n"
		"Reads in the stellar evolution tracks and generates a file storing "
		"the interpolation information based on those tracks. This file can "
		"later be read in to re-create stellar evolution interpolation much "
		"faster than generating it from the tracks again.\n"
		"\n"
		"Args:\n"
		"    filename: the name of the output file (default "
		DEFAULT_STELLAR_INTERP ").\n"
		"\n"
		"Raises:\n"
		"    TypeError: the input argument was not a string."},
	{"evolve_orbit", (PyCFunction)evolve_orbit, METH_VARARGS|METH_KEYWORDS, 
		"Calculates the orbital evolution for a planet-star system.\n"
		"\n"
		"Args:\n"
		"    star_mass: stellar mass in solar masses "
		"(default " STRING(DEFAULT_MSTAR) ")\n"
		"    planet_mass: planet mass in Jupiter masses "
		"(default " STRING(DEFAULT_MPLANET) ")\n"
		"    planet_radius: planet radius in Jupiter radii "
		"(default " STRING(DEFAULT_RPLANET) ")\n"
		"    Qstar: tidal quality factor of the star "
		"(default " STRING(DEFAULT_QSTAR) ")\n"
		"    P0: initial orbital period in days "
		"(default " STRING(DEFAULT_P0) ")\n"
		"    a0: initial semimajor axis in AU (no default)\n"
		"    wind_K: with strength in solar units "
		"(default " STRING(DEFAULT_WINDK) ")\n"
		"    wind_wsat: wind saturation frequency "
		"(default " STRING(DEFAULT_WIND_WSAT) ")\n"
		"    core_env_coupling: the timescale for coupling of the "
		"stellar core and envelope in Gyr "
		"(default " STRING(DEFAULT_CORE_ENV_COUPLING) ")\n"
		"    star_P0: initial spin period of the star in days (the disk "
		"locking period) default " STRING(DEFAULT_STAR_P0) ")\n"
		"    star_W0: initial spin frequency of the star in rad/day "
		"(no default)\n"
		"    disk_lifetime: the age at which the disk dissipates "
		"releasing the surface rotation of the star (default "
		STRING(DEFAULT_DISK_DISSIPATION) ").\n"
		"    stellar_evolution_interp: the filename of a previously stored "
		"stellar evolution interpolation (default "
		DEFAULT_STELLAR_INTERP ").\n"
		"Any argument (other than stellar_evolution_interp) can be either a "
		"single real value or a list of values. If multiple arguments are "
		"lists they must all be of the same length.\n"
		"Returns:\n"
		"    A dictionary indexed by variable name containing a list of the "
		"values the corresponding variable gets at each step. The keys "
		"are:\n"
		"        t: the age.\n"
		"        evol_mode: the evolution mode.\n"
		"        a: the semimajor axis.\n"
		"        a_dot: the age derivative of the semimajor axis.\n"
		"        Lconv: the convective envelope angular momentum.\n"
		"        Lconv_dot: the age derivative of Lconv.\n"
		"        Lrad: the radiative core angular momentum.\n"
		"        Lrad_dot: the age derivative of Lrad.\n"
		"If any of the input arguments are lists, a list of dictionaries is "
		"returned of the same length as the input argument lists.\n"
		"Raises:\n"
		"    AssertionError: if both (P0 and a0) or (star_P0 and star_w0) "
		"are specified.\n"
		"    TypeError: if the input arguments were not of the correct "
		"types.\n"
		"    MemoryError: if allocating some result object fails.\n"
		"    IndexError: if setting some entry in an output list fails.\n"
		"    KeyError: if setting some entry in the output dictionary fails."
	},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpoet(void)
{
	(void) Py_InitModule("poet", poetMethods);
}
