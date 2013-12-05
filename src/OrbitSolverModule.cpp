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
#define DEFAULT_WIND_WSAT			2.2
#define DEFAULT_CORE_ENV_COUPLING	0.03
#define DEFAULT_STAR_P0				10.0
#define DEFAULT_DISK_DISSIPATION	0.005
#define DEFAULT_STELLAR_INTERP		"interp_state_data"

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
	if(!PyArg_ParseTupleAndKeywords(pos_args, kw_args, "|ddddddddddds",
				kwlist, &star_mass, &planet_mass, &planet_radius, &Qstar,
				&P0, &a0, &wind_K, &wind_wsat, &core_env_coupling, &star_P0,
				&star_W0, &disk_lifetime, &stellar_evolution_interp))
		return NULL;
	if(!std::isnan(P0) && !std::isnan(a0)) {
		PyErr_SetString(PyExc_AssertionError,
				"Only one of P0 or a0 arguments should be specified.");
		return NULL;
	} else if(std::isnan(P0)) P0=DEFAULT_P0;
	if(std::isnan(a0))
		a0=std::pow(AstroConst::G*star_mass*AstroConst::solar_mass*
				std::pow(P0*AstroConst::day, 2)/4.0/M_PI/M_PI, 1.0/3.0)/
			AstroConst::AU;
	if(!std::isnan(star_P0) && !std::isnan(star_W0)) {
		PyErr_SetString(PyExc_AssertionError,
				"Only one of star_P0 or star_W0 arguments should be "
				"specified.");
		return NULL;
	} else if(std::isnan(star_P0)) star_P0=DEFAULT_STAR_P0;
	if(std::isnan(star_W0)) star_W0=2.0*M_PI/star_P0;
	try {
		YRECEvolution evol;
		evol.load_state(stellar_evolution_interp);
		Star star(star_mass, Qstar, wind_K, wind_wsat, core_env_coupling, 0,
				star_W0, disk_lifetime, evol);
		Planet planet(&star, planet_mass, planet_radius, a0);
		StellarSystem system(&star, &planet);
		OrbitSolver solver(std::min(star.core_formation_age(), disk_lifetime),
				star.get_lifetime(), 1e-5, Inf, Inf);
		solver(system, Inf, disk_lifetime, a0);
		return pythonize_tabulated_orbit(solver);
	} catch (Error::General err) {
		std::ostringstream msg;
		msg << err.what() << ": " << err.get_message();
		PyErr_SetString(PyExc_RuntimeError, msg.str().c_str());
		return NULL;
	}
}

static PyMethodDef OrbitSolverMethods[] = {
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
		"Any argument can be either a single real value or a list of values."
	    " If multiple arguments are lists they must all be of the same "
		"length.\n"
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

PyMODINIT_FUNC initOrbitSolver(void)
{
	(void) Py_InitModule("OrbitSolver", OrbitSolverMethods);
}
