#include "YRECIO.h"
#include <Python.h>
#include <OrbitSolver.h>

#define STR(x) #x
#define STRING(str) STR(str)
#define DEFAULT_MSTAR				1.0
#define DEFAULT_MPLANET				1.0
#define DEFAULT_QSTAR				1e7
#define DEFAULT_P0					3.0
#define DEFAULT_WINDK				0.17
#define DEFAULT_WIND_WSAT			2.2
#define DEFAULT_CORE_ENV_COUPLING	0.03
#define DEFAULT_STAR_P0				10.0
#define DEFAULT_DISK_DISSIPATION	0.005

#ifdef __cplusplus
extern "C" {
#endif 

	static PyObject *serialize_stellar_evolution(PyObject *self,
			PyObject *args)
	{
		const char *filename;
		if (!PyArg_ParseTuple(args, "s", &filename)) {
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

	static PyObject *evolve_orbit(PyObject *self, PyObject *kw_args)
	{
		Py_RETURN_NONE;
	}

	static PyMethodDef OrbitSolverMethods[] = {
		{"serialize_stellar_evolution", serialize_stellar_evolution,
			METH_VARARGS, "Reads in the stellar evolution tracks and "
				"generates a file storing the interpolation information "
				"based on those tracks.\n\nThis file can later be read in to"
				" re-create stellar evolution interpolation much faster than"
				" generating in from the tracks again. The only argument is "
				"the name of the output file."},
		{"evolve", evolve_orbit, METH_KEYWORDS, "Calculates the orbital "
			"evolution for a given planet-star system from given initial "
			"conditions.\n\nOnly keyword arguments are accepted (no "
			"positional) with the following names:\n."
			"\t* star_mass: stellar mass in solar masses "
			"(default " STRING(DEFAULT_MSTAR) ")\n"
			"\t* planet_mass: planet mass in Jupiter masses "
			"(default " STRING(DEFAULT_MPLANET) ")\n"
			"\t* Qstar: tidal quality factor of the star "
			"(default " STRING(DEFAULT_QSTAR) ")\n"
			"\t* P0: initial orbital period in days "
			"(default " STRING(DEFAULT_P0) ")\n"
			"\t* a0: initial semimajor axis in AU (no default)\n"
			"\t* wind_K: with strength in solar units "
			"(default " STRING(DEFAULT_WINDK) ")\n"
			"\t* wind_wsat: wind saturation frequency "
			"(default " STRING(DEFAULT_WIND_WSAT) ")\n"
			"\t* core_env_coupling: the timescale for coupling of the "
			"stellar core and envelope in Gyr "
			"(default " STRING(DEFAULT_CORE_ENV_COUPLING) ")\n"
			"\t* star_P0: initial spin period of the star in days "
			"(default " STRING(DEFAULT_START_P0) ")\n"
		    "\t* star_W0: initial spin frequency of the star in rad/day "
			"(no default)\n"
			"\t* disk_lifetime: the age at which the disk dissipates "
			"releasing the surface rotation of the star (default "
			STRING(DEFAULT_DISK_DISSIPATION) ").\n"
			"An exception will be raised if both (P0 and a0) or (star_P0 and"
			" star_w0) are specified."}
	};

	PyMODINIT_FUNC initOrbitSolver(void)
	{
		(void) Py_InitModule("OrbitSolver", OrbitSolverMethods);
	}

#ifdef __cplusplus
}
#endif
