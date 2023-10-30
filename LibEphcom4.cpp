/* 
 * LibEphcom4.cpp
 * Dynamically loads the ephcom4 library and provides geocentric
 * positions used in the Mars opposition calculator.
 * Uses Alan Irwin's 'timeephem' library as at 1st July 2020.
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <dlfcn.h>
#include "LibEphcom4.h"


/* Paths to TimeEphem library by Alan Irwin and JPL Digital Ephemeris model.
 *   https://sourceforge.net/projects/timeephem/"
 *   https://ssd.jpl.nasa.gov/planets/eph_export.html"
 */
const char EphcomBaseDir[] =
	"/home/thomas/astro/timeephem/ephcom";
const char EphcomDLL[] = "build/src/libs/libephcomc.so";
const char JPLDataFilename[] =
	"/home/thomas/astro/timeephem/ephcom/de432/eph1980_2100.432";

/*** REFERENCE ROUTINES ***/
void convert_ICRF_eq_to_ecliptic(const double *pos, double *pos_ecl);
void convert_ICRF_eq_to_mean_ecliptic(double JDE, const double *pos, double *pos_ecl);
void xyz_to_sph(double *pos, double *sph);

int LibEphcom4::loadLib()
{
	size_t len = strlen(EphcomBaseDir) + strlen(EphcomDLL) + 3;
	char *EphcomFilename = (char*) malloc(len);

	sprintf(EphcomFilename, "%s/%s", EphcomBaseDir, EphcomDLL);
	libDLL = dlopen(EphcomFilename, RTLD_LAZY);
	if (!libDLL)
	{
		std::cerr << "Failed to open library. " << dlerror() << std::endl;
		exit (1);
	}
	free (EphcomFilename);

	/* Register API functions */
	fnEphcom4Open = reinterpret_cast<ephcom_context* (*)(const char*)> (dlsym(libDLL, "ephcom_open"));
	fnEphcom4Close = reinterpret_cast<int (*)( ephcom_context * )> (dlsym(libDLL, "ephcom_close"));
	fnEphcom4GetCoords = reinterpret_cast<int (*)( ephcom_context *, ephcom_coords * )> (dlsym(libDLL, "ephcom_get_coords"));
	fnEphcom4Interpolate = reinterpret_cast<int (*)( ephcom_context *, const double *, unsigned, int, int, double * )> (dlsym(libDLL, "ephcom_interpolate_relative"));
	fnEphcom4Alloc2dChar = reinterpret_cast<int (*)( char ***, unsigned, unsigned )> (dlsym(libDLL, "ephcom_Alloc2dChar"));
	fnEphcom4Free2dChar = reinterpret_cast<int (*)( char **, unsigned )> (dlsym(libDLL, "ephcom_Free2dChar"));
	fnEphcom4ReadConstants = reinterpret_cast<int (*)( ephcom_context *, char **, unsigned, unsigned, double *, unsigned *, double *, double *, double *)>  (dlsym(libDLL, "ephcom_read_constants"));

	assert (fnEphcom4Open != NULL);
	assert (fnEphcom4Close != NULL);
	assert (fnEphcom4GetCoords != NULL);
	assert (fnEphcom4Interpolate != NULL);
	assert (fnEphcom4Alloc2dChar != NULL);
	assert (fnEphcom4Free2dChar != NULL);
	assert (fnEphcom4ReadConstants != NULL);

	return 0;
}

void LibEphcom4::initContext(ephcom_object ephcomID)
{
	std::unique_lock<std::mutex> lck {ctxm};

	loadLib();

	if (fnEphcom4Open != NULL)
	{
		ctx = (*fnEphcom4Open)(JPLDataFilename);
	}
	
	if (ctx == NULL)
	{
		std::cerr << "ERROR: Context cannot be initialised." << std::endl;
		return;
	}

	int            actual_nvalues;
	unsigned       nvalues = EPHCOM_MAXCON;
	char           **cnames;
	double         *values, ss[3], et2[2];
	double         au, emrat;
	unsigned       numde;

	if ( ( values = (double *) calloc( 1, (size_t) nvalues * sizeof ( double ) ) ) == NULL )
	{
		std::cerr << "ERROR: Can't allocate memory" << std::endl;
		return;
	}
	fnEphcom4Alloc2dChar( &cnames, nvalues, EPHCOM_MAX_CNAM_LENGTH + 1 );
	actual_nvalues = fnEphcom4ReadConstants( ctx, cnames, nvalues, EPHCOM_MAX_CNAM_LENGTH, values, &numde, &au, &emrat, ss );
	if ( actual_nvalues < 0 )
	{
		std::cerr << "ERROR: Call to ephcom_read_constants failed." << std::endl;
		return;
	}

	fnEphcom4Free2dChar( cnames, nvalues );
	free( values );

	/* Populate the ephemeris properties */
	tJDStart = ss[0];
	tJDStop =  ss[1];
	tJDStep =  ss[2];

	/* use J2000.0 ecliptic/equinoxes */
	equinoxJ2000 = true;
	targetID = ephcomID;

	std::cout << "JPL DE" << numde << " ephemeris session started." << std::endl;
}

void LibEphcom4::deinitContext()
{
	std::unique_lock<std::mutex> lck {ctxm};
	if (fnEphcom4Close != NULL)
	{
		(*fnEphcom4Close)(ctx);
		ctx = NULL;
	}

	dlclose(libDLL);
	std::cout << "JPL ephemeris session ended." << std::endl;
}

int LibEphcom4::setObject(ephcom_object ephcomID)
{
	if (ephcomID > EPHCOM_ORDERED_INDICES)
	{
		std::cerr << "setObject: Target object out of range. Expecting index from "
		          << EPHCOM_MERCURY << " to "
		          << EPHCOM_ORDERED_INDICES << ")." << std::endl;
		return JPLEPH_ERR_TARGET_OUTOFRANGE;
	}
	targetID = ephcomID;
	return JPLEPH_ERR_OK;
}

ephcom_object LibEphcom4::getObject()
{
	return targetID;
}

void LibEphcom4::setMeanEquinox()
{
	equinoxJ2000 = false;
}

int LibEphcom4::computeDE(int target, const double time, jpl_PosData *data)
{
	if (time < tJDStart || time > tJDStop)
	{
		return JPLEPH_ERR_TIME_OUTOFRANGE;
	}
	if (targetID >= EPHCOM_MAXOBJECTS)
	{
		std::cerr << "computeDE: new target: " << target << std::endl;
		targetID = static_cast<ephcom_object> (target);
	}
	ephcom_header header1;
	ephcom_coords coords;
	double  r_Geo[6], r_Helio[6];

	/* Set up coords */
 	coords.et2[0] = time; coords.et2[1] = 0.0;
	coords.bary = EPHCOM_BARYCENTRIC_COORDS;
	coords.km = EPHCOM_UNITS_AU;
	coords.seconds = EPHCOM_UNITS_DAYS;
	for (int i = 0; i < EPHCOM_MAXOBJECTS; i++ )
	{
		coords.list[i] = 0;
	}
	coords.list[EPHCOM_EARTH] = 1;
	coords.list[targetID] = 1;

	int err = fnEphcom4GetCoords(ctx, &coords);
	if (err != 0)
	{
		std::cerr << "fnEphcom4GetCoords: An error occurred: " << err << std::endl;
		return JPLEPH_ERR_NODATA;
	}

	err = fnEphcom4Interpolate(ctx, coords.et2, 0, target, EPHCOM_EARTH, r_Geo); /* Geocentric position */
	if (err != 0)
	{
		std::cerr << "fnEphcom4Interpolate: An error occurred: " << err << std::endl;
		return JPLEPH_ERR_NODATA;
	}
	err = fnEphcom4Interpolate(ctx, coords.et2, 0, target, EPHCOM_SUN,   r_Helio); /* Heliocentric position */
	if (err != 0)
	{
		std::cerr << "fnEphcom4Interpolate: An error occurred: " << err << std::endl;
		return JPLEPH_ERR_NODATA;
	}

	/*
	if (equinoxJ2000)
	{
		double r_Geo_mean_equinox[3];
		new_equinox(testjd, r_Geo, r_Geo_mean_equinox);
		for (i=0; i<3; i++)
			r_Geo[i] = r_Geo_mean_equinox[i];
	}
	*/

	/* compute geo coords to RA and dec */
	double  sph_coords_Geo[3], sph_coords_Helio[3];
	xyz_to_sph(r_Geo, sph_coords_Geo);
	//deg_to_hms(sph_coords_Geo[1], &RightAsc);
	//deg_to_dms(sph_coords_Geo[2], &Declination);

	xyz_to_sph(r_Helio, sph_coords_Helio);
	
	/* copy back the results */
	data->dist_Geo = sph_coords_Geo[0];
	data->dist_Helio = sph_coords_Helio[0];
	return JPLEPH_ERR_OK;
}

int LibEphcom4::computeGeocentricPos(int target, const double time, jpl_PosDataVer2 *data)
{
	if (time < tJDStart || time > tJDStop)
	{
		return JPLEPH_ERR_TIME_OUTOFRANGE;
	}

	const double et2[2] = {time, 0.0};

	/* To hold x, xdot, y, ydot, z, zdot for all bodies */
	double  r_Geo[6], r_Geo_ecliptic[6];
	int err = fnEphcom4Interpolate(ctx, et2, 0, target, EPHCOM_EARTH, r_Geo);

	if (err != 0)
	{
		std::cerr << "fnEphcom4Interpolate: An error occurred: " << err << std::endl;
		return JPLEPH_ERR_NODATA;
	}

	if (equinoxJ2000)
	{
		/* convert ICRF equatorial to ecliptic using J2000.0 equinox */
		convert_ICRF_eq_to_ecliptic(r_Geo, r_Geo_ecliptic);
	}
	else
	{
		/* Convert ICRF equatorial frame to mean equinox */
		convert_ICRF_eq_to_mean_ecliptic(time, r_Geo, r_Geo_ecliptic);
	}

	/* copy rectangular coords */
	for (int j=0; j<3; ++j)
	{
		data->posXYZ[j] = r_Geo_ecliptic[j];
	}

	/* convert Geocentric coords to spherical position */
	double  sph_coords_Geo[3];
	xyz_to_sph(r_Geo_ecliptic, sph_coords_Geo);

	data->radius    = sph_coords_Geo[0];
	data->longitude = sph_coords_Geo[1]; /* [0,360] deg */
	data->latitude  = sph_coords_Geo[2]; /* [-90,90] deg */
	return JPLEPH_ERR_OK;
}

int LibEphcom4::computeHeliocentricPos(int target, const double time, jpl_PosDataVer2 *data)
{
	if (time < tJDStart || time > tJDStop)
	{
		return JPLEPH_ERR_TIME_OUTOFRANGE;
	}
	const double et2[2] = {time, 0.0};

	/* To hold x, xdot, y, ydot, z, zdot for all bodies */
	double  r_Helio[6], r_Helio_ecliptic[6];
	int err = fnEphcom4Interpolate(ctx, et2, 0, target, EPHCOM_SUN, r_Helio);

	if (err != 0)
	{
		std::cerr << "fnEphcom4Interpolate: An error occurred." << std::endl;
		return JPLEPH_ERR_NODATA;
	}

	/* convert ICRF equatorial to ecliptic using J2000.0 mean equinox */
	convert_ICRF_eq_to_ecliptic(r_Helio, r_Helio_ecliptic);

	/* copy rectangular coords */
	for (int j=0; j<3; ++j)
	{
		data->posXYZ[j] = r_Helio_ecliptic[j];
	}

	/* convert Heliocentric coords to spherical position */
	double  sph_coords_Helio[3], sph_coords_Helio_mean_equinox[3];
	xyz_to_sph(r_Helio_ecliptic, sph_coords_Helio);
	data->radius = sph_coords_Helio[0];
	data->longitude = sph_coords_Helio[1]; /* [0,360] deg */
	data->latitude = sph_coords_Helio[2];  /* [-90,90] deg */
	return JPLEPH_ERR_OK;
}

/*********************************************
 * REFERENCE ROUTINES (see jpl_frontend.c)   *
 * Applicable to any ephemeris               *
 *********************************************/
#define PI  3.14159265358979323846
#define DEG_PER_RAD  (180.0/PI)
#define RAD_PER_DEG  (PI/180.0)

void convert_ICRF_eq_to_ecliptic(const double *pos, double *pos_ecl)
{
	/* Mean obliquity of ecliptic at J2000.0 epoch */
	double epsilon0 = (23.0 + (26.0 + 21.4059/60) / 60) * RAD_PER_DEG;

	/* Seidelmann & Kovalevsky (2002) A&A 392, 341 */
	pos_ecl[0] =  pos[0]; // X coord unchanged
	pos_ecl[1] =  pos[1]*cos(epsilon0) + pos[2]*sin(epsilon0);
	pos_ecl[2] = -pos[1]*sin(epsilon0) + pos[2]*cos(epsilon0);
}

void convert_ICRF_eq_to_mean_ecliptic(double JDE, const double *pos, double *pos_ecl)
{
	/* FIXME: doesn't account for precession */

	double T = ((JDE - 2451545.0) / 36525); /* Julian centuries before/after 2000 */

	/* Mean obliquity of ecliptic at epoch of the date.
	   See Meeus, eqn. 22.2 */
	double epsilon0 = (23.0 + (26.0 + 21.448/60) / 60)
		- 46.8150*T / 3600
		- 0.00059*T*T / 3000
		+ 0.001813*T*T*T / 3600;

	epsilon0 *= RAD_PER_DEG;

	/* Seidelmann & Kovalevsky (2002) A&A 392, 341 */
	pos_ecl[0] =  pos[0]; // X coord unchanged
	pos_ecl[1] =  pos[1]*cos(epsilon0) + pos[2]*sin(epsilon0);
	pos_ecl[2] = -pos[1]*sin(epsilon0) + pos[2]*cos(epsilon0);
}

/* convert a position in 3d Cartesian space to spherical coordinates */
void  xyz_to_sph(double *pos, double *sph)
{
  double  x=pos[0], y=pos[1], z=pos[2];
  double  rho, theta, phi;

  rho = sqrt(x*x + y*y);
  if(rho > 0)
    phi = atan(z/rho);
  else
    phi = PI/2;

  theta = atan2(y, x);  /* take care of the sign ambiguity */
  if (theta < 0)
    theta += 2*PI;   /* theta is in the range [0, 2*PI] */

  /* convert to degrees */
  sph[0] = sqrt(rho*rho + z*z);
  sph[1] = theta * 180/PI;
  sph[2] = phi * 180/PI;
}


