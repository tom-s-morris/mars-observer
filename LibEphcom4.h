/* 
 * LibEphcom4.h
 * Dynamically loads the ephcom4 library and provides geocentric
 * positions used in the Mars opposition calculator.
 * Uses Alan Irwin's 'timeephem' library as at 1st July 2020.
 */

#ifndef LIB_EPHCOM4_H
#define LIB_EPHCOM4_H
#include <mutex>

extern "C" {
#include "ephcomc.h"
}


struct  DegMinSec
{
	int  deg, min;
	double  sec;
};
struct  HourMinSec
{
	int  hour, min;
	double  sec;
};

struct jpl_PosData
{
	struct HourMinSec ra;	/* right ascension */
	struct DegMinSec dec;	/* declination */
	double dist_Helio, dist_Geo;
};

struct jpl_PosDataVer2
{
	double posXYZ[3];	/* rectangular coords (AU) */
	double radius;		/* radius (AU) */
	double longitude, latitude; /* long, lat (deg) */
};

#define JPLEPH_ERR_OK				(0)
#define JPLEPH_ERR_NODATA			(1)
#define JPLEPH_ERR_TIME_OUTOFRANGE	(2)
#define JPLEPH_ERR_TARGET_OUTOFRANGE (3)


class LibEphcom4
{
private:
	ephcom_context *ctx;
	std::mutex ctxm;
	void *libDLL;
	double tJDStart, tJDStop, tJDStep;
	ephcom_object targetID;
	ephcom_object secondaryID;
	bool equinoxJ2000;

	/* Internal API */
	ephcom_context* (*fnEphcom4Open)( const char * );
	int (*fnEphcom4Close)( ephcom_context * );
	int (*fnEphcom4GetCoords)( ephcom_context *, ephcom_coords * );
	int (*fnEphcom4Interpolate)( ephcom_context *, const double *, unsigned, int, int, double * );
	int (*fnEphcom4Alloc2dChar)( char ***, unsigned, unsigned );
	int (*fnEphcom4Free2dChar)( char **, unsigned );
	int (*fnEphcom4ReadConstants)( ephcom_context *, char **, unsigned, unsigned, double *values, unsigned *numde, double *au, double *emrat, double *sss );
public:
	LibEphcom4(): ctxm() {}	// Ensure the mutex is initialised to an unlocked state
	int loadLib(void);
	void initContext(ephcom_object ephcomID);
	void deinitContext(void);

	/* Set target planet for subsequent computations */
	int setObject(ephcom_object ephcomID);
	ephcom_object getObject(void);
	/* Set secondary planet for subsequent computations of e.g. conjunctions */
	int setSecondaryObject(ephcom_object ephcomID);
	ephcom_object getSecondaryObject(void);
	void setMeanEquinox();

	/* Combined geo+helio ephemeris data, as for original jpleph_compute_de405() */
	int computeDE(int target, const double time, jpl_PosData *);

	/* Geocentric position */
	int computeGeocentricPos(int target, const double time, jpl_PosDataVer2 *);

	/* Heliocentric position, for comparison to VSOP87 theory */
	int computeHeliocentricPos(int target, const double time, jpl_PosDataVer2 *);
};

#endif
