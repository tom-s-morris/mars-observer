/**********************************************************************

 EventsFinder.cpp
 January 2024 Thomas Morris

 Algorithms to find time of an event in an interval. Supports:
 - oppositions
 - solar longitude passages (seasons)
 - conjunctions between 2 planets
**********************************************************************/
#include <iostream>
#include <iomanip>
#include <cmath>
#include "LibEphcom4.h"
#include "EventsFinder.h"
#include "MarsObserver.h"


	/** Conjunctions between 2 planets/objects using ecliptic coords.
	 */
double EventsFinder::conjunction_in_longitude(LibEphcom4& ephcom, double time)
{
	struct jpl_PosDataVer2 position1, position2;
	double delta;
	ephcom.computeGeocentricPos(ephcom.getObject(),          time, &position1);
	ephcom.computeGeocentricPos(ephcom.getSecondaryObject(), time, &position2);
	delta = position1.longitude - position2.longitude;
	if (delta >  180) delta -= 360;
	if (delta < -180) delta += 360;
	return delta;
}
	/** Method 1: Compute the longitude difference
	 *	At opposition the planetary longitude matches the Earth's.
	 *	For an eccentric orbit this is not necessarily the closest approach.
	 */
double EventsFinder::delta_longitude(LibEphcom4& ephcom, double time)
{
	struct jpl_PosDataVer2 position_Earth, position_Mars;

	ephcom.computeHeliocentricPos(EPHCOM_EARTH,       time, &position_Earth);
	ephcom.computeHeliocentricPos(ephcom.getObject(), time, &position_Mars);
	return position_Mars.longitude - position_Earth.longitude;
}

	/** Method 2: Compute the longitude difference using geocentric coords.
	 *	Mars's longitude is 180 degrees from the Sun's.
	 *	NOTE: Uses the mean equinox of the date.
	 */
double EventsFinder::delta_longitude_mean_equinox(LibEphcom4& ephcom, const double time)
{
	struct jpl_PosDataVer2 position_Sun, position_Mars;
	double delta;

	ephcom.computeGeocentricPos(EPHCOM_SUN,         time, &position_Sun);
	ephcom.computeGeocentricPos(ephcom.getObject(), time, &position_Mars);
	delta = position_Mars.longitude - position_Sun.longitude + 180;
	if (delta >  180) delta -= 360;
	if (delta < -180) delta += 360;
	return delta;
}

double EventsFinder::solar_longitude(LibEphcom4& ephcom, const double time)
{
	struct jpl_PosDataVer2 position_Mars;
	double longitude;

	ephcom.computeHeliocentricPos(EPHCOM_MARS, time, &position_Mars);
	longitude = solarLongitude(position_Mars, time);  /* Ls */
	return longitude;
}

	/** Compute the geocentric distance (in AU).
	 */
double EventsFinder::geo_distance(LibEphcom4& ephcom, double time)
{
	struct jpl_PosDataVer2 position_Mars;

	ephcom.computeGeocentricPos(ephcom.getObject(), time, &position_Mars);
	return position_Mars.radius;
}

	/**	Utility function. Finds root of f(t) = f0 in specified
	 *	interval [t1, t2], y1=f(t1) changes sign to y2=f(t2)
	 */
inline
double EventsFinder::find_zero(double t1, double t2, double y1, double y2, const double f0,
				 LibEphcom4& ephcom,
				 double func(LibEphcom4&, const double) )
{
	double t3 = 0;

	y1 -= f0; y2 -= f0;
	for (int i=0; i<10; i++)
	{
		/* estimate the slope */
		double m = (y2 - y1) / (t2 - t1);
		t3 = t2 - (y2 / m);  /* new estimate of time of opposition, t1 < t3 < t2 */

		if (t1 > t3 || t3 > t2)
		{
			/* Linear approx is no longer valid */
			std::cerr << "Error: Convergence problems, giving up!" << std::endl;
			std::cerr << std::setprecision(12) << t3 << " is not contained in [" << t1 << "," << t2 << "]" << std::endl;
			break;
		}

		double delta = func(ephcom, t3) - f0;
		if (fabs(delta) < 1e-4 || fabs(t2 - t1) < 1e-3)
		{
			//std::cout << "Converged after " << i << " iterations." << std::endl;
			break;
		}

		if (((delta > 0) && (y1 > 0))
			|| ((delta < 0) && (y1 < 0)))
			{ t1 = t3; y1 = delta; } /* replace lower end-point */
		else
			{ t2 = t3; y2 = delta; } /* replace upper end-point */
	}

	return t3;
}

	/**	Utility function. Finds a minima/maxima of f'(t) = 0 in specified
	 *	interval [t1, t2] using the Golden Section.
	 */
inline
double EventsFinder::find_max_min(double t1, double t2,
				    LibEphcom4& ephcom,
				    double func(LibEphcom4&, double) )
{
	double t3 = 0, t4 = 0;
	const double gr = (sqrt(5) + 1) / 2;

	for (int i=0; i<50; i++)
	{
		// sample points based on Golden section
		t3 = t2 - (t2 - t1) / gr;
		t4 = t1 + (t2 - t1) / gr;

		double y3 = func(ephcom, t3);
		double y4 = func(ephcom, t4);

		if (y3 < y4) {
			t2 = t4;
		} else {
			t1 = t3;
		}

		if ((t2 - t1) < 1e-3) {
			//std::cout << "Converged after " << i << " iterations." << std::endl;
			break;
		}
	}

	return (t1 + t2) / 2.0;
}

double EventsFinder::find_conjunction(LibEphcom4& ephcom, double t1, double t2, double y1, double y2)
{
	return EventsFinder::find_zero(t1, t2, y1, y2, 0, ephcom, EventsFinder::conjunction_in_longitude);
}

double EventsFinder::find_opposition(LibEphcom4& ephcom, double t1, double t2, double y1, double y2)
{
	return EventsFinder::find_zero(t1, t2, y1, y2, 0, ephcom, EventsFinder::delta_longitude_mean_equinox);
}

double EventsFinder::find_equinox(LibEphcom4& ephcom, double t1, double t2, double y1, double y2, const double L0)
{
	/* TODO: check discontinuity around Ls = 360 degrees */
	return EventsFinder::find_zero(t1, t2, y1, y2, L0, ephcom, EventsFinder::solar_longitude);
}

double EventsFinder::find_closest_approach(LibEphcom4& ephcom, double t1, double t2)
{
	return EventsFinder::find_max_min(t1, t2, ephcom, EventsFinder::geo_distance);
}


