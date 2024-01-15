/**********************************************************************

 EventsFinder.h
 January 2024 Thomas Morris

 Algorithms to find time of an event in an interval. Supports:
 - oppositions
 - solar longitude passages (seasons)
 - conjunctions between 2 planets
**********************************************************************/
#include "LibEphcom4.h"

namespace EventsFinder
{

/** Conjunctions between 2 planets/objects using ecliptic coords.
	 */
double conjunction_in_longitude(LibEphcom4& ephcom, double time);
	/** Method 1: Compute the longitude difference
	 *	At opposition the planetary longitude matches the Earth's.
	 *	For an eccentric orbit this is not necessarily the closest approach.
	 */
double delta_longitude(LibEphcom4& ephcom, double time);
	/** Method 2: Compute the longitude difference using geocentric coords.
	 *	Mars's longitude is 180 degrees from the Sun's.
	 *	NOTE: Uses the mean equinox of the date.
	 */
double delta_longitude_mean_equinox(LibEphcom4& ephcom, const double time);
double solar_longitude(LibEphcom4& ephcom, const double time);
	/** Compute the geocentric distance (in AU).
	 */
double geo_distance(LibEphcom4& ephcom, double time);

inline
double find_zero(double t1, double t2, double y1, double y2, const double f0,
				 LibEphcom4& ephcom,
				 double func(LibEphcom4&, const double) );
inline
double find_max_min(double t1, double t2,
				    LibEphcom4& ephcom,
				    double func(LibEphcom4&, double) );

double find_conjunction(LibEphcom4& ephcom, double t1, double t2, double y1, double y2);
	/** Bracketed time of opposition = [t1, t2]
	 * y1 > 0, y2 < 0, t2 > t1
	 */
double find_opposition(LibEphcom4& ephcom, double t1, double t2, double y1, double y2);
double find_equinox(LibEphcom4& ephcom, double t1, double t2, double y1, double y2, const double L0);
double find_closest_approach(LibEphcom4& ephcom, double t1, double t2);
} // namespace EventsFinder

