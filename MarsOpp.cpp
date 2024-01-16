/**********************************************************************
 
 MarsOpp.cpp
 April 2020 Thomas Morris
 
Compute oppositions of Mars from 2010 to 2030 using two passes.
1. scan the timespan in large steps to find events.
2. a gradient descent method to find time with high accuracy.
**********************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cassert>
#include <string>
#include <vector>

#include "EventsFinder.h"
#include "LibEphcom4.h"
#include "LReader.h"
#include "MarsDataLogger.h"

// Mars observables, Meeus "Astro Algorithms" ch.42
#include "MarsObserver.h"


using namespace std;
using namespace EventsFinder;

struct EphemDatum
{
	string			name, date, RA, Dec;
	double			rHelio, rGeo;
	double			deltaLong;
	string			tRise, tCulminate, tSet;
	/* Heliocentric ecliptic long/lat and solar longitude/declination */
	double			longHelio, latHelio, L_S, D_S;
	double          longitudeCM; /*!< longitude of central meridian */

	void setTime(double time);
};

class MarsOpp;

#include "MarsOppSeries.hpp"

string::size_type getWidth(const vector<string>& v);
vector<string> addFrame(const vector<string>& v);

void calcSolsticeDates2020(MarsOpp& mars2020);
void calculateDataFromConfig(MarsOpp& planetOpp, LConfig& conf);

void calendar_date_from_JD(double jdate, int *dmy, double *frac);
void test_calendar_date_from_JD(void);
double JD_from_calendar_date(double day, int month, int year);

double arc_length(double a1, double a2, double d1, double d2);

// Sets the date string from the ephemeris JD
void EphemDatum::setTime(double time)
{
	int dmy[3];
	double timeofday;
	/* Find calendar date/time from JD */
	calendar_date_from_JD(time, dmy, &timeofday);
	ostringstream oss;
	oss << dmy[0] << "/" << dmy[1] << "/" << dmy[2];
	date = oss.str();
}

#if defined(JPLEPH_VER1)
void test_VSOP87_comparison()
{
	/* See Meeus ch.32
	 * Note that the VSOP87 uses a dynamical (mean) ecliptic while DE405 uses
	 * the ICRF = J2000.0.
	 */
	struct jpl_PosDataVer2 position_Venus;
	struct jpl_PosDataVer2 position_Saturn;

	double JD1 = 2448976.5; /* 20th Dec 1992 0h DT */
	double JD2 = 2451178.5 + 207; /* 26th July 1999 0h DT */

	jpleph_init("/data/thomas/astro/astro_sara_pp/data/64");
	jpleph_compute_de405_heliopos(EPHCOM_VENUS, JD1, &position_Venus);
	printf("Longitude:   %12lf | 26.11428 deg (VSOP87)\n", position_Venus.longitude);
	printf("Latitude:    %12lf | -2.62070 deg (VSOP87)\n", position_Venus.latitude);
	printf("Radius (AU): %12lf |  0.724603 AU (VSOP87)\n",  position_Venus.radius);

	jpleph_compute_de405_heliopos(EPHCOM_SATURN, JD2, &position_Saturn);
	printf("Longitude:   %12lf | 39.9723901 deg (VSOP87, high precision)\n", position_Saturn.longitude);
}
#endif

void test_VSOP87_comparison_ephcom4()
{
	/* See Meeus ch.32
	 * Note that the VSOP87 uses a dynamical (mean) ecliptic while DE405 uses
	 * the ICRF version 2.0 = J2000.0.
	 */
	struct jpl_PosDataVer2 position_Venus;
	struct jpl_PosDataVer2 position_Saturn;

	double JD1 = 2448976.5; /* 20th Dec 1992 0h DT */
	double JD2 = 2451178.5 + 207; /* 26th July 1999 0h DT */

	// Create an EPHCOM3/4 context
	LibEphcom4 Ephcom4;
	Ephcom4.initContext(EPHCOM_VENUS);

	Ephcom4.computeHeliocentricPos(EPHCOM_VENUS, JD1, &position_Venus);
	cout << "** VENUS **" << endl;
	printf("Longitude:   %12lf | 26.11428 deg (VSOP87)\n", position_Venus.longitude);
	printf("Latitude:    %12lf | -2.62070 deg (VSOP87)\n", position_Venus.latitude);
	printf("Radius (AU): %12lf |  0.724603 AU (VSOP87)\n",  position_Venus.radius);

	Ephcom4.setObject(EPHCOM_SATURN);
	Ephcom4.computeHeliocentricPos(EPHCOM_SATURN, JD2, &position_Saturn);
	cout << "** SATURN **" << endl;
	printf("Longitude:   %12lf | 39.9723901 deg (VSOP87, high precision)\n", position_Saturn.longitude);

	Ephcom4.deinitContext();
}


/////////////////////////////////////////////////////////////////////
class MarsOpp
{
public:
	MarsOpp()
	{
		Ephcom.initContext(EPHCOM_MARS);
		MyPrint.Flags = Printing::PRINT_FLAG_LIST;
		ephcomErrStat = 0;
	}

	~MarsOpp()
	{
		Ephcom.deinitContext();	
	}

	void setObject(ephcom_object object, ephcom_object object2) {
		EphcomID = object;
		Ephcom.setObject(object);
		SecondaryID = object2;
		Ephcom.setSecondaryObject(object2);
	}

	void generateEphem(double timeStart, double timeInterval, int numTimes, bool checkMethod);
	void computeEphem(double timeOpp, bool prettyPrint = false);
	double angularSize(double delta);
	void boxPrint(EphemDatum& data);
	void listPrint(double time, double distanceHelio, double delta);

private:
	LibEphcom4 Ephcom;
	ephcom_object EphcomID;
	ephcom_object SecondaryID;
	MarsDataLogger dataLog;
	int ephcomErrStat;

	/* modal flags */
	static const int MODE_OPPOSITION = 1; /*!< Find oppositions in date range */
	static const int MODE_CLOSEST = 2;    /*!< TODO: Find closest approaches (perihelia) in date range */
	static const int MODE_EQUINOX = 4;    /*!< TODO: Find equinoxes in date range */
	static const int MODE_CONJUNCTION = 8;/*!< Find conjunctions between two planets in date range */
	int mode = MODE_CONJUNCTION;

	struct Printing
	{
		static const int PRINT_FLAG_LIST = 1;
		static const int PRINT_FLAG_BOX = 2;
		int Flags;
		bool hasFlag(int flag) {
			return Flags & flag;
		}
		string Prefix;
	} MyPrint;
};
/////////////////////////////////////////////////////////////////////


void MarsOpp::generateEphem(double timeStart, double timeInterval, int numTimes, bool checkMethod)
{
	int			idxList = 0, i;
	double		time;

	/* iterate over number of time samples */
	double lastDelta = 0, lastDeltaConj = 0;
	double lastLong = 0, solarLong;
	for(i=0; i<numTimes; i++)
	{
		/* calculate the current time */
		time = timeStart + i*timeInterval;

		if (ephcomErrStat != 0) {
		  cerr << "An error occurred during the JPL ephemeris calculation: " 
		    << "iter = " << i << ", time = " << time
		    << endl;
		  break;
		}

	if (mode & MODE_OPPOSITION)
	{
		/* Look for opposition */
		double deltaLong = EventsFinder::delta_longitude(Ephcom, time);
		//cout << std::setprecision (12) << time << ", " << deltaLong << endl;
		if (-30 < deltaLong && deltaLong < 0 && lastDelta > 0)
		{
			/* Earth has surpassed Mars => opposition occurred */
			double timeOpp = find_opposition(Ephcom, time - timeInterval, time, lastDelta, deltaLong);
			MyPrint.Prefix = string("[OPP] "); // mark 'opposition'
			computeEphem(timeOpp);
			MyPrint.Prefix = string(""); // reset label

			if (checkMethod)
			{
				cout << endl;

				/* compare with method in Meeus ch.36 */
				MarsOppSeries mySeries;
				//JupiterOppSeries mySeries;
				mySeries.find_nearest_opposition(timeOpp);
			}
		}
		lastDelta = deltaLong;
	}

	if (mode & MODE_CONJUNCTION)
	{
		/* Look for conjunctions.
		 * Case (2) below finds conjunctions during retrograde motion of the
		 * interior planet, i.e. triple conjunctions.
		 */
		double deltaLong = EventsFinder::conjunction_in_longitude(Ephcom, time);
		if ((deltaLong < 15 && deltaLong > 0 && lastDeltaConj < 0) ||
		    (deltaLong > -15 && deltaLong < 0 && lastDeltaConj > 0))
		{
			double timeOpp = find_conjunction(Ephcom, time - timeInterval, time, lastDeltaConj, deltaLong);
			MyPrint.Prefix = string("[CONJ]"); // mark 'conjunction'
			computeEphem(timeOpp);
			MyPrint.Prefix = string(""); // reset label
		}
		lastDeltaConj = deltaLong;
	}

	if (mode & MODE_CLOSEST)
	{
		/* Look for close approach. Needs two historical values */
		double delta0 = EventsFinder::geo_distance(Ephcom, time - timeInterval);
		double delta1 = EventsFinder::geo_distance(Ephcom, time);
		double delta2 = EventsFinder::geo_distance(Ephcom, time + timeInterval);
		if (delta0 > delta1 && delta1 < delta2)
		{
			/* Minima found. Search interval. */
			double timeCloseApproach = find_closest_approach(Ephcom, time - timeInterval/2, time + timeInterval/2);
			MyPrint.Prefix = string("[CLO] "); // mark 'close approach'
			computeEphem(timeCloseApproach);
			MyPrint.Prefix = string(""); // reset label
		}
	}

	if (mode & MODE_EQUINOX)
	{
		double solarLong = EventsFinder::solar_longitude(Ephcom, time);
		//cout << std::setprecision (12) << time << ", " << solarLong << endl;
		if ((lastLong > 345 && solarLong > 0 && solarLong < 15) ||
			(lastLong < 90 && solarLong > 90 && solarLong < 105) ||
			(lastLong < 180 && solarLong > 180 && solarLong < 195) ||
			(lastLong < 270 && solarLong > 270 && solarLong < 285))
		{
			double timeEquinox;
			if (lastLong > 345 && solarLong > 0 && solarLong < 15)
			{
				/* Spring equinox occurred */
				lastLong -= 360;
				timeEquinox = find_equinox(Ephcom, time - timeInterval, time, lastLong, solarLong, 0);
				MyPrint.Prefix = string("[SPR] ");
			}
			else if (lastLong < 90 && solarLong > 90 && solarLong < 105)
			{
				/* Summer solstice occurred */
				timeEquinox = find_equinox(Ephcom, time - timeInterval, time, lastLong, solarLong, 90);
				MyPrint.Prefix = string("[SUM] ");
			}
			else if (lastLong < 180 && solarLong > 180 && solarLong < 195)
			{
				/* Autumn equinox occurred */
				timeEquinox = find_equinox(Ephcom, time - timeInterval, time, lastLong, solarLong, 180);
				MyPrint.Prefix = string("[AUT] ");
			}
			else if (lastLong < 270 && solarLong > 270 && solarLong < 285)
			{
				/* Winter solstice occurred */
				timeEquinox = find_equinox(Ephcom, time - timeInterval, time, lastLong, solarLong, 270);
				MyPrint.Prefix = string("[WIN] ");
			}
			computeEphem(timeEquinox, false);
			MyPrint.Prefix = string(""); // reset label
		}
		lastLong = solarLong;
	}
	} /* next time */

	dataLog.writeToCSVFile("data_logger.csv");
}

void MarsOpp::computeEphem(double timeOpp, bool prettyPrint)
{
	/*
	 * Main 'Compute ephemeris' function.
	 * Calculates heliocentric and geocentric equatorial coords
	 * at the specified time.
	 * Based on 'jpleph_compute_de405' in the libJPLEPH ver1.
	 */
	jpl_PosData observables;
	if (Ephcom.computeDE(EphcomID, timeOpp, &observables) != 0)
	{
		ephcomErrStat = 1;
	}

	/*
	 * Calculates geocentric ecliptic coords at the specified
	 * time and J2000.0 equinox.
	 */
	jpl_PosDataVer2 geoData;
	if (Ephcom.computeGeocentricPos(EphcomID, timeOpp, &geoData) != 0)
	{
		ephcomErrStat = 1;
	}
	jpl_PosDataVer2 geoData2;
	Ephcom.computeGeocentricPos(Ephcom.getSecondaryObject(), timeOpp, &geoData2);

	if (MyPrint.hasFlag(Printing::PRINT_FLAG_LIST))
	{
		listPrint(timeOpp, observables.dist_Helio, geoData.radius);
	}

	/*
	 * Calculates heliocentric ecliptic coords at the specified
	 * time and J2000.0 equinox.
	 */
	jpl_PosDataVer2 helioData;
	jpl_PosDataVer2 posSun;
	if (Ephcom.computeHeliocentricPos(EphcomID, timeOpp, &helioData) != 0)
	{
		ephcomErrStat = 1;
	}

	if (EphcomID == EPHCOM_MARS)
	{
		calculateAreographicSolarCoords(helioData, timeOpp, posSun);
	}
	else
	{
		// TODO: solar coords for other planets
		posSun.longitude = 0;
		posSun.latitude = 0;
	}

	if (prettyPrint || MyPrint.hasFlag(Printing::PRINT_FLAG_BOX))
	{
		EphemDatum	data;

		/* Populate 'datum' structure and print nicely */
		data.name = string("Mars");
		data.setTime(timeOpp);
		data.rHelio = observables.dist_Helio;
		data.rGeo = geoData.radius;
		if (EphcomID == EPHCOM_MARS)
		{
			data.longitudeCM = longitudeCM(geoData, timeOpp, calcLightDelay(geoData.radius));
		}

		/* Heliocentric ecliptic long/lat and solar longitude/declination */
		if (EphcomID != EPHCOM_SUN)
		{
			data.longHelio = helioData.longitude;
			data.latHelio = helioData.latitude;
			data.L_S = posSun.longitude;
			data.D_S = posSun.latitude;
		}
		else
		{
			data.longHelio = geoData.longitude;
			data.latHelio = geoData.latitude;
		}

		boxPrint(data);
	}

	{
		/* CSV file */
		MarsData item;
		item.JDay   = timeOpp;
		calendar_date_from_JD(timeOpp, item.dmy, &item.timeofday);
		item.Delta  = geoData.radius;
		item.radius = observables.dist_Helio;
		item.diameter = angularSize(item.Delta);
		item.separation =
		  60 * arc_length(geoData.longitude, geoData2.longitude, geoData.latitude, geoData2.latitude);
		item.omega = longitudeCM(geoData, timeOpp, calcLightDelay(geoData.radius));
		item.L_S = posSun.longitude;
		item.D_S = posSun.latitude;
		dataLog.add(item);
	}
}

double MarsOpp::angularSize(double delta)
{
	switch (EphcomID)
	{
		case EPHCOM_MERCURY:
			return 6.73/delta;
		case EPHCOM_VENUS:
			return 16.69/delta;
		case EPHCOM_MARS:
			return 9.36/delta;
		case EPHCOM_JUPITER:
			return 197.15/delta;
		case EPHCOM_SATURN:
			return 166.19/delta;
		case EPHCOM_URANUS:
			return 70.48/delta;
		case EPHCOM_NEPTUNE:
			return 68.29/delta;
		case EPHCOM_PLUTO:
			return 3.295/delta;
		default:
			break;
	}
	return 0;
}

void MarsOpp::listPrint(double time, double distanceHelio, double delta)
{
	int			dmy[3];
	double		timeofday;

	/* Find calendar date/time from JD */
	calendar_date_from_JD(time, dmy, &timeofday);

	char		buffer[64], bufDateTime[64];
	sprintf(bufDateTime, "%02d/%02d/%04d %02d:%02d",
			dmy[0], dmy[1], dmy[2],
			(int) floor(timeofday*24),
			(int) round(60*(timeofday*24 - floor(timeofday*24))));

	sprintf(buffer, "%14.5lf", time);
	cout << MyPrint.Prefix << buffer << ", " <<
		bufDateTime << ", " <<
		std::setprecision (12) << distanceHelio << " AU, " <<
		std::setprecision (12) << delta << " AU, " <<
		std::setprecision (6) << angularSize(delta) << " arcsec" <<
		endl;
}

void MarsOpp::boxPrint(EphemDatum& data)
{
	vector<string> info;
	const int data_size = 9;
	ostringstream oss[data_size];

	oss[0] << data.name << " on " << data.date;
	oss[1] << "Solar distance = " << std::fixed << std::setprecision (6) << data.rHelio << " AU.";
	oss[2] << "Helio long. L =  " << std::fixed << std::setprecision (6) << data.longHelio << " deg.";
	oss[3] << "Helio lat. b =   " << std::fixed << std::setprecision (6) << data.latHelio << " deg.";
	oss[4] << "Geo distance =   " << std::fixed << std::setprecision (6) << data.rGeo << " AU.";
	oss[5] << "Disk diameter =  " << std::fixed << std::setprecision (3) << angularSize(data.rGeo) << " arcsec.";
	oss[6] << "Solar long. L_S= " << std::fixed << std::setprecision (6) << data.L_S << " deg.";
	oss[7] << "Solar dec. D_S = " << std::fixed << std::setprecision (6) << data.D_S << " deg.";
	oss[8] << "CM long.       = " << std::fixed << std::setprecision (4) << data.longitudeCM << " deg.";

	for(int i=0; i<data_size; ++i)
	{
		info.push_back(oss[i].str());
	}

	vector<string> framed = addFrame(info);
	for (vector<string>::const_iterator it = framed.begin(); it != framed.end(); ++it)
	{
		cout << *it << endl;
	}
}

/*
 * Pretty print/frame helper functions
 */
string::size_type getWidth(const vector<string>& v)
{
	string::size_type maxlen = 0;
	for (vector<string>::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		maxlen = max(maxlen, it->size() );
	}
	return maxlen;
}

vector<string> addFrame(const vector<string>& v)
{
	vector<string> result;
	string::size_type width = max(32UL, getWidth(v));
	string border(width + 4, '*');

	// top border
	result.push_back(border);

	// title row
	result.push_back("***" + v[0] + string(width - v[0].size(), '*') + "*");

	// interior rows, bordered by an asterisk and a space
	for (vector<string>::size_type i = 1; i != v.size(); ++i)
	{
		result.push_back("* " + v[i] + string(width - v[i].size(), ' ') + " *");
	}

	// bottom border
	result.push_back(border);
	return result;
}
/*************************************/


void calcSolsticeDates2020(MarsOpp& mars2020)
{
	// BAA Handbook 2020 opposition (interpolated from nearest data)
	mars2020.computeEphem(JD_from_calendar_date( 8.0, 4, 2020), true); // L_S = 180
	mars2020.computeEphem(JD_from_calendar_date( 1.9, 9, 2020), true); // L_S = 270
	mars2020.computeEphem(JD_from_calendar_date( 6.8, 2, 2021), true); // L_S = 0
	mars2020.computeEphem(JD_from_calendar_date(24.2, 8, 2021), true); // L_S = 90
}

void calcOppData2014()
{
	// BAA 2014 opposition data, see JBAA
	MarsOpp mars2014;

	mars2014.computeEphem(JD_from_calendar_date(18, 4, 2013), true);
	mars2014.computeEphem(JD_from_calendar_date(31, 7, 2013), true);
	mars2014.computeEphem(JD_from_calendar_date( 1, 1, 2014), true);
	mars2014.computeEphem(JD_from_calendar_date(15, 2, 2014), true);
	mars2014.computeEphem(JD_from_calendar_date(17, 8, 2014), true);
	mars2014.computeEphem(JD_from_calendar_date(11, 1, 2015), true);
}

void calculateDataFromConfig(MarsOpp& planetOpp, LConfig& conf)
{
	double JDstart;
	int numSteps  = 180;
	const int stepJDays = 7;

	// Set the ephemeris/opposition search parameters
	planetOpp.setObject( (ephcom_object) conf.ephcomID,
	                     (ephcom_object) conf.secondaryID );

	if (conf.snapshotMode)
	{
		JDstart = JD_from_calendar_date(conf.Tsnapshot[0] + conf.TFsnapshot,
				conf.Tsnapshot[1], conf.Tsnapshot[2]);
		numSteps = 1;
	}
	else
	{
		JDstart = JD_from_calendar_date(conf.startDate[0], conf.startDate[1], conf.startDate[2]);
		double JDend = JD_from_calendar_date(conf.endDate[0], conf.endDate[1], conf.endDate[2]);
		numSteps = (JDend - JDstart) / stepJDays;
	}

	if (numSteps > 1)
	{
		// Opposition search function
		planetOpp.generateEphem(JDstart, stepJDays, numSteps, false);
	}
	else
	{
		planetOpp.computeEphem(JDstart, true);
	}
}

int main(int argc, char **argv)
{
	MarsOpp mars2020;
	LReader luaReader;
	LConfig conf;
	double JDstart;

	if (argc == 2)
	{
		luaReader.open(argv[1]);
		luaReader.read(conf);
		luaReader.close();
		calculateDataFromConfig(mars2020, conf);
	}
	else
	{
		// Find next Mars opposition after today's date
		mars2020.setObject(EPHCOM_MARS, EPHCOM_MERCURY /* none */ );
		time_t today = time(NULL);
		tm today_d = *gmtime(&today);
		JDstart = JD_from_calendar_date(today_d.tm_mday, today_d.tm_mon + 1, today_d.tm_year + 1900);
		mars2020.computeEphem(JDstart, true);
	}

	//calcOppData2014();
	//calcSolsticeDates2020(mars2020);

	//compute_planet_longitudes(2458939.5, 7, 10000, 0);

	//MarsOppSeries mars2020;
	//for(int j=0; j<10; j++)
	//	(void)mars2020.compute_next_opposition(j);
	//mars2020.compute_next_opposition(2729.0);

	//test_VSOP87_comparison_ephcom4();

	return 0;
}

void calendar_date_from_JD(double jdate, int *dmy, double *frac)
{
	int a, b, c, d, e, z, alpha;
	assert(jdate > 0);

	/* Meeus p.63 */
	z = (int) floor(jdate + 0.5);
	*frac = jdate + 0.5 - z;

	if (z < 2299161)
	{
		a = z;
	}
	else
	{
		alpha = (int) floor((z - 1867216.25) / 36524.25);
		a = z + 1 + alpha - alpha/4;
	}

	b = a + 1524;
	c = (int) floor((b - 122.1) / 365.25);
	d = (int) floor(365.25*c);
	e = (int) floor((b - d) / 30.6001);

	/* day, month, year */
	dmy[0] = b - d - (int) floor(30.6001*e);
	dmy[1] = (e < 14) ? (e - 1) : (e - 13);
	dmy[2] = (dmy[1] > 2) ? (c - 4716) : (c - 4715);
}

void test_calendar_date_from_JD(void)
{
	int dmy[3];
	double f;

	calendar_date_from_JD(2436116.31, dmy, &f);
	printf("%4d-%02d-%6lf | 1957-10-4.81\n", dmy[2], dmy[1], dmy[0] + f);

	calendar_date_from_JD(1842713.0, dmy, &f);
	printf("%4d-%02d-%6lf |  333-01-27.5\n", dmy[2], dmy[1], dmy[0] + f);

	calendar_date_from_JD(1507900.13, dmy, &f);
	printf("%4d-%02d-%6lf | -584-05-28.63\n", dmy[2], dmy[1], dmy[0] + f);

	double Sputnik_launch = JD_from_calendar_date(4.81, 10, 1957);
	printf("%10.2lf |  2436116.31\n", Sputnik_launch);
}

double JD_from_calendar_date(double day, int month, int year)
{
	if (month == 1 || month == 2)
	{
		year--; month += 12;
	}

	/* Assume Gregorian calendar */
	int a = year/100;
	int b = 2 - a + (a/4);
	double JD = floor(365.25*(year+4716)) + floor(30.601*(month+1)) + day + b - 1524.5;

	return JD;
}

static inline
double rad2deg(double radians)
{
	return 180.0 / M_PI * radians;
}

static inline
double deg2rad(double degrees)
{
	return M_PI / 180.0 * degrees;
}

/* 'Haversine' function */
static inline
double hav(double theta)
{
	double s = sin(theta/2);
	return s*s;
}

/* Inverse 'Haversine' function */
static inline
double ahav(double h)
{
	return 2 * asin(sqrt(h));
}

double arc_length(double a1, double a2, double d1, double d2)
{
	/* See Meeus, equation 17.5 */
	double x = hav( deg2rad(d1-d2) ) + cos(deg2rad(d1)) * cos(deg2rad(d2)) * hav(deg2rad(a1 - a2));
	return rad2deg(ahav(x));
}

