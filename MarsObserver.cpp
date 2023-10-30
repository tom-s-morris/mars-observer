/* 
 * MarsObserver.cpp
 * Calculates observable parameters for Mars:
 * - Areocentric longitude/declination of Sun (seasons)
 *
 * Uses Alan Irwin's 'timeephem' library as at 1st July 2020.
 * Formulae based on Meeus, Ch. 42
 *
 * Test values, interpolated from BAA Handbook 2020
 * Ls = 180, 8.0 Apr 2020  -- N autumn equinox
 * Ls = 270, 1.9 Sept 2020 -- N winter solstice, S summer solstice
 * Ls = 0,   6.8 Feb 2021  -- N spring equinox
 * Ls = 90,  24.2 Aug 2021 -- N summer solstice, S winter equinox
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include "MarsObserver.h"

using namespace std;


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

/*
 * Light travel time in days per AU
 */
double calcLightDelay(double delta)
{
	const double millionKmInAU = 149.597870700; /* 1AU in million km */
	const double cInMillionKmPerSec = 0.299792458;

	return delta * millionKmInAU / cInMillionKmPerSec / 86400;
}

double declination(double lambda, double beta, double lambda_N, double beta_N, double T)
{
	double x= -sin(beta_N)*sin(beta) -cos(beta_N)*cos(beta)*cos(lambda_N-lambda);
	return rad2deg(asin(x));
}

double solarLongitude(jpl_PosDataVer2& helioData, double time)
{
	double T = (time - 2451545.0) / 36525.0; // Modified JD in centuries

	/* Longitude of the ascending node/Spring equinox */
	double L_node   = deg2rad( 49.5581 + 0.7721*T);
	double L_Sp     = deg2rad(251);
	double L_S;

	/* Transformation matrix for
	 * heliocentric ecliptic --> Martian orbital plane in J2000.0 system */
	const double Q[3][3] =
		{ {    0.913454,   -0.405761,   -0.030970 },
		  {    0.406198,    0.913739,    0.009152 },
		  {    0.024585,   -0.020939,    0.999478 } };

	/* Radius not used in angle calculation */
	const double x[3] = {cos(deg2rad(helioData.longitude)) * cos(deg2rad(helioData.latitude)),
	                     sin(deg2rad(helioData.longitude)) * cos(deg2rad(helioData.latitude)),
	                     sin(deg2rad(helioData.latitude)) };
	double y[3];

	for (int i=0; i<3; i++)
	{
		y[i] = 0;
		for (int j=0; j<3; j++)
		{
			y[i] += Q[i][j] * x[j];
		}
	}

	L_S = atan2(y[1], y[0]) + L_Sp;
	if (L_S > 2*M_PI) L_S -= 2*M_PI;

	return rad2deg(L_S);
}

double longitudeCM(jpl_PosDataVer2& geoData, double time, double tLight)
{
	/* Rotational phase of Mars = hour angle (in degrees), p.289 */
	double W = 11.504 + 350.89200025 * (time - tLight - 2433282.5);

	/* Geocentric longitude of Mars in its orbit, in Mars's equatorial plane. 
	 * This is the correction due to the changing viewing angle from Earth.
	 * For maximum accuracy, the calculation would need to account for the
	 * position of the observer on Earth's surface.
	 * See Ch. 42 pp.289-290
	 */
	/* Mean obliquity of ecliptic at J2000.0 epoch [should use eqn (22.2)] */
	double epsilon0 = deg2rad( (23.0 + (26.0 + 21.4059/60) / 60) );

	double u, v, x = geoData.posXYZ[0], y = geoData.posXYZ[1], z = geoData.posXYZ[2];
	u = y * cos(epsilon0) - z * sin(epsilon0);
	v = y * sin(epsilon0) + z * cos(epsilon0);

	/* R.A. and declination of Martian North pole */
	double alpha0 = deg2rad( 317.681 );
	double delta0 = deg2rad(  52.886 );
	/* Convert equatorial coords to local equatorial coords on Mars */
	double alpha = atan2(u, x);
	double delta = atan2(v, sqrt(x*x + u*u));
	double zeta =  atan2(sin(delta0)*cos(delta)*cos(alpha0 - alpha) - sin(delta)*cos(delta0),
		cos(delta)*sin(alpha0 - alpha) );

	printf("(alpha, delta, zeta) = (%f, %f, %f)\n",
		rad2deg(alpha), rad2deg(delta), rad2deg(zeta));
	double omega = fmod(W - rad2deg(zeta), 360.0);

	return omega;
}

/*
 * Computes the areographic (Mars-centred) coordinates of the Sun.
 * This information indicates the season on Mars.
 */
void calculateAreographicSolarCoords(jpl_PosDataVer2& helioData, double time, struct jpl_PosDataVer2& helioData2)
{
	/* Position of Mars's north pole in ecliptic coordinates */
	double T = (time - 2451545.0) / 36525.0; // Modified JD in centuries
	double lambda_N = deg2rad(352.9065 + 1.17330*T);
	double beta_N   = deg2rad( 63.2818 - 0.00394*T);

	/* Longitude of the ascending node/Spring equinox */
	double L_node   = deg2rad( 49.5581 + 0.7721*T);

	/* Planetocentric declination of the Sun */
	double D_S = declination(deg2rad(helioData.longitude),
	                         deg2rad(helioData.latitude),
	                         lambda_N, beta_N, T);

	/* Planetocentric longitude of the Sun
	 * See https://www.giss.nasa.gov/tools/mars24/help/notes.html
	 */
	//double L_S= asin(sin(D_S)/cos(beta_N));

	double L_S = solarLongitude(helioData, time);

	helioData2.longitude = L_S;
	helioData2.latitude  = D_S;
	helioData2.radius    = helioData.radius;

	string marsSeason;
	if (helioData2.longitude > 0 && helioData2.longitude <= 90)
	{
		marsSeason = string("Northern hemisphere spring.");
	}
	else if (helioData2.longitude <= 180 )
	{
		marsSeason = string("Northern hemisphere summer.");
	}
	else if (helioData2.longitude <= 270)
	{
		marsSeason = string("Northern hemisphere autumn.");
	}
	else
	{
		marsSeason = string("Northern hemisphere winter.");
	}

	if (helioData2.longitude > 180 && helioData2.longitude < 360)
	{
		marsSeason += string(" Dust season active.");
	}

	cout << "\t" << "L_S = " << helioData2.longitude
	          << ", " << marsSeason << endl;
	cout << "\t" << "D_S = " << helioData2.latitude << endl;
}

