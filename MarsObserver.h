/* 
 * MarsObserver.h
 * Calculates observable parameters for Mars:
 * - Areocentric longitude/declination of Sun (seasons)
 *
 * Uses Alan Irwin's 'timeephem' library as at 1st July 2020.
 * Formulae based on Meeus, Ch. 42
 */

#include "LibEphcom4.h"

/*
 * Light travel time in days per AU
 */
double calcLightDelay(double delta);

double longitudeCM(jpl_PosDataVer2& geoData, double time, double tLight);

double declination(double lambda, double beta, double lambda_N, double beta_N, double T);

double solarLongitude(jpl_PosDataVer2& helioData, double time);

/*
 * Computes the areographic (Mars-centred) coordinates of the Sun.
 * This information indicates the season on Mars.
 */
void calculateAreographicSolarCoords(jpl_PosDataVer2& helioData, double time, struct jpl_PosDataVer2& helioData2);

