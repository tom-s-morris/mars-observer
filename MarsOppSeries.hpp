/**********************************************************
 * Method from Meeus, astronomical algorithms, ch.36
 *********************************************************/
class MarsOppSeries
{
private:
	const double millionKmInAU = 149.597870700; /* 1AU in million km */

	double Mars_C0[11] = {
		 -0.3088,
		-17.6965,
		+18.3131,
		 -0.2162,
		 -4.5028,
		 +0.8797,
		 +0.7666,
		 -0.3636,
		 +0.0402,
		 +0.0737,
		 -0.0980};

	double Mars_C1[11] = {
		0,
		+0.0363,
		+0.0467,
		-0.0198,
		-0.0019,
		+0.0058,
		-0.0050,
		-0.0001,
		+0.0032,
		-0.0008,
		-0.0011};

	double Mars_C2[11] = {
		+2e-5,
		+5e-5,
		-6e-5,
		-1e-5,
		+7e-5,
		-2e-5,
		-3e-5,
		+2e-5,
		0,
		0,
		0};

int find_k(double Year)
{
	// next 'k' to give opposition in the current year
	return floor ((365.2425*Year + 1721060 - A) / B) + 1;
}

protected:
	ephcom_object targetID;
	int numTerms;
	double A;
	double B;
	double M0;
	double M1;

	double *C0;
	double *C1;
	double *C2;

virtual double series_correction(double M, double t)
{
	double JDE1 = 0;

	/* convert M to radians */
	M *= M_PI/180.0;

	//cout << "series_correction(): numTerms = " << numTerms << ", C0[0] = " << C0[0] << endl;

	/* const terms */
	for(int a=1; a<numTerms; a+=2)
		JDE1 += C0[a] * sin((a+1)/2 * M);
	for(int a=0; a<numTerms; a+=2)
		JDE1 += C0[a] * cos(a/2 * M);

	/* linear terms */
	for(int a=1; a<numTerms; a+=2)
		JDE1 += C1[a] * sin((a+1)/2 * M) * t;
	for(int a=0; a<numTerms; a+=2)
		JDE1 += C1[a] * cos(a/2 * M) * t;

	if (C2)
	{
	/* quadratic terms */
	for(int a=1; a<numTerms; a+=2)
		JDE1 += C2[a] * sin((a+1)/2 * M) * t*t;
	for(int a=0; a<numTerms; a+=2)
		JDE1 += C2[a] * cos(a/2 * M) * t*t;
	}
	return JDE1;
}

virtual void set_coeffs()
{
	targetID = EPHCOM_MARS;
	C0 = Mars_C0;
	C1 = Mars_C1;
	C2 = Mars_C2;
	numTerms = 11;
	A = 2452097.382;
	B = 779.936104; // Mars's synodic period (d)
	M0 = 181.9573;  // mean anomaly of Earth at opposition
	M1 = 48.705244;
}

public:
	double JDE_mean, JDE_true;
	jpl_PosData observables;

MarsOppSeries()
{
	set_coeffs();
}

double compute_next_opposition(double time)
{
	int k = find_k(time);
	cout << "Year = " << time << ", k = " << k << endl;
	return compute_next_opposition(k);
}

double compute_next_opposition(int k)
{
	double JDE = A + k*B; /* time of mean opposition */
	double M = M0 + k*M1; /* mean anomaly of the Earth */

	JDE_mean = JDE;

	/* Series for correction to mean JDE (Table 36.B) */
	double t = ((JDE - 2451545.0) / 36525); /* Julian centuries before/after 2000 */
	JDE += series_correction(M, t);

	JDE_true = JDE;
#if defined(JPLEPH_VER1)
	// Legacy EPHCOM v1.0 library and JPL DE405
	jpleph_compute_de405(target, JDE, &observables);
#else
	// Create an EPHCOM3/4 context
	LibEphcom4 Ephcom4;
	Ephcom4.initContext(targetID);
	Ephcom4.setMeanEquinox();
	Ephcom4.computeDE(static_cast<int> (targetID), JDE, &observables);
	Ephcom4.deinitContext();
#endif
	return JDE;
}

void find_nearest_opposition(double timeOpp)
{
	double JD2000 = 2451545.0;
	double year = 2000 + floor((timeOpp - JD2000) / 365.25);
	compute_next_opposition(year);
	pretty_print();

	cout << std::setprecision (12) <<
		"Meeus: " << JDE_true << ", " <<
		"JPL = " << timeOpp << ", " <<
		"delta = " << timeOpp - JDE_true << endl;
}

void pretty_print()
{
	cout << "+------------------------------+" << endl;
	cout << std::setprecision (12) << "| Mean JDE = " << JDE_mean << "     |" << endl;
	cout << std::setprecision (12) << "| True JDE = " << JDE_true << "     |" << endl;
	cout << std::setprecision (12) << "| Delta    = " << observables.dist_Geo << " AU |" << endl;
	cout << std::setprecision (12) << "|          = " << 
		observables.dist_Geo * millionKmInAU << " mkm |" << endl;
	cout << std::setprecision (6)  << "| Diameter = " <<
		9.36/observables.dist_Geo << " \"        |" << endl;
	cout << "+------------------------------+" << endl << endl;
}
};


class JupiterOppSeries : public MarsOppSeries
{
private:
	/* Identical to Mars class, with coeffs for Jupiter. */
	double Jupiter_C0[7] = {
		-0.1029,
		-1.9658,
		+6.1537,
		-0.2081,
		-0.1116,
		+0.0074,
		-0.0097};

	double Jupiter_C1[7] = {
		0,
		-0.0056,
		+0.0210,
		-0.0013,
		-0.0010,
		+0.0001,
		-0.0001};

	double Jupiter_C2[7] = {
		-9e-5,
		+7e-5,
		-6e-5,
		0,
		0,
		0,
		0};

	// auxiliary angle, see Meeus p.251
	double Jupiter_C0_a[2] = {
		0,
		+0.3642};

	double Jupiter_C1_a[2] = {
		+0.0144,
		-0.0019};

	double Jupiter_C2_a[2] = {
		-8e-5,
		-2.9e-4};

virtual void set_coeffs()
{
	targetID = EPHCOM_JUPITER;
	C0 = Jupiter_C0;
	C1 = Jupiter_C1;
	C2 = Jupiter_C2;
	numTerms = 7; // 7 terms in Fourier series expansion in Earth's anomaly (M)
	A = 2451870.628;
	B = 398.884046; // Jupiter's synodic period (d)
	M0 = 318.4681;  // mean anomaly of Earth at opposition
	M1 = 33.140229; // advance of mean anomaly between successive oppositions
}

virtual double series_correction(double M, double t)
{
	double JDE1 = MarsOppSeries::series_correction(M, t);

	/* add correction due to auxiliary angle */
	double a = 82.74 + 40.76 * t;
	a *= M_PI/180.0;

	/* const terms */
	JDE1 += Jupiter_C0_a[0] * sin(a) +
		Jupiter_C0_a[1] * cos(a);

	/* linear terms */
	JDE1 += Jupiter_C1_a[0] * sin(a) * t +
		Jupiter_C1_a[1] * cos(a) * t;

	/* quadratic terms */
	JDE1 += Jupiter_C2_a[0] * sin(a) * t*t +
		Jupiter_C2_a[1] * cos(a) * t*t;
	return JDE1;
}

public:
JupiterOppSeries()
{
	set_coeffs();
}


};
/*******************************************************************************/

