#include <string>
#include <list>

struct MarsData
{
	double   JDay;
	int      dmy[3];
	double   timeofday;
	double   Delta;  // distance to Earth (AU)
	double   radius; // distance to Sun (AU)
	double   RA;     // right ascension
	double   Dec;    // declination
	double   diameter; // angular size (arcsec)
	double   separation; // separation (arcmin) for e.g. conjunctions 
	double   omega;  // longitude of central meridian (CM) (degrees)
	double   L_S;    // solar longitude (degrees)
	double   D_E;    // planetocentric declination of Earth
	double   D_S;    // planetocentric declination of Sun

	std::string dateToString() const;
};

class MarsDataLogger
{
public:
	void add(MarsData& item);
	void writeToCSVFile(const char *filename);

private:
	std::list<MarsData> data;
};

