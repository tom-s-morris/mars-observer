#include "MarsDataLogger.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

void MarsDataLogger::add(MarsData& item)
{
	data.push_back(item);
}

void MarsDataLogger::writeToCSVFile(const char *filename)
{
	std::ofstream out_file;
	out_file.open(filename);

	for(std::list<MarsData>::const_iterator it = data.begin();
		it != data.end();
		++it)
	{
		std::ostringstream oss;
		oss << std::fixed << std::setprecision (5) << it->JDay;
		oss << ",";
		oss << it->dateToString();
		oss << ",";
		oss << std::fixed << std::setprecision (6) << it->Delta;
		oss << ",";
		oss << std::fixed << std::setprecision (6) << it->radius;
		oss << ",";
		oss << std::fixed << std::setprecision (4) << it->diameter;
		oss << ",";
		//oss << std::fixed << std::setprecision (3) << it->omega;
		//oss << ",";
		oss << std::fixed << std::setprecision (2) << it->L_S;
		oss << ",";
		oss << std::fixed << std::setprecision (2) << it->D_S;
		out_file << oss.str() << std::endl;
	}
	out_file.close();
}

std::string MarsData::dateToString() const
{
	char bufDateTime[64];
	sprintf(bufDateTime, "%02d/%02d/%04d %02d:%02d",
			dmy[0], dmy[1], dmy[2],
			(int) floor(timeofday*24),
			(int) round(60*(timeofday*24 - floor(timeofday*24))));
	return std::string(bufDateTime);
}

