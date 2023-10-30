/**********************************************************************
 
 LReader.h
 March 2022 Thomas Morris

Reads a configuration file using the Lua language.
Example:

	-- configuration file for program `MarsOpp'
	-- define ephemeris time range + body
	planet = Mars
	start = {d=5, m=3, y=2022}
	end =   {d=4, m=3, y=2032}
**********************************************************************/


struct LConfig
{
	bool snapshotMode; // calculate ephemeris for a single point in time
	int startDate[3];
	int endDate[3];
	int Tsnapshot[3];
	double TFsnapshot;
	int ephcomID;
};

struct lua_State;

class LReader
{
public:
	void open(char *configFile);
	void read(LConfig& config);
	void close();
private:
	void setPlanetsTable();
	int getField (const char *key);
	bool getFieldDMY (const char *datekey, int *dmy);
	lua_State *L;
};

