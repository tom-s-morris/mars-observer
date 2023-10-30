/**********************************************************************
 
 LReader.cpp
 March 2022 Thomas Morris

Reads a configuration file using the Lua language.
Example:

	-- configuration file for program `MarsOpp'
	-- define ephemeris time range + body
	planet = Mars
	Tstart = {d=5, m=3, y=2022}
	Tend =   {d=4, m=3, y=2032}
**********************************************************************/

#include <iostream>
#include <string>
#include <cstdio>
#include <iomanip>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}
#include "LReader.h"

extern "C" {
#include "ephcomc.h" /* for planets listed in ephemeris */
}

struct PlanetList
{
	const char *name;
	ephcom_object ephcomID;
};

const PlanetList planets[] =
{
	{"Sun",      EPHCOM_SUN},
	{"Mercury",  EPHCOM_MERCURY},
	{"Venus",    EPHCOM_VENUS},
	{"Mars",     EPHCOM_MARS},
	{"Jupiter",  EPHCOM_JUPITER},
	{"Saturn",   EPHCOM_SATURN},
	{"Uranus",   EPHCOM_URANUS},
	{"Neptune",  EPHCOM_NEPTUNE},
	{"Pluto",    EPHCOM_PLUTO},
	{NULL,       EPHCOM_MAXOBJECTS}
};

void LReader::open(char *configFile)
{
	L = luaL_newstate();
	luaopen_base(L);
	luaopen_io(L);
	luaopen_string(L);
	luaopen_math(L);

	/* Set planets as Lua global variables before running file */
	setPlanetsTable();

	if (luaL_loadfile(L, configFile) || lua_pcall(L, 0, 0, 0))
	{
		std::cerr << "Cannot run configuration file: "
		          << lua_tostring(L, -1) << std::endl;
	}
}

void LReader::read(LConfig& config)
{
	config.snapshotMode = false;

	if ( !getFieldDMY("Tstart", &config.startDate[0]) ||
		 !getFieldDMY("Tend",   &config.endDate[0]) )
	{
		/* fallback to snapshot mode */
		getFieldDMY("Tsnapshot", &config.Tsnapshot[0]);
		config.TFsnapshot = 0.0; // fractional part of Julian Day
		config.snapshotMode = true;
	}

	/* read planet/body */
	lua_getglobal(L, "planet");
	if (!lua_isinteger(L, -1) ||
		(lua_tointeger(L, -1) > EPHCOM_PLUTO &&
		 lua_tointeger(L, -1) != EPHCOM_SUN))
	{
		std::cerr << "'planet' is not a valid object (expecting index from "
		          << EPHCOM_MERCURY << " to "
		          << EPHCOM_PLUTO << ", or "
		          << EPHCOM_SUN << " for the Sun)." << std::endl;
		return;
	}
	config.ephcomID = lua_tointeger(L, -1);
}

void LReader::close()
{
	lua_close(L);
}

bool LReader::getFieldDMY (const char *datekey, int *dmy)
{
	lua_getglobal(L, datekey);
	if (!lua_istable(L, -1))
	{
		std::cerr << "'" << datekey << "' is not a valid Date table" << std::endl;
		return false;
	}

	dmy[0] = getField("d");
	dmy[1] = getField("m");
	dmy[2] = getField("y");
	return true;
}

/* Reads a table element from the stack */
int LReader::getField (const char *key)
{
	int result;
	/* Assume that table is on the stack top */
	lua_pushstring(L, key);
	lua_gettable(L, -2);
	if (!lua_isnumber(L, -1))
	{
		std::cerr << "Invalid component in Tstart/Tend date" << std::endl;
	}
	result = (int)lua_tonumber(L, -1);
	lua_pop(L, 1);  /* remove number */
	return result;
}

void LReader::setPlanetsTable()
{
	int i = 0;
	while (planets[i].name != NULL)
	{
		/* push ID to top of stack and assign to global variable */
		lua_pushinteger(L, planets[i].ephcomID);
		lua_setglobal(L, planets[i++].name);
	}
}

