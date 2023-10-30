// Copyright (C) 1994-2004 Paul Hardy
// Copyright (C) 2011-2015 Alan W. Irwin
// Copyright (C) 2013 David Howells <dhowells@redhat.com>
//
// This file is part of the timeephem software project.
//
// timeephem is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; version 2.1 of the License.
//
// timeephem is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with timeephem; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

//   Note that this software is not a product of the Jet Propulsion
//   Laboratory; it just uses and supports their ASCII and binary
//   ephemeris files.  Please don't mail JPL concerning any bugs.
//   Send bug reports or suggestions to airwin@users.sourceforge.net instead.

//! @file
//! Header information for the ephcom library.
//!

#ifndef __EPHCOM_H__
#define __EPHCOM_H__

#include <stdio.h>


#if defined ( WIN32 ) && !defined ( __CYGWIN__ )
  #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __WATCOMC__ )
    #define inline
  #elif defined ( __GNUC__ )
    #define inline    inline
  #else
    #define inline
  #endif
#endif

// Set up for dealing with function visibility issues.
#include "ephcomdll.h"

#define EPHCOM_SECONDS_PER_DAY    86400.
#define EPHCOM_J2000_EPOCH        2451545.0

// IAU definitions of time constants from
// IAU Resolution B1.9, 2000, <http://www.iau.org/static/resolutions/IAU2000 French.pdf>
// IAU Resolution B3, 2006, <http://www.iau.org/static/resolutions/IAU2006 Resol3.pdf>
// that are required for converting between TDB and TCB, between
// TT and TCG, and ultimately between TT and TDB.
// These first two IAU constants are defined in units of days.
// uncrustify cannot handle floating point constants in macros.

#define EPHCOM_IAU_T_0               2443144.5003725
#define EPHCOM_IAU_TDB_0             ( 6.55e-5 / EPHCOM_SECONDS_PER_DAY )
// These constants are dimensionless.
#define EPHCOM_IAU_LB                1.550519768e-8
#define EPHCOM_IAU_LG                6.969290134e-10

#define EPHCOM_PI                    3.14159265358979323846

#define EPHCOM_VERSION               "1.0"

// Maximum JPL ascii format version for which this software is designed.
#define EPHCOM_MAX_ASCII_VERSION     1

// Maximum JPL binary format version for which this software is designed.
#define EPHCOM_MAX_BINARY_VERSION    1

// These kind of objects are typically used for compact indices which
// are only of concern if you are using the libephcom API.  They can
// be ignored if you restrict your calls to the libephcomc API (a
// higher level API than that of libephcom).

enum ephcom_compact_object
{
    EPHCOM_COMPACT_MERCURY         = 0,
    EPHCOM_COMPACT_VENUS           = 1,
    EPHCOM_COMPACT_EMBARY          = 2,       // Earth-Moon Barycenter
    EPHCOM_COMPACT_MARS            = 3,
    EPHCOM_COMPACT_JUPITER         = 4,
    EPHCOM_COMPACT_SATURN          = 5,
    EPHCOM_COMPACT_URANUS          = 6,
    EPHCOM_COMPACT_NEPTUNE         = 7,
    EPHCOM_COMPACT_PLUTO           = 8,
    EPHCOM_COMPACT_GEOMOON         = 9,       // Geocentric Lunar coordinates
    EPHCOM_COMPACT_SUN             = 10,
    EPHCOM_COMPACT_NUTATION        = 11,
    EPHCOM_COMPACT_LIBRATION       = 12,
    EPHCOM_COMPACT_LUNAR_CORE      = 13,    // Lunar core angles
    EPHCOM_COMPACT_TEI             = 14,    // Time-ephemeris integral
    EPHCOM_COMPACT_TEV             = 15,    // Time-ephemeris vector
    EPHCOM_COMPACT_ORDERED_INDICES = 16,    // Maximum number of ordered compact body indices
    EPHCOM_COMPACT_MAXOBJECTS      = 400    // Maximum number of compact body indices including arbitrarily ordered ones such as asteroids.
};

// Offset between non-compact and compact indices for large indices.
#define EPHCOM_NONCOMPACT_OFFSET    3

//   These kind of objects are typically used for non-compact indices.
enum ephcom_object
{
    EPHCOM_MERCURY         = 0,
    EPHCOM_VENUS           = 1,
    EPHCOM_EARTH           = 2,
    EPHCOM_MARS            = 3,
    EPHCOM_JUPITER         = 4,
    EPHCOM_SATURN          = 5,
    EPHCOM_URANUS          = 6,
    EPHCOM_NEPTUNE         = 7,
    EPHCOM_PLUTO           = 8,
    EPHCOM_MOON            = 9,         // Moon relative to SSBARY
    EPHCOM_SUN             = 10,
    EPHCOM_SSBARY          = 11,        // Solar System Barycenter
    EPHCOM_EMBARY          = 12,        // Earth-Moon Barycenter
    EPHCOM_NUTATION        = 13,
    EPHCOM_LIBRATION       = 14,
    EPHCOM_GEOMOON         = 15,        // Geocentric Lunar coordinates
    EPHCOM_LUNAR_CORE      = 16,        // Lunar core angles
    EPHCOM_TEI             = 17,        // Time-ephemeris integral
    EPHCOM_TEV             = 18,        // Time-ephemeris vector
    // Maintenance: should equal EPHCOM_COMPACT_ORDERED_INDICES + EPHCOM_NONCOMPACT_OFFSET
    EPHCOM_ORDERED_INDICES = 19,        // Maximum number of ordered body indices
    // Maintenance: should equal EPHCOM_COMPACT_MAXOBJECTS + EPHCOM_NONCOMPACT_OFFSET
    EPHCOM_MAXOBJECTS      = 403        // Maximum number of body indices including arbitrarily ordered ones such as asteroids
};

// Maximum # characters to allow in ephcom ascii input line
// Maintenance: should equal 6*EPHCOM_COMPACT_MAXOBJECTS to allow
// input of ipt values for asteroids.
#define EPHCOM_MAXLINE            2400

// Maximum number of characters in header names.
#define EPHCOM_MAX_CNAM_LENGTH    6

// Maximum # characters of information in ttl[0-3]
// Maintenance: should be equal to 14*EPHCOM_MAX_CNAM_LENGTH
#define EPHCOM_MAXTTL    84

// Maximum # constant names and values in header.
#define EPHCOM_MAXCON    4000

// Maximum # of coordinates and time derivatives of those coordinates in pv.
// I doubt this will ever change for future ephemerides since we live
// in a universe with 3 spatial dimensions and one time dimension!  :-)
#define EPHCOM_MAXPVCOORD                 6
#define EPHCOM_MINJD                      -999999999.5
#define EPHCOM_MAXJD                      999999999.5

// Documented offset between the NAIF ID code for asteroids and corresponding IAU number for asteroids.
#define EPHCOM_NAIF_ID_ASTEROID_OFFSET    2000000

EPHCOMDLLIMPEXP extern const char *const ephcom_object_names[EPHCOM_ORDERED_INDICES];

//
// Names of the objects in the compact Chebyshev coefficient arrays.
//
// Index corresponds to type enum ephcom_compact_object
//

static const char *const   ephcom_compact_coeff_name[EPHCOM_COMPACT_ORDERED_INDICES] = {
    "Mercury",
    "Venus",
    "EMBary",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
    "Geocentric Moon",
    "Sun",
    "Nutation",
    "Libration",
    "Lunar Core",
    "Time-ephemeris Integral",
    "Time-ephemeris Vector"
};

EPHCOMDLLIMPEXP extern const char ephcom_month_names[12][4];

enum ephcom_calendar_type
{
    EPHCOM_JULIAN_CALENDAR    = -1,
    EPHCOM_AUTOMATIC_CALENDAR = 0,
    EPHCOM_GREGORIAN_CALENDAR = 1,
};

//! This struct holds all the information contained in a JPLEPH header.
//! When an ASCII or binary header is read, this structure is populated.
//! Fill out this structure before writing an ASCII or binary header, and
//! before performing any interpolations.
//!
struct ephcom_header
{
    //! Binary format version.
    unsigned binary_version;

    //! Data record size in 4-byte words, i.e., the actual record size in
    //! bytes is always 4 times this value.
    unsigned ksize;

    //! Number of double-precision coefficients in data blocks, i.e.,
    //! the starting and stopping Julian day numbers (preferably
    //! integral or half-integral for maximum numerical precision) +
    //! ncoeff-2 Chebyshev coefficients.  In other words, the record
    //! size in bytes is 8*ncoeff.
    unsigned ncoeff;

    //! Hold up to 14*6=EPHCOM_MAXTTL characters + "\n" + null.
    char     ttl[3][EPHCOM_MAXTTL + 2];

    //! Number of defined values for cnam and cval.
    unsigned ncon;

    //! Maximum number of characters in a given cnam string for this particular
    //! ephemeris (excluding trailing NULL).  Must be less than or equal
    //! to EPHCOM_MAX_CNAM_LENGTH
    unsigned cnam_length;

    //! Hold up to EPHCOM_MAXCON names each with up to
    //! EPHCOM_MAX_CNAM_LENGTH characters (excluding trailing NULL).
    char     cnam[EPHCOM_MAXCON][EPHCOM_MAX_CNAM_LENGTH + 1];

    //! Number of values for cval, to compare with ncon.
    unsigned nval;

    //! Constant values, corresponding to cnam names.
    double   cval[EPHCOM_MAXCON];

    //! Astronomical unit in km.
    double   au;

    //! Earth-Moon mass ratio.
    double   emrat;

    //! Speed of light, km/sec.
    double   clight;

    //! Ephemeris number.
    unsigned numde;

    //! Lunar ephemeris number (can be same # as numde).
    unsigned numle;

    //! Time ephemeris number (can be zero, but this may be set to the
    //! equivalent denum used to generate the time ephemeris data if
    //! those data are included in the ephemeris).
    unsigned numte;

    //! Start Julian day number, stop Julian day number, and step size (in days).
    double   ss[3];

    //! Maximum number of possible compact indices for this particular ephemeris.
    //! Note 13 <= npt <= EPHCOM_COMPACT_MAXOBJECT where the lower limit
    //! is required to be compatible with traditional JPL ephemerides.
    unsigned npt;

    //! Index pointers into Chebyshev coefficients.  If index is the
    //! compact index referring to a particular ephemeris quantity, then ipt[index][0],
    //! ipt[index][1], and ipt[index][2] refer to the first index
    //! within the Chebyshev coefficients, the number of Chebyshev
    //! coefficients, and the number of epoch sub-intervals within the
    //! epoch interval for that particular ephemeris quantity.
    //! ipt is only defined if its first index is less than npt.
    unsigned ipt[EPHCOM_COMPACT_MAXOBJECTS][3];

    //! Number of coordinates for each ephemeris quantity defined by a
    //! compact index.  These values are 1 for TEI, 2 for nutation,
    //! and 3 for all other kinds of ephemeris quantities.
    //! ncoords is only defined if its index is less than npt.
    unsigned ncoords[EPHCOM_COMPACT_MAXOBJECTS];

    //! Calculated maximum number of Chebyshev coefficients for all
    //! kinds of ephemeris data for this particular ephemeris.
    unsigned maxcheby;

    //! When the Chebyshev coefficients corresponding to the time ephemeris
    //! integral are read (either in ascii or binary form) the following
    //! linear transformation is applied when tei_transform is true.
    //! tei = tei_offset + tei_scale*tei_input.
    unsigned tei_transform;
    double   tei_offset;
    double   tei_scale;

    //! When the Chebyshev coefficients corresponding to the time ephemeris
    //! vector are read (either in ascii or binary form) the following
    //! scaling transformation is applied when tev_transform is true.
    //! tev = tev_scale*tev_input.
    unsigned tev_transform;
    double   tev_scale;

    //! Derived quantities associated with "asteroids" where that term means any
    //! solar system object not handled by the "ordered" indices, i.e.,
    //! any solar system object other than the major planets, Pluto, the Sun
    //! and the Moon.

    //! Number of asteroids.  This should be set equal to the maximum of
    //! 0 and npt - EPHCOM_COMPACT_ORDERED_INDICES.
    unsigned nasteroids;

    //! Index within _ordered_ constant names that corresponds to the
    //! first "MA####" header constant name encountered.
    unsigned asteroid_start;

    //! This array keeps track of the naif identification code of the asteroid indices
    //! which those indices are defined by either (1)
    //! the ephcom compact index of the object minus EPHCOM_COMPACT_ORDERED_INDICES
    //! or (2)
    //! the ephcom index of the object minus EPHCOM_ORDERED_INDICES.
    //!
    //! This array of NAIF ID values is calculated to correspond to
    //! the mass constants in the header named "MA####", where the
    //! number "####" corresponds to the IAU designation of the
    //! asteroid, and 2 000 000 + "####" corresponds to the NAIF ID of
    //! the object.  This array is sorted by NAIF ID.

    unsigned asteroid_id[EPHCOM_MAXOBJECTS - EPHCOM_ORDERED_INDICES];

    //! Number of defined indices in asteroid_id
    unsigned nasteroid_id;

    //! Compact index of "center" index of asteroid ephemerides.  A
    //! negative value means the asteroid motions are referred to SSB
    //! which would correspond to NAIF ID code of 0.  Experience so
    //! far with the only publicly distributed asteroid ephemeris I
    //! have access to (i.e., the JPL asteroid ephemeris for de430) is
    //! the asteroid motions are referred to the Sun which corresponds to NAIF ID of 10
    //! and also compact index of 10.
    int asteroid_center;
};
typedef struct ephcom_header ephcom_header_t;

// struct needed for sorting of names and values by name.
struct ephcom_sort
{
    char   name[EPHCOM_MAX_CNAM_LENGTH + 1];
    double value;
};
typedef struct ephcom_sort ephcom_sort_t;

enum ephcom_distance_units
{
    EPHCOM_UNITS_AU = 0,
    EPHCOM_UNITS_KM = 1,
};

enum ephcom_time_units
{
    EPHCOM_UNITS_DAYS    = 0,
    EPHCOM_UNITS_SECONDS = 1,
};

enum ephcom_coordinates
{
    EPHCOM_ADJUST_FOR_CENTRE  = 0,
    EPHCOM_BARYCENTRIC_COORDS = 1,
};

// Specification of endianness which must agree with numbers in
// ../test_host_endianness.c (used by the CMake-based build system to
// set the EPHCOM_HOST_ENDIANNESS macro for file.c).
enum ephcom_endianness
{
    EPHCOM_LITTLE_ENDIAN = 1,
    EPHCOM_BIG_ENDIAN    = 2,
    EPHCOM_DEC_ENDIAN    = 3,
};

//! This struct contains information that controls the determination
//! of results consisting of interpolated positions and velocities of
//! planets, Sun, and Moon; interpolated nutation and libration angles
//! and their time derivatives; and possibly other kinds of
//! user-defined data and time derivatives.  When the desired results
//! are determined, they are also stored in this struct.
//!
//! To populate this structure, have the ephemeris file open, and set:
//!
//! km, seconds, bary, et2[0], and et2[1], and list in the ephcom_coords struct whose
//! pointer is an argument to ephcom_get_coords().
//!
//! Then call ephcom_get_coords() to determine the coordinates
//! specified by list, see the ephcomc library routine
//! ephcom_interpolate_list(_deprecate) for an example.
//!
//! The ephcom indices values in list below can be one of the following values.
//! <table border>
//!   <tr> <td><b>Numerical Index</b></td> <td><b>Symbolic Index</b></td><td><b>Identification</b></td> </tr>
//!   <tr> <td>0</td> <td>EPHCOM_MERCURY</td> <td>Mercury</td> </tr>
//!   <tr> <td>1</td> <td>EPHCOM_VENUS</td> <td>Venus</td> </tr>
//!   <tr> <td>2</td> <td>EPHCOM_EARTH</td> <td>Earth</td> </tr>
//!   <tr> <td>3</td> <td>EPHCOM_MARS</td> <td>Mars</td> </tr>
//!   <tr> <td>4</td> <td>EPHCOM_JUPITER</td> <td>Jupiter</td> </tr>
//!   <tr> <td>5</td> <td>EPHCOM_SATURN</td> <td>Saturn</td> </tr>
//!   <tr> <td>6</td> <td>EPHCOM_URANUS</td> <td>Uranus</td> </tr>
//!   <tr> <td>7</td> <td>EPHCOM_NEPTUNE</td> <td>Neptune</td> </tr>
//!   <tr> <td>8</td> <td>EPHCOM_PLUTO</td> <td>Pluto</td> </tr>
//!   <tr> <td>9</td> <td>EPHCOM_MOON</td> <td>Moon</td> </tr>
//!   <tr> <td>10</td> <td>EPHCOM_SUN</td> <td>Sun</td> </tr>
//!   <tr> <td>11</td> <td>EPHCOM_SSBARY</td> <td>Solar System Barycenter</td> </tr>
//!   <tr> <td>12</td> <td>EPHCOM_EMBARY</td> <td>Earth-Moon Barycenter</td> </tr>
//!   <tr> <td>13</td> <td>EPHCOM_NUTATION</td> <td>Nutation Angles</td> </tr>
//!   <tr> <td>14</td> <td>EPHCOM_LIBRATION</td> <td>Libration Angles</td> </tr>
//!   <tr> <td>15</td> <td>EPHCOM_GEOMOON</td> <td>Moon (Geocentric)</td> </tr>
//!   <tr> <td>16</td> <td>EPHCOM_LUNAR_CORE</td> <td>Lunar Core Angles</td> </tr>
//!   <tr> <td>17</td> <td>EPHCOM_TEI</td> <td>Time-ephemeris Integral</td> </tr>
//!   <tr> <td>18</td> <td>EPHCOM_TEV</td> <td>Time-ephemeris Vector</td> </tr>
//!   <tr> <td>19+</td> <td>EPHCOM_ORDERED_INDICES+</td> <td>"Asteroid" indices</td> </tr>
//!   <tr> <td>2000000+</td> <td>EPHCOM_NAIF_ID_ASTEROID_OFFSET+</td> <td>NAIF ID's to be transformed to "Asteroid" indices</td> </tr>
//! </table>
//!
struct ephcom_coords
{
    //! 1 = positions in km; 0 = positions in AU.
    enum ephcom_distance_units km;
    //! 1 = timescale is seconds; 0 = timescale is days.
    enum ephcom_time_units     seconds;
    //! 1 = Barycentric coordinates; 0 = adjust for center.
    enum ephcom_coordinates    bary;
    //! object to use as center (instead of SSBARY).
    unsigned center;
    //! Julian Day of interpolation. For best precision, et2[0] should be an exact integral or exact half-integral number of days while et2[1] should be a correction to et2[0] between 0. and 1.
    double   et2[2];
    //! Maximum number of possible non-compact indices for this particular ephemeris
    //! (= header.npt + EPHCOM_NONCOMPACT_OFFSET).
    unsigned maxobjects;

    //! list of indices (interpreted as in above table) to search for in ephemeris.
    unsigned list[EPHCOM_MAXOBJECTS];
    //! number of indices in list.
    unsigned nlist;
    //! x, y, z Position & Velocity.  The first index of pv is interpreted as in the above table.
    double   pv[EPHCOM_MAXOBJECTS][EPHCOM_MAXPVCOORD];
};
typedef struct ephcom_coords ephcom_coords_t;

//
// Ephemeris dataset file context
//
struct ephcom_context
{
    // Copy of the filename for error printing purposes
    char          *filename;

    // Length of filename string
    size_t        len_filename;

    // FILE pointer for the file.
    FILE          *file;

    // Useful data describing the contents
    ephcom_header header;

    // Endianness of binary ephemeris file
    // that is being processed.
    enum ephcom_endianness binary_file_endianness;

    // Number of last block read (if not 0).
    unsigned last_blocknum;

    // Pointer to an array that will contain the block of Chebyshev
    // coefficient data that has been read that has a Julian date range
    // that contains the specified Julian date in coords where the last
    // interpolation occurred.
    //
    double *datablock;

    // Pointers to arrays and stuff that will by used by ephcom_cheby().
    // One array is for position components and the other for velocity
    // components.
    //
    double *pc;
    double *vc;
    double twox;
    double lastx;
};
typedef struct ephcom_context ephcom_context_t;

// The following 0-suffix routines should only be used by utilities
// that use the low-level libephcom API.  User-level code should
// be using the high-level libephcomc API instead which includes
// wrappers for these without the 0 suffix on the name.

// Open an ephcom file for reading.  This will read the header and constants
// which can then be accessed from the ephcom struct.  A NULL return indicates
// an error occurred.
EPHCOMDLLIMPEXP ephcom_context *
ephcom_open0( const char *filename );

// Close an ephcom file.
EPHCOMDLLIMPEXP int
ephcom_close0( ephcom_context *ctx );

// Useful macro for Julian date comparison.
// JPL ephemeris Julian dates (so far) range from 6.e5 (3000 BC) to 2.8e6 (3000 AD).  So
// double precision (64-bit floating point) should be able to represent a date to within
// ~ 1.e-15*2.8e6 = 2.8e-9 days (or 0.24 ms).  Make criterion ~100 times that value since
// all JPL ephemerides up to now have dates which are offset by 0.5 from an integer so
// that a criterion of 0.25 days would probably work.

#define JULIAN_DATE_CRITERION    ( 3.e-7 )

// Need this logic to compare dates rather than exact equalities because
// (for DE422 at least) rounding errors have crept into the dates.
#define IF_SAME_DATE( a, b )    ( fabs( a - b ) <= JULIAN_DATE_CRITERION )

// Here is where the public API for libephcom is described.  This API
// has been changed considerably in recent releases, but that should
// not matter to users because libephcom is a low-level library that
// user routines should not be linking directly to.  Instead, user
// routines should be linking to the ephcomc library which has a more
// stable (and higher-level) API.

// Read a JPL Ephemeris ASCII header from the FILE pointed to by the
// infp argument and store and return values in the ephcom_header
// struct pointed to by the header argument.  Write any errors to
// stderr.
EPHCOMDLLIMPEXP int
ephcom_readascii_header( const char *filename, FILE *infp, ephcom_header *header );

// Read a foreign (i.e., written in an ascii format that
// ephcom_readascii_header would not understand) ephemeris ASCII
// header from the FILEs pointed to by the infp_header and infps array
// argument and store and return values in the ephcom_header struct
// pointed to by the header argument.  Write any errors to stderr.
EPHCOMDLLIMPEXP int
foreign_ephcom_readascii_header( const char *filename_header, const char * const * filenames, FILE *infp_header, FILE **infps, unsigned ascii_kind, ephcom_header *header );

// Read a block of data coefficients from a JPL ASCII ephemeris file.
// Returns number of coefficients read, 0 at EOF.
EPHCOMDLLIMPEXP int
ephcom_readascii_block( const char *filename, FILE *infp,
                        ephcom_header *header, double *datablock );

// Read a block of data coefficients from a collection of foreign
// (i.e., written in an ascii format that ephcom_readascii_block
// would not understand) ASCII ephemeris files.  Returns number of
// coefficients read, 0 at EOF.
EPHCOMDLLIMPEXP int
foreign_ephcom_readascii_block( const char * const * filenames, FILE **infps,
                                unsigned ascii_kind, ephcom_header *header, double *datablock );

// Calculate additional header information consistent with npt, ipt[<index>][1], and ipt[<index>][2].
EPHCOMDLLIMPEXP int
ephcom_ipt_fixup( const char *filename, unsigned if_check, ephcom_header *header );

// Read a JPL Ephemeris header in binary format.  Store values in
// an ephcom_header struct.
EPHCOMDLLIMPEXP int
ephcom_readbinary_header( ephcom_context *ctx );

// Read a JPL Ephemeris data block in binary format.
EPHCOMDLLIMPEXP int
ephcom_readbinary_block( ephcom_context *ctx,
                         unsigned blocknum   // Data block number, starting with 0
                         );

// Write header information in ASCII format.
EPHCOMDLLIMPEXP int
ephcom_writeascii_header( const char *filename, FILE *outfp, ephcom_header *header );

// Write coefficient block information in ASCII format.
EPHCOMDLLIMPEXP int
ephcom_writeascii_block( const char *filename, FILE *outfp, ephcom_header *header,
                         unsigned blocknum, double *datablock );

// Write a JPL Ephemeris header in binary format.
EPHCOMDLLIMPEXP int
ephcom_writebinary_header( const char *filename, FILE *outfp, ephcom_header *header );

// Write a block of data coefficients in JPL binary file format.
EPHCOMDLLIMPEXP int
ephcom_writebinary_block( const char *filename, FILE *outfp, ephcom_header *header,
                          unsigned blocknum, double *datablock );

// ephcom_parse_block() - Parse a binary block of data.  Warning: verbose!
//                        Writes parsed output to file pointer outfp.
EPHCOMDLLIMPEXP int
ephcom_parse_block( const char *filename, FILE *outfp,
                    ephcom_header *header, double *datablock );

// Read in a double precision value from the given file with bytes swapped
// if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
// the bytes will get swapped, on Motorola or SPARC they won't.
EPHCOMDLLIMPEXP int
ephcom_indouble( ephcom_context *ctx, double *_value );

// Read in an unsigned (4--byte) value to the given file with bytes swapped
// if necessary to match host order (Little- or DEC- Endian).  On Intel 80x86
// the bytes will get swapped, on Motorola or SPARC they won't.
EPHCOMDLLIMPEXP int
ephcom_inunsigned( ephcom_context *ctx, unsigned *_value );

// ephcom_get_coords() - Interpolate positions and velocities at given time.
EPHCOMDLLIMPEXP int
ephcom_get_coords( ephcom_context *ctx, ephcom_coords *coords );

// ephcom_cal2jd0() - convert calendar date and time to JD.
EPHCOMDLLIMPEXP void
ephcom_cal2jd0( double *tjd2, const int *idate, enum ephcom_calendar_type calendar_type );

// ephcom_jd2cal0() - convert JD to calendar date and time.
EPHCOMDLLIMPEXP void
ephcom_jd2cal0( int *idate, const double *tjd2, enum ephcom_calendar_type calendar_type );

// Sort the cnam and cvals arrays of an existing header into ascending order in
// cnam.

EPHCOMDLLIMPEXP int
ephcom_sort_header_constants( ephcom_header *header );

// ephcom_update_asteroids - update asteroid-related quantities in the header to be
// consistent with the remaining quantities in the header.
EPHCOMDLLIMPEXP int
ephcom_update_asteroids( ephcom_header * header );

// Search an ordered table for a key.
EPHCOMDLLIMPEXP int
ephcom_search( const void *key, const void *base, int n, size_t size, int *low, unsigned ( *ge )( const void *keyval, const void *datum ) );

#endif  // __EPHCOM_H__
