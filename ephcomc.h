// Copyright (C) 2011-2018 Alan W. Irwin
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

// ephcomc.h - header information for the ephcomc
// C library for ephcom which contains high-level API that is suitable for
// the ephcom_testeph and ephcom_vtransit applications as well as
// for creating the ephcom language bindings.

#ifndef __EPHCOMC_H__
#define __EPHCOMC_H__

#include "ephcomdll.h"
#include "ephcom.h"

// Bit fields in the control variable:

#define EPHCOM_CONTROL_KM         1
#define EPHCOM_CONTROL_SECONDS    2

// C Utility routines that help to interface ephcom language bindings
// but which are only meant to be available at the C level.

EPHCOMCDLLIMPEXP int
ephcom_Alloc2dChar( char ***fp, unsigned nx, unsigned ny );

EPHCOMCDLLIMPEXP int
ephcom_Free2dChar( char **f, unsigned nx );

EPHCOMCDLLIMPEXP int
ephcom_Alloc2dDouble( double ***fp, unsigned nx, unsigned ny );

EPHCOMCDLLIMPEXP int
ephcom_Free2dDouble( double **f, unsigned nx );

// Actual C API to be interfaced to all ephcom language bindings.

// The next few routines are just wrappers (named without a 0 suffix)
// for 0-suffixed routines that are declared in ephcom.h for
// the libephcom low-level API.  These wrappers make these low-level
// API routines also available to the high-level API for libephcomc which
// allows high-level API users and bindings code
// to link only to the libephcomc library.

//
// Open an ephcom file for reading.  This will read the header and constants
// which can then be accessed from the ephcom struct.  A NULL return indicates
// an error occurred.
//
EPHCOMCDLLIMPEXP
ephcom_context *ephcom_open( const char *filename );

// Close an ephcom file.
EPHCOMCDLLIMPEXP int
ephcom_close( ephcom_context *ctx );

// ephcom_cal2jd() - convert calendar date and time to JD.
EPHCOMCDLLIMPEXP void
ephcom_cal2jd( double *tjd2, const int *idate, enum ephcom_calendar_type calendar_type );

// ephcom_jd2cal() - convert JD to calendar date and time.
EPHCOMCDLLIMPEXP void
ephcom_jd2cal( int *idate, const double *tjd2, enum ephcom_calendar_type calendar_type );

// End of the wrappers for 0-suffixed routines declared in ephcom.h

// Recommended, thread-safe way to handle the context-sensitive part of our API.
// These are also bound by our Python and Fortran bindings.

EPHCOMCDLLIMPEXP int
ephcom_read_constants( ephcom_context *ctx, char **cnames, unsigned nvalues, unsigned cname_length, double *values, unsigned *numde, double *au, double *emrat, double *sss );

EPHCOMCDLLIMPEXP int
ephcom_interpolate_relative( ephcom_context *ctx, const double *et2, unsigned control, int ntarget, int ncenter, double *pv );

EPHCOMCDLLIMPEXP int
ephcom_interpolate_list( ephcom_context *ctx, const double *et2, unsigned control, unsigned nlist, const unsigned *list, double **pv );

// Deprecated, thread-unsafe way to handle the context-sensitive part of our API.  We will
// probably keep these indefinitely because all but ephcom_interpolate_list_deprecated
// are required for the
// deprecated f77 part of our Fortran binding, and they could also prove useful for
// bindings for other languages which might not have a way to pass C pointers such
// as an ephcom_context * as language arguments.

EPHCOMCDLLIMPEXP int
ephcom_read_constants_deprecated( const char *filename, char **cnames, unsigned nvalues, unsigned cname_length, double *values, unsigned *numde, double *au, double *emrat, double *sss );

EPHCOMCDLLIMPEXP int
ephcom_interpolate_relative_deprecated( const char *filename, const double *et2, unsigned control, int ntarget, int ncenter, double *pv );

EPHCOMCDLLIMPEXP int
ephcom_interpolate_list_deprecated( const char *filename, const double *et2, unsigned control, unsigned nlist, const unsigned *list, double **pv );

#endif  // __EPHCOMC_H__
