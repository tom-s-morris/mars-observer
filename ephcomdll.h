#ifndef __EPHCOM_DLL_H
#define __EPHCOM_DLL_H

#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define EPHCOMDLLEXPORT    __declspec( dllexport )
      #define EPHCOMDLLIMPORT    __declspec( dllimport )
    #else
      #define EPHCOMDLLEXPORT
      #define EPHCOMDLLIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define EPHCOMDLLEXPORT    __declspec( dllexport )
    #define EPHCOMDLLIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// EPHCOMDLLEXPORT to be visible.
    #define EPHCOMDLLEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define EPHCOMDLLIMPORT
  #endif
#endif

// For an unknown compiler or static build we clear the macros
#ifndef EPHCOMDLLEXPORT
  #define EPHCOMDLLEXPORT
  #define EPHCOMDLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case), if
// cmake didn't set the corresponding macro xxxx_EXPORTS when the
// corresponding library is built (DLLIMPEXP is set to DLLEXPORT
// then)
#if defined ( gnulliver_EXPORTS )
  #define GNULLIVERDLLIMPEXP    EPHCOMDLLEXPORT
  #define GNULLIVERDLLIMPEXP_DATA( type )    EPHCOMDLLEXPORT type
#else
  #define GNULLIVERDLLIMPEXP    EPHCOMDLLIMPORT
  #define GNULLIVERDLLIMPEXP_DATA( type )    EPHCOMDLLIMPORT type
#endif

#if defined ( ephcom_EXPORTS )
  #define EPHCOMDLLIMPEXP    EPHCOMDLLEXPORT
  #define EPHCOMDLLIMPEXP_DATA( type )    EPHCOMDLLEXPORT type
#else
  #define EPHCOMDLLIMPEXP    EPHCOMDLLIMPORT
  #define EPHCOMDLLIMPEXP_DATA( type )    EPHCOMDLLIMPORT type
#endif

#if defined ( ephcom_spice_EXPORTS )
  #define EPHCOM_SPICEDLLIMPEXP    EPHCOMDLLEXPORT
  #define EPHCOM_SPICEDLLIMPEXP_DATA( type )    EPHCOMDLLEXPORT type
#else
  #define EPHCOM_SPICEDLLIMPEXP    EPHCOMDLLIMPORT
  #define EPHCOM_SPICEDLLIMPEXP_DATA( type )    EPHCOMDLLIMPORT type
#endif

#if defined ( ephcomc_EXPORTS )
  #define EPHCOMCDLLIMPEXP    EPHCOMDLLEXPORT
  #define EPHCOMCDLLIMPEXP_DATA( type )    EPHCOMDLLEXPORT type
#else
  #define EPHCOMCDLLIMPEXP    EPHCOMDLLIMPORT
  #define EPHCOMCDLLIMPEXP_DATA( type )    EPHCOMDLLIMPORT type
#endif

#endif // __EPHCOM_DLL_H
