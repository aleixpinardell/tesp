/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TESP_INPUT_OUTPUT
#define TESP_INPUT_OUTPUT

#include <string>

namespace tesp
{

//! Get filename for path.
static inline std::string getFilenameForPath( const std::string path, const double includeExtension = true )
{
    const int firstIndex = path.find_last_of("\\/") + 1;
    int lastIndex = path.find_last_of(".") - 1;
    if ( includeExtension || lastIndex <= firstIndex ) {
        lastIndex = path.length() - 1;
    }
    return path.substr( firstIndex, lastIndex - firstIndex + 1 );
}
static inline std::string getFilenameForPath( const char* path, const double includeExtension = true )
{
    std::string strPath( path );
    return getFilenameForPath( strPath, includeExtension );
}

//! Get directory for path.
static inline std::string getDirectoryForPath( const std::string path )
{
    return path.substr( 0, path.find_last_of("\\/") + 1 );
}
static inline std::string getDirectoryForPath( const char* path )
{
    std::string strPath( path );
    return getDirectoryForPath( strPath );
}

//! Get directory of the file from which the function is called.
#define __DIR__ getDirectoryForPath( __FILE__ )

} // namespace tesp

#endif // TESP_INPUT_OUTPUT
