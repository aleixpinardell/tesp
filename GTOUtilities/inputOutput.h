/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef GTO_TOOLS_INPUT_OUTPUT
#define GTO_TOOLS_INPUT_OUTPUT

#include <string>

namespace gto_utilities
{
namespace input_output
{

//! Get path for current directory.
static inline std::string getLastDirectoryForPath( const char *path )
{
    std::string filePath( path );
    return filePath.substr( 0, filePath.find_last_of("\\/") + 1 );
}
#define __DIR__ getLastDirectoryForPath( __FILE__ )

//! Path for GTO output directory.
const static std::string GTO_OUTPUT_PATH = "/Volumes/TarDisk/Documents/Clase/MSc/Thesis/Thesis/gto_output/";

} // namespace input_output
} // namespace gto_utilities

#endif // GTO_TOOLS_INPUT_OUTPUT
