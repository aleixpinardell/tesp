/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TESP_TABULATED_DRAG_COEFFICIENTS_H
#define TESP_TABULATED_DRAG_COEFFICIENTS_H

#include <string>
#include <vector>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tesp {

//! Function to read the altitude and drag coefficients from a file.
/*!
 *  Function to read the altitude and drag coefficients from a file.
 *  \param filePath The path of the file containing the altitude [km] and drag coefficients (separated by tabs).
 *  \return A pair containing a vector of altitudes and a vector of Vector3d with aerodynamic coefficients.
 */
std::pair< std::vector< double >, std::vector< Eigen::Vector3d > > readAltitudesAndDragCoefficientsFromFile(
        const std::string& filePath );


} // namespace tesp

#endif // TESP_TABULATED_DRAG_COEFFICIENTS_H
