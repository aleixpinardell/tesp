/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tabulatedDragCoefficients.h"

namespace tesp {

//! Function to read the altitude and drag coefficients from a file.
std::pair< std::vector< double >, std::vector< Eigen::Vector3d > > readAltitudesAndDragCoefficientsFromFile(
        const std::string& filePath )
{
    std::pair< std::vector< double >, std::vector< Eigen::Vector3d > > altitudesAndAerodynamicCoefficients;

    Eigen::MatrixXd CDh = tudat::input_output::readMatrixFromFile( filePath );
    Eigen::VectorXd altitudes_km = CDh.col(0);
    std::vector< double > altitudes;
    for ( unsigned int i = 0; i < altitudes_km.size( ); i++ )
    {
        altitudes.push_back( altitudes_km(i) * 1.0E3 );
    }
    altitudesAndAerodynamicCoefficients.first = altitudes;

    Eigen::VectorXd dragCoefficients = CDh.col(1);
    std::vector< Eigen::Vector3d > aerodynamicCoefficients;
    for ( unsigned int i = 0; i < dragCoefficients.size( ); i++ )
    {
        aerodynamicCoefficients.push_back( dragCoefficients(i) * Eigen::Vector3d::UnitX( ) );
    }
    altitudesAndAerodynamicCoefficients.second = aerodynamicCoefficients;

    return altitudesAndAerodynamicCoefficients;
}

} // namespace tesp
