/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "graceGravityModels.h"

#include "inputOutput.h"

namespace tesp {

//! Function to create a gravity field model based on a GRACE model.
tudat::simulation_setup::SphericalHarmonicsGravityFieldSettings
GraceGravityFieldSettings::createGravityFieldSettingsForGraceModel( const GraceGravityModel graceGravityModel )
{
    std::string filename;
    switch ( graceGravityModel ) {
        case ggm02c: filename = "GGM02C.txt"; break;
        case ggm02s: filename = "GGM02S.txt"; break;
        default: throw std::runtime_error( "GRACE gravity model not recognised" );
    }
    Eigen::MatrixXd coefficients = tudat::input_output::readMatrixFromFile( __DIR__ + filename );

    Eigen::VectorXd leftColumn = coefficients.col(0);
    Eigen::VectorXd rightColumn = coefficients.col(1);

    double gravitationalParameter = leftColumn(0);
    double referenceRadius = rightColumn(0);
    int degreeSquared = leftColumn.size() - 1;
    int degree = sqrt(degreeSquared);
    int order = degree;

    Eigen::VectorXd cosines = leftColumn.segment(1,degreeSquared);
    Eigen::VectorXd sines = rightColumn.segment(1,degreeSquared);

    Eigen::MatrixXd cosineCoefficients = Eigen::Map<Eigen::MatrixXd>(cosines.data(), degree, order).transpose();
    Eigen::MatrixXd sineCoefficients = Eigen::Map<Eigen::MatrixXd>(sines.data(), degree, order).transpose();

    return SphericalHarmonicsGravityFieldSettings(
                gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, "IAU_Earth" );
}

} // namespace tesp
