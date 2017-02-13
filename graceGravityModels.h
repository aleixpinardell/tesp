/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TESP_GRACE_GRAVITY_MODELS_H
#define TESP_GRACE_GRAVITY_MODELS_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
// #include "tudat/Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"

namespace tesp {

//! List of GRACE gravity field models available
/*!
 *  List of GRACE gravity field models available. Gravity field models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum GraceGravityModel
{
    ggm02c,
    ggm02s
};

//! Derived class of SphericalHarmonicsGravityFieldSettings defining settings for Earth's gravity field according
//! to the a GRACE gravity model (GGM02C up to degree and order 200 or GGM02S up to degree and order 160).
class GraceGravityFieldSettings: public tudat::simulation_setup::SphericalHarmonicsGravityFieldSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param graceGravityModel GRACE gravity model to be used. `ggm02s` for a model based only on GRACE data
     *  or `ggm02c` for a model based on GRACE data and previous data.
     */
    GraceGravityFieldSettings( const GraceGravityModel graceGravityModel ):
        SphericalHarmonicsGravityFieldSettings( createGravityFieldSettingsForGraceModel( graceGravityModel ) )
    {  }

private:
    //! Function to create a gravity field model based on a GRACE model.
    /*!
     *  Function to create a gravity field model based on a GRACE model.
     *  \param graceGravityModel GRACE gravity model to be used. `ggm02s` for a model based only on GRACE data
     *  or `ggm02c` for a model based on GRACE data and previous data.
     *  \return Spherical harmonics gravity field model created according to the specified GRACE model.
     */
    SphericalHarmonicsGravityFieldSettings createGravityFieldSettingsForGraceModel(
            const GraceGravityModel graceGravityModel );

};

} // namespace tesp

#endif // TESP_GRACE_GRAVITY_MODELS_H
