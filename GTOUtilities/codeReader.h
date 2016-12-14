/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef GTO_TOOLS_CODE_READER
#define GTO_TOOLS_CODE_READER

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace gto_utilities
{
namespace code_reader
{

//! List of integer orbit variables that can be changed.
/*!
 *  List of orbital variables that can be changed and are stored as an `int`.
 */
enum IntOrbitVariable
{
    degreeGeopotential          = 0,
    orderGeopotential           = 1,
    propagatorType              = 17,
    integratorType              = 18,
    saveFrequency               = 22
};

//! List of bool orbit variables that can be changed.
/*!
 *  List of orbital variables that can be changed and are stored as a `bool`.
 */
enum BoolOrbitVariable
{
    sunGravity                  = 2,
    moonGravity                 = 3,
    atmosphericDrag             = 4,
    solarRadiationPressure      = 5
};

//! List of double orbit variables that can be changed.
/*!
 *  List of orbital variables that can be changed and are stored as a `double`.
 */
enum DoubleOrbitVariable
{
    bodyMass                    = 6,    // kg
    bodyReferenceArea           = 7,    // m^2
    bodyDragCoefficient         = 8,    // -
    bodySRPDragCoefficient      = 9,    // -
    initialEpoch                = 10,   // years since J2000
    perigeeAltitude             = 11,   // km
    apogeeAltitude              = 12,   // km
    inclination                 = 13,   // deg
    rightAscensionAscendingNode = 14,   // deg
    argumentPerigee             = 15,   // deg
    trueAnomaly                 = 16,   // deg
    stepsize                    = 19,   // s
    maximumPropagationPeriod    = 20,   // years
    integratorErrorTolerance    = 21    // -
};

//! Function that returns the value for an integer variable of a GTO defined by a code identifier.
/*!
 *  Function that returns the value for an integer variable of a GTO defined by a code identifier.
 *  \param codeIdentifier Hexadecimal code identifying the orbit.
 *  \param variable The orbit variable of which the value is desired.
 *  \return The value of the variable for the specified orbit.
 */
int intValueForOrbitWithCode( const std::string& codeIdentifier, const IntOrbitVariable variable );

//! Function that returns the value for a bool variable of a GTO defined by a code identifier.
/*!
 *  Function that returns the value for a bool variable of a GTO defined by a code identifier.
 *  \param codeIdentifier Hexadecimal code identifying the orbit.
 *  \param variable The orbit variable of which the value is desired.
 *  \return The value of the variable for the specified orbit.
 */
bool boolValueForOrbitWithCode(const std::string& codeIdentifier, const BoolOrbitVariable variable);

//! Function that returns the value for a double variable of a GTO defined by a code identifier.
/*!
 *  Function that returns the value for a double variable of a GTO defined by a code identifier.
 *  \param codeIdentifier Hexadecimal code identifying the orbit.
 *  \param variable The orbit variable of which the value is desired.
 *  \return The value of the variable for the specified orbit.
 */
double doubleValueForOrbitWithCode(const std::string& codeIdentifier, const DoubleOrbitVariable variable);

} // namespace code_reader
} // namespace gto_utilities

#endif // GTO_TOOLS_CODE_READER
