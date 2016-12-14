/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "codeReader.h"

namespace gto_utilities
{
namespace code_reader
{

//! Function that returns the value for an integer variable of a GTO defined by a code identifier.
int intValueForOrbitWithCode( const std::string& codeIdentifier, const IntOrbitVariable variable )
{
    int intValue = 0;
    switch ( variable ) {

    case degreeGeopotential: {
        switch ( codeIdentifier[variable] ) {
        case '0':   intValue = 1;      break;
        case '1':   intValue = 2;      break;
        case '2':   intValue = 3;      break;
        case '3':   intValue = 4;      break;
        case '4':   intValue = 5;      break;
        case '5':   intValue = 6;      break;
        case '6':   intValue = 7;      break;
        case '7':   intValue = 8;      break;
        case '8':   intValue = 9;      break;
        case '9':   intValue = 10;     break;
        case 'A':   intValue = 12;     break;
        case 'B':   intValue = 15;     break;
        case 'C':   intValue = 20;     break;
        case 'D':   intValue = 50;     break;
        case 'E':   intValue = 100;    break;
        case 'F':   intValue = 200;    break;
        default:    std::runtime_error( "Did not recognize value for DEGREE OF GEOPOTENTIAL" );
        }
    } break;

    case orderGeopotential: {
        switch ( codeIdentifier[variable] ) {
        case '0':   intValue = intValueForOrbitWithCode(codeIdentifier, degreeGeopotential);      break;
        case '1':   intValue = 0;      break;
        case '2':   intValue = 1;      break;
        case '3':   intValue = 2;      break;
        case '4':   intValue = 3;      break;
        case '5':   intValue = 4;      break;
        case '6':   intValue = 5;      break;
        default:    std::runtime_error( "Did not recognize value for ORDER OF GEOPOTENTIAL" );
        }
    } break;

    case propagatorType: {
        switch ( codeIdentifier[variable] ) {
        case '0':   intValue = 0;       break;
        case '1':   intValue = 1;       break;
        default:    std::runtime_error( "Did not recognize value for PROPAGATOR TYPE" );
        }
    } break;

    case integratorType: {
        switch ( codeIdentifier[variable] ) {
        case '0':   intValue = -1;       break;
        case '1':   intValue = 0;       break;
        case '2':   intValue = 1;       break;
        case '3':   intValue = 2;       break;
        case '4':   intValue = 3;       break;
        default:    std::runtime_error( "Did not recognize value for INTEGRATOR TYPE" );
        }
    } break;

    case saveFrequency: {
        switch ( codeIdentifier[variable] ) {
        case '0':   intValue = 1     ;    break;
        case '1':   intValue = 2     ;    break;
        case '2':   intValue = 5     ;    break;
        case '3':   intValue = 10    ;    break;
        case '4':   intValue = 20    ;    break;
        case '5':   intValue = 50    ;    break;
        case '6':   intValue = 100   ;    break;
        case '7':   intValue = 200   ;    break;
        case '8':   intValue = 500   ;    break;
        case '9':   intValue = 1000  ;    break;
        case 'A':   intValue = 2000  ;    break;
        case 'B':   intValue = 5000  ;    break;
        case 'C':   intValue = 10000 ;    break;
        case 'D':   intValue = 20000 ;    break;
        case 'E':   intValue = 50000 ;    break;
        case 'F':   intValue = 100000;    break;
        default:    std::runtime_error( "Did not recognize value for SAVE FREQUENCY" );
        }
    } break;


    default: std::runtime_error( "Did not recognize variable for int GTO code" );

    }
    return intValue;
}


//! Function that returns the value for a bool variable of a GTO defined by a code identifier.
bool boolValueForOrbitWithCode( const std::string& codeIdentifier, const BoolOrbitVariable variable )
{
    bool boolValue = false;

    switch ( codeIdentifier[variable] ) {
    case '0':   boolValue = false;   break;
    case '1':   boolValue = true;    break;
    default:    std::runtime_error( "Did not recognize variable for bool GTO code" );
    }

    return boolValue;
}


//! Function that returns the value for a double variable of a GTO defined by a code identifier.
double doubleValueForOrbitWithCode( const std::string& codeIdentifier, const DoubleOrbitVariable variable )
{
    double doubleValue = 0.0;
    switch ( variable ) {

    case bodyMass: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 1;            break;
        case '1':   doubleValue = 10;           break;
        case '2':   doubleValue = 100;          break;
        case '3':   doubleValue = 200;          break;
        case '4':   doubleValue = 300;          break;
        case '5':   doubleValue = 400;          break;
        case '6':   doubleValue = 500;          break;
        case '7':   doubleValue = 1E+3;         break;
        case '8':   doubleValue = 2E+3;         break;
        case '9':   doubleValue = 3E+3;         break;
        case 'A':   doubleValue = 4E+3;         break;
        case 'B':   doubleValue = 5E+3;         break;
        case 'C':   doubleValue = 6E+3;         break;
        case 'D':   doubleValue = 8E+3;         break;
        case 'E':   doubleValue = 10E+3;        break;
        case 'F':   doubleValue = 20E+3;        break;
        default:    std::runtime_error( "Did not recognize value for MASS OF ROCKET" );
        }
    } break;

    case bodyReferenceArea: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.1;      break;
        case '1':   doubleValue = 1;        break;
        case '2':   doubleValue = 10;       break;
        case '3':   doubleValue = 15;       break;
        case '4':   doubleValue = 20;       break;
        case '5':   doubleValue = 25;       break;
        case '6':   doubleValue = 50;       break;
        case '7':   doubleValue = 100;      break;
        default:    std::runtime_error( "Did not recognize value for AREA OF ROCKET" );
        }
    } break;

    case bodyDragCoefficient: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = -1.0;     break;
        case '1':   doubleValue = 1.0;      break;
        case '2':   doubleValue = 1.5;      break;
        case '3':   doubleValue = 2.0;      break;
        case '4':   doubleValue = 2.1;      break;
        case '5':   doubleValue = 2.2;      break;
        case '6':   doubleValue = 2.5;      break;
        case '7':   doubleValue = 3.0;      break;
        case '8':   doubleValue = 5.0;      break;
        case '9':   doubleValue = 10.0;     break;
        case 'A':   doubleValue = 1.8;      break;
        case 'B':   doubleValue = 1.9;      break;
        default:    std::runtime_error( "Did not recognize value for CD (DRAG)" );
        }
    } break;

    case bodySRPDragCoefficient: {
        switch ( codeIdentifier[variable] ) {
        case '0': {
            double cd_drag = doubleValueForOrbitWithCode(codeIdentifier, bodyDragCoefficient);
            if ( cd_drag < 0 ) {
                doubleValue = 2.2;
            } else {
                doubleValue = cd_drag;
            }
            break;
        }
        case '1':   doubleValue = 1.0;       break;
        case '2':   doubleValue = 1.5;       break;
        case '3':   doubleValue = 2.0;       break;
        case '4':   doubleValue = 2.1;       break;
        case '5':   doubleValue = 2.2;       break;
        case '6':   doubleValue = 2.5;       break;
        case '7':   doubleValue = 3.0;       break;
        case '8':   doubleValue = 5.0;       break;
        case '9':   doubleValue = 10.0;      break;
        case 'A':   doubleValue = 1.8;       break;
        case 'B':   doubleValue = 1.9;       break;
        default:    std::runtime_error( "Did not recognize value for CD (SRP)" );
        }
    } break;

    case initialEpoch: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.00;     break;
        case '1':   doubleValue = 0.05;     break;
        case '2':   doubleValue = 0.10;     break;
        case '3':   doubleValue = 0.15;     break;
        case '4':   doubleValue = 0.20;     break;
        case '5':   doubleValue = 0.25;     break;
        case '6':   doubleValue = 0.30;     break;
        case '7':   doubleValue = 0.40;     break;
        case '8':   doubleValue = 0.50;     break;
        case '9':   doubleValue = 0.75;     break;
        case 'A':   doubleValue = 1.00;     break;
        case 'B':   doubleValue = 1.25;     break;
        case 'C':   doubleValue = 1.50;     break;
        case 'D':   doubleValue = 1.75;     break;
        case 'E':   doubleValue = 2.00;     break;
        case 'F':   doubleValue = 10.00;    break;
        default:    std::runtime_error( "Did not recognize value for START" );
        }
    } break;

    case perigeeAltitude: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 100;       break;
        case '1':   doubleValue = 125;       break;
        case '2':   doubleValue = 150;       break;
        case '3':   doubleValue = 175;       break;
        case '4':   doubleValue = 200;       break;
        case '5':   doubleValue = 250;       break;
        case '6':   doubleValue = 300;       break;
        case '7':   doubleValue = 400;       break;
        case '8':   doubleValue = 500;       break;
        case '9':   doubleValue = 750;       break;
        case 'A':   doubleValue = 1000;      break;
        case 'B':   doubleValue = 1500;      break;
        case 'C':   doubleValue = 2000;      break;
        case 'D':   doubleValue = 3000;      break;
        case 'E':   doubleValue = 4000;      break;
        case 'F':   doubleValue = 50000;     break;
        default:    std::runtime_error( "Did not recognize value for PERIGEE ALTITUDE" );
        }
    } break;

    case apogeeAltitude: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 35780;     break;
        case '1':   doubleValue = 15000;     break;
        case '2':   doubleValue = 16000;     break;
        case '3':   doubleValue = 17000;     break;
        case '4':   doubleValue = 18000;     break;
        case '5':   doubleValue = 19000;     break;
        case '6':   doubleValue = 20000;     break;
        case '7':   doubleValue = 22000;     break;
        case '8':   doubleValue = 25000;     break;
        case '9':   doubleValue = 30000;     break;
        case 'A':   doubleValue = 32000;     break;
        case 'B':   doubleValue = 33000;     break;
        case 'C':   doubleValue = 34000;     break;
        case 'D':   doubleValue = 35000;     break;
        case 'E':   doubleValue = 36000;     break;
        case 'F':   doubleValue = 37000;     break;
        default:    std::runtime_error( "Did not recognize value for APOGEE ALTITUDE" );
        }
    } break;

    case inclination: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.0;      break;
        case '1':   doubleValue = 1.0;      break;
        case '2':   doubleValue = 2.0;      break;
        case '3':   doubleValue = 5.0;      break;
        case '4':   doubleValue = 7.0;      break;
        case '5':   doubleValue = 10.0;     break;
        case '6':   doubleValue = 15.0;     break;
        case '7':   doubleValue = 20.0;     break;
        case '8':   doubleValue = 23.4;     break;
        case '9':   doubleValue = 30.0;     break;
        case 'A':   doubleValue = 40.0;     break;
        case 'B':   doubleValue = 50.0;     break;
        case 'C':   doubleValue = 60.0;     break;
        case 'D':   doubleValue = 70.0;     break;
        case 'E':   doubleValue = 80.0;     break;
        case 'F':   doubleValue = 90.0;     break;
        default:    std::runtime_error( "Did not recognize value for INCLINATION" );
        }
    } break;

    case rightAscensionAscendingNode: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.0;         break;
        case '1':   doubleValue = 2.0;         break;
        case '2':   doubleValue = 5.0;         break;
        case '3':   doubleValue = 10.0;        break;
        case '4':   doubleValue = 25.0;        break;
        case '5':   doubleValue = 45.0;        break;
        case '6':   doubleValue = 75.0;        break;
        case '7':   doubleValue = 90.0;        break;
        case '8':   doubleValue = 120.0;       break;
        case '9':   doubleValue = 150.0;       break;
        case 'A':   doubleValue = 180.0;       break;
        case 'B':   doubleValue = 210.0;       break;
        case 'C':   doubleValue = 240.0;       break;
        case 'D':   doubleValue = 270.0;       break;
        case 'E':   doubleValue = 300.0;       break;
        case 'F':   doubleValue = 330.0;       break;
        default:    std::runtime_error( "Did not recognize value for RAAN" );
        }
    } break;

    case argumentPerigee: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.0;         break;
        case '1':   doubleValue = 2.0;         break;
        case '2':   doubleValue = 5.0;         break;
        case '3':   doubleValue = 10.0;        break;
        case '4':   doubleValue = 25.0;        break;
        case '5':   doubleValue = 45.0;        break;
        case '6':   doubleValue = 75.0;        break;
        case '7':   doubleValue = 90.0;        break;
        case '8':   doubleValue = 120.0;       break;
        case '9':   doubleValue = 150.0;       break;
        case 'A':   doubleValue = 180.0;       break;
        case 'B':   doubleValue = 210.0;       break;
        case 'C':   doubleValue = 240.0;       break;
        case 'D':   doubleValue = 270.0;       break;
        case 'E':   doubleValue = 300.0;       break;
        case 'F':   doubleValue = 330.0;       break;
        default:    std::runtime_error( "Did not recognize value for ARGUMENT OF PERIGEE" );
        }
    } break;

    case trueAnomaly: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.0;         break;
        case '1':   doubleValue = 2.0;         break;
        case '2':   doubleValue = 5.0;         break;
        case '3':   doubleValue = 10.0;        break;
        case '4':   doubleValue = 25.0;        break;
        case '5':   doubleValue = 45.0;        break;
        case '6':   doubleValue = 75.0;        break;
        case '7':   doubleValue = 90.0;        break;
        case '8':   doubleValue = 120.0;       break;
        case '9':   doubleValue = 150.0;       break;
        case 'A':   doubleValue = 180.0;       break;
        case 'B':   doubleValue = 210.0;       break;
        case 'C':   doubleValue = 240.0;       break;
        case 'D':   doubleValue = 270.0;       break;
        case 'E':   doubleValue = 300.0;       break;
        case 'F':   doubleValue = 330.0;       break;
        default:    std::runtime_error( "Did not recognize value for TRUE ANOMALY" );
        }
    } break;

    case stepsize: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 1;         break;
        case '1':   doubleValue = 10;        break;
        case '2':   doubleValue = 20;        break;
        case '3':   doubleValue = 30;        break;
        case '4':   doubleValue = 40;        break;
        case '5':   doubleValue = 50;        break;
        case '6':   doubleValue = 60;        break;
        case '7':   doubleValue = 80;        break;
        case '8':   doubleValue = 100;       break;
        case '9':   doubleValue = 150;       break;
        case 'A':   doubleValue = 200;       break;
        case 'B':   doubleValue = 500;       break;
        case 'C':   doubleValue = 1000;      break;
        case 'D':   doubleValue = 3600;      break;
        case 'E':   doubleValue = 14400;     break;
        case 'F':   doubleValue = 86400;     break;
        default:    std::runtime_error( "Did not recognize value for STEPSIZE" );
        }
    } break;

    case maximumPropagationPeriod: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 0.5;     break;
        case '1':   doubleValue = 1;       break;
        case '2':   doubleValue = 2;       break;
        case '3':   doubleValue = 3;       break;
        case '4':   doubleValue = 4;       break;
        case '5':   doubleValue = 5;       break;
        case '6':   doubleValue = 6;       break;
        case '7':   doubleValue = 7;       break;
        case '8':   doubleValue = 8;       break;
        case '9':   doubleValue = 9;       break;
        case 'A':   doubleValue = 10;      break;
        case 'B':   doubleValue = 12;      break;
        case 'C':   doubleValue = 15;      break;
        case 'D':   doubleValue = 25;      break;
        case 'E':   doubleValue = 50;      break;
        case 'F':   doubleValue = 100;     break;
        default:    std::runtime_error( "Did not recognize value for MAX PROP TIME" );
        }
    } break;

    case integratorErrorTolerance: {
        switch ( codeIdentifier[variable] ) {
        case '0':   doubleValue = 1.00E-12;     break;
        case '1':   doubleValue = 1.00E-01;     break;
        case '2':   doubleValue = 1.00E-02;     break;
        case '3':   doubleValue = 1.00E-03;     break;
        case '4':   doubleValue = 1.00E-04;     break;
        case '5':   doubleValue = 1.00E-05;     break;
        case '6':   doubleValue = 1.00E-06;     break;
        case '7':   doubleValue = 1.00E-07;     break;
        case '8':   doubleValue = 1.00E-08;     break;
        case '9':   doubleValue = 1.00E-09;     break;
        case 'A':   doubleValue = 1.00E-10;     break;
        case 'B':   doubleValue = 1.00E-11;     break;
        case 'C':   doubleValue = 1.00E-13;     break;
        case 'D':   doubleValue = 1.00E-14;     break;
        case 'E':   doubleValue = 1.00E-15;     break;
        case 'F':   doubleValue = 1.00E-20;     break;
        default:    std::runtime_error( "Did not recognize value for INTEGRATOR ERROR TOLERANCE" );
        }
    } break;


    default: std::runtime_error( "Did not recognize variable for double GTO code" );

    }
    return doubleValue;
}


} // namespace code_reader
} // namespace gto_utilities
