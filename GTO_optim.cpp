/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "GTOUtilities/inputOutput.h"
#include "GTOUtilities/codeReader.h"
#include "GTOUtilities/graceGravityModels.h"
#include "GTOUtilities/tabulatedDragCoefficients.h"


//! Execute propagation of orbit of Rocket around the Earth.
int main()
{

    ////////////////////////////////////////////////////////////////////////////
    ///////////            USING STATEMENTS              ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace gravitation;
    using namespace interpolators;
    using namespace input_output;

    using namespace gto_utilities::input_output;
    using namespace gto_utilities::code_reader;
    using namespace gto_utilities::environment;


    /////////////////////////////////////////////////////////////////////////////////
    ////////     CREATE ENVIRONMENT AND VEHICLE       ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    // Read and translate settings from TXT file
    std::string codeFilename = "currentOptimisation";

    std::ifstream file;
    std::string line;
    file.open( __DIR__ + codeFilename + ".txt" );
    std::string originalSettingsCode;
    std::string settingsCode;
    std::map< char, std::vector< double > > xyzValues;
    std::map< char, std::vector< int > > xyzIndexes;
    while ( std::getline(file,line) )
    {
        if ( line[0] != '%' )
        {
            std::vector< std::string > parts;
            boost::split(parts, line, boost::is_any_of("\t"));
            if ( line[0] == '#' )
            {
                // The letter to be replaced by values. X, Y or Z (maximum three-variable optimisation)
                char hashLetter = line[1];
                if ( hashLetter != 'X' || hashLetter != 'Y' || hashLetter != 'Z' )
                {
                    std::runtime_error( "Unrecognised hash letter for optimisation (only X Y Z accepted)" );
                }

                // How are the values specified
                //  :   -->     first   last    spacing
                //  !   -->     first   last    count
                //  =   -->     array of values separated by tabs
                std::string specification = parts[1];
                std::vector< double > values;
                if ( specification == ":" || specification == "!" )
                {
                    double first = std::stod( parts[2] );
                    double last  = std::stod( parts[3] );

                    double spacing = std::stod( parts[4] );
                    if ( specification == "!" )
                    {
                        spacing = ( last - first ) / ( spacing - 1 );
                    }
                    double currentValue = first;
                    while ( currentValue < last )
                    {
                        values.push_back( currentValue );
                        currentValue += spacing;
                    }
                    values.push_back( last );
                }
                else if ( specification == "=" )
                {
                    for ( int i = 2; i < parts.size(); i++ )
                    {
                        values.push_back( std::stod( parts[i]) );
                    }
                }
                else
                {
                    std::runtime_error( "Unrecognised optimisation specification character (only : ! = accepted)" );
                }

                xyzValues[ hashLetter ] = values;
            }
            else
            {
                originalSettingsCode = parts[0];
                settingsCode = originalSettingsCode;
                for ( int i = 0; i < settingsCode.size(); i++ )
                {
                    if ( settingsCode.at(i) == 'X' || settingsCode.at(i) == 'Y' || settingsCode.at(i) == 'Z' )
                    {
                        char key = settingsCode.at(i);
                        std::vector < int > currentIndexes;
                        if ( xyzIndexes.count( key ) )
                        {
                            currentIndexes = xyzIndexes[ key ];
                        }
                        currentIndexes.push_back( i );
                        xyzIndexes[ key ] = currentIndexes;

                        settingsCode.replace( i, 1, "0" );
                    }
                }
            }
        }
    }

    if ( xyzValues.count( 'X' ) == 0 ) { xyzValues[ 'X' ] = { 0.0 }; }
    if ( xyzValues.count( 'Y' ) == 0 ) { xyzValues[ 'Y' ] = { 0.0 }; }
    if ( xyzValues.count( 'Z' ) == 0 ) { xyzValues[ 'Z' ] = { 0.0 }; }

    for ( unsigned int xIndex = 0; xIndex < xyzValues['X'].size(); xIndex++ )
    {
        double xValue = xyzValues['X'][xIndex];

        for ( unsigned int yIndex = 0; yIndex < xyzValues['Y'].size(); yIndex++ )
        {
            double yValue = xyzValues['Y'][yIndex];

            for ( unsigned int zIndex = 0; zIndex < xyzValues['Z'].size(); zIndex++ )
            {
                double zValue = xyzValues['Z'][zIndex];

                std::cout << std::endl << "Propagating orbit: " << originalSettingsCode << std::endl;
                std::cout << "X: " << xValue << ",  Y: " << yValue << ",  Z: " << zValue << std::endl;

                // Read all double variables using code_reader first
                std::map< DoubleOrbitVariable, double > variablesMap;
                std::vector< DoubleOrbitVariable > doubleVars = { bodyMass, bodyReferenceArea, bodyDragCoefficient, bodySRPDragCoefficient, initialEpoch, perigeeAltitude, apogeeAltitude, inclination, rightAscensionAscendingNode, argumentPerigee, trueAnomaly, stepsize, maximumPropagationPeriod, integratorErrorTolerance };
                for ( int i = 0; i < doubleVars.size(); i++ )
                {
                    DoubleOrbitVariable doubleVar = doubleVars[i];
                    variablesMap[ doubleVar ] = doubleValueForOrbitWithCode(settingsCode, doubleVar);
                }

                // Then replace those which where specified originally by X, Y or Z by current values of loops
                for( std::map< char, std::vector< int > >::iterator iter = xyzIndexes.begin(); iter != xyzIndexes.end(); ++iter )
                {
                    char hashLetter = iter->first;
                    std::vector< int > indexes = iter->second;

                    double currentValue;
                    switch ( hashLetter ) {
                    case 'X': currentValue = xValue; break;
                    case 'Y': currentValue = yValue; break;
                    case 'Z': currentValue = zValue; break;
                    default: break;
                    }

                    for ( int j = 0; j < indexes.size(); j++ )
                    {
                        variablesMap[ static_cast<DoubleOrbitVariable>( indexes[j] ) ] = currentValue;
                    }
                }


                // Perturbations
                int n = intValueForOrbitWithCode(settingsCode, degreeGeopotential);
                int m = intValueForOrbitWithCode(settingsCode, orderGeopotential);
                bool sunG = boolValueForOrbitWithCode(settingsCode, sunGravity);
                bool moonG = boolValueForOrbitWithCode(settingsCode, moonGravity);
                bool drag = boolValueForOrbitWithCode(settingsCode, atmosphericDrag);
                bool srp = boolValueForOrbitWithCode(settingsCode, solarRadiationPressure);

                // Body
                double m_rocket = variablesMap[ bodyMass ];
                double A_rocket = variablesMap[ bodyReferenceArea ];
                double CD_drag = variablesMap[ bodyDragCoefficient ];
                double CD_srp = variablesMap[ bodySRPDragCoefficient ];

                // Initial state
                double delay_years = variablesMap[ initialEpoch ];
                double h_p_km = variablesMap[ perigeeAltitude ];
                double h_a_km = variablesMap[ apogeeAltitude ];
                double inc = variablesMap[ inclination ];
                double raan = variablesMap[ rightAscensionAscendingNode ];
                double arg_perigee = variablesMap[ argumentPerigee ];
                double true_anom = variablesMap[ trueAnomaly ];

                // Propagation / integration settings
                int ptype = intValueForOrbitWithCode(settingsCode, propagatorType);
                TranslationalPropagatorType propagator_type = static_cast<TranslationalPropagatorType>(ptype);

                int itype = intValueForOrbitWithCode(settingsCode, integratorType);
                AvailableIntegrators integrator_type;
                RungeKuttaCoefficients::CoefficientSets rungekutta_type;
                if ( itype == -1 ) {
                    integrator_type = rungeKutta4;
                } else {
                    integrator_type = rungeKuttaVariableStepSize;
                    rungekutta_type = static_cast<RungeKuttaCoefficients::CoefficientSets>(itype);
                }

                double stepSize = variablesMap[ stepsize ];
                double integratorTolerance = variablesMap[ integratorErrorTolerance ];
                double maxPropPeriod = variablesMap[ maximumPropagationPeriod ];
                int exportFrequency = intValueForOrbitWithCode(settingsCode, saveFrequency);

//                std::cout << n << std::endl;
//                std::cout << m << std::endl;
//                std::cout << sunG << std::endl;
//                std::cout << moonG << std::endl;
//                std::cout << drag << std::endl;
//                std::cout << srp << std::endl;
//                std::cout << m_rocket << std::endl;
//                std::cout << A_rocket << std::endl;
//                std::cout << CD_drag << std::endl;
//                std::cout << CD_srp << std::endl;
//                std::cout << delay_years << std::endl;
//                std::cout << h_p_km << std::endl;
//                std::cout << h_a_km << std::endl;
//                std::cout << inc << std::endl;
//                std::cout << raan << std::endl;
//                std::cout << arg_perigee << std::endl;
//                std::cout << true_anom << std::endl;
//                std::cout << stepSize << std::endl;
//                std::cout << maxPropPeriod << std::endl;
//                std::cout << integratorTolerance << std::endl;
//                std::cout << exportFrequency << std::endl;


                // Set simulation time settings.
                const double year = tudat::physical_constants::JULIAN_YEAR;
                const double delay = delay_years * year;
                const double simulationStartEpoch = delay;
                const double simulationEndEpoch = delay + maxPropPeriod*year;

                // Define body settings for simulation.
                std::vector< std::string > bodiesToCreate;
                bodiesToCreate.push_back( "Sun" );
                bodiesToCreate.push_back( "Earth" );
                bodiesToCreate.push_back( "Moon" );

                // Create body objects.
                std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings(
                            bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

                for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
                {
                    bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
                    bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
                }

                // EARTH
                bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GraceGravityFieldSettings >( ggm02c );
                bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

                // MOON
                bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

                NamedBodyMap bodyMap = createBodies( bodySettings );


                ////////////////////////////////////////////////////////////////////////
                //////             CREATE VEHICLE            ///////////////////////////
                ////////////////////////////////////////////////////////////////////////

                // Create spacecraft object.
                bodyMap[ "Rocket" ] = boost::make_shared< simulation_setup::Body >( );
                bodyMap[ "Rocket" ]->setConstantBodyMass( m_rocket );
                double referenceArea = A_rocket;

                boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;
                if ( CD_drag < 0.0 )
                {
                    // Read CD as a function of altitude
                    auto hCD = readAltitudesAndDragCoefficientsFromFile( __DIR__ + "tabulatedDragCoefficients.txt" );
                    std::vector< double > altitudes = hCD.first;
                    std::vector< Eigen::Vector3d > aerodynamicCoefficients = hCD.second;

                    // Create interpolator
                    boost::shared_ptr< InterpolatorSettings > interpolatorSettings = boost::make_shared< InterpolatorSettings >( OneDimensionalInterpolatorTypes::linear_interpolator );

                    // Tabulated aerodynamic settings
                    aerodynamicCoefficientSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >( altitudes, aerodynamicCoefficients, referenceArea, aerodynamics::AerodynamicCoefficientsIndependentVariables::altitude_dependent, interpolatorSettings, 1, 1 );
                }
                else
                {
                    // Constant aerodynamic settings
                    aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >( referenceArea, CD_drag * Eigen::Vector3d::UnitX( ), 1, 1 );
                }
                // Aerodynamics interface
                bodyMap[ "Rocket" ]->setAerodynamicCoefficientInterface( createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Rocket" ) );


                // SRP interface
                boost::shared_ptr< RadiationPressureInterfaceSettings > rocketRadiationPressureSettings = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", referenceArea, CD_srp, boost::assign::list_of( "Earth" ) );
                bodyMap[ "Rocket" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface( rocketRadiationPressureSettings, "Rocket", bodyMap ) );


                // Finalize body creation.
                setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

                ////////////////////////////////////////////////////////////////////////////////
                //////            CREATE ACCELERATIONS          ////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                // Define propagator settings variables.
                SelectedAccelerationMap accelerationMap;
                std::vector< std::string > bodiesToPropagate;
                std::vector< std::string > centralBodies;

                // Define propagation settings.
                std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfRocket;

                accelerationsOfRocket[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( n, m ) );

                if ( sunG )
                {
                    accelerationsOfRocket[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                }

                if ( moonG )
                {
                    accelerationsOfRocket[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                }

                if ( srp )
                {
                    accelerationsOfRocket[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );
                }

                if ( drag )
                {
                    accelerationsOfRocket[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic ) );
                }

                accelerationMap[ "Rocket" ] = accelerationsOfRocket;
                bodiesToPropagate.push_back( "Rocket" );
                centralBodies.push_back( "Earth" );

                basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


                ////////////////////////////////////////////////////////////////////////////
                //////             CREATE PROPAGATION SETTINGS            //////////////////
                ////////////////////////////////////////////////////////////////////////////

                // Set Keplerian elements for Rocket.
                double meanEarthRadius = 6371.0E+3;
                double h_p = h_p_km*1.0E3;
                double h_a = h_a_km*1.0E3;
                double a = (h_a + h_p)/2.0 + meanEarthRadius;
                double e = (h_a - h_p)/(h_a + h_p + 2.0 * meanEarthRadius);

                Vector6d rocketInitialStateInKeplerianElements;
                rocketInitialStateInKeplerianElements( semiMajorAxisIndex ) = a;
                rocketInitialStateInKeplerianElements( eccentricityIndex ) = e;
                rocketInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( inc );
                rocketInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( arg_perigee );
                rocketInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( raan );
                rocketInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( true_anom );

                double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
                const Vector6d rocketInitialState = convertKeplerianToCartesianElements( rocketInitialStateInKeplerianElements, earthGravitationalParameter );


                // Hybrid termination conditions
                std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;

                // Time limit
                constituentSettings.push_back( boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

                // Altitude limit
                boost::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings = boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >( boost::make_shared< propagators::SingleDependentVariableSaveSettings >( propagators::altitude_dependent_variable, "Rocket" ), 100.0E3, 1 );
                constituentSettings.push_back( altitudeTerminationSettings );

                // Stop if ANY of the two is met
                boost::shared_ptr< PropagationTerminationSettings > terminationSettings = boost::make_shared< propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

                // Save dependent variables
                std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;
                // dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Rocket" ) );
                // dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( aerodynamic_moment_coefficients_dependent_variable, "Rocket" ) );
                // dependentVariablesToSave.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >( spherical_harmonic_gravity, "Rocket", "Earth", 0 ) );
                // dependentVariablesToSave.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >( third_body_central_gravity, "Rocket", "Sun", 0 ) );
                // dependentVariablesToSave.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >( third_body_central_gravity, "Rocket", "Moon", 0 ) );

                // dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable, "Rocket", "Earth" ) );
                // dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable, "Sun", "Earth" ) );

                boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings = boost::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave, 0 ) ;

                // Translational propagator settings
                boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > > ( centralBodies, accelerationModelMap, bodiesToPropagate, rocketInitialState, terminationSettings, propagator_type );

                // Create multiple propagator settings
                std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsList;
                propagatorSettingsList.push_back( translationalPropagatorSettings );

                boost::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings = boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsList, terminationSettings, dependentVariableSaveSettings );

                boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                if ( integrator_type == rungeKutta4 )
                {
                    integratorSettings = boost::make_shared< IntegratorSettings< > > ( rungeKutta4, simulationStartEpoch, stepSize, exportFrequency );
                }
                else
                {
                    integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >( rungeKuttaVariableStepSize, simulationStartEpoch, stepSize, rungekutta_type, 1.0E-3, 1.0E+4, integratorTolerance, integratorTolerance, exportFrequency );
                }


                /////////////////////////////////////////////////////////////////////////////////
                //////////             PROPAGATE ORBIT            ///////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////

                auto t_ini = std::chrono::steady_clock::now();

                // Create simulation object and propagate dynamics.
                SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings, true, false, false );
                std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput = dynamicsSimulator.getDependentVariableHistory( );

                auto t_end = std::chrono::steady_clock::now();
                auto t_comp = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_ini).count() * 1.0e-3;
                auto t_prop = ( (--integrationResult.end( ) )->first - simulationStartEpoch ) / tudat::physical_constants::JULIAN_DAY;
                auto prop_rate = t_prop / t_comp;
                std::cout << "Computation time: " << t_comp << " seconds" << std::endl;
                std::cout << "Propagation time: " << t_prop << " days" << std::endl;
                std::cout << "Propagation rate: " << prop_rate << " days/second" << std::endl;


                ////////////////////////////////////////////////////////////////////////////////
                ////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                // Convert history to Keplerian
                std::map< double, Eigen::VectorXd > keplerianResult = integrationResult;
                for(std::map< double, Eigen::VectorXd >::iterator iter = integrationResult.begin(); iter != integrationResult.end(); ++iter)
                {
                    double key = iter->first;
                    Vector6d value = iter->second;
                    keplerianResult[key] = convertCartesianToKeplerianElements( value, earthGravitationalParameter );
                }
                Eigen::VectorXd finalKeplerianState = ( --keplerianResult.end( ) )->second;
                double a_f = finalKeplerianState(0);
                double e_f = finalKeplerianState(1);
                double h_p_f = (a_f * ( 1.0 - e_f ) - meanEarthRadius)*1.0E-3;

                auto nlim = std::numeric_limits< double >::digits10;

                // Create directory for all the results of the orbits used for this optimisation
                std::string folder = "propagated/optim/" + originalSettingsCode + "/";
                boost::filesystem::path dir( GTO_OUTPUT_PATH + folder );
                boost::filesystem::create_directory( dir );

                // Filename with values of X Y Z
                std::stringstream xStream;
                xStream << std::fixed << std::setprecision(8) << xValue;
                std::string xString = xStream.str();

                std::stringstream yStream;
                yStream << std::fixed << std::setprecision(8) << yValue;
                std::string yString = yStream.str();

                std::stringstream zStream;
                zStream << std::fixed << std::setprecision(8) << zValue;
                std::string zString = zStream.str();

                std::string filename = "X" + xString + "Y" + yString + "Z" + zString;

                // Write propagation and computation times to file.
                std::string resultFilename = folder + filename + "_result.txt";
                std::map< int, double > resultMap;
                resultMap[0] = t_prop;
                resultMap[1] = t_comp;
                resultMap[2] = h_p_f;
                writeDataMapToTextFile( resultMap, resultFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

                // Write perturbed satellite propagation history to file.
                std::string historyFilename = folder + filename + "_history.txt";
                writeDataMapToTextFile( keplerianResult, historyFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

                // Write additional variables to file.
                // std::string variablesFilename = folder + filename + "_variables.txt";
                // writeDataMapToTextFile( dependentVariableOutput, variablesFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

                std::cout << std::endl;

            }   // end Z for
        }       // end Y for
    }           // end X for

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}


