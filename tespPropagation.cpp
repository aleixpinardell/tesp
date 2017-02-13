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
#include <chrono>

#include "tespSettings.h"
#include "inputOutput.h"
#include "graceGravityModels.h"
#include "tabulatedDragCoefficients.h"


int main( int argc, char* argv[] )
{
    if ( argc != 2 )
    {
        std::cerr << "Usage: tesp \"relative or absolute path to a 'tespin' file\"" << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        ////////////////////////////////////////////////////////////////////////////
        ///////////            USING STATEMENTS              ///////////////////////
        ////////////////////////////////////////////////////////////////////////////

        using namespace tudat;
        using namespace tudat::simulation_setup;
        using namespace tudat::propagators;
        using namespace tudat::numerical_integrators;
        using namespace tudat::orbital_element_conversions;
        using namespace tudat::basic_mathematics;
        using namespace tudat::basic_astrodynamics;
        using namespace tudat::gravitation;
        using namespace tudat::interpolators;
        using namespace tudat::input_output;

        using namespace tesp;


        ////////////////////////////////////////////////////////////////////////////
        ///////////      READ SETTINGS FROM INPUT FILE       ///////////////////////
        ////////////////////////////////////////////////////////////////////////////

        std::string inputFilePath = argv[1];
        TespSettings settings( inputFilePath );

        if ( ! settings.checkRequestedOutputExists() )
        {
            ////////////////////////////////////////////////////////////////////////
            //////           CREATE ENVIRONMENT          ///////////////////////////
            ////////////////////////////////////////////////////////////////////////

            auto t_ini = std::chrono::steady_clock::now();

            // Load Spice kernels.
            spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
            spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
            spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

            // Set simulation time settings.
            double simulationStartEpoch = settings.initialEpoch;
            const double simulationEndEpoch = settings.maximumFinalEpoch();

            // Define body settings for simulation.
            std::vector< std::string > bodiesToCreate;
            bodiesToCreate.push_back( "Sun" );
            bodiesToCreate.push_back( "Earth" );
            bodiesToCreate.push_back( "Moon" );

            // Create body objects.
            std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
            if ( settings.preloadCelestialBodiesData )
            {
                bodySettings = getDefaultBodySettings(
                            bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
            }
            else
            {
                bodySettings = getDefaultBodySettings( bodiesToCreate );
            }

            for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
            {
                std::string body = bodiesToCreate.at( i );
                bodySettings[ body ]->ephemerisSettings->resetFrameOrientation( "J2000" );
                bodySettings[ body ]->rotationModelSettings->resetOriginalFrame( "J2000" );
            }

            // EARTH
            bodySettings[ "Earth" ]->gravityFieldSettings =
                    boost::make_shared< GraceGravityFieldSettings >( settings.geopotentialModel );

            if ( settings.atmosphericModel == exponential_atmosphere )
            {
                const double densityScaleHeight = 7.99E3;
                const double constantTemperature = 273.0;
                const double densityAtZeroAltitude = 1.225;
                const double specificGasConstant = 287.058;
                bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >(
                            densityScaleHeight, constantTemperature, densityAtZeroAltitude, specificGasConstant );
            }
            else
            {
                bodySettings[ "Earth" ]->atmosphereSettings =
                        boost::make_shared< AtmosphereSettings >( nrlmsise00 );
            }

            // MOON
            bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

            // Create bodies
            NamedBodyMap bodyMap = createBodies( bodySettings );


            ////////////////////////////////////////////////////////////////////////
            //////             CREATE VEHICLE            ///////////////////////////
            ////////////////////////////////////////////////////////////////////////

            // Create spacecraft object.
            bodyMap[ "Body" ] = boost::make_shared< simulation_setup::Body >( );
            bodyMap[ "Body" ]->setConstantBodyMass( settings.bodyMass );

            boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;
            if ( settings.bodyDragCoefficient < 0.0 )
            {
                // Read CD as a function of altitude
                auto hCD = readAltitudesAndDragCoefficientsFromFile( __DIR__ + "tabulatedDragCoefficients.txt" );
                std::vector< double > altitudes = hCD.first;
                std::vector< Eigen::Vector3d > aerodynamicCoefficients = hCD.second;

                // Create interpolator
                boost::shared_ptr< InterpolatorSettings > interpolatorSettings = boost::make_shared<
                        InterpolatorSettings >( OneDimensionalInterpolatorTypes::linear_interpolator );

                // Tabulated aerodynamic settings
                aerodynamicCoefficientSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                            altitudes, aerodynamicCoefficients, settings.bodyCrossSectionalArea,
                            aerodynamics::AerodynamicCoefficientsIndependentVariables::altitude_dependent,
                            interpolatorSettings, 1, 1 );
            }
            else
            {
                // Constant aerodynamic settings
                aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                            settings.bodyCrossSectionalArea,
                            settings.bodyDragCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );
            }
            // Aerodynamics interface
            bodyMap[ "Body" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Body" ) );


            // SRP interface
            std::vector< std::string > occultingBodies = {};
            if ( ! settings.ignoreEclipses )
            {
                occultingBodies.push_back( "Earth" );
            }
            boost::shared_ptr< RadiationPressureInterfaceSettings > bodyRadiationPressureSettings =
                    boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                        "Sun", settings.bodyCrossSectionalArea, settings.bodyRadiationPressureCoefficient,
                        occultingBodies );

            bodyMap[ "Body" ]->setRadiationPressureInterface(
                        "Sun", createRadiationPressureInterface( bodyRadiationPressureSettings, "Body", bodyMap ) );


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

            accelerationsOfRocket[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >(
                                                            settings.geopotentialDegree, settings.geopotentialOrder ) );

            if ( settings.thirdBodyAttractionSun )
            {
                accelerationsOfRocket[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                              basic_astrodynamics::central_gravity ) );
            }

            if ( settings.thirdBodyAttractionMoon )
            {
                accelerationsOfRocket[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                               basic_astrodynamics::central_gravity ) );
            }

            if ( settings.solarRadiationPressure )
            {
                accelerationsOfRocket[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                              basic_astrodynamics::cannon_ball_radiation_pressure ) );
            }

            if ( settings.atmosphericDrag )
            {
                accelerationsOfRocket[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                basic_astrodynamics::aerodynamic ) );
            }

            accelerationMap[ "Body" ] = accelerationsOfRocket;
            bodiesToPropagate.push_back( "Body" );
            centralBodies.push_back( "Earth" );

            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


            ////////////////////////////////////////////////////////////////////////////
            //////             CREATE PROPAGATION SETTINGS            //////////////////
            ////////////////////////////////////////////////////////////////////////////

            // Set Keplerian elements for Rocket.

            Vector6d bodyInitialState;
            Vector6d bodyInitialStateKep;

            double earthGravitationalParameter =
                    bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

            // Check if the initial epoch is read from inut file or the propagation will be resumed from old output
            if ( settings.resuming )
            {
                simulationStartEpoch = settings.resumingEpoch;

                if ( settings.resumingStateIsCartesian )
                {
                    bodyInitialState = settings.resumingState;
                }
                else
                {
                    bodyInitialStateKep = settings.resumingState;
                }
            }
            else
            {
                bodyInitialStateKep( semiMajorAxisIndex ) = settings.initialSemimajorAxis * 1.0E3;
                bodyInitialStateKep( eccentricityIndex ) = settings.initialEccentricity;
                bodyInitialStateKep( inclinationIndex ) =
                        unit_conversions::convertDegreesToRadians( settings.initialInclination );
                bodyInitialStateKep( argumentOfPeriapsisIndex ) =
                        unit_conversions::convertDegreesToRadians( settings.initialArgumentPerigee );
                bodyInitialStateKep( longitudeOfAscendingNodeIndex ) =
                        unit_conversions::convertDegreesToRadians( settings.initialLongitudeAscendingNode );
                bodyInitialStateKep( trueAnomalyIndex ) =
                        unit_conversions::convertDegreesToRadians( settings.initialTrueAnomaly );
            }

            if ( ( settings.resuming && ! settings.resumingStateIsCartesian ) || ! settings.resuming )
            {
                bodyInitialState = convertKeplerianToCartesianElements(
                            bodyInitialStateKep, earthGravitationalParameter );
            }


            // Hybrid termination conditions
            std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;

            // Time limit
            constituentSettings.push_back(
                        boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

            if ( settings.reentryAltitude > 0 )
            {
                // Altitude limit
                boost::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings = boost::make_shared<
                        propagators::PropagationDependentVariableTerminationSettings >(
                            boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                propagators::altitude_dependent_variable, "Body" ),
                            settings.reentryAltitude * 1.0E3, 1 );

                constituentSettings.push_back( altitudeTerminationSettings );
            }

            // Stop if ANY of the two is met
            boost::shared_ptr< PropagationTerminationSettings > terminationSettings = boost::make_shared<
                    propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

            // Save dependent variables
            std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;

            if ( settings.outputSunPosition )
            {
                dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                        relative_position_dependent_variable, "Sun", "Earth" ) );
            }

            if ( settings.outputMoonPosition )
            {
                dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                        relative_position_dependent_variable, "Moon", "Earth" ) );
            }

            boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave, 0 ) ;

            // Translational propagator settings
            boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > > (
                        centralBodies, accelerationModelMap, bodiesToPropagate, bodyInitialState,
                        terminationSettings, settings.propagatorType );

            // Create multiple propagator settings
            std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsList;
            propagatorSettingsList.push_back( translationalPropagatorSettings );

            boost::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< MultiTypePropagatorSettings< double > >(
                        propagatorSettingsList, terminationSettings, dependentVariableSaveSettings );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings;
            if ( settings.integratorType == rungeKutta4 )
            {
                integratorSettings = boost::make_shared< IntegratorSettings< > > (
                            settings.integratorType, simulationStartEpoch,
                            settings.integratorFixedStepsize, settings.outputOneInEveryIntegrationSteps );
            }
            else
            {
                integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                            settings.integratorType, simulationStartEpoch, settings.integratorInitialStepsize,
                            settings.integratorSet, 1.0E-10, 1.0E+10, settings.integratorErrorTolerance,
                            settings.integratorErrorTolerance, settings.outputOneInEveryIntegrationSteps );
            }


            /////////////////////////////////////////////////////////////////////////////////
            //////////             PROPAGATE ORBIT            ///////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );

            std::map< double, Eigen::VectorXd > integrationResult =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
                    dynamicsSimulator.getDependentVariableHistory( );

            // Determine total computation time
            auto t_end = std::chrono::steady_clock::now();

            double computationTime =
                    std::chrono::duration_cast< std::chrono::milliseconds >( t_end - t_ini ).count() * 1.0e-3;

            // double lastEpoch = ( --integrationResult.end() )->first;


            /////////////////////////////////////////////////////////////////////////////////
            //////////              EXPORT RESULTS            ///////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////

            // Combine all the results in a map

            const int kepCols = settings.outputBodyKeplerianState * 6;
            const int carCols = settings.outputBodyCartesianState * 6;
            int varCols = 0;
            if ( dependentVariableOutput.size() > 0 )
            {
                varCols = ( ( dependentVariableOutput.begin() )->second ).rows();
            }

            std::map< double, Eigen::VectorXd > resultsMap;
            for( std::map< double, Eigen::VectorXd >::iterator iter = integrationResult.begin();
                 iter != integrationResult.end(); ++iter )
            {
                if ( ! settings.outputOnlyLastIntegrationStep || iter == --integrationResult.end() )
                {
                    double epoch = iter->first;

                    Eigen::VectorXd kepState( kepCols );
                    if ( settings.outputBodyKeplerianState )
                    {
                        Vector6d carState = iter->second;
                        kepState = convertCartesianToKeplerianElements( carState, earthGravitationalParameter );
                    }

                    Eigen::VectorXd carState( carCols );
                    if ( settings.outputBodyCartesianState )
                    {
                        carState = iter->second;
                    }

                    Eigen::VectorXd varState( varCols );
                    if ( dependentVariableOutput.size() > 0 )
                    {
                        Eigen::VectorXd dependentVariableRow = dependentVariableOutput[ epoch ];
                        varState = dependentVariableRow;
                    }

                    Eigen::VectorXd resultRow( kepCols + carCols + varCols );
                    resultRow << kepState, carState, varState;

                    resultsMap[ epoch ] = resultRow;
                }
            }

            // If propagation was resumed from output file, include also the old steps in the results map

            if ( settings.resuming )
            {
                for( std::map< double, Eigen::VectorXd >::iterator iter = settings.resultsMap.begin();
                     iter != settings.resultsMap.end(); ++iter )
                {
                    resultsMap[ iter->first ] = iter->second;
                }
            }

            settings.computationTime = computationTime;
            settings.resultsMap = resultsMap;
            settings.exportResults();
            settings.propagationSucceeded();
        }
    }
    return EXIT_SUCCESS;
}

