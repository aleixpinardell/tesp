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

#include "Tudat/Astrodynamics/Propagators/DSST/forces/zonalSphericalHarmonicGravity.h"


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
        using namespace tudat::propagators::sst::force_models;

        using namespace tesp;


        ////////////////////////////////////////////////////////////////////////////
        ///////////      READ SETTINGS FROM INPUT FILE       ///////////////////////
        ////////////////////////////////////////////////////////////////////////////

        auto t_ini = std::chrono::steady_clock::now();

        std::string inputFilePath = argv[1];
        TespSettings settings( inputFilePath );

        if ( ! settings.checkRequestedOutputExists() )
        {
            ////////////////////////////////////////////////////////////////////////
            //////           CREATE ENVIRONMENT          ///////////////////////////
            ////////////////////////////////////////////////////////////////////////

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
                if ( settings.spaceWeatherFileRelativePath.length() > 0 )
                {
                    bodySettings[ "Earth" ]->atmosphereSettings =
                            boost::make_shared< NRLMSISE00AtmosphereSettings >(
                                settings.spaceWeatherFilePath.string() );
                }
                else
                {
                    bodySettings[ "Earth" ]->atmosphereSettings =
                            boost::make_shared< AtmosphereSettings >( nrlmsise00 );
                }
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
            std::vector< std::string > occultingBodies;
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

            if ( settings.geopotentialDegree >= 2 )  // spherical harmonics
            {
                accelerationsOfRocket[ "Earth" ].push_back(
                            boost::make_shared< SphericalHarmonicAccelerationSettings >(
                                settings.geopotentialDegree, settings.geopotentialOrder ) );
            }
            else  // only central gravity
            {
                accelerationsOfRocket[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
            }

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

            Eigen::Vector6d bodyInitialState;
            Eigen::Vector6d bodyInitialStateKep;

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

                // Transform to mean elements
                /* FIXME
                if ( settings.propagatorType == dsst )
                {
                    const double mu = 3.986004415e14;
                    const double J2 = 0.00108264;
                    const double R  = 6378136.30;
                    Eigen::Vector6d shortPeriodTermsDueToJ2 =
                            sst::force_models::transformToOsculatingUsingJ2Contribution(
                                bodyInitialStateKep, mu, J2, R ) - bodyInitialStateKep;

                    // std::cout << shortPeriodTermsDueToJ2.transpose() << std::endl;

                    bodyInitialStateKep -= shortPeriodTermsDueToJ2;
                }
                */
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
                /*if ( settings.atmosphericDrag )
                {*/
                // Altitude limit
                boost::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings =
                        boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                            boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                propagators::periapsis_altitude_dependent_variable, "Body", "Earth" ),
                            settings.reentryAltitude * 1.0E3, 1 );

                constituentSettings.push_back( altitudeTerminationSettings );
                /*}
                else
                {
                    std::cout << "Could not create a terminating condition based on an altitude limit because "
                                 "atmospheric drag is turned off." << std::endl;
                }*/
            }

            // Stop if ANY of the two is met
            boost::shared_ptr< PropagationTerminationSettings > terminationSettings = boost::make_shared<
                    propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

            // Determine dependent variables to save
            std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

            if ( settings.outputSunPosition )
            {
                dependentVariables.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                  relative_position_dependent_variable, "Sun", "Earth" ) );
            }

            if ( settings.outputMoonPosition )
            {
                dependentVariables.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                  relative_position_dependent_variable, "Moon", "Earth" ) );
            }

            // DSST output
            for ( auto ent: settings.dsstOutputSettings )
            {
                PropagationDependentVariables dependentVariable = ent.first;
                std::vector< ForceIdentifier > forceIDs = ent.second;
                for ( ForceIdentifier forceID: forceIDs )
                {
                    if ( forceID == ForceIdentifier() )
                    {
                        dependentVariables.push_back(
                                    boost::make_shared< SingleDependentVariableSaveSettings >(
                                        dependentVariable, "Body" ) );
                    }
                    else
                    {
                        dependentVariables.push_back(
                                    boost::make_shared< ForceIdentifiableDependentVariableSaveSettings >(
                                        dependentVariable, "Body", forceID.body, forceID.type ) );
                    }
                }
            }

            // Create dependent variable save settings
            boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariables, 0 );

            // Translational propagator settings
            boost::shared_ptr< DSSTTranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                    boost::make_shared< DSSTTranslationalStatePropagatorSettings< double > > (
                        centralBodies, accelerationModelMap, bodiesToPropagate, bodyInitialState,
                        terminationSettings, settings.propagatorType );

            // DSST specific settings
            if ( settings.propagatorType == dsst )
            {
                // Settings for thrid-body attraction Sun
                translationalPropagatorSettings->setForceModelSettings(
                            ForceIdentifier( "Sun", third_body_central_gravity ),
                            boost::make_shared< ConservativeSettings >(
                                settings.SThirdBodyAttractionSun, settings.NThirdBodyAttractionSun ) );

                // Settings for thrid-body attraction Sun
                translationalPropagatorSettings->setForceModelSettings(
                            ForceIdentifier( "Moon", third_body_central_gravity ),
                            boost::make_shared< ConservativeSettings >(
                                settings.SThirdBodyAttractionMoon, settings.NThirdBodyAttractionMoon ) );

                // Settings for atmospheric drag caused by Earth
                translationalPropagatorSettings->setForceModelSettings(
                            ForceIdentifier( "Earth", aerodynamic ),
                            boost::make_shared< AtmosphericDragSettings >(
                                settings.altitudeLimitEarthAtmosphericDrag * 1.0E+3,
                                settings.numberOfQuadratureNodesEarthAtmosphericDrag,
                                settings.scalableNumberOfQuadratureNodesEarthAtmosphericDrag ) );

                if ( settings.ignoreEclipses ) {
                    // Settings for solar radiation pressure when eclipses ARE ignored (conservative)
                    translationalPropagatorSettings->setForceModelSettings(
                                ForceIdentifier( "Sun", cannon_ball_radiation_pressure ),
                                boost::make_shared< ConservativeSettings >(
                                    settings.SSolarRadiationPressure, settings.NSolarRadiationPressure ) );
                } else {
                    // Settings for solar radiation pressure when eclipses are NOT ignored (non-conservative)
                    translationalPropagatorSettings->setForceModelSettings(
                                ForceIdentifier( "Sun", cannon_ball_radiation_pressure ),
                                boost::make_shared< RadiationPressureSettings >(
                                    settings.numberOfQuadratureNodesSolarRadiationPressure,
                                    settings.scalableNumberOfQuadratureNodesSolarRadiationPressure ) );
                }
            }

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
                            settings.integratorErrorTolerance, settings.outputOneInEveryIntegrationSteps,
                            0.8, 100.0, 0.01 );
            }


            /////////////////////////////////////////////////////////////////////////////////
            //////////             PROPAGATE ORBIT            ///////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////

            auto t_prop = std::chrono::steady_clock::now();

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );

            std::map< double, Eigen::VectorXd > integrationResult =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            // Determine propagation time
            auto t_end = std::chrono::steady_clock::now();
            double propagationTime =
                    std::chrono::duration_cast< std::chrono::milliseconds >( t_end - t_prop ).count() * 1.0e-3;
            double computationTime =
                    std::chrono::duration_cast< std::chrono::milliseconds >( t_end - t_ini ).count() * 1.0e-3;

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
                    dynamicsSimulator.getDependentVariableHistory( );

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
                        Eigen::Vector6d carState = iter->second;
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

            settings.resultsMap = resultsMap;
            settings.propagationTime = propagationTime;
            settings.computationTime = computationTime;
            settings.propagationTerminationReason = dynamicsSimulator.getPropagationTerminationReason();
            settings.exportResults();
            settings.propagationSucceeded();
        }
    }
    return EXIT_SUCCESS;
}

