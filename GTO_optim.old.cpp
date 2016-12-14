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

#include "GTOUtilities/graceGravityModels.h"
#include "GTOUtilities/inputOutput.h"
#include "GTOUtilities/tabulatedDragCoefficients.h"


//! Execute propagation of orbit of Rocket around the Earth.
int main()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace interpolators;
    using namespace input_output;

    using namespace gto_utilities::input_output;
    using namespace gto_utilities::environment;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    // Vectors of independent variables
    for ( double delay_years = 0.0; delay_years <= 1.0; delay_years += 0.1 )
    {
        for ( double arg_perigee = 0.0; arg_perigee < 360.0; arg_perigee += 36.0 )
        {

            std::cout << std::endl << "Propagating orbit with delay " << delay_years << " years, argument of perigee " << arg_perigee << " deg" << std::endl;

            // Perturbations
            int n = 7;
            int m = 3;
            bool sunG = true;
            bool moonG = true;
            bool drag = true;
            bool srp = true;

            // Body
            double m_rocket = 400.0;
            double A_rocket = 10.0;
            double CD_drag = -1.0;
            double CD_srp = 2.2;

            // Initial state
            double h_p_km = 175.0;
            double h_a_km = 35780.0;
            double inc = 23.4;
            double raan = 0.0;
            double true_anom = 0.0;

            // Propagation / integration settings
            TranslationalPropagatorType propagator_type = cowell;
            AvailableIntegrators integrator_type = rungeKuttaVariableStepSize;
            RungeKuttaCoefficients::CoefficientSets rungekutta_type = RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78;
            double stepSize = 60.0;
            double maxPropPeriod = 1.0;

            //    std::cout << n << std::endl;
            //    std::cout << m << std::endl;
            //    std::cout << sunG << std::endl;
            //    std::cout << moonG << std::endl;
            //    std::cout << drag << std::endl;
            //    std::cout << srp << std::endl;
            //    std::cout << m_rocket << std::endl;
            //    std::cout << A_rocket << std::endl;
            //    std::cout << CD_drag << std::endl;
            //    std::cout << CD_srp << std::endl;
            //    std::cout << delay_years << std::endl;
            //    std::cout << h_p_km << std::endl;
            //    std::cout << h_a_km << std::endl;
            //    std::cout << inclination << std::endl;
            //    std::cout << raan << std::endl;
            //    std::cout << argumentOfPerigee << std::endl;
            //    std::cout << trueAnomaly << std::endl;
            //    std::cout << stepSize << std::endl;
            //    std::cout << maxPropPeriod << std::endl;


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
            std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

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

            /////////////////////////////////////////////////////////////////////////////////////////////////////
            /////             CREATE VEHICLE            /////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create spacecraft object.
            bodyMap[ "Rocket" ] = boost::make_shared< simulation_setup::Body >( );
            bodyMap[ "Rocket" ]->setConstantBodyMass( m_rocket );
            double referenceArea = A_rocket;

            if ( drag ) {
                boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;
                if ( CD_drag < 0 )
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
            }

            // SRP interface
            boost::shared_ptr< RadiationPressureInterfaceSettings > rocketRadiationPressureSettings = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", referenceArea, CD_srp, boost::assign::list_of( "Earth" ) );
            bodyMap[ "Rocket" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface( rocketRadiationPressureSettings, "Rocket", bodyMap ) );


            // Finalize body creation.
            setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

            ///////////////////////////////////////////////////////////////////////////////////////////////////
            ///            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////

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

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set Keplerian elements for Rocket.
            double meanEarthRadius = 6371.0E+3;
            double h_p = h_p_km * 1.0E3;
            double h_a = h_a_km * 1.0E3;
            double a = ( h_a + h_p ) / 2.0 + meanEarthRadius;
            double e = ( h_a - h_p ) / ( h_a + h_p + 2.0 * meanEarthRadius );

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
            // dependentVariablesToSave.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >( third_body_central_gravity, "Rocket", "Sun", 1 ) );
            // dependentVariablesToSave.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >( third_body_central_gravity, "Rocket", "Moon", 1 ) );

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
                integratorSettings = boost::make_shared< IntegratorSettings< > > ( rungeKutta4, simulationStartEpoch, stepSize );
            }
            else
            {
                integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >( rungeKuttaVariableStepSize, simulationStartEpoch, stepSize, rungekutta_type, stepSize/100.0, stepSize*100.0 );
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////
            /////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////

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


            ///////////////////////////////////////////////////////////////////////////////////////////////////
            ///        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////

            // Convert history to Keplerian
            std::map< double, Eigen::VectorXd > keplerianResult = integrationResult;
            for( std::map< double, Eigen::VectorXd >::iterator iter = integrationResult.begin(); iter != integrationResult.end(); ++iter )
            {
                double key = iter->first;
                Vector6d value = iter->second;
                keplerianResult[key] = convertCartesianToKeplerianElements( value, earthGravitationalParameter );
            }
            Eigen::VectorXd finalKeplerianState = ( --keplerianResult.end( ) )->second;
            double a_f = finalKeplerianState(0);
            double e_f = finalKeplerianState(1);
            double h_p_f = (a_f * ( 1 - e_f ) - meanEarthRadius)*1.0E-3;

            auto nlim = std::numeric_limits< double >::digits10;

            // Write propagation and computation times to file.
            std::string settingsCode = std::to_string( (int)round( delay_years * 100.0 ) ) + "_" + std::to_string( (int)( arg_perigee ) );
            std::string resultFilename = "propagated/optim/" + settingsCode + "_result.txt";
            std::map< int, double > resultMap;
            resultMap[0] = t_prop;
            resultMap[1] = t_comp;
            resultMap[2] = h_p_f;
            writeDataMapToTextFile( resultMap, resultFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

            // Write perturbed satellite propagation history to file.
            std::string historyFilename = "propagated/optim/" + settingsCode + "_history.txt";
            writeDataMapToTextFile( keplerianResult, historyFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

//            // Write additional variables to file.
//            std::string variablesFilename = "propagated/optim/" + settingsCode + "_variables.txt";
//            writeDataMapToTextFile( dependentVariableOutput, variablesFilename, GTO_OUTPUT_PATH, "", nlim, nlim, "," );

            std::cout << std::endl;

        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

