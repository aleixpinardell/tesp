/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TESP_SETTINGS_H
#define TESP_SETTINGS_H

#include <string>
#include <boost/filesystem.hpp>

#include "graceGravityModels.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tesp {

typedef GraceGravityModel                                                       GeopotentialModel;

typedef tudat::simulation_setup::AtmosphereTypes                                AtmosphericModel;

typedef tudat::propagators::TranslationalPropagatorType                         PropagatorType;

typedef tudat::numerical_integrators::AvailableIntegrators                      IntegratorType;
typedef tudat::numerical_integrators::RungeKuttaCoefficients::CoefficientSets   IntegratorSet;

enum SkippedPropagationOutput {
    none,
    reference,
    duplicate
};

enum EventResponse {
    allow,
    warn,
    error
};

typedef boost::filesystem::path path;

const double EARTH_AVERAGE_RADIUS = 6371010;  // m

class TespSettings
{
protected:

    double initialPerigeeAltitude =                         -1;         // km, <0 = not defined
    double initialApogeeAltitude =                          -1;         // km, <0 = not defined

    path inputDirectory;
    path inputFileName;
    path inputPath;

    path spaceWeatherFileRootDirectory;

    std::vector< path > includedPaths;
    path outputPath;

    path outputRootDirectory = "";
    path searchOutputRootDirectory = "";

    path duplicateSourcePath;

    std::string outputContentHeader( bool printIncludedFilePaths = true );
    std::string outputContentInputSettings();
    std::string outputContentOutputDescription();
    std::string propagationTerminationCause();
    std::string outputContentOutputResults();

    std::string resultsText;
    void parseResultsFromOutputFile();
    void canonizePaths();


public:
    TespSettings() { }

    TespSettings( const std::string inputFilePath, const bool outputFile = false );

    // const boost::filesystem::path& getInputPath() const { return inputPath; }
    // const boost::filesystem::path& getOutputPath() const { return outputPath; }

    bool checkRequestedOutputExists();
    void exportResults();

    void propagationSucceeded();


    // To be defined after the propagation

    std::map< double, Eigen::VectorXd > resultsMap;
    double propagationTime;
    double computationTime;
    tudat::propagators::PropagationTerminationReason propagationTerminationReason;


    // In case the initial state has to be read from an output file in which the final Cartesian state is provided

    double maximumFinalEpoch() {
        return std::min( endEpoch, initialEpoch + propagationPeriod );
    }
    bool resuming = false;
    double resumingEpoch = 0.0;
    Eigen::Matrix< double, 6, 1 > resumingState;
    bool resumingStateIsCartesian = false;


    // ACCELERATION MODEL

    GeopotentialModel geopotentialModel =                   GeopotentialModel::ggm02c;
    int geopotentialDegree =                                2;
    int geopotentialOrder =                                 0;

    bool thirdBodyAttractionSun =                           true;
    bool thirdBodyAttractionMoon =                          true;

    AtmosphericModel atmosphericModel =                     AtmosphericModel::nrlmsise00;
    std::string spaceWeatherFileRelativePath =              "";
    path spaceWeatherFilePath;
    bool atmosphericDrag =                                  true;

    bool solarRadiationPressure =                           true;
    bool ignoreEclipses =                                   false;


    // BODY PROPERTIES

    double bodyMass =                                       100.0;      // kg
    double bodyCrossSectionalArea =                         10.0;       // m^2

    bool altitudeDependentDragCoefficient =                 false;
    double bodyDragCoefficient =                            2.5;        // -
    double bodyRadiationPressureCoefficient =               1.5;        // -


    // INITIAL STATE

    double initialEpoch =                                   0.0;        // seconds since J2000  ( = 1-Jan-2000 @ 12 )

    double initialSemimajorAxis =                           7371.0;     // km
    double initialEccentricity =                            0.0;        // -
    double initialInclination =                             0.0;        // deg
    double initialArgumentPerigee =                         0.0;        // deg
    double initialLongitudeAscendingNode =                  0.0;        // deg
    double initialTrueAnomaly =                             0.0;        // deg


    // PROPAGATION SETTINGS

    bool preloadCelestialBodiesData =                       true;

    double endEpoch =                                       504835200;  // seconds since J2000  ( = 31-Dec-2015 @ 12 )
    double propagationPeriod = tudat::physical_constants::SIDEREAL_YEAR * 10.0;       // seconds
    double reentryAltitude =                                0.0;        // km

    PropagatorType propagatorType =                         PropagatorType::cowell;

    IntegratorType integratorType =                         IntegratorType::rungeKutta4;
    IntegratorSet integratorSet =                           IntegratorSet::rungeKuttaFehlberg78;
    double integratorFixedStepsize =                        60.0;       // s
    double integratorInitialStepsize =                      60.0;       // s
    double integratorErrorTolerance =                       1.0E-12;    // -

    // DSST specific settings
    double altitudeLimitEarthAtmosphericDrag                     = 600.0;  // km
    int numberOfQuadratureNodesEarthAtmosphericDrag              = 40;
    bool scalableNumberOfQuadratureNodesEarthAtmosphericDrag     = true;
    int numberOfQuadratureNodesSolarRadiationPressure            = 10;
    bool scalableNumberOfQuadratureNodesSolarRadiationPressure   = false;
    int SThirdBodyAttractionSun                                  = 2;
    int NThirdBodyAttractionSun                                  = 2;
    int SThirdBodyAttractionMoon                                 = 2;
    int NThirdBodyAttractionMoon                                 = 2;
    int SSolarRadiationPressure                                  = 2;
    int NSolarRadiationPressure                                  = 2;


    // OUTPUT SETTINGS - CURRENT RUN

    std::string outputDirectoryPath =                       "";
    std::string outputFileName =                            "$FILENAME";

    int outputNumericPrecision =                            8;

    bool outputInputSettings =                              true;
    bool outputResultsColumnsDescriptions =                 true;
    bool outputPropagationTime =                            true;
    bool outputComputationTime =                            true;
    bool outputPropagationTerminationCause =                true;

    bool outputBodyKeplerianState =                         true;
    bool outputBodyCartesianState =                         false;
    bool outputSunPosition =                                false;
    bool outputMoonPosition =                               false;

    bool outputDSSTMeanElementRates =                       false;
    bool outputDSSTMeanElementRatesZonalTerms =             false;
    bool outputDSSTMeanElementRatesSunGravity =             false;
    bool outputDSSTMeanElementRatesMoonGravity =            false;
    bool outputDSSTMeanElementRatesAtmosphericDrag =        false;
    bool outputDSSTMeanElementRatesSolarRadiationPressure = false;

    bool outputDSSTShortPeriodTerms =                       false;
    bool outputDSSTShortPeriodTermsZonalTerms =             false;
    bool outputDSSTShortPeriodTermsSunGravity =             false;
    bool outputDSSTShortPeriodTermsMoonGravity =            false;
    bool outputDSSTShortPeriodTermsAtmosphericDrag =        false;
    bool outputDSSTShortPeriodTermsSolarRadiationPressure = false;

    bool outputDSSTComputationTimes =                       false;
    bool outputDSSTComputationTimesZonalTerms =             false;
    bool outputDSSTComputationTimesSunGravity =             false;
    bool outputDSSTComputationTimesMoonGravity =            false;
    bool outputDSSTComputationTimesAtmosphericDrag =        false;
    bool outputDSSTComputationTimesSolarRadiationPressure = false;

    std::map< tudat::propagators::PropagationDependentVariables,
    std::vector< tudat::propagators::sst::force_models::ForceIdentifier > > dsstOutputSettings;

    int outputOneInEveryIntegrationSteps =                  1;
    bool outputOnlyLastIntegrationStep =                    false;

    bool printSuccessMessageAfterPropagation =              true;
    bool deleteInputFileAfterPropagation =                  false;


    // OUTPUT SETTINGS - PREVIOUS RUNS

    bool skipPropagationIfOutputFileExists =                false;
    bool skipPropagationIfResultsFound =                    true;
    bool allowResumingPropagationFromResults =              true;

    std::string searchOutputDirectoryPath =                 "";
    std::string searchOutputFilenamePattern =               R"(.+)";
    bool searchInSubdirectories =                           false;
    bool searchIgnoreFileContents =                         false;

    SkippedPropagationOutput skippedPropagationOutput =     SkippedPropagationOutput::duplicate;


    // WARNINGS AND ERRORS

    EventResponse nonDefinedVariable =                      allow;
    EventResponse unrecognizedVariableName =                warn;

    EventResponse multiplyDefinedVariableDifFiles =         allow;
    EventResponse multiplyDefinedVariableSameFile =         warn;

    EventResponse nonValidValue =                           error;

    EventResponse skippingPropagation =                     warn;
    EventResponse resumingPropagationFromOutput =           warn;

};


} // namespace tesp

#endif // TESP_SETTINGS_H
