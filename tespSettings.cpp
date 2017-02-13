/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <fstream>
#include <iostream>
#include <iterator>
// #include <regex>
#include <map>
#include <boost/filesystem.hpp>

#include "tespSettings.h"

namespace tesp {

std::vector< int > occurencesOfElement( std::string element, std::vector< std::string > vector )
{
    std::vector< int > indexes;

    for ( unsigned int i = 0; i < vector.size(); i++ )
    {
        if ( vector[ i ] == element )
        {
            indexes.push_back( i );
        }
    }

    return indexes;
}

void printWarningOrError( std::string variableName, std::string message,
                          EventResponse response, std::string continueMessage = "Using default value" )
{
    if ( response == error )
    {
        std::cerr << variableName << " " << message << "." << std::endl;
        std::terminate();
    }
    else
    {
        std::cout << variableName << " " << message << ". " << continueMessage << "." << std::endl;
    }
}

void checkMultipleDefinitions( std::string name, std::vector< std::string > names, std::vector< int > files,
                               EventResponse difFiles, EventResponse sameFile, int index )
{
    if ( difFiles != allow || sameFile != allow )
    {
        std::vector< int > occ = occurencesOfElement( name, names );
        if ( occ.size() > 1 )
        {
            if ( occ.back() == index )
            {
                bool multiDefInSameFile = false;
                bool multiDefInDifFiles = false;
                for ( unsigned int o = 0; o < occ.size(); o++ )
                {
                    int f = files[ occ[ o ] ];
                    for ( unsigned int oo = 0; oo < occ.size(); oo++ )
                    {
                        if ( oo != o )
                        {
                            int ff = files[ occ[ oo ] ];
                            if ( ff == f )
                            {
                                multiDefInSameFile = true;
                            }
                            if ( ff != f )
                            {
                                multiDefInDifFiles = true;
                            }
                        }
                    }
                }
                if ( difFiles != allow && multiDefInDifFiles )
                {
                    printWarningOrError( name, "multiply defined in different input files",
                                         difFiles, "Using last definition" );
                }
                if ( sameFile != allow && multiDefInSameFile )
                {
                    printWarningOrError( name, "multiply defined in the same input file",
                                         sameFile, "Using last definition" );
                }
            }
        }
    }
}

void assignInsideRange( double &variable, std::string name, std::string value, double min, double max,
                        bool minIncluded, bool maxIncluded, EventResponse response,
                        bool ignoreLowerLimit = false, bool ignoreUpperLimit = false )
{
    double val = std::stod( value );
    bool aboveLower = ignoreLowerLimit || ( minIncluded ? val >= min : val > min );
    bool belowUpper = ignoreUpperLimit || ( maxIncluded ? val <= max : val < max );
    if ( aboveLower && belowUpper )
    {
        variable = val;
    }
    else if ( response != allow )
    {
        printWarningOrError( name, "value not valid", response );
    }
}

void assignInsideRange( int &variable, std::string name, std::string value, int min, int max,
                        bool minIncluded, bool maxIncluded, EventResponse response,
                        bool ignoreLowerLimit = false, bool ignoreUpperLimit = false )
{
    int val = std::stoi( value );
    bool aboveLower = ignoreLowerLimit || ( minIncluded ? val >= min : val > min );
    bool belowUpper = ignoreUpperLimit || ( maxIncluded ? val <= max : val < max );
    if ( aboveLower && belowUpper )
    {
        variable = val;
    }
    else if ( response != allow )
    {
        printWarningOrError( name, "value not valid", response );
    }
}

void assignBoolean( bool &variable, std::string name, std::string value, EventResponse response )
{
    if ( value == "NO" )
    {
        variable = false;
    }
    else if ( value == "YES" )
    {
        variable = true;
    }
    else if ( response != allow )
    {
        printWarningOrError( name, "value not valid", response );
    }
}

void assignDate( double &variable, std::string name, std::string value, EventResponse response )
{
    // Check if the user provided a formatted date
    double formattedDate = false;
    boost::regex dateFormat(
        R"((\d{4})-(\d{1,2})-(\d{1,2})(?: (\d{2}))?(?::(\d{2}))?(?::(\d{2}(?:.\d+)?))?)" );
    boost::smatch dateMatch;
    std::string tempValue = value;
    std::vector< double > dateComponents = { 2000, 1, 1, 12, 0, 0 };
    while( boost::regex_search( tempValue, dateMatch, dateFormat, boost::match_not_dot_newline ) )
    {
        for ( unsigned int d = 1; d < dateMatch.size(); d++ )
        {
            std::string match = dateMatch[ d ].str();
            if ( match.length() > 0 )
            {
                formattedDate = true;
                dateComponents[ d - 1 ] = std::stod( match );
            }
        }
        tempValue = dateMatch.suffix();
    }

    // If user provided a formatted date, transform it to seconds since J2000
    // https://en.wikipedia.org/wiki/Julian_day#Converting_Julian_or_Gregorian_calendar_date_to_Julian_day_number
    if ( formattedDate )
    {
        double year =   dateComponents[ 0 ];
        double month =  dateComponents[ 1 ];
        double day =    dateComponents[ 2 ];
        double hour =   dateComponents[ 3 ];
        double minute = dateComponents[ 4 ];
        double second = dateComponents[ 5 ];

        // Check valid date (non-exhaustive)
        if ( month <= 12 || day <= 31 || hour <= 23 || minute <= 59 || second < 60 )
        {
            double a = floor( ( 14.0 - month )/12.0 );
            double y = year + 4800.0 - a;
            double m = month + 12.0 * a - 3.0;

            // Julina day number
            double JDN = day + floor( ( 153.0 * m + 2.0 ) / 5.0 ) + 365.0 * y + floor( y / 4.0 )
                             - floor( y / 100.0 ) + floor( y / 400.0 ) - 32045.0;

            // Julian day
            double JD = JDN + ( hour - 12.0 ) / 24.0 + minute / 1440.0 + second / 86400.0;

            // Substract Julian day for J2000 --> then transform from julian days to seconds
            variable = ( JD - 2451545.0 ) * 86400.0;
        }
        else if ( response != allow )
        {
            printWarningOrError( name, "value not valid", response );
        }
    }
    else // If user did not provide a formatted date, then understand it directly as seconds since J2000
    {
        variable = std::stod( value );
    }
}



void TespSettings::parseResultsFromOutputFile()
{
    // Process results
    boost::regex newLineEx( R"((.+)\n)" );
    boost::smatch newLineMatch;
    std::vector< std::string > lines;
    std::string newLineSuffix = resultsText;
    while( boost::regex_search( newLineSuffix, newLineMatch, newLineEx, boost::match_not_dot_newline ) )
    {
        lines.push_back( newLineMatch[ 1 ] );
        newLineSuffix = newLineMatch.suffix();
    }

    for ( std::string line : lines )
    {
        boost::regex valuesEx( R"((\S+))" );
        boost::smatch valuesMatch;
        std::string valuesSuffix = line;
        std::vector< double > values;
        while( boost::regex_search( valuesSuffix, valuesMatch, valuesEx, boost::match_not_dot_newline ) )
        {
            values.push_back( std::stod( valuesMatch[ 1 ] ) );
            valuesSuffix = valuesMatch.suffix();
        }

        double epoch = values[ 0 ];
        Eigen::VectorXd vector( values.size() - 1 );
        for ( unsigned int v = 1; v < values.size(); v++ )
        {
            vector[ v - 1 ] = values[ v ];
        }

        resultsMap[ epoch ] = vector;
    }
}


TespSettings::TespSettings( const std::string inputFilePath, const bool outputFile )
{
    // Add .tespin extension if not present
    std::string inputFilePathWithExtension = inputFilePath;
    if ( ! outputFile && inputFilePathWithExtension.find(".tespin") == std::string::npos )
    {
        inputFilePathWithExtension += ".tespin";
    }

    // Determine the directory and filename (without file extension) of the input file
    inputDirectory = boost::filesystem::canonical( inputFilePathWithExtension ).parent_path();
    inputFileName = path( inputFilePathWithExtension ).stem();

    inputPath = inputDirectory / path( inputFilePathWithExtension ).filename();


    // std::string inputDirectory = getDirectoryForPath( inputFilePath );
    // std::string inputFileName = getFilenameForPath( inputFilePath, false );


    // PARSE

    std::vector< std::string > includedFilesPaths = { inputPath.string() };
    std::vector< std::string > includedFilesContents;

    // Look for included files recursively
    std::vector< bool > includeStatementsProcessed = { false };

    bool allFilesHaveBeenIncluded = false;
    while ( ! allFilesHaveBeenIncluded )
    {
        // Determine the index of the first of includedFiles whose include statements have not been processed yet
        int currentIndex = -1;
        for ( unsigned int i = 0; i < includeStatementsProcessed.size(); i++ )
        {
            if ( ! includeStatementsProcessed[ i ] )
            {
                currentIndex = i;
                break;
            }
        }

        // If currentIndex has not been updated, then all files have already been included
        if ( currentIndex == -1 )
        {
            allFilesHaveBeenIncluded = true;
        }
        else
        {
            // Flag current file's #INCLUDE statements as processed
            includeStatementsProcessed[ currentIndex ] = true;

            std::string currentFile = includedFilesPaths[ currentIndex ];
            // Add .tespin extension if not present
            if ( ! outputFile && currentFile.find(".tespin") == std::string::npos )
            {
                currentFile += ".tespin";
            }

            // Replace current file path by its canonical path
            path currentFilePath = boost::filesystem::canonical( currentFile, inputDirectory );
            includedFilesPaths[ currentIndex ] = currentFilePath.string();

            // Open current file
            if ( ! boost::filesystem::exists( currentFilePath ) )
            {
                std::cerr << "The" << ( includedFilesPaths.size() > 1 ? " included " : " " )
                          << "input file \"" << currentFilePath << "\" does not exist." << std::endl;
                std::terminate();
            }
            std::ifstream ifs( currentFilePath.string() );
            std::string contents( ( std::istreambuf_iterator<char>(ifs) ),
                                  ( std::istreambuf_iterator<char>() ) );

            // Remove comments
            boost::regex commentEx( R"(\%.*)" );
            contents = boost::regex_replace( contents, commentEx, std::string( "" ), boost::match_not_dot_newline );

            // Look for #INCLUDE statements in current file contents
            boost::regex includeEx( R"(\#INCLUDE\s+(["'])(.*)\1)" );
            boost::smatch includeMatch;
            std::string tempContents = contents;
            std::vector< std::string > currentIncludedPaths;
            std::vector< bool > pathsProcessed;

            // If this is not an output file, then remove its comments and store the found included files
            if ( ! outputFile )
            {
                while ( boost::regex_search( tempContents, includeMatch, includeEx, boost::match_not_dot_newline ) )
                {
                    boost::ssub_match submatch = includeMatch[ 2 ];
                    currentIncludedPaths.push_back( submatch.str() );
                    pathsProcessed.push_back( false );
                    tempContents = includeMatch.suffix();
                }

                // Remove #INCLUDE statements
                contents = boost::regex_replace( contents, includeEx, std::string( "" ),
                                                 boost::match_not_dot_newline );

                // Add #INCLUDE'd files to vector of includedFilesPaths and mark them as not processed
                includedFilesPaths.insert( includedFilesPaths.begin() + currentIndex,
                                           currentIncludedPaths.begin(), currentIncludedPaths.end() );

                includeStatementsProcessed.insert( includeStatementsProcessed.begin() + currentIndex,
                                                   pathsProcessed.begin(), pathsProcessed.end() );
            }

            // Store current file contents for later processing
            includedFilesContents.insert( includedFilesContents.begin() + currentIndex, contents );
        }
    }

    // Vector of defined variable names and values (to deal with multiple definitions later)
    std::vector< std::string > variableNames;
    std::vector< std::string > variableValues;
    std::vector< int > variableDefinitionFiles;

    // Loop through all the included files
    for ( unsigned int f = 0; f < includedFilesContents.size(); f++ )
    {
        std::string input = includedFilesContents[ f ];

        // Add included file to member variable so that it can be accessed later
        if ( f < includedFilesContents.size() - 1 )
        {
            includedPaths.push_back( includedFilesPaths[ f ] );
        }


        // If the file is an output file, we need to treat it differently
        if ( outputFile )
        {
            // Read computation time
            boost::regex computationTimeEx( R"(COMPUTATION TIME = (.*?)\s+)" );
            boost::smatch computationTimeMatch;
            std::string computationTimeSuffix = input;
            while( boost::regex_search( computationTimeSuffix, computationTimeMatch, computationTimeEx,
                                        boost::match_not_dot_newline) )
            {
                computationTime = std::stod( computationTimeMatch[ 1 ] );
                computationTimeSuffix = computationTimeMatch.suffix();
                break;
            }

            // Read results (store them in memeber variable for later processing, if needed)
            boost::regex resultsEx( R"(RESULTS = \n([\s\S]*))" );
            boost::smatch resultsMatch;
            std::string resultsSuffix = input;
            while ( boost::regex_search( resultsSuffix, resultsMatch, resultsEx, boost::match_not_dot_newline ) )
            {
                resultsText = resultsMatch[ 1 ];
                resultsSuffix = resultsMatch.suffix();
                input = resultsMatch.prefix();
                break;
            }

            // Read input settings, assign it to input variable, and ignore all the rest
            boost::regex inputSettingsEx( R"(\% INPUT SETTINGS([\s\S]*)\% OUTPUT DESCRIPTION)" );
            boost::smatch inputSettingsMatch;
            std::string tempinputSettingsText = input;
            while ( boost::regex_search( tempinputSettingsText, inputSettingsMatch, inputSettingsEx,
                                       boost::match_not_dot_newline ) )
            {
                std::cout << "a match!" << std::endl;
                input = inputSettingsMatch[ 1 ];
                tempinputSettingsText = inputSettingsMatch.suffix();
                break;
            }

            // Remove comments
            boost::regex commentEx( R"(\%.*)" );
            input = boost::regex_replace( input, commentEx, std::string( "" ), boost::match_not_dot_newline );
        }


        // Look for variables "between quotes" or 'single quotes'
        boost::regex variableQuotesEx( R"((\S+)\s+(["'])(.*)\2.*)" );
        boost::smatch variableQuotesMatch;
        std::string tempQuotesInput = input;
        while( boost::regex_search( tempQuotesInput, variableQuotesMatch, variableQuotesEx,
                                    boost::match_not_dot_newline ) )
        {
            boost::ssub_match submatchName = variableQuotesMatch[ 1 ];
            boost::ssub_match submatchValue = variableQuotesMatch[ 3 ];
            variableNames.push_back( submatchName.str() );
            variableValues.push_back( submatchValue.str() );
            variableDefinitionFiles.push_back( f );
            tempQuotesInput = variableQuotesMatch.suffix();
        }

        // Remove lines matching variables "between quotes"
        input = boost::regex_replace( input, variableQuotesEx, std::string( "" ), boost::match_not_dot_newline );

        // Look for the rest of variables
        boost::regex variableEx( R"((\S+)\s+(\S+).*)" );
        boost::smatch variableMatch;
        std::string tempVarInput = input;
        while( boost::regex_search( tempVarInput, variableMatch, variableEx, boost::match_not_dot_newline ) )
        {
            boost::ssub_match submatchName = variableMatch[ 1 ];
            boost::ssub_match submatchValue = variableMatch[ 2 ];
            variableNames.push_back( submatchName.str() );
            variableValues.push_back( submatchValue.str() );
            variableDefinitionFiles.push_back( f );
            tempVarInput = variableMatch.suffix();
        }

    }


    // Determine warnings and errors

    std::vector< std::string > warningErrorsVariablesNames = { "NON_DEFINED_VARIABLE",
                                                               "UNRECOGNIZED_VARIABLE_NAME",
                                                               "MULTIPLY_DEFINED_VARIABLE_DIF_FILES",
                                                               "MULTIPLY_DEFINED_VARIABLE_SAME_FILE",
                                                               "NON_VALID_VALUE",
                                                               "SKIPPING_PROPAGATION",
                                                               "RESUMING_PROPAGATION_FROM_OUTPUT" };

    if ( outputFile )
    {
        nonDefinedVariable = allow;
        unrecognizedVariableName = allow;
        multiplyDefinedVariableDifFiles = allow;
        multiplyDefinedVariableSameFile = allow;
        nonValidValue = allow;
        skippingPropagation = allow;
    }
    else
    {
        for ( unsigned int v = 0; v < warningErrorsVariablesNames.size(); v++ )
        {
            std::string name = warningErrorsVariablesNames[ v ];
            std::vector< int > occ = occurencesOfElement( name, variableNames );
            if ( occ.size() > 0 )
            {
                int pos = occ.back();
                std::string value = variableValues[ pos ];

                if ( name == "NON_DEFINED_VARIABLE" )
                {
                    if ( value == "USE_DEFAULT" ) { nonDefinedVariable = allow;                 }
                    else if ( value == "WARN" ) {   nonDefinedVariable = warn;                  }
                    else if ( value == "ERROR" ) {  nonDefinedVariable = error;                 }
                }
                else if ( name == "UNRECOGNIZED_VARIABLE_NAME" )
                {
                    if ( value == "IGNORE" ) {      unrecognizedVariableName = allow;           }
                    else if ( value == "WARN" ) {   unrecognizedVariableName = warn;            }
                    else if ( value == "ERROR" ) {  unrecognizedVariableName = error;           }
                }
                else if ( name == "MULTIPLY_DEFINED_VARIABLE_DIF_FILES" )
                {
                    if ( value == "USE_LAST" ) {    multiplyDefinedVariableDifFiles = allow;    }
                    else if ( value == "WARN" ) {   multiplyDefinedVariableDifFiles = warn;     }
                    else if ( value == "ERROR" ) {  multiplyDefinedVariableDifFiles = error;    }
                }
                else if ( name == "MULTIPLY_DEFINED_VARIABLE_SAME_FILE" )
                {
                    if ( value == "USE_LAST" ) {    multiplyDefinedVariableSameFile = allow;    }
                    else if ( value == "WARN" ) {   multiplyDefinedVariableSameFile = warn;     }
                    else if ( value == "ERROR" ) {  multiplyDefinedVariableSameFile = error;    }
                }
                else if ( name == "NON_VALID_VALUE" )
                {
                    if ( value == "USE_DEFAULT" ) { nonValidValue = allow;                      }
                    else if ( value == "WARN" ) {   nonValidValue = warn;                       }
                    else if ( value == "ERROR" ) {  nonValidValue = error;                      }
                }
                else if ( name == "SKIPPING_PROPAGATION" )
                {
                    if ( value == "SILENT" ) {      skippingPropagation = allow;                }
                    else if ( value == "WARN" ) {   skippingPropagation = warn;                 }
                }
                else if ( name == "RESUMING_PROPAGATION_FROM_OUTPUT" )
                {
                    if ( value == "SILENT" ) {      resumingPropagationFromOutput = allow;      }
                    else if ( value == "WARN" ) {   resumingPropagationFromOutput = warn;       }
                }
            }
        }
    }

    // Handle initial perigee/apogee altitude

    double Re = 6371.0; // km
    double sma = initialSemimajorAxis;
    double ecc = initialEccentricity;
    double hp = sma * ( 1 - ecc ) - Re;
    double ha = sma * ( 1 + ecc ) - Re;

    // Read initial perigee/apogee altitudes
    std::vector< int > occp = occurencesOfElement( "INITIAL_PERIGEE_ALTITUDE", variableNames );
    if ( occp.size() > 0 )
    {
        checkMultipleDefinitions( "INITIAL_PERIGEE_ALTITUDE", variableNames, variableDefinitionFiles,
                                  multiplyDefinedVariableDifFiles, multiplyDefinedVariableSameFile, occp.back() );
    }

    std::vector< int > occa = occurencesOfElement( "INITIAL_APOGEE_ALTITUDE", variableNames );
    if ( occa.size() > 0 )
    {
        checkMultipleDefinitions( "INITIAL_APOGEE_ALTITUDE", variableNames, variableDefinitionFiles,
                                  multiplyDefinedVariableDifFiles, multiplyDefinedVariableSameFile, occa.back() );
    }

    if ( occp.size() == occa.size() )
    {
        for ( unsigned int o = 0; o < occp.size(); o++ )
        {
            int op = occp[ o ];
            int oa = occa[ o ];

            assignInsideRange( initialPerigeeAltitude, "INITIAL_PERIGEE_ALTITUDE", variableValues[ op ],
                               0.0, 0.0, false, false, nonValidValue, false, true );

            assignInsideRange( initialApogeeAltitude, "INITIAL_APOGEE_ALTITUDE", variableValues[ oa ],
                               0.0, 0.0, false, false, nonValidValue, false, true );

            hp = initialPerigeeAltitude;
            ha = initialApogeeAltitude;

            // Replace semi-major axis and eccentricity
            sma = ( ha + hp + 2.0 * Re ) / 2.0;
            variableNames[ op ] = "INITIAL_SEMIMAJOR_AXIS";
            std::stringstream smaStream;
            smaStream << std::fixed << std::setprecision(8) << sma;
            variableValues[ op ] = smaStream.str();

            ecc = ( ha - hp ) / ( ha + hp + 2.0 * Re );
            variableNames[ oa ] = "INITIAL_ECCENTRICITY";
            std::stringstream eccStream;
            eccStream << std::fixed << std::setprecision(8) << ecc;
            variableValues[ oa ] = eccStream.str();
        }
    }
    else
    {
        std::cerr << "Could not determine initial state:\n"
                  << "INITIAL_PERIGEE_ALTITUDE and INITIAL_APOGEE_ALTITUDE defined different number of times.\n\n"
                  << std::endl;
        std::terminate();
    }

    // Determine the rest of variables
    std::vector< std::string > recognizedVariablesNames = { "GEOPOTENTIAL_MODEL",
                                                            "GEOPOTENTIAL_DEGREE",
                                                            "GEOPOTENTIAL_ORDER",
                                                            "THIRD_BODY_ATTRACTION_SUN",
                                                            "THIRD_BODY_ATTRACTION_MOON",
                                                            "ATMOSPHERIC_MODEL",
                                                            "ATMOSPHERIC_DRAG",
                                                            "SOLAR_RADIATION_PRESSURE",
                                                            "IGNORE_ECLIPSES",
                                                            "BODY_MASS",
                                                            "BODY_CROSS_SECTIONAL_AREA",
                                                            "ALTITUDE_DEPENDENT_DRAG_COEFFICIENT",
                                                            "BODY_DRAG_COEFFICIENT",
                                                            "BODY_RADIATION_PRESSURE_COEFFICIENT",
                                                            "INITIAL_EPOCH",
                                                            "INITIAL_SEMIMAJOR_AXIS",
                                                            "INITIAL_ECCENTRICITY",
                                                            "INITIAL_INCLINATION",
                                                            "INITIAL_ARGUMENT_PERIGEE",
                                                            "INITIAL_LONGITUDE_ASCENDING_NODE",
                                                            "INITIAL_TRUE_ANOMALY",
                                                            "PRELOAD_CELESTIAL_BODIES_DATA",
                                                            "END_EPOCH",
                                                            "PROPAGATION_PERIOD",
                                                            "REENTRY_ALTITUDE",
                                                            "PROPAGATOR_NAME",
                                                            "INTEGRATOR_NAME",
                                                            "INTEGRATOR_FIXED_STEPSIZE",
                                                            "INTEGRATOR_INITIAL_STEPSIZE",
                                                            "INTEGRATOR_ERROR_TOLERANCE",
                                                            "OUTPUT_NUMERIC_PRECISION",
                                                            "OUTPUT_INPUT_SETTINGS",
                                                            "OUTPUT_RESULTS_COLUMNS_DESCRIPTIONS",
                                                            "OUTPUT_COMPUTATION_TIME",
                                                            "OUTPUT_BODY_KEPLERIAN_STATE",
                                                            "OUTPUT_BODY_CARTESIAN_STATE",
                                                            "OUTPUT_SUN_POSITION",
                                                            "OUTPUT_MOON_POSITION",
                                                            "OUTPUT_ONE_IN_EVERY_INTEGRATION_STEPS",
                                                            "OUTPUT_ONLY_LAST_INTEGRATION_STEP",
                                                            "OUTPUT_DIRECTORY_PATH",
                                                            "OUTPUT_FILE_NAME",
                                                            "PRINT_SUCCESS_MESSAGE_AFTER_PROPAGATION",
                                                            "DELETE_INPUT_FILE_AFTER_PROPAGATION",
                                                            "SEARCH_OUTPUT_DIRECTORY_PATH",
                                                            "SEARCH_OUTPUT_FILENAME_PATTERN",
                                                            "SEARCH_IN_SUBDIRECTORIES",
                                                            "SKIP_PROPAGATION_IF_OUTPUT_FILE_EXISTS",
                                                            "SKIP_PROPAGATION_IF_RESULTS_FOUND",
                                                            "ALLOW_RESUMING_PROPAGATION_FROM_RESULTS",
                                                            "OUTPUT_FOR_SKIPPED_PROPAGATION",
                                                            "SEARCH_IGNORE_FILE_CONTENTS" };

    // Check if all the recognized variables have been provided by the user
    if ( nonDefinedVariable != allow )
    {
        for ( unsigned int i = 0; i < recognizedVariablesNames.size(); i++ )
        {
            std::string recognizedName = recognizedVariablesNames[ i ];
            if ( occurencesOfElement( recognizedName, variableNames ).size() == 0 )
            {
                printWarningOrError( recognizedName, "not defined", nonDefinedVariable );
            }
        }
    }

    // Check all variables provided by user
    for ( unsigned int i = 0; i < variableNames.size(); i++ )
    {
        std::string name = variableNames[ i ];

        // Check if recognized variable or warn/error variable
        std::vector< int > occre = occurencesOfElement( name, recognizedVariablesNames );
        std::vector< int > occwe = occurencesOfElement( name, warningErrorsVariablesNames );
        if ( unrecognizedVariableName != allow )
        {
            if ( occre.size() + occwe.size() == 0 )
            {
                printWarningOrError( name, "is not a recognized variable name",
                                     unrecognizedVariableName, "Variable ignored" );
            }
        }

        // Only consider recognized variables (warn/error variables were parsed earlier)
        if ( occre.size() > 0 )
        {
            checkMultipleDefinitions( name, variableNames, variableDefinitionFiles,
                                      multiplyDefinedVariableDifFiles, multiplyDefinedVariableSameFile, i );

            int pos = occurencesOfElement( name, variableNames ).back();
            std::string value = variableValues[ pos ];
            path definitionFilePath( includedFilesPaths[ variableDefinitionFiles[ pos ] ] );
            EventResponse res = nonValidValue;

            if ( name == "GEOPOTENTIAL_MODEL" )
            {
                if ( value == "GGM02C" )
                {
                    geopotentialModel = GraceGravityModel::ggm02c;
                }
                else if ( value == "GGM02S" )
                {
                    geopotentialModel = GraceGravityModel::ggm02s;
                }
                else
                {
                    printWarningOrError( name, "value is not valid", nonValidValue );
                }
            }
            else if ( name == "GEOPOTENTIAL_DEGREE" )
            {
                assignInsideRange( geopotentialDegree, name, value, 1, 160, true, true, res );
            }
            else if ( name == "GEOPOTENTIAL_ORDER" )
            {
                assignInsideRange( geopotentialOrder, name, value, 0, 160, true, true, res );
            }
            else if ( name == "THIRD_BODY_ATTRACTION_SUN" )
            {
                assignBoolean( thirdBodyAttractionSun, name, value, res );
            }
            else if ( name == "THIRD_BODY_ATTRACTION_MOON" )
            {
                assignBoolean( thirdBodyAttractionMoon, name, value, res );
            }
            else if ( name == "ATMOSPHERIC_MODEL" )
            {
                if ( value == "EXPONENTIAL" )
                {
                    atmosphericModel = AtmosphericModel::exponential_atmosphere;
                }
                else if ( value == "NRLMSISE00" )
                {
                    atmosphericModel = AtmosphericModel::nrlmsise00;
                }
                else
                {
                    printWarningOrError( name, "value is not valid", nonValidValue );
                }
            }
            else if ( name == "ATMOSPHERIC_DRAG" )
            {
                assignBoolean( atmosphericDrag, name, value, res );
            }
            else if ( name == "SOLAR_RADIATION_PRESSURE" )
            {
                assignBoolean( solarRadiationPressure, name, value, res );
            }
            else if ( name == "IGNORE_ECLIPSES" )
            {
                assignBoolean( ignoreEclipses, name, value, res );
            }
            else if ( name == "BODY_MASS" )
            {
                assignInsideRange( bodyMass, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "BODY_CROSS_SECTIONAL_AREA" )
            {
                assignInsideRange( bodyCrossSectionalArea, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "ALTITUDE_DEPENDENT_DRAG_COEFFICIENT" )
            {
                assignBoolean( altitudeDependentDragCoefficient, name, value, res );
            }
            else if ( name == "BODY_DRAG_COEFFICIENT" )
            {
                assignInsideRange( bodyDragCoefficient, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "BODY_RADIATION_PRESSURE_COEFFICIENT" )
            {
                assignInsideRange( bodyRadiationPressureCoefficient, name, value, 1.0, 2.0, true, true, res );
            }
            else if ( name == "INITIAL_EPOCH" )
            {
                assignDate( initialEpoch, name, value, res );
            }
            else if ( name == "INITIAL_SEMIMAJOR_AXIS" )
            {
                assignInsideRange( initialSemimajorAxis, name, value, 0.0, 0.0, false, false, res, false, true );
                // initialPerigeeAltitude = -1;
            }
            else if ( name == "INITIAL_ECCENTRICITY" )
            {
                assignInsideRange( initialEccentricity, name, value, 0.0, 1.0, true, false, res );
                // initialApogeeAltitude = -1;
            }
            else if ( name == "INITIAL_INCLINATION" )
            {
                assignInsideRange( initialInclination, name, value, 0.0, 180.0, true, true, res );
            }
            else if ( name == "INITIAL_ARGUMENT_PERIGEE" )
            {
                assignInsideRange( initialArgumentPerigee, name, value, 0.0, 360.0, true, true, res );
            }
            else if ( name == "INITIAL_LONGITUDE_ASCENDING_NODE" )
            {
                assignInsideRange( initialLongitudeAscendingNode, name, value, 0.0, 360.0, true, true, res );
            }
            else if ( name == "INITIAL_TRUE_ANOMALY" )
            {
                assignInsideRange( initialTrueAnomaly, name, value, 0.0, 360.0, true, true, res );
            }
            else if ( name == "PRELOAD_CELESTIAL_BODIES_DATA" )
            {
                assignBoolean( preloadCelestialBodiesData, name, value, res );
            }
            else if ( name == "END_EPOCH" )
            {
                assignDate( endEpoch, name, value, res );
            }
            else if ( name == "PROPAGATION_PERIOD" )
            {
                assignInsideRange( propagationPeriod, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "REENTRY_ALTITUDE" )
            {
                assignInsideRange( reentryAltitude, name, value, 0.0, 0.0, true, false, res, false, true );
            }
            else if ( name == "PROPAGATOR_NAME" )
            {
                if ( value == "COWELL" )
                {
                    propagatorType = PropagatorType::cowell;
                }
                else if ( value == "ENCKE" )
                {
                    propagatorType = PropagatorType::encke;
                }
                else
                {
                    printWarningOrError( name, "value is not valid", nonValidValue );
                }
            }
            else if ( name == "INTEGRATOR_NAME" )
            {
                if ( value == "RK4" )
                {
                    integratorType = IntegratorType::rungeKutta4;
                }
                else if ( value == "RK45" )
                {
                    integratorType = IntegratorType::rungeKuttaVariableStepSize;
                    integratorSet = IntegratorSet::rungeKuttaFehlberg45;
                }
                else if ( value == "RK56" )
                {
                    integratorType = IntegratorType::rungeKuttaVariableStepSize;
                    integratorSet = IntegratorSet::rungeKuttaFehlberg56;
                }
                else if ( value == "RK78" )
                {
                    integratorType = IntegratorType::rungeKuttaVariableStepSize;
                    integratorSet = IntegratorSet::rungeKuttaFehlberg78;
                }
                else if ( value == "DP78" )
                {
                    integratorType = IntegratorType::rungeKuttaVariableStepSize;
                    integratorSet = IntegratorSet::rungeKutta87DormandPrince;
                }
                else
                {
                    printWarningOrError( name, "value is not valid", nonValidValue );
                }
            }
            else if ( name == "INTEGRATOR_FIXED_STEPSIZE" )
            {
                assignInsideRange( integratorFixedStepsize, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "INTEGRATOR_INITIAL_STEPSIZE" )
            {
                assignInsideRange( integratorInitialStepsize, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "INTEGRATOR_ERROR_TOLERANCE" )
            {
                assignInsideRange( integratorErrorTolerance, name, value, 0.0, 0.0, false, false, res, false, true );
            }
            else if ( name == "OUTPUT_NUMERIC_PRECISION" )
            {
                assignInsideRange( outputNumericPrecision, name, value, 1, 0, true, false, res, false, true );
            }
            else if ( name == "OUTPUT_INPUT_SETTINGS" )
            {
                assignBoolean( outputInputSettings, name, value, res );
            }
            else if ( name == "OUTPUT_RESULTS_COLUMNS_DESCRIPTIONS" )
            {
                assignBoolean( outputResultsColumnsDescriptions, name, value, res );
            }
            else if ( name == "OUTPUT_COMPUTATION_TIME" )
            {
                assignBoolean( outputComputationTime, name, value, res );
            }
            else if ( name == "OUTPUT_BODY_KEPLERIAN_STATE" )
            {
                assignBoolean( outputBodyKeplerianState, name, value, res );
            }
            else if ( name == "OUTPUT_BODY_CARTESIAN_STATE" )
            {
                assignBoolean( outputBodyCartesianState, name, value, res );
            }
            else if ( name == "OUTPUT_SUN_POSITION" )
            {
                assignBoolean( outputSunPosition, name, value, res );
            }
            else if ( name == "OUTPUT_MOON_POSITION" )
            {
                assignBoolean( outputMoonPosition, name, value, res );
            }
            else if ( name == "OUTPUT_ONE_IN_EVERY_INTEGRATION_STEPS" )
            {
                assignInsideRange( outputOneInEveryIntegrationSteps, name, value,
                                   1, 0, true, false, res, false, true );
            }
            else if ( name == "OUTPUT_ONLY_LAST_INTEGRATION_STEP" )
            {
                assignBoolean( outputOnlyLastIntegrationStep, name, value, res );
            }
            else if ( name == "OUTPUT_DIRECTORY_PATH" )
            {
                outputDirectoryPath = value;
                outputRootDirectory = definitionFilePath.parent_path();
            }
            else if ( name == "OUTPUT_FILE_NAME" )
            {
                outputFileName = value;
            }
            else if ( name == "PRINT_SUCCESS_MESSAGE_AFTER_PROPAGATION" )
            {
                assignBoolean( printSuccessMessageAfterPropagation, name, value, res );
            }
            else if ( name == "DELETE_INPUT_FILE_AFTER_PROPAGATION" )
            {
                assignBoolean( deleteInputFileAfterPropagation, name, value, res );
            }
            else if ( name == "SEARCH_OUTPUT_DIRECTORY_PATH" )
            {
                searchOutputDirectoryPath = value;
                searchOutputRootDirectory = definitionFilePath.parent_path();
            }
            else if ( name == "SEARCH_OUTPUT_FILENAME_PATTERN" )
            {
                searchOutputFilenamePattern = value;
            }
            else if ( name == "SEARCH_IN_SUBDIRECTORIES" )
            {
                assignBoolean( searchInSubdirectories, name, value, res );
            }
            else if ( name == "SKIP_PROPAGATION_IF_OUTPUT_FILE_EXISTS" )
            {
                assignBoolean( skipPropagationIfOutputFileExists, name, value, res );
            }
            else if ( name == "SKIP_PROPAGATION_IF_RESULTS_FOUND" )
            {
                assignBoolean( skipPropagationIfResultsFound, name, value, res );
            }
            else if ( name == "ALLOW_RESUMING_PROPAGATION_FROM_RESULTS" )
            {
                assignBoolean( allowResumingPropagationFromResults, name, value, res );
            }
            else if ( name == "OUTPUT_FOR_SKIPPED_PROPAGATION" )
            {
                if ( value == "NONE" )
                {
                    skippedPropagationOutput = SkippedPropagationOutput::none;
                }
                else if ( value == "REFERENCE" )
                {
                    skippedPropagationOutput = SkippedPropagationOutput::reference;
                }
                else if ( value == "DUPLICATE" )
                {
                    skippedPropagationOutput = SkippedPropagationOutput::duplicate;
                }
                else
                {
                    printWarningOrError( name, "value is not valid", nonValidValue );
                }
            }
            else if ( name == "SEARCH_IGNORE_FILE_CONTENTS" )
            {
                assignBoolean( searchIgnoreFileContents, name, value, res );
            }
        }
    }

    canonizePaths();
}


void TespSettings::canonizePaths()
{
    // Update output file name
    outputFileName = boost::regex_replace( outputFileName, boost::regex( R"(\$FILENAME)" ),
                                         inputFileName.string() );
    if ( outputFileName.find(".tespout") == std::string::npos )
    {
        outputFileName += ".tespout";
    }

    // If the root directories haven't been defined, they are the same as the input directory (given as full path)
    if ( outputRootDirectory.string().length() == 0 )
    {
        outputRootDirectory = inputDirectory;
    }
    if ( searchOutputRootDirectory.string().length() == 0 )
    {
        searchOutputRootDirectory = inputDirectory;
    }

    // Check if output directory exists; create it if it doesn't.
    path fullOutputDirectoryPath = outputRootDirectory / path( outputDirectoryPath );
    if ( ! boost::filesystem::exists( fullOutputDirectoryPath ) )
    {
        boost::filesystem::create_directories( fullOutputDirectoryPath );
    }

    // Determine the output full path
    path canonicalOutputPath = boost::filesystem::canonical( path( outputDirectoryPath ), outputRootDirectory );
    outputPath = canonicalOutputPath / path( outputFileName );

    // Check if search output directory exists; if it doesn't, do not allow skipping propagations.
    path fullSearchOutputDirectoryPath = searchOutputRootDirectory / path( searchOutputDirectoryPath );
    if ( ! boost::filesystem::exists( fullSearchOutputDirectoryPath ) )
    {
        skipPropagationIfResultsFound = false;
        skipPropagationIfOutputFileExists = false;
    }
    else
    {
        // Determine the search full path
        searchOutputDirectoryPath = boost::filesystem::canonical(
                    path( searchOutputDirectoryPath ), searchOutputRootDirectory ).string();
    }
}


void addSectionTitleToStream( std::stringstream &stream, std::string title )
{
    stream << std::setfill( '%' );
    stream  << "\n"
            << "\n"
            << std::right << std::setw(80) << "\n"
            << std::left  << std::setw(79) << "% " + title + " " << "\n"
            << std::right << std::setw(80) << "\n";
    stream << std::setfill( ' ' ) << std::left;
}

void addSubsectionTitleToStream( std::stringstream &stream, std::string title )
{
    stream  << "\n"
            << "% " << title << "\n";
}

void addNameToStream( std::stringstream &stream, std::string name )
{
    stream << std::setw(50) << name;
}

void addCommentToStream( std::stringstream &stream, std::string comment )
{
    if ( comment.length() > 0 )
    {
        stream << "% " << comment;
    }
    stream << "\n";
}

void addColumnDescriptionsToStream( std::stringstream &stream, int &currentColumn,
                                    std::vector< std::string> descriptions, std::vector< std::string> units )
{
    stream << "\n";
    for ( unsigned int i = 0; i < descriptions.size(); i++ )
    {
        stream << "COLUMN " << std::setw(2) << currentColumn++ << "  =  " <<
                  std::setw(40) << descriptions[ i ] << "% " << units[ i ] << "\n";
    }
}

void addToStream( std::stringstream &stream, std::string name, bool value, std::string comment = "" )
{
    addNameToStream( stream, name );
    stream << std::setw(20) << ( value ? "YES" : "NO" );
    addCommentToStream( stream, comment );
}

void addToStream( std::stringstream &stream, std::string name, int value, std::string comment = "" )
{
    addNameToStream( stream, name );
    stream << std::setw(20) << value;
    addCommentToStream( stream, comment );
}

void addToStream( std::stringstream &stream, std::string name, double value, std::string comment = "" )
{
    addNameToStream( stream, name );
    stream << std::setprecision(16) << std::setw(20) << value;
    addCommentToStream( stream, comment );
}

void addToStream( std::stringstream &stream, std::string name, std::string value,
                  bool quoted = false, std::string comment = "" )
{
    addNameToStream( stream, name );
    if ( quoted )
    {
        stream << "\"" << value << "\" ";
    }
    else
    {
        stream << std::setw(20) << value;
    }
    addCommentToStream( stream, comment );
}

void addToStream( std::stringstream &stream, std::string name, int value,
                  std::vector< std::string > enumValues, std::string comment = "" )
{
    addToStream( stream, name, enumValues[ value ], false, comment );
}


std::string TespSettings::outputContentHeader( bool printIncludedFilePaths )
{
    std::stringstream stream;

    stream << std::left
           << "% \n"
           << "% Output file generated by TESP\n"
           << "% \n";

    if ( printIncludedFilePaths && includedPaths.size() > 0 )
    {
        stream << "% Included files:\n";
        for ( path &includedPath : includedPaths )
        {
            stream << "% " << includedPath.string() << "\n";
        }
        stream << "% \n";
    }

    stream << "% Original input file:\n"
           << "% " << inputPath.string() << "\n"
           << "% \n";

    if ( duplicateSourcePath.string().size() > 0 )
    {
        stream << "% Propagation was skipped. The requested results were copied from:\n"
               << "% " << duplicateSourcePath.string() << "\n"
               << "% \n";
    }

    return stream.str();
}


std::string TespSettings::outputContentInputSettings() {
    std::stringstream stream;
    stream << std::left;

    addSectionTitleToStream( stream, "INPUT SETTINGS" );

    addSubsectionTitleToStream( stream, "ACCELERATION MODEL" );
    std::vector< std::string > geopotentialModels = { "GGM02C", "GGM02S" };
    addToStream( stream, "GEOPOTENTIAL_MODEL",                      geopotentialModel, geopotentialModels );
    addToStream( stream, "GEOPOTENTIAL_DEGREE",                     geopotentialDegree );
    addToStream( stream, "GEOPOTENTIAL_ORDER",                      geopotentialOrder );
    addToStream( stream, "THIRD_BODY_ATTRACTION_SUN",               thirdBodyAttractionSun );
    addToStream( stream, "THIRD_BODY_ATTRACTION_MOON",              thirdBodyAttractionMoon );
    std::vector< std::string > atmosphericModels = { "EXPONENTIAL", "", "NRLMSISE00" };
    addToStream( stream, "ATMOSPHERIC_MODEL",                       atmosphericModel, atmosphericModels );
    addToStream( stream, "ATMOSPHERIC_DRAG",                        atmosphericDrag );
    addToStream( stream, "SOLAR_RADIATION_PRESSURE",                solarRadiationPressure );
    addToStream( stream, "IGNORE_ECLIPSES",                         ignoreEclipses );

    addSubsectionTitleToStream( stream, "BODY PROPERTIES" );
    addToStream( stream, "BODY_MASS",                               bodyMass );
    addToStream( stream, "BODY_CROSS_SECTIONAL_AREA",               bodyCrossSectionalArea );
    addToStream( stream, "ALTITUDE_DEPENDENT_DRAG_COEFFICIENT",     altitudeDependentDragCoefficient );
    if ( ! altitudeDependentDragCoefficient ) {
        addToStream( stream, "BODY_DRAG_COEFFICIENT",               bodyDragCoefficient );
    }
    addToStream( stream, "BODY_RADIATION_PRESSURE_COEFFICIENT",     bodyRadiationPressureCoefficient );

    addSubsectionTitleToStream( stream, "INITIAL STATE" );
    addToStream( stream, "INITIAL_EPOCH",                           initialEpoch, "seconds since J2000" );
    if ( initialPerigeeAltitude == -1 ) {
        addToStream( stream, "INITIAL_SEMIMAJOR_AXIS",              initialSemimajorAxis, "km" );
    } else {
        addToStream( stream, "INITIAL_PERIGEE_ALTITUDE",            initialPerigeeAltitude, "km" );
    }
    if ( initialApogeeAltitude == -1 ) {
        addToStream( stream, "INITIAL_ECCENTRICITY",                initialEccentricity );
    } else {
        addToStream( stream, "INITIAL_APOGEE_ALTITUDE",             initialApogeeAltitude, "km" );
    }
    addToStream( stream, "INITIAL_INCLINATION",                     initialInclination, "deg" );
    addToStream( stream, "INITIAL_ARGUMENT_PERIGEE",                initialArgumentPerigee, "deg" );
    addToStream( stream, "INITIAL_LONGITUDE_ASCENDING_NODE",        initialLongitudeAscendingNode, "deg" );
    addToStream( stream, "INITIAL_TRUE_ANOMALY",                    initialTrueAnomaly, "deg" );

    addSubsectionTitleToStream( stream, "PROPAGATION SETTINGS" );
    // addToStream( stream, "PRELOAD_CELESTIAL_BODIES_DATA",        preloadCelestialBodiesData );
    addToStream( stream, "END_EPOCH",                               endEpoch, "seconds since J2000" );
    addToStream( stream, "PROPAGATION_PERIOD",                      propagationPeriod, "sidereal years" );
    addToStream( stream, "REENTRY_ALTITUDE",                        reentryAltitude, "km" );
    std::vector< std::string > propagatorTypes = { "COWELL", "ENCKE" };
    addToStream( stream, "PROPAGATOR_NAME",                         propagatorType, propagatorTypes );
    if ( integratorType == IntegratorType::rungeKutta4 ) {
        addToStream( stream, "INTEGRATOR_NAME",                     std::string( "RK4" ) );
        addToStream( stream, "INTEGRATOR_FIXED_STEPSIZE",           integratorFixedStepsize, "s" );
    } else {
        std::vector< std::string > integratorSets = { "RK45", "RK56", "RK78", "DP78" };
        addToStream( stream, "INTEGRATOR_NAME",                     integratorSet, integratorSets );
        addToStream( stream, "INTEGRATOR_INITIAL_STEPSIZE",         integratorInitialStepsize, "s" );
        addToStream( stream, "INTEGRATOR_ERROR_TOLERANCE",          integratorErrorTolerance );
    }

    addSubsectionTitleToStream( stream, "OUTPUT SETTINGS" );
    addToStream( stream, "OUTPUT_NUMERIC_PRECISION",                outputNumericPrecision );
    addToStream( stream, "OUTPUT_COMPUTATION_TIME",                 outputComputationTime );
    addToStream( stream, "OUTPUT_BODY_KEPLERIAN_STATE",             outputBodyKeplerianState );
    addToStream( stream, "OUTPUT_BODY_CARTESIAN_STATE",             outputBodyCartesianState );
    addToStream( stream, "OUTPUT_SUN_POSITION",                     outputSunPosition );
    addToStream( stream, "OUTPUT_MOON_POSITION",                    outputMoonPosition );
    addToStream( stream, "OUTPUT_ONE_IN_EVERY_INTEGRATION_STEPS",   outputOneInEveryIntegrationSteps );
    addToStream( stream, "OUTPUT_ONLY_LAST_INTEGRATION_STEP",      outputOnlyLastIntegrationStep );
    // addToStream( stream, "OUTPUT_DIRECTORY_PATH",                   outputDirectoryPath, true );
    // addToStream( stream, "OUTPUT_FILE_NAME",                        outputFileName, true );

    return stream.str();
}


std::string TespSettings::outputContentOutputDescription() {
    std::stringstream stream;
    stream << std::left;

    addSectionTitleToStream( stream, "OUTPUT DESCRIPTION" );

    int currentCol = 0;

    std::vector< std::string > descriptions = { "EPOCH" };
    std::vector< std::string > units = { "seconds since J2000" };
    addColumnDescriptionsToStream( stream, currentCol, descriptions, units );

    if ( outputBodyKeplerianState )
    {
        std::vector< std::string > descriptions = { "BODY SEMIMAJOR AXIS",
                                                    "BODY ECCENTRICITY",
                                                    "BODY INCLINATION",
                                                    "BODY ARGUMENT OF PERIGEE",
                                                    "BODY LONGITUDE OF THE ASCENDING NODE",
                                                    "BODY TRUE ANOMALY" };
        std::vector< std::string > units = { "m", "-", "rad", "rad", "rad", "rad" };
        addColumnDescriptionsToStream( stream, currentCol, descriptions, units );
    }

    if ( outputBodyCartesianState )
    {
        std::vector< std::string > descriptions = { "BODY POSITION X", "BODY POSITION Y", "BODY POSITION Z",
                                                    "BODY VELOCITY X", "BODY VELOCITY Y", "BODY VELOCITY Z" };
        std::vector< std::string > units = { "m", "m", "m", "m/s", "m/s", "m/s" };
        addColumnDescriptionsToStream( stream, currentCol, descriptions, units );
    }

    if ( outputSunPosition )
    {
        std::vector< std::string > descriptions = { "SUN POSITION X", "SUN POSITION Y", "SUN POSITION Z" };
        std::vector< std::string > units = { "m", "m", "m" };
        addColumnDescriptionsToStream( stream, currentCol, descriptions, units );
    }

    if ( outputMoonPosition )
    {
        std::vector< std::string > descriptions = { "MOON POSITION X", "MOON POSITION Y", "MOON POSITION Z" };
        std::vector< std::string > units = { "m", "m", "m" };
        addColumnDescriptionsToStream( stream, currentCol, descriptions, units );
    }

    return stream.str();
}


std::string TespSettings::outputContentOutputResults() {
    std::stringstream stream;
    stream << std::left;

    addSectionTitleToStream( stream, "OUTPUT RESULTS" );
    stream << "\n";

    int width = outputNumericPrecision + 7;
    stream << std::left << std::setprecision( outputNumericPrecision );

    // Get settings description and add computation time
    if ( outputComputationTime )
    {
        stream << "COMPUTATION_TIME = " << computationTime << "  % s\n\n";
    }

    // Iterate through map and add values to stream
    stream << "RESULTS = \n";
    for( std::map< double, Eigen::VectorXd >::iterator iter = resultsMap.begin();
         iter != resultsMap.end(); ++iter )
    {
        stream << std::setw( width ) << iter->first;
        Eigen::VectorXd vector = iter->second;
        for ( unsigned int j = 0; j < vector.rows(); j++ )
        {
            stream << std::setw( width ) << vector[ j ];
        }
        stream << "\n";
    }

    return stream.str();
}


bool TespSettings::checkRequestedOutputExists()
{
    // If the user does not allow to skip propagations, return directly with false (as if the output did not exist)
    if ( ! skipPropagationIfResultsFound && ! skipPropagationIfOutputFileExists )
    {
        return false;
    }

    // If the user specifies to skip propagation if a file already exists in the output path, then return with true
    if ( skipPropagationIfOutputFileExists && boost::filesystem::exists( outputPath ) )
    {
        std::cout << "Skipped propagation of: \"" << inputPath.string() << "\"" << std::endl;
        return true;
    }

    // Try to find previous output files with the results that have been requested.
    // If any error occurs, then return false, i.e. propagate the orbit.
    // We don't want a corrupt output file to prevent the propagation from being carried out...
    /*try
    {*/
        bool requestedOutputExists = false;

        // Determine whether the requested output already exists
        std::vector< path > candidatePaths;

        if ( searchOutputFilenamePattern.find(".tespout") == std::string::npos )
        {
            searchOutputFilenamePattern += ".tespout";
        }
        boost::regex filenamePattern( searchOutputFilenamePattern );

        if ( searchInSubdirectories )
        {
            boost::filesystem::recursive_directory_iterator iter( searchOutputDirectoryPath );
            boost::filesystem::recursive_directory_iterator iter_end;
            for ( ; iter != iter_end; iter++ )
            {
                // Skip if not a file
                if( ! boost::filesystem::is_regular_file( iter->status() ) ) continue;

                // Skip if no match
                std::string filename = iter->path().filename().string();
                boost::smatch match;
                if ( ! boost::regex_search( filename, match, filenamePattern,
                                            boost::match_not_dot_newline ) ) continue;

                candidatePaths.push_back( iter->path() );
            }
        }
        else
        {
            boost::filesystem::directory_iterator iter( searchOutputDirectoryPath );
            boost::filesystem::directory_iterator iter_end;
            for ( ; iter != iter_end; iter++ )
            {
                // Skip if not a file
                if( ! boost::filesystem::is_regular_file( iter->status() ) ) continue;

                // Skip if no match
                std::string filename = iter->path().filename().string();
                boost::smatch match;
                if ( ! boost::regex_search( filename, match, filenamePattern,
                                          boost::match_not_dot_newline ) ) continue;

                candidatePaths.push_back( iter->path() );
            }
        }

        // Determine if any of the candidates matches `this` (same input settings and at least same output requested)
        path matchingCandidatePath;
        TespSettings matchingCandidate;
        for ( path &candidate : candidatePaths )
        {
            // Initialize specifying that it is an output file what is being read ( second argument = true )
            TespSettings candidateSettings( candidate.string(), true );

            // Continue to the next candidate if any of the comparisons fails
            if ( ! searchIgnoreFileContents )
            {
                // ACCELERATION MODEL
                if ( candidateSettings.geopotentialModel != geopotentialModel ) continue;
                if ( candidateSettings.geopotentialDegree != geopotentialDegree ) continue;
                if ( candidateSettings.geopotentialOrder != geopotentialOrder ) continue;
                if ( candidateSettings.thirdBodyAttractionSun != thirdBodyAttractionSun ) continue;
                if ( candidateSettings.thirdBodyAttractionMoon != thirdBodyAttractionMoon ) continue;
                if ( candidateSettings.atmosphericModel != atmosphericModel ) continue;
                if ( candidateSettings.atmosphericDrag != atmosphericDrag ) continue;
                if ( candidateSettings.solarRadiationPressure != solarRadiationPressure ) continue;
                if ( candidateSettings.ignoreEclipses != ignoreEclipses ) continue;

                // BODY PROPERTIES
                if ( candidateSettings.bodyMass != bodyMass ) continue;
                if ( candidateSettings.bodyCrossSectionalArea != bodyCrossSectionalArea ) continue;
                if ( candidateSettings.altitudeDependentDragCoefficient != altitudeDependentDragCoefficient )
                    continue;
                if ( candidateSettings.bodyDragCoefficient != bodyDragCoefficient ) continue;
                if ( candidateSettings.bodyRadiationPressureCoefficient != bodyRadiationPressureCoefficient )
                    continue;

                // INITIAL STATE
                if ( candidateSettings.initialEpoch != initialEpoch ) continue;
                if ( candidateSettings.initialSemimajorAxis != initialSemimajorAxis ) continue;
                if ( candidateSettings.initialEccentricity != initialEccentricity ) continue;
                if ( candidateSettings.initialInclination != initialInclination ) continue;
                if ( candidateSettings.initialArgumentPerigee != initialArgumentPerigee ) continue;
                if ( candidateSettings.initialLongitudeAscendingNode != initialLongitudeAscendingNode ) continue;
                if ( candidateSettings.initialTrueAnomaly != initialTrueAnomaly ) continue;
                if ( candidateSettings.initialPerigeeAltitude != initialPerigeeAltitude ) continue;
                if ( candidateSettings.initialApogeeAltitude != initialApogeeAltitude ) continue;

                // PROPAGATION SETTINGS
                if ( candidateSettings.maximumFinalEpoch() < maximumFinalEpoch()
                     && ! allowResumingPropagationFromResults ) continue;
                if ( candidateSettings.reentryAltitude != reentryAltitude ) continue;
                if ( candidateSettings.propagatorType != propagatorType ) continue;
                if ( candidateSettings.integratorType != integratorType ) continue;
                if ( integratorType == IntegratorType::rungeKutta4 )
                {
                    if ( candidateSettings.integratorFixedStepsize != integratorFixedStepsize ) continue;
                }
                else
                {
                    if ( candidateSettings.integratorSet != integratorSet ) continue;
                    if ( candidateSettings.integratorInitialStepsize != integratorInitialStepsize ) continue;
                    if ( candidateSettings.integratorErrorTolerance != integratorErrorTolerance ) continue;
                }

                // OUTPUT SETTINGS
                if ( candidateSettings.outputNumericPrecision < outputNumericPrecision ) continue;
                if ( candidateSettings.outputComputationTime < outputComputationTime ) continue;
                if ( candidateSettings.outputBodyKeplerianState < outputBodyKeplerianState ) continue;
                if ( candidateSettings.outputBodyCartesianState < outputBodyCartesianState ) continue;
                if ( candidateSettings.outputSunPosition < outputSunPosition ) continue;
                if ( candidateSettings.outputMoonPosition < outputMoonPosition ) continue;
                if ( candidateSettings.outputOneInEveryIntegrationSteps != outputOneInEveryIntegrationSteps )
                        continue;
                if ( candidateSettings.outputOnlyLastIntegrationStep > outputOnlyLastIntegrationStep ) continue;
            }

            // If no conditions has failed, then we've found a matching candidate
            matchingCandidatePath = candidate;
            matchingCandidate = candidateSettings;
            requestedOutputExists = true;
            break;
        }

        if ( requestedOutputExists )
        {
            // If can propagate further, but no states found in matching file, then start propagation from zero
            if ( maximumFinalEpoch() > matchingCandidate.maximumFinalEpoch() &&
                 ! matchingCandidate.outputBodyCartesianState && ! matchingCandidate.outputBodyKeplerianState )
            {
                return false;
            }

            // Determine if propagation is going to be resumed from output
            resuming = allowResumingPropagationFromResults &&
                    matchingCandidate.maximumFinalEpoch() < maximumFinalEpoch() &&
                    ( matchingCandidate.outputBodyCartesianState || matchingCandidate.outputBodyKeplerianState );

            // Warn if necessary
            if ( skippingPropagation == warn && ! resuming )
            {
                std::cout << "Skipped propagation of: \"" << inputPath.string() << "\"" << std::endl;
            }

            // If output exists, the propagation is not going to be resumed and the file is the same, don't propagate
            // Return true, i.e. tell the propagator that the output already exists, and it will do nothing
            if ( ! resuming && ( matchingCandidatePath.string() == outputPath.string() ) )
            {
                return true;
            }

            // Generate output with a reference to file
            if ( skippedPropagationOutput == reference && matchingCandidatePath != outputPath && ! resuming )
            {
                std::ofstream outputFile( outputPath.string() );
                outputFile << outputContentHeader( false )
                           << "% Requested results were found in an output file from a previous propagation:\n"
                           << "% " << matchingCandidatePath.string() << "\n"
                           << "% \n";
                outputFile.close();
            }
            else if ( skippedPropagationOutput != none ) // or duplicate file, resume, etc.
            {
                // Change the matchingCandidate paths to the ones for the current run
                matchingCandidate.inputPath = inputPath;
                matchingCandidate.inputDirectory = inputDirectory;
                matchingCandidate.inputFileName = inputFileName;
                matchingCandidate.outputDirectoryPath = outputDirectoryPath;
                matchingCandidate.outputFileName = outputFileName;
                matchingCandidate.outputRootDirectory = outputRootDirectory;
                matchingCandidate.searchOutputRootDirectory = searchOutputRootDirectory;
                matchingCandidate.canonizePaths();

                // Define duplicate original source path so that it can be indicated in the new file header
                matchingCandidate.duplicateSourcePath = matchingCandidatePath;

                // Update results: only duplicate those requested for this run
                matchingCandidate.parseResultsFromOutputFile();

                std::vector< bool > candidateOutputs = {
                    matchingCandidate.outputBodyKeplerianState, matchingCandidate.outputBodyCartesianState,
                    matchingCandidate.outputSunPosition, matchingCandidate.outputMoonPosition };

                std::vector< bool > outputs = { outputBodyKeplerianState, outputBodyCartesianState,
                                                outputSunPosition, outputMoonPosition };

                std::vector< int > outputsLengths = { 6, 6, 3, 3 };

                std::vector< int > candidateStartIndexes;
                std::vector< int > startIndexes;

                int currentCandidateStartIndex = 0;
                int currentStartIndex = 0;
                for ( unsigned int o = 0; o < outputsLengths.size(); o++ )
                {
                    if ( candidateOutputs[ o ] )
                    {
                        candidateStartIndexes.push_back( currentCandidateStartIndex );
                        currentCandidateStartIndex += outputsLengths[ o ];
                    }
                    else
                    {
                        candidateStartIndexes.push_back( -1 );
                    }
                    if ( outputs[ o ] )
                    {
                        startIndexes.push_back( currentStartIndex );
                        currentStartIndex += outputsLengths[ o ];
                    }
                    else
                    {
                        startIndexes.push_back( -1 );
                    }
                }


                // Populate new results map (only those requested by the user)

                std::map< double, Eigen::VectorXd > requestedResultsMap;
                for( std::map< double, Eigen::VectorXd >::iterator iter = matchingCandidate.resultsMap.begin();
                     iter != matchingCandidate.resultsMap.end(); ++iter )
                {
                    if ( ! outputOnlyLastIntegrationStep /*|| iter == --matchingCandidate.resultsMap.end()*/ )
                    {
                        double epoch = iter->first;

                        Eigen::VectorXd oldResult = iter->second;
                        Eigen::VectorXd requestedResult( currentStartIndex );
                        for ( unsigned int j = 0; j < outputs.size(); j++ )
                        {
                            if ( outputs[ j ] ) // if this specific output has been requested, get if from old output
                            {
                                requestedResult << oldResult.segment ( candidateStartIndexes[ j ],
                                                                       outputsLengths[ j ] );
                            }
                        }
                        requestedResultsMap[ epoch ] = requestedResult;

                        if ( epoch > maximumFinalEpoch() ) break;
                    }
                }
                resultsMap = requestedResultsMap;


                // Handle resume propagation from output
                if ( resuming )
                {
                    if ( resumingPropagationFromOutput == warn )
                    {
                        std::cout << "Resuming propagation defined in input file \""
                                  << inputPath.string() << "\" from the results found in output file \""
                                  << matchingCandidatePath.string() << "\"" << std::endl;
                    }

                    double lastEpoch = (--matchingCandidate.resultsMap.end())->first;
                    Eigen::VectorXd lastVector = (--matchingCandidate.resultsMap.end())->second;

                    resumingEpoch = lastEpoch;

                    resumingState = lastVector.segment< 6 >(
                                candidateStartIndexes[ ! matchingCandidate.outputBodyKeplerianState ] );
                    resumingStateIsCartesian = ! matchingCandidate.outputBodyKeplerianState;
                }


                // Duplicate file
                if ( ! resuming )
                {
                    // Update requested output variables

                    matchingCandidate.outputNumericPrecision = outputNumericPrecision;
                    matchingCandidate.outputComputationTime = outputComputationTime;
                    matchingCandidate.outputBodyKeplerianState = outputBodyKeplerianState;
                    matchingCandidate.outputBodyCartesianState = outputBodyCartesianState;
                    matchingCandidate.outputSunPosition = outputSunPosition;
                    matchingCandidate.outputMoonPosition = outputMoonPosition;
                    matchingCandidate.outputOnlyLastIntegrationStep = outputOnlyLastIntegrationStep;

                    matchingCandidate.exportResults();
                }
            }
        }

        return requestedOutputExists && ! resuming;
    /*}
    catch ( ... )
    {
        return false;
    }*/
}


void TespSettings::exportResults()
{
    std::ofstream outputFile( outputPath.string() );
    outputFile << outputContentHeader();
    if ( outputInputSettings )
    {
        outputFile << outputContentInputSettings();
    }
    if ( outputResultsColumnsDescriptions )
    {
        outputFile << outputContentOutputDescription();
    }
    outputFile << outputContentOutputResults();
    outputFile.close();
}


void TespSettings::propagationSucceeded()
{
    if ( printSuccessMessageAfterPropagation )
    {
        std::cout << "Propagation of " << inputPath << " succeeded." << std::endl;
    }

    if ( deleteInputFileAfterPropagation )
    {
        boost::filesystem::remove( inputPath );
    }
}


} // namespace tesp
