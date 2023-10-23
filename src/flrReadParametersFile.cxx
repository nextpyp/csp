/*=========================================================================

 Name:        flrReadParametersFile.cxx
 Category:    Implementation file
 Description: Implementation of the parameters structure and related functions.
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrReadParametersFile.cxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/
#include "flrReadParametersFile.hxx"
#include "vtkMath.h"
#include <refine3d_cspt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <math.h>
#include <cstdlib>

/** not allow to use std::to_string using old compiler, so write an own one ...
    used to convert any type of objects to string  
*/
template <typename T>
std::string to_string(T value) {

    std::ostringstream os ;
    os << value ;
    return os.str() ;
}
/** Helper function - check if the file exists
 */
bool checkExistence(const char* filename) {

    ifstream Infield(filename);
    return Infield.good();
}


std::string RemoveQuotations( char value[] ) {
    std::string str = to_string(value);
    str.erase(remove(str.begin(), str.end(), '"'), str.end());
    return str;
}

/**
 * Read PYP parameters supporting serial iterations (i.e. refine_rhref = "1:2:3:4")
 */
double GetParamInListByIter( char value[], int iter, std::string delimiter ) {
    
    // we don't know if it is a sequence of boolean ("T:F:T:F") or numbers ("20:14:10:7")
    vector< double > params;
    
    // remove the quotation marks 
    std::string str = RemoveQuotations( value );
    //str.erase(remove(str.begin(), str.end(), '"'), str.end());
    
    double param;
    std::string substr;

    int start = 0;
    int end = str.find(delimiter);

    while (end!=-1) {
        substr = str.substr(start, end - start);
        
        if ( tolower(substr[0]) == 't' ) param = 1.0;
        else if ( tolower(substr[0]) == 'f' ) param = 0.0;
        else param = atof( substr.c_str() );

        params.push_back(param);
        start = end + delimiter.size();
        end = str.find(delimiter, start);
    }
    substr = str.substr(start, end - start).c_str();
        
    if ( tolower(substr[0]) == 't' ) param = 1.0;
    else if ( tolower(substr[0]) == 'f' ) param = 0.0;
    else param = atof( substr.c_str() );

    params.push_back(param);


    try {
        param = params.at( iter-2 );
    }
    catch ( const std::out_of_range& oor ) {
        param = params.at( params.size()-1 );
    }

    return param;
}




void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] );

void ProcessCommandLineOptions( int argc, char * argv[], SParameters *params ) {
}

/**
 * int LoadConfiguration( SParameters *params, char filename[80] )
 * */
int LoadConfiguration( SParameters *params, char filename[80] ) {
  char line[120], str1[80], str2[80];
  FILE *file;

  file = fopen( filename, "r" );
  if (!file) {
    printf( "File %s does not exist.\n", filename );
    return 1;
  }
  while (!feof( file )) {
    fgets( line, sizeof line, file ); // Get one line
    sscanf( line, "%s %s", str1, str2 ); // Scan the line

    if (!strcmp( str1, "BFactor" )) params->mBFactor = atof( str2 );
    if (!strcmp( str1, "InnerRadius" )) params->mInnerRadius = atof( str2 );
    if (!strcmp( str1, "OuterRadius" )) params->mOuterRadius = atof( str2 );
    //if (!strcmp( str1, "DataFolder" )) strcpy( params->mDataFolder, str2 );
    if (!strcmp( str1, "WorkingFolder" )) strcpy( params->mWorkingFolder, str2 );
    //if (!strcmp( str1, "FrealignFolder" )) strcpy( params->mFrealignFolder, str2 );
    if (!strcmp( str1, "NumberOfRings" )) params->mNumberOfRings = atoi( str2 );
    if (!strcmp( str1, "ParFile" )) strcpy( params->mParFile, str2 );
    if (!strcmp( str1, "OptimizerMaxIter" )) params->mOptimizerMaxIter = atof( str2 );
    if (!strcmp( str1, "OptimizerValueTolerance" )) params->mOptimizerValueTolerance = atof( str2 );
    if (!strcmp( str1, "NumberOfIterations" )) params->mNumberOfIterations = atoi( str2 );
    if (!strcmp( str1, "NumberOfRandomIterations" )) params->mNumberOfRandomIterations = atoi( str2 );
    if (!strcmp( str1, "OptimizerStepTolerance" )) params->mOptimizerStepTolerance = atof( str2 );
    if (!strcmp( str1, "OptimizerStepLength" )) params->mOptimizerStepLength = atof( str2 );
    if (!strcmp( str1, "NumberOfThreads" )) params->mNumberOfThreads = atoi( str2 );

    if (!strcmp( str1, "RefineMicrographTiltAngles" )) params->mRefineMicrographTiltAngles = atoi( str2 );
    if (!strcmp( str1, "RefineMicrographTiltAxisAngles" )) params->mRefineMicrographTiltAxisAngles = atoi( str2 );
    if (!strcmp( str1, "RefineParticlesTheta" )) params->mRefineParticlesTheta = atoi( str2 );
    if (!strcmp( str1, "RefineParticlesPsi" )) params->mRefineParticlesPsi = atoi( str2 );
    if (!strcmp( str1, "RefineParticlesPhi" )) params->mRefineParticlesPhi = atoi( str2 );

    if (!strcmp( str1, "ToleranceMicrographTiltAngles" )) params->mToleranceMicrographTiltAngles = atof( str2 );
    if (!strcmp( str1, "ToleranceMicrographTiltAxisAngles" )) params->mToleranceMicrographTiltAxisAngles = atof( str2 );
    if (!strcmp( str1, "ToleranceMicrographShifts" )) params->mToleranceMicrographShifts = atof( str2 );
    
    if (!strcmp( str1, "ToleranceMicrographDefocus1" )) params->mToleranceMicrographDefocus1 = atof( str2 );
    if (!strcmp( str1, "ToleranceMicrographDefocus2" )) params->mToleranceMicrographDefocus2 = atof( str2 );
    if (!strcmp( str1, "ToleranceMicrographAstigmatism" )) params->mToleranceMicrographAstigmatism = atof( str2 );

    if (!strcmp( str1, "ToleranceParticlesTheta" )) params->mToleranceParticlesTheta = atof( str2 );
    if (!strcmp( str1, "ToleranceParticlesPsi" )) params->mToleranceParticlesPsi = atof( str2 );
    if (!strcmp( str1, "ToleranceParticlesPhi" )) params->mToleranceParticlesPhi = atof( str2 );
    if (!strcmp( str1, "ToleranceParticlesShifts" )) params->mToleranceParticlesShifts = atof( str2 );

    if (!strcmp( str1, "UseImagesForRefinementMax" )) params->mUseImagesForRefinementMax = atoi( str2 );
    if (!strcmp( str1, "UseImagesForRefinementMin" )) params->mUseImagesForRefinementMin = atoi( str2 );
    
    if (!strcmp( str1, "UseImagesForReconstructionMax" )) params->mUseImagesForReconstructionMax = atoi( str2 );
    if (!strcmp( str1, "UseImagesForReconstructionMin" )) params->mUseImagesForReconstructionMin = atoi( str2 );

    if (!strcmp( str1, "DebugFull" )) params->mDebugFull = atoi( str2 );
    if (!strcmp( str1, "DebugData" )) params->mDebugData = atoi( str2 );
    if (!strcmp( str1, "DebugInfo" )) params->mDebugInfo = atoi( str2 );
    if (!strcmp( str1, "DebugNone" )) params->mDebugNone = atoi( str2 );
    if (!strcmp( str1, "DebugBasic" )) params->mDebugBasic = atoi( str2 );

  }

  // Generate other variables.
  strncpy( params->mBaseName, params->mParFile, strcspn( params->mParFile, "." ) );
  sprintf( params->mFrealignFolder, "%s/frealign/", params->mWorkingFolder );
  sprintf( params->mDataFolder, "%s/data/", params->mWorkingFolder );

  fclose( file );
  return 0;
}

int LoadParameters( Refine3dParameters *refine_params, SParameters *csp_params, char pyp_config[80] ) {

    char line[1000], key[200], equal[2], value[800];
    bool use_fboost = false, param_bool;
    FILE *file;
    double param_float;
    

    // get where the scratch folder is by PYP_SCRATCH environmental variable
    char* scratch_path = std::getenv("PYP_SCRATCH");
    if (scratch_path != NULL) {
        refine_params->scratch = scratch_path;
    }

    // First parse pyp_config file 
    file = fopen( pyp_config, "r" );
    
    if (!file) {
        printf( "%s file does not exist.\n", pyp_config );
        return EXIT_FAILURE;
    }
    // First get the current iteration number
    while( !feof( file ) ) {
        fgets( line, sizeof line, file );
        sscanf( line, "%s %s %s", key, equal, value );
        if (!strcmp( key, "refine_iter" )) refine_params->currentIter = atoi( value );
    }

    rewind( file );

    while ( !feof( file ) ) {
        fgets( line, sizeof line, file );
        sscanf( line, "%s %s %s", key, equal, value );
        
        /**
         * standard parameters
         */
        if (!strcmp( key, "data_mode" )) {
            if (RemoveQuotations(value).c_str()[0] == 's') {
                refine_params->isSPR = true;
            }
            else if (RemoveQuotations(value).c_str()[0] == 't') {
                refine_params->isTOMO = true;
            }
            else {
                printf("Data mode %s unrecognized. \n", RemoveQuotations(value).c_str());
                exit(1);
            }
        }
        if (!strcmp( key, "scope_pixel" )) refine_params->scopePixel = atof( value );
        if (!strcmp( key, "data_bin" )) refine_params->dataBin = atoi( value );
        if (!strcmp( key, "extract_bin" )) refine_params->particleBin = atoi( value );
        if (!strcmp( key, "scope_voltage" )) refine_params->beamEnergy = atof( value );
        if (!strcmp( key, "scope_cs" )) refine_params->sphericalAbberation = atof( value );
        if (!strcmp( key, "scope_wgh" )) refine_params->amplitutdeContrast = atof( value );
        if (!strcmp( key, "particle_mw" )) refine_params->molecularMass = atof( RemoveQuotations(value).c_str() );
        if (!strcmp( key, "extract_box" )) refine_params->boxsize = atoi( value );
        if (!strcmp( key, "particle_rad" )) { 
            refine_params->outerMask = atof( RemoveQuotations(value).c_str() );
            refine_params->maskRadiusGlobalSearch = refine_params->outerMask * 1.5;
        }
        if (!strcmp( key, "particle_sym" )) strcpy( refine_params->symmetry, RemoveQuotations(value).c_str() );
         
        /**
         * refine3d-related parameters
         */
        if (!strcmp( key, "data_set" )) strcpy( csp_params->mBaseName, RemoveQuotations( value ).c_str() );
        if ( (!strcmp( key, "refine_fssnr" )) && ( strlen(value) > 0 ) ) refine_params->useStatistics = GetParamInListByIter(value, refine_params->currentIter);
        if ( (!strcmp( key, "refine_priors" )) && ( strlen(value) > 0 ) ) refine_params->usePriors = GetParamInListByIter(value, refine_params->currentIter);
        if ( (!strcmp( key, "refine_fboost" )) && ( strlen(value) > 0) ) use_fboost = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "refine_fboostlim" )) refine_params->resLimitSignCC = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "refine_iblow" )) refine_params->paddingFactor = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "refine_rhref" )) refine_params->highResLimit = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "refine_rlref" )) refine_params->lowResLimit = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_particle_blur_sigma" )) refine_params->particle_blur_sigma = GetParamInListByIter(value, refine_params->currentIter); 
        /**
         * CSP-related parameters
         */
        if (!strcmp( key, "csp_OptimizerMaxIter" )) csp_params->mOptimizerMaxIter = atoi( value );
        if (!strcmp( key, "csp_OptimizerValueTolerance" )) csp_params->mOptimizerValueTolerance = atof( value );
        if (!strcmp( key, "csp_NumberOfRandomIterations" )) csp_params->mNumberOfRandomIterations = atoi( value );
        if (!strcmp( key, "csp_OptimizerValueTolerance" )) csp_params->mOptimizerStepTolerance = atof( value );
        if (!strcmp( key, "csp_OptimizerStepLength" )) csp_params->mOptimizerStepLength = atof( value );

        if (!strcmp( key, "csp_RefineProjectionCutoff" )) csp_params->mRefineProjectionCutoff = atoi( value );

        if (!strcmp( key, "csp_ToleranceMicrographTiltAngles" )) csp_params->mToleranceMicrographTiltAngles = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceMicrographTiltAxisAngles" )) csp_params->mToleranceMicrographTiltAxisAngles = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceMicrographShifts" )) csp_params->mToleranceMicrographShifts = GetParamInListByIter(value, refine_params->currentIter);
    
        if (!strcmp( key, "csp_ToleranceMicrographDefocus1" )) csp_params->mToleranceMicrographDefocus1 = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceMicrographDefocus2" )) csp_params->mToleranceMicrographDefocus2 = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceMicrographAstigmatism" )) csp_params->mToleranceMicrographAstigmatism = GetParamInListByIter(value, refine_params->currentIter);

        if (!strcmp( key, "csp_ToleranceParticlesTheta" )) csp_params->mToleranceParticlesTheta = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceParticlesPsi" )) csp_params->mToleranceParticlesPsi = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceParticlesPhi" )) csp_params->mToleranceParticlesPhi = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_ToleranceParticlesShifts" )) csp_params->mToleranceParticlesShifts = GetParamInListByIter(value, refine_params->currentIter);

        if (!strcmp( key, "csp_UseImagesForRefinementMax" )) csp_params->mUseImagesForRefinementMax = (int) GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_UseImagesForRefinementMin" )) csp_params->mUseImagesForRefinementMin = (int) GetParamInListByIter(value, refine_params->currentIter);
    
        if (!strcmp( key, "csp_UseImagesForReconstructionMax" )) csp_params->mUseImagesForReconstructionMax = (int) GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_UseImagesForReconstructionMin" )) csp_params->mUseImagesForReconstructionMin = (int) GetParamInListByIter(value, refine_params->currentIter);

        if (!strcmp( key, "csp_DebugFull" )) csp_params->mDebugFull = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_DebugData" )) csp_params->mDebugData = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_DebugInfo" )) csp_params->mDebugInfo = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_DebugNone" )) csp_params->mDebugNone = GetParamInListByIter(value, refine_params->currentIter);
        if (!strcmp( key, "csp_DebugBasic" )) csp_params->mDebugBasic = GetParamInListByIter(value, refine_params->currentIter);
        
    }
    refine_params->scopePixel = refine_params->scopePixel * refine_params->dataBin * refine_params->particleBin;
    if ( !use_fboost ) refine_params->resLimitSignCC = 30.0;
    if ( refine_params->resLimitSignCC < 2.0 * refine_params->scopePixel ) refine_params->resLimitSignCC = refine_params->scopePixel;
    if ( refine_params->highResLimit < 2.0 * refine_params->scopePixel ) refine_params->highResLimit = 2.0 * refine_params->scopePixel;
     
    if (csp_params->mToleranceParticlesPsi > 0) csp_params->mRefineParticlesPsi = true;
    if (csp_params->mToleranceParticlesPhi > 0) csp_params->mRefineParticlesPhi = true;
    if (csp_params->mToleranceParticlesTheta > 0) csp_params->mRefineParticlesTheta = true;
    
    if (csp_params->mUseImagesForRefinementMin < 0) csp_params->mUseImagesForRefinementMin = 0; 

    // statistic file from previous iteration 
    sprintf(refine_params->statisticFile, "frealign/scratch/%s_r01_%02d_statistics.txt", csp_params->mBaseName, refine_params->currentIter-1);

    fclose( file );
    
    return 0;
}



/********************************************************
 * void ShowConfiguration( SParameters *params )
 *******************************************************/
void ShowConfiguration( SParameters *params ) {
  printf( "[Parameters]\n" );
  printf( "B-factor: %f\n", params->mBFactor );
  printf( "InnerRadius: %f\n", params->mInnerRadius );
  printf( "Number of Rings: %d\n", params->mNumberOfRings );

  printf( "Debug Data: %d\n", params->mDebugData );
  printf( "Debug Full: %d\n", params->mDebugFull );
  printf( "Debug None: %d\n", params->mDebugNone );
  printf( "Debug Info: %d\n", params->mDebugInfo );

  printf( "Debug Basic: %d\n", params->mDebugBasic );
  printf( "OptimizerMaxIter: %d\n", params->mOptimizerMaxIter );
  printf( "OptimizerValueTolerance: %e\n", params->mOptimizerValueTolerance );
  printf( "OuterRadius: %f\n", params->mOuterRadius );
  printf( "Number of Threads: %d\n", params->mNumberOfThreads );

  printf( "ParFile: %s\n", params->mParFile );
  printf( "Basename: %s\n", params->mBaseName );
  printf( "Working Folder: %s\n", params->mWorkingFolder );
  printf( "Data Folder: %s\n", params->mDataFolder );
  printf( "Frealign Folder: %s\n", params->mFrealignFolder );
  printf( "Image size: %dx%d\n", params->mImageSizeX, params->mImageSizeY );
  printf( "Number of iterations: %d\n", params->mNumberOfIterations );
  printf( "Number of random iterations: %d\n", params->mNumberOfRandomIterations );
  printf( "OptimizerStepTolerance: %f\n", params->mOptimizerStepTolerance );
  printf( "OptimizerStepLength: %f\n", params->mOptimizerStepLength );

  printf( "Refine Micrograph Tilt Angles: %d\n", params->mRefineMicrographTiltAngles );
  printf( "Refine Micrograph Tilt Axis Angles: %d\n", params->mRefineMicrographTiltAxisAngles );
  printf( "Refine Particles Theta Euler Angle: %d\n", params->mRefineParticlesTheta );
  printf( "Refine Particles Psi Euler Angle: %d\n", params->mRefineParticlesPsi );
  printf( "Refine Particles Phi Euler Angle: %d\n", params->mRefineParticlesPhi );

  printf( "Tolerance Micrograph Tilt Angles: %f\n", params->mToleranceMicrographTiltAngles );
  printf( "Tolerance Micrograph Tilt Axis Angles: %f\n", params->mToleranceMicrographTiltAxisAngles );
  printf( "Tolerance Micrograph Shifts: %f\n", params->mToleranceMicrographShifts );

  printf( "Tolerance Micrograph Defocus 1: %f\n", params->mToleranceMicrographDefocus1 );
  printf( "Tolerance Micrograph Defocus 2: %f\n", params->mToleranceMicrographDefocus2 );
  printf( "Tolerance Micrograph Astigmatism: %f\n", params->mToleranceMicrographAstigmatism );
  
  printf( "Tolerance Particles Theta: %f\n", params->mToleranceParticlesTheta );
  printf( "Tolerance Particles Psi: %f\n", params->mToleranceParticlesPsi );
  printf( "Tolerance Particles Phi: %f\n", params->mToleranceParticlesPhi );
  printf( "Tolerance Particles Shifts: %f\n", params->mToleranceParticlesShifts );

  printf( "\n" );
}



/********************************************************
 * int ProcessConfiguration( char filename[80], SParameters *params )
 *******************************************************/
int ProcessConfiguration( SParameters *params, char parxfile[300], vector < CTiltSeries > & pTiltSeries, int mode, int first_index, int last_index, int frame_index, vector<double> & weights, vector<int> & images_to_extract, vector< pair<double, double> > & frame_shifts, map<int, float*> & priors_average, map<int, float*> & priors_variance ) {
    FILE *file;
    char filename[200];
    char line[2000];
    int ind, film, ptlind, scanor, counter, mag, logp;
    EulerAngleType psi, theta, phi, tiltan, tiltaxisan;
    CoordinateType shx, shy;
    float df1, df2, angast, ppindex, occ, sigma, score, change, dosexx, cnfdnc, ptlccx;
    CEulerAngles eulerAngles;
    CEulerAngles zeroEulerAngles( 0, 0, 0 );
    CPosition2D particleProjectionShift, micrographTiltAxisShift, frameShift, projectionFrameShift;
    EulerAngleType micrographTiltAxisAngle;
    EulerAngleType normal[3];
    EulerAngleType matrix[16];
    CMicrograph* pMicrograph;
    CParticle* pParticle;
    CParticleProjection* pPP;
    CPosition3D position3D( 0, 0, 0 );
    int particleIndex;
    int micrographNumberID;
    int tiltSeriesIndex;
    int ind_in_stack = 0;

    vector <int> countPerFrame;
    map<int, float*>::iterator it;
    map<int, int> countPerTilt;

    EulerAngleType ppsi, ptheta, pphi;
    ppsi = ptheta = pphi = 0;

    
    file = fopen( parxfile, "r" );
    if (!file) {
        printf( "File %s does not exist.\nExiting...\n", filename );
        exit( EXIT_FAILURE );
    }

    params->mTotalNumberOfParticlesProjections = 0;
       
    int input_actual_counter = 0;

    while (!feof( file )) {
        ind = -1;
        line[0] = 'C';
        fgets( line, sizeof line, file );
        if (line[0] != 'C') {
        
        sscanf(
                line,
                "%d %lf %lf %lf %lf %lf %d %d %f %f %f %f %d %f %f %f %d %lf %lf %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &ind, &psi, &theta, &phi, &shx, &shy, &mag, &film, &df1, &df2, &angast, &occ, &logp, &sigma, &score, &change, &ptlind, &tiltan,
                &dosexx, &scanor, &cnfdnc, &ptlccx, &tiltaxisan, &normal[0], &normal[1], &normal[2], &matrix[0],
                &matrix[1], &matrix[2], &matrix[3], &matrix[4], &matrix[5], &matrix[6], &matrix[7], &matrix[8],
                &matrix[9], &matrix[10], &matrix[11], &matrix[12], &matrix[13], &matrix[14], &matrix[15],
                &ppsi, &ptheta, &pphi	);
        
        
        if (priors_average.find(scanor) == priors_average.end()) {
            // not found 
            float* average = new float[NUM_COL_REFINE3D];
            float* variance = new float[NUM_COL_REFINE3D];
            for (int j = 0; j < NUM_COL_REFINE3D; j++) {
                average[j] = 0.0;
                variance[j] = 0.0; 
            }
            priors_average.insert( {scanor, average} );
            priors_variance.insert( {scanor, variance} );
            countPerTilt.insert( {scanor, 0} );
        }
        
        countPerTilt.at(scanor)++;

        priors_average.at(scanor)[0] += 0.0; 
        priors_average.at(scanor)[1] += phi; priors_average.at(scanor)[2] += theta; priors_average.at(scanor)[3] += psi; priors_average.at(scanor)[4] += shx; priors_average.at(scanor)[5] += shy;
        priors_average.at(scanor)[6] += mag; priors_average.at(scanor)[7] += film; priors_average.at(scanor)[8] += df1; priors_average.at(scanor)[9] += df2; priors_average.at(scanor)[10] += angast; priors_average.at(scanor)[11] += ppindex;
        priors_average.at(scanor)[12] += occ; priors_average.at(scanor)[13] += logp; priors_average.at(scanor)[14] += sigma; priors_average.at(scanor)[15] += score; priors_average.at(scanor)[16] += change; priors_average.at(scanor)[17] += 0.0;
        
        priors_variance.at(scanor)[0] += powf(0.0, 2); 
        priors_variance.at(scanor)[1] += powf(phi, 2); priors_variance.at(scanor)[2] += powf(theta, 2); priors_variance.at(scanor)[3] += powf(psi, 2); priors_variance.at(scanor)[4] += powf(shx, 2); priors_variance.at(scanor)[5] += powf(shy, 2);
        priors_variance.at(scanor)[6] += powf(mag, 2); priors_variance.at(scanor)[7] += powf(film, 2); priors_variance.at(scanor)[8] += powf(df1, 2); priors_variance.at(scanor)[9] += powf(df2, 2); priors_variance.at(scanor)[10] += powf(angast, 2); priors_variance.at(scanor)[11] += powf(ppindex, 2);
        priors_variance.at(scanor)[12] += powf(occ, 2); priors_variance.at(scanor)[13] += powf(logp, 2); priors_variance.at(scanor)[14] += powf(sigma, 2); priors_variance.at(scanor)[15] += powf(score, 2); priors_variance.at(scanor)[16] += powf(change, 2); priors_variance.at(scanor)[17] += powf(0.0, 2);
        
        input_actual_counter++;
        
        // compute frame weight 
        while ( weights.size() < cnfdnc + 1 ) {
            weights.push_back( 0.0 );
            countPerFrame.push_back( 0 );
            refine3dParams.numberFrames = weights.size();
        }
        weights[cnfdnc] += score;
        countPerFrame[cnfdnc] ++;

        // skip particle/micrograph that are not in the user-specified range 
        if ( (mode == 0) || (mode == 4) || (mode == 6) || (mode == -1)) {
            if ( ( scanor < first_index ) || ( ( scanor > last_index ) && ( last_index != -1 ) ) ) {
                continue;
            }
        }
        else {
            if ( ( ptlind < first_index ) || ( ( ptlind > last_index ) && ( last_index != -1 ) ) ) {
                continue;
            }
        }
        

        // images to be extracted 
        images_to_extract.push_back( ind );
        // frame shifts to be applied to extracted images to generate weighted running averages
        frame_shifts.push_back( { matrix[14], matrix[15] } );

        // if this is not the frame for refinement, DO NOT process the metadata as followed
        // but we need to get the index in stack correct
        /**
        if ( cnfdnc != frame_index ) {
            ind_in_stack++;
            continue;
        }*/

        // frealignx metric
        //sscanf(
        //        line,
        //        "%d %lf %lf %lf %lf %lf %d %d %f %f %f %f %f %d %f %f %f %d %lf %lf %d %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        //        &ind, &psi, &theta, &phi, &shx, &shy, &mag, &film, &df1, &df2, &angast, &ppindex, &occ, &logp, &sigma, &score, &change, &ptlind, &tiltan,
        //        &dosexx, &scanor, &cnfdnc, &ptlccx, &tiltaxisan, &normal[0], &normal[1], &normal[2], &matrix[0],
        //        &matrix[1], &matrix[2], &matrix[3], &matrix[4], &matrix[5], &matrix[6], &matrix[7], &matrix[8],
        //        &matrix[9], &matrix[10], &matrix[11], &matrix[12], &matrix[13], &matrix[14], &matrix[15],
        //		  &ppsi, &ptheta, &pphi );
        

        // Create the particle projection, add it to the micrograph and to
        // the tilt series. Also create the micrograph when it is
        // necessary.
        eulerAngles.SetAngles( psi, theta, phi );
        particleProjectionShift.SetPosition( shx, shy );
        particleIndex = ptlind;
        // assume micrograph X, Y shifts are stored in matrix 12 and 13 repectively
        micrographTiltAxisShift.SetPosition(  matrix[12], matrix[13]);
        micrographTiltAxisAngle = tiltaxisan;
        micrographNumberID = scanor;
        tiltSeriesIndex = film;

        frameShift.SetPosition( 0.0, 0.0 );
        projectionFrameShift.SetPosition( matrix[14], matrix[15] );

        CEulerAngles initialEulerAngles( ppsi, ptheta, pphi );

        while ( pTiltSeries.size() < tiltSeriesIndex + 1 ) {
            CTiltSeries series;
            pTiltSeries.push_back( series );
        }
        pTiltSeries[tiltSeriesIndex].SetIndex( tiltSeriesIndex );

        pPP = new CParticleProjection( ind, eulerAngles, particleProjectionShift, projectionFrameShift, mag, tiltSeriesIndex, df1, df2, angast,
                                        ppindex, occ, logp, sigma, score, change, particleIndex, tiltan, micrographTiltAxisAngle, dosexx, micrographNumberID, (int) cnfdnc, ptlccx, micrographTiltAxisShift, ind_in_stack );
        
        // If the micrograph where belongs the particle projection is not
        // present in the tilt series, create it and add it to the tilt
        // series.
        if (!pTiltSeries[ tiltSeriesIndex ].HasMicrograph( micrographNumberID, (int) cnfdnc )) {
            pMicrograph = new CMicrograph( tiltan, micrographTiltAxisAngle, micrographTiltAxisShift, frameShift, micrographNumberID, (int) cnfdnc,
                                        tiltSeriesIndex, 0.0, 0.0, 0.0 ); 
            pTiltSeries[ tiltSeriesIndex ].AddMicrograph( pMicrograph );
        }
        // If the particle is not present in the tilt series, create and
        // add it. The orientation of the particle is (0,0,0) around the
        // refined normal computed by tomography. This normal refined is
        // read in the .parx.
        if (!pTiltSeries[tiltSeriesIndex].HasParticle( particleIndex, (int) cnfdnc )) {
            pParticle = new CParticle( initialEulerAngles, position3D, particleIndex, (int) cnfdnc, pTiltSeries[tiltSeriesIndex].mTiltSeriesIndex, normal,
                                    matrix, occ );
            pTiltSeries[tiltSeriesIndex].AddParticle( pParticle );
        }

        // update the max frame index in this dataset 
        pTiltSeries[tiltSeriesIndex].SetMaxFrameIndex( (int) cnfdnc );

        // Add the particle projection to the micrograph in the tilt series.
        counter = pTiltSeries[tiltSeriesIndex].AddParticleProjection( pPP );
        if (counter != -1) {
            params->mTotalNumberOfParticlesProjections++;
        } else {
            printf( "[ERROR] AddParticleProjection: a particle projection was added but its counter is -1.\n" );
            exit( EXIT_FAILURE );
        } // end if else

        ind_in_stack ++;  
        } // end if
    } // end while
    fclose( file );
    
    // calculating the weights
    double maxMeanScore = 0.0;
    const int weight_power = 1;
    double total = 0.0; 
    for ( int i = 0; i < weights.size(); i++) {
        weights[i] = pow( (weights[i] / countPerFrame[i]), weight_power);
        total += weights[i];
    }
    //double max_weight = *max_element(weights.begin(), weights.end());
    //double min_weight = *min_element(weights.begin(), weights.end());

    for ( int i = 0; i < weights.size(); i++ ) {
        //weights[i] = ( weights[i] - min_weight ) / ( max_weight - min_weight );
        weights[i] /= total;
    }
    for (it = priors_average.begin(); it != priors_average.end(); it++) {
        // printf("SCANORD %d:\n", it->first);
        for (int i = 0; i < NUM_COL_REFINE3D; i++) {
            priors_average.at(it->first)[i] /= countPerTilt.at(it->first); 
            priors_variance.at(it->first)[i] /= countPerTilt.at(it->first);
            priors_variance.at(it->first)[i] -= powf(priors_average.at(it->first)[i],2);
            // printf("Col = %d, avg = %f, variance = %f\n", i, priors_average.at(it->first)[i], priors_variance.at(it->first)[i]);
        }
    }
    return params->mTotalNumberOfParticlesProjections;
}

/** Write the extended .parx file. */
int WriteParFile( const flrCostFunction::ParametersType & parameters,
                  LoopModeType mode,
                  int objectIndex,
                  int frameIndex,
                  int minScanOrderToUseForRefinement,
                  int maxScanOrderToUseForRefinement,
                  bool frame_refine,
                  ReconstructedVolume* ref_3d,
                  CImage* image_set, 
                  ImageProjectionComparison** comparison_objects,
                  bool useCistemOptimizor,
                  bool evaluateAllProjections,
                  map<int, float*> priors_average, 
                  map<int, float*> priors_variance ) {

    int ind, film, ptlind, scanor, mag, logp, ind_in_stack;
    EulerAngleType psi, theta, phi, tiltan;
    CoordinateType shx, shy;
    
    float df1, df2, angast, ppindex, occ, sigma, score, change, dosexx, cnfdnc, ptlccx;
    float avgScore_before = 0.0, avgScore_after = 0.0;
    float fx = 0.0, fy = 0.0, mx = 0.0, my = 0.0;
    CEulerAngles orientation;
    CPosition2D shift, micrographShift, frameShift;

    CEulerAngles eA;
    EulerAngleType newEulerAngles[3];
    CoordinateType vshift[2];
    EulerAngleType tiltAngle;
    EulerAngleType tiltAxisAngle;
    EulerAngleType normal[3];
    EulerAngleType matrix[16];
    int mindex;
    int pindex;
    int PPCounter;

    int prior_tilt = 0;
    
    int numberOfParticleProjections, numberOfLocalProjections;
    
    CTiltSeries* pTiltSeries;
    CParticleProjection* pPP;
    CMicrograph* pMicrograph;
    CParticle* pParticle;
    vector<CParticleProjection*> pParticleProjectionSet;
    
    // Get the tilt series
    pTiltSeries = CContainer::Instance()->mTiltSeriesPointer;

    vector<int> positionsList; 
    int projCounter;

    switch (mode) {
        case OnMicrographRotation: 
        case OnMicrographTranslation: 
        case OnMicrographDefocus: 
        case OnMicrograph: {
            prior_tilt = objectIndex;
            if ( !frame_refine ) {
                // we just wanna refine micrograph i frame 0
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByMicrograph( objectIndex, frameIndex, positionsList );
            }
            else {    
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByMicrograph( objectIndex, positionsList );
            }
            //pMicrograph = pTiltSeries->mMicrograph.at( pTiltSeries->GetMicrographPositionByIndex( objectIndex, 0 ) );
            //numberOfParticleProjections = pMicrograph->mParticleProjectionSet.size();
            //pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
            break;
        }
        case OnParticleRotation:
        case OnParticleTranslation: 
        case OnParticle: {
            prior_tilt = priors_average.begin()->first;
            if ( !frame_refine ) {
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByParticle( objectIndex, 0, positionsList );
            }
            else {
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByParticle( objectIndex, positionsList );
            }
            //pParticle = pTiltSeries->mParticle.at( pTiltSeries->GetParticlePositionByIndex( objectIndex, 0 ) );
            //numberOfParticleProjections = pParticle->mParticleProjectionSet.size();
            //pParticleProjectionSet = pParticle->mParticleProjectionSet;
            break;
        }
        default:
        break;
    }
    if ( numberOfParticleProjections > 0 ) {
        
        float** p = new float*[numberOfParticleProjections];
        CPosition2D* shifts_before_refine = new CPosition2D[numberOfParticleProjections];
        
        projCounter = 0;
        
        for ( int pos = 0; pos < positionsList.size(); pos++ ) {
            
            switch (mode) {
                
                case OnMicrographRotation: 
                case OnMicrographDefocus:
                case OnMicrographTranslation: 
                case OnMicrograph: {
                    pMicrograph = pTiltSeries->mMicrograph.at(positionsList[pos]);
                    pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
                    numberOfLocalProjections = pMicrograph->mParticleProjectionSet.size(); 
                    break;
                } 

                case OnParticleRotation: 
                case OnParticleTranslation: 
                case OnParticle: {
                    pParticle = pTiltSeries->mParticle.at(positionsList[pos]);
                    pParticleProjectionSet = pParticle->mParticleProjectionSet;
                    numberOfLocalProjections = pParticle->mParticleProjectionSet.size();
                    break;
                }

                default: {
                    printf( "[ERROR] ComputeAnglesInSPManifold: no correct mode selected (mode:%d)", mode );
                    exit( EXIT_FAILURE );
                }
                

            }        
            for ( PPCounter = 0; PPCounter < numberOfLocalProjections; PPCounter++ ) {
                
                switch (mode) {
                    case OnMicrographRotation: 
                    case OnMicrograph: {

                        // Get the particle projection
                        pPP = pParticleProjectionSet.at( PPCounter );

                        // Get the orientation of the particle.
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        // Get the normal and matrix refinement of the particle.
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );

                        // Get the angles of the micrograph.
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle() + parameters[2];
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle() + parameters[3];

                        break;
                    }
                    case OnMicrographDefocus: 
                    case OnMicrographTranslation: {
                        
                        pPP = pParticleProjectionSet.at( PPCounter );

                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        break;
                    }


                    case OnParticleRotation: {

                        pPP = pParticleProjectionSet.at( PPCounter );

                        // Get the angles of the micrograph.
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();

                        // Get the orientation of the particle.
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA.SetTheta( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetTheta() + parameters[3] );
                        eA.SetPsi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPsi() + parameters[4] );
                        eA.SetPhi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPhi() + parameters[5] );

                        // Get the normal and matrix refinement of the particle.
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );

                        break;
                    } 
                    case OnParticleTranslation: {
                        
                        pPP = pParticleProjectionSet.at( PPCounter );

                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        matrix[3] = matrix[3] + parameters[0]; // pX 
                        matrix[7] = matrix[7] + parameters[1]; // pY
                        matrix[11] = matrix[11] + parameters[2]; // pZ
                        break;
                    }
                    
                    case OnParticle: {
                        pPP = pParticleProjectionSet.at( PPCounter );

                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        
                        eA.SetTheta( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetTheta() + parameters[3] );
                        eA.SetPsi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPsi() + parameters[4] );
                        eA.SetPhi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPhi() + parameters[5] );

                        matrix[3] = matrix[3] + parameters[0]; // pX 
                        matrix[7] = matrix[7] + parameters[1]; // pY
                        matrix[11] = matrix[11] + parameters[2]; // pZ
                        break;
                    }

                    

                    default: {
                        printf( "[ERROR] ComputeAnglesInSPManifold: no correct mode selected (mode:%d)", mode );
                        exit( EXIT_FAILURE );
                    }
                } 

                EulerAngleType ppsi, ptheta, pphi;
                ppsi = eA.GetPsi();
                ptheta = eA.GetTheta();
                pphi = eA.GetPhi();
                    
                // use local microgrpah geometry instead of global ones
                if ( frame_refine && mode == OnMicrographTranslation ) {
                    micrographShift = pPP->mMicrographTiltAxisShift;
                    tiltAngle = pPP->mMicrographTiltAngle;
                    tiltAxisAngle = pPP->mMicrographTiltAxisAngle;
                }   
                else {
                    micrographShift = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisShift();
                }

                // Map the angles for the particle projection
                SPAEulerAngles( eA.GetPsi(), eA.GetTheta(), eA.GetPhi(), tiltAngle, tiltAxisAngle, normal, matrix, newEulerAngles, vshift );

                // Save the angles into the data structure
                eA.SetAngles( newEulerAngles[0], newEulerAngles[1], newEulerAngles[2] );
                    
                // when refining defocus/astig, do NOT want to change current rot/shifts
                if ( mode != OnMicrographDefocus ) {
                    pPP->SetOrientation( eA );
                    pPP->SetShift( vshift[0], vshift[1] );
                }

                // Write the line
                ind = pPP->mIndexInStack; // ind
                orientation = pPP->mOrientation; // psi, theta, phi
                psi = orientation.GetPsi();
                theta = orientation.GetTheta();
                phi = orientation.GetPhi();
                shift = pPP->mShift; // shx, shy
                frameShift = pPP->mFrameShift; 

                // add micrograph shifts if OnMicrographShift mode
                if ( mode == OnMicrographTranslation || mode == OnMicrograph ) {
                    shx = shift.mCoordX + micrographShift.mCoordX + parameters[0]; //+ frameShift.mCoordX;
                    shy = shift.mCoordY + micrographShift.mCoordY + parameters[1]; // + frameShift.mCoordY;
                    if ( image_set->IsRaw() ) {
                        // frame shift should be already applied to (weighted) averages
                        // NOTE: frame shifts should be zero if no frames at all
                        shx += frameShift.mCoordX;
                        shy += frameShift.mCoordY;
                    }
                }
                else if ( mode == OnMicrographDefocus ) {
                    shx = shift.mCoordX;
                    shy = shift.mCoordY;
                }
                else {
                    shx = shift.mCoordX + micrographShift.mCoordX; // + frameShift.mCoordX;
                    shy = shift.mCoordY + micrographShift.mCoordY; // + frameShift.mCoordY;
                    if ( image_set->IsRaw() ) {
                        shx += frameShift.mCoordX;
                        shy += frameShift.mCoordY;
                    }
                }
                
                mag = pPP->mMagnification; // mag
                film = pPP->mTiltSeriesIndex; // film
                    
                // add defocus offset and astig if OnMicrographDefocus mode
                if ( mode == OnMicrographDefocus ) {
                    df1 = pPP->mDefocus01 + parameters[3];
                    df2 = pPP->mDefocus02 + parameters[4];
                    angast = pPP->mAstigmatism + parameters[5];
                }
                else {
                    df1 = pPP->mDefocus01;
                    df2 = pPP->mDefocus02;
                    angast = pPP->mAstigmatism;
                }
        
                ppindex = pPP->mPpindex; // ppindex
                occ = pPP->mOcc; // occ
                logp = pPP->mLogp; // logp
                sigma = pPP->mSigma; // sigma
                score = pPP->mScore; // score
                change = pPP->mChange; // change
                ptlind = pPP->mParticleIndex; //ptlind
                tiltan = pPP->mMicrographTiltAxisAngle; // tiltan;
                dosexx = pPP->mDose; // dosexx
                scanor = pPP->mMicrographsIndex; // scanor
                cnfdnc = pPP->mFrameIndex; // cnfdnc
                ptlccx = pPP->mCCX; // ptlccx
                ind_in_stack = pPP->mIndexInLocalStack;

                p[projCounter] = new float[NUM_COL_REFINE3D];
                p[projCounter][0] = ind;
                p[projCounter][1] =  psi; p[projCounter][2] = theta; p[projCounter][3] = phi;  p[projCounter][4] = shx; p[projCounter][5] = shy;
                p[projCounter][6] = mag; p[projCounter][7] = film;  p[projCounter][8] = df1; p[projCounter][9] = df2; p[projCounter][10] = angast;  
                p[projCounter][11] = ppindex; p[projCounter][12] = occ; p[projCounter][13] = logp; p[projCounter][14] = sigma; p[projCounter][15] = score; p[projCounter][16] = change; 
                
                if (evaluateAllProjections || (scanor >=  minScanOrderToUseForRefinement) && ((maxScanOrderToUseForRefinement == -1) || (scanor <= maxScanOrderToUseForRefinement))){
                    p[projCounter][17] = ind_in_stack;
                }
                else { p[projCounter][17] = -1; }
                
                shifts_before_refine[projCounter].SetPosition( shx, shy );
                avgScore_before += score;
                projCounter ++;
            }
        }

        avgScore_before /= projCounter;
        RunRefine3d( p, 
                    image_set->imageStack, 
                    ref_3d, 
                    comparison_objects,
                    useCistemOptimizor, 
                    evaluateAllProjections,
                    image_set->IsRaw(), 
                    frame_refine, 
                    refine3dParams.usePriors,
                    priors_average.at(prior_tilt), 
                    priors_variance.at(prior_tilt), 
                    numberOfParticleProjections, 
                    refine3dParams.symmetry, 
                    refine3dParams.scopePixel, 
                    refine3dParams.boxsize, 
                    refine3dParams.beamEnergy, 
                    refine3dParams.sphericalAbberation, 
                    refine3dParams.amplitutdeContrast, 
                    refine3dParams.molecularMass, 
                    refine3dParams.outerMask, 
                    refine3dParams.lowResLimit, 
                    refine3dParams.highResLimit, 
                    refine3dParams.resLimitSignCC, 
                    refine3dParams.ResLimitClassification );
        projCounter = 0;
        for ( int pos = 0; pos < positionsList.size(); pos++ ) {
            switch (mode) {              
                case OnMicrographRotation: 
                case OnMicrographDefocus:
                case OnMicrographTranslation: 
                case OnMicrograph: {
                    pMicrograph = pTiltSeries->mMicrograph.at(positionsList[pos]);
                    pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
                    numberOfLocalProjections = pMicrograph->mParticleProjectionSet.size(); 
                    break;
                } 
                case OnParticleRotation: 
                case OnParticleTranslation: 
                case OnParticle : {
                    pParticle = pTiltSeries->mParticle.at(positionsList[pos]);
                    pParticleProjectionSet = pParticle->mParticleProjectionSet;
                    numberOfLocalProjections = pParticle->mParticleProjectionSet.size();
                    break;
                }
            }

            for ( PPCounter = 0; PPCounter < numberOfLocalProjections; PPCounter++ ) {
                pPP = pParticleProjectionSet.at( PPCounter );
                
                // only if we use cistem optimizor (only call CostFunction->evaluate once)
                if (useCistemOptimizor && mode == 3) {
                    if (frame_refine) {
                        pPP->mScore = p[projCounter][15];
                        fx = pPP->mFrameShift.mCoordX + ( p[projCounter][4] - shifts_before_refine[projCounter].mCoordX );
                        fy = pPP->mFrameShift.mCoordY + ( p[projCounter][5] - shifts_before_refine[projCounter].mCoordY );
                        pPP->mFrameShift.SetPosition( fx, fy ); 
                        // printf("%d %f\n", projCounter, p[projCounter][15] );
                    } 
                    else {
                        mx = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisShift().mCoordX + ( p[projCounter][4] - shifts_before_refine[projCounter].mCoordX ); 
                        my = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisShift().mCoordY + ( p[projCounter][5] - shifts_before_refine[projCounter].mCoordY );
                        micrographShift.SetPosition( mx, my );
                        pTiltSeries->mMicrograph.at( mindex )->SetTiltAxisShift( micrographShift );
                        pTiltSeries->mParticle.at( pindex )->SetOrientation( -p[projCounter][2], -p[projCounter][1], -p[projCounter][3] );
                        // printf("%d %f\n", projCounter, p[projCounter][15] );
                    }
                }
                else { 
                    pPP->mScore = p[projCounter][15];
                    avgScore_after += pPP->mScore; 
                }
                avgScore_after += pPP->mScore;
                projCounter++;
            }
        }
        avgScore_after /= projCounter;
        if (useCistemOptimizor) {
            // print only if we use refine3d internal optimizor
            // printf("\n### Particle frame refinement at tilt %d ###\n", objectIndex);
            // printf("    Start from average score : %f \n", avgScore_before);
            // printf("            To average score : %f \n", avgScore_after);
        }
            
        // delete array to avoid mem leak 
        for ( int i = 0; i < numberOfParticleProjections; i++ ) {
            delete[] p[i];
        }
        delete[] p;
        delete[] shifts_before_refine;
    }
    return numberOfParticleProjections;
}


/** Write the extended .parx file. */
int WriteMultipleParFileNew( vector< CTiltSeries > & dataset, char filename[200], char outfilename[200], LoopModeType mode, bool frame_refine, bool extract_particle_frames, int minScanOrderUsedForReconstruction, int maxScanOrderUsedForReconstruction ) {
  FILE *file;
  int ind, film, ptlind, scanor, mag, logp;
  EulerAngleType psi, theta, phi, tiltan;
  CoordinateType shx, shy;
  
  float df1, df2, angast, occ, sigma, score, change, dosexx, cnfdnc, ptlccx;
  CEulerAngles orientation;
  CPosition2D shift, micrographShift, frameShift;

  CEulerAngles eA;
  EulerAngleType newEulerAngles[3];
  CoordinateType vshift[2];
  EulerAngleType tiltAngle;
  EulerAngleType tiltAxisAngle;
  EulerAngleType normal[3];
  EulerAngleType matrix[16];
  int mindex;
  int findex;
  int pindex;
  unsigned int PPCounter;
  unsigned int numberOfParticleProjections;

  CTiltSeries* pTiltSeries;
  CParticleProjection* pPP;
  CMicrograph* pMicrograph;
  CParticle* pParticle;

  // Get the tilt series
  pTiltSeries = CContainer::Instance()->mTiltSeriesPointer;

  std :: ifstream src( filename );
  stringstream tiltseries;
  src >> tiltseries.rdbuf();

  // Open the .parx file to write
  file = fopen( outfilename, "w" );
  if (!file) {
	  printf( "[ERROR] WriteMultipleParFileNew: Could not open %s for writing.\n", outfilename );
	  return -1;
  }

  // Write captions
  //fprintf( file, "C FREALIGN parameter file\n");
  //fprintf( file, "C Created by WriteMultipleParFileNew.cxx\n");
  fprintf( file, "C     1       2       3       4         5         6       7     8        9       10      11      12        13         14      15      16       17      18        19       20       21       22        23        24        25        26        27        28        29        30        31        32        33        34        35        36        37        38        39        40        41        42        43        44        45\n");
  fprintf( file, "C    NO     PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1      DF2  ANGAST     OCC      LOGP      SIGMA   SCORE  CHANGE   PTLIND   TILTAN   DOSEXX   SCANOR   CNFDNC   PTLCCX      AXIS     NORM0     NORM1     NORM2  MATRIX00  MATRIX01  MATRIX02  MATRIX03  MATRIX04  MATRIX05  MATRIX06  MATRIX07  MATRIX08  MATRIX09  MATRIX10  MATRIX11  MATRIX12  MATRIX13  MATRIX14  MATRIX15      PPSI    PTHETA      PPHI\n");

  
  
  // write one file for each tilt series
  for ( int i = 0; i < dataset.size(); i++ ){

	int index;
	char currentSeries[200];
	tiltseries >> index;
	tiltseries >> currentSeries;
		
	pTiltSeries = &dataset[i];

	// Determine the number of particle projections to write given the
	// LoopMode and the object.
	numberOfParticleProjections = pTiltSeries->mParticleProjectionSet.size();

	// Go over the set of particle projections mapping the CSP and USP
	// manifold and writting the parx file.
	for ( PPCounter = 0; PPCounter < numberOfParticleProjections; PPCounter++ ) {

		// Get the particle projection
		pPP = pTiltSeries->mParticleProjectionSet.at( PPCounter );

		// Get the angles of the micrograph.
		mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
        findex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), pPP->GetFrameIndex() );
		
        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
		tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();

		// Get the orientation, normal and matrix refinement of the particle.
		pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
		eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();
		pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
		pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );

		EulerAngleType ppsi, ptheta, pphi;
		ppsi = eA.GetPsi();
		ptheta = eA.GetTheta();
		pphi = eA.GetPhi();


        if ( frame_refine ) {
            tiltAngle = pPP->mMicrographTiltAngle;
            tiltAxisAngle = pPP->mMicrographTiltAxisAngle;
            micrographShift = pPP->mMicrographTiltAxisShift;
        }
        else {
            micrographShift = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisShift();
        }

		// Map the angles for the particle projection
		SPAEulerAngles( eA.GetPsi(), eA.GetTheta(), eA.GetPhi(), tiltAngle, tiltAxisAngle, normal, matrix, newEulerAngles, vshift );
        
		// Save the angles into the data structure
		eA.SetAngles( newEulerAngles[0], newEulerAngles[1], newEulerAngles[2] );
		
        if ( mode != OnMicrographDefocus && mode != OnSkip ){ 
            pPP->SetOrientation( eA );
            pPP->SetShift( vshift[0], vshift[1] );
        }

		// Write the line
		ind = pPP->mIndexInStack; // ind
		orientation = pPP->mOrientation; // psi, theta, phi
		psi = orientation.GetPsi();
		theta = orientation.GetTheta();
		phi = orientation.GetPhi();
		shift = pPP->mShift; // shx, shy
        frameShift = pPP->mFrameShift;
        
        if ( mode == OnMicrographDefocus || mode == OnSkip ) {
            shx = shift.mCoordX;
            shy = shift.mCoordY;
        }    
        else {
            shx = shift.mCoordX + micrographShift.mCoordX + frameShift.mCoordX; //pTiltSeries->mMicrograph.at( findex )->GetFrameShift().mCoordX;//+ frameShift.mCoordX;
            shy = shift.mCoordY + micrographShift.mCoordY + frameShift.mCoordY; // pTiltSeries->mMicrograph.at( findex )->GetFrameShift().mCoordY;// + frameShift.mCoordY;
        }
        
        matrix[12] = micrographShift.mCoordX;
        matrix[13] = micrographShift.mCoordY;

        mag = pPP->mMagnification; // mag
		film = pPP->mTiltSeriesIndex; // film
		df1 = pPP->mDefocus01 + pTiltSeries->mMicrograph.at( mindex )->GetTiltDefocus1Offset(); // df1
		df2 = pPP->mDefocus02 + pTiltSeries->mMicrograph.at( mindex )->GetTiltDefocus2Offset(); // df2
		angast = pPP->mAstigmatism + pTiltSeries->mMicrograph.at( mindex )->GetTiltAstigmatism(); // angast
        matrix[14] = pTiltSeries->mMicrograph.at( findex )->GetFrameShift().mCoordX + frameShift.mCoordX;
        matrix[15] = pTiltSeries->mMicrograph.at( findex )->GetFrameShift().mCoordY + frameShift.mCoordY; 

        occ = pPP->mOcc; // occ
        logp = pPP->mLogp; // LogP
		sigma = pPP->mSigma; // sigma
		score = (extract_particle_frames && !frame_refine) ? pPP->mPreviousScore : pPP->mScore; // score
		change = pPP->mChange; // change
		ptlind = pPP->mParticleIndex; //ptlind
		tiltan = pPP->mMicrographTiltAxisAngle; // tiltan;
		dosexx = pPP->mDose; // dosexx
		scanor = pPP->mMicrographsIndex; // scanor
		cnfdnc = pPP->mFrameIndex; // cnfdnc
		ptlccx = pPP->mCCX; // ptlccx

        
        // we are not truncate negative scores to zeros in objective function, but we do so here for the parfile
        score = (score <= 0.0) ? 0.01 : score;
        
		/**
         * Metirc new
         */
        
		fprintf( file,
		"%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f%9d%9.2f%9.2f%9d%9.2f%9.2f",
		ind, psi, theta, phi, shx, shy, mag, film, df1, df2, angast, occ, logp, sigma, score, change, ptlind, tiltAngle,
		dosexx, scanor, cnfdnc, ptlccx );
		fprintf( file,
			"%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
			tiltAxisAngle, normal[0], normal[1], normal[2], matrix[0],
			matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6], matrix[7], matrix[8],
			matrix[9], matrix[10], matrix[11], matrix[12], matrix[13], matrix[14], matrix[15],
			ppsi, ptheta, pphi	);
        
    
        /**
         * Metric frealignx
         */
        /**
        fprintf( file,
		"%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%8.2f%10d%11.4f%8.2f%8.2f%9d%9.2f%9.2f%9d%9.2f%9.2f",
		ind, psi, theta, phi, shx, shy, mag, film, df1, df2, angast, 0.0, occ, logp, sigma, score, change, ptlind, tiltAngle,
		dosexx, scanor, cnfdnc, ptlccx );
		fprintf( file,
			"%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
			tiltAxisAngle, normal[0], normal[1], normal[2], matrix[0],
			matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6], matrix[7], matrix[8],
			matrix[9], matrix[10], matrix[11], matrix[12], matrix[13], matrix[14], matrix[15],
			ppsi, ptheta, pphi	);
        */
	}

	}

	// Close the .par file
	fclose( file );

	src.close();
	
	return PPCounter;
}


void GetRefinedFrameShiftsAtTilt(vector< CTiltSeries > & dataset, int tilt_index, double ** refined_frame_shifts) {

    CTiltSeries* pTiltSeries;
    CParticleProjection* pPP;
    CMicrograph* pMicrograph;
    CParticle* pParticle;

    CPosition2D frameShift;

    vector<int> positionsList; 

    vector<CParticleProjection*> pParticleProjectionSet;
    
    int projCounter = 0;

    for ( int i = 0; i < dataset.size(); i++ ){
            
        pTiltSeries = &dataset[i];

        // int numberOfParticleProjections = pTiltSeries->mParticleProjectionSet.size();
        int numberOfParticleProjections = pTiltSeries->GetAllPositionsByMicrograph( tilt_index, positionsList );
        
        for ( int pos = 0; pos < positionsList.size(); pos++ ) {
            pMicrograph = pTiltSeries->mMicrograph.at(positionsList[pos]);
            pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
            int numberOfLocalProjections = pMicrograph->mParticleProjectionSet.size(); 
            
            for ( int PPCounter = 0; PPCounter < numberOfLocalProjections; PPCounter++ ) {
                pPP = pParticleProjectionSet.at( PPCounter );
                frameShift = pPP->mFrameShift;
                
                refined_frame_shifts[projCounter][0] = pPP->mIndexInLocalStack;
                refined_frame_shifts[projCounter][1] = frameShift.mCoordX;
                refined_frame_shifts[projCounter][2] = frameShift.mCoordY;
                // printf("%f %f %f \n", refined_frame_shifts[projCounter][0], refined_frame_shifts[projCounter][1], refined_frame_shifts[projCounter][2]);
                projCounter++;
            }
        }
    }
}



void GetRefinedFrameShifts(vector< CTiltSeries > & dataset, double ** refined_frame_shifts) {

    CTiltSeries* pTiltSeries;
    CParticleProjection* pPP;
    CMicrograph* pMicrograph;
    CParticle* pParticle;

    CPosition2D frameShift;

    vector<int> positionsList; 

    vector<CParticleProjection*> pParticleProjectionSet;
    
    for ( int i = 0; i < dataset.size(); i++ ){
            
        pTiltSeries = &dataset[i];

        int numberOfParticleProjections = pTiltSeries->mParticleProjectionSet.size();

            
        for ( int PPCounter = 0; PPCounter < numberOfParticleProjections; PPCounter++ ) {

            // Get the particle projection
            pPP = pTiltSeries->mParticleProjectionSet.at( PPCounter );

            frameShift = pPP->mFrameShift;
            
            refined_frame_shifts[PPCounter][0] = pPP->mIndexInLocalStack;
            refined_frame_shifts[PPCounter][1] = frameShift.mCoordX;
            refined_frame_shifts[PPCounter][2] = frameShift.mCoordY;
            
            //printf("%d %f %f\n", PPCounter, frameShift.mCoordX, frameShift.mCoordY);
        }
    }
}



Image** ExtractParticles( char mrc_file[1000], bool use_frames, bool stack_avg, char allboxes[1000], vector<int> images_to_extract, vector< pair<double, double> > frame_shifts, int write_to_disk, char output_stack[1000] ) {
    /** Read the coordinates from allboxes file, and extract particles from the mrc file(s)
    */

    char line[200];
    int coord_x, coord_y, idx_tilt, idx_frame=0;
    vector<int> tmp;
    vector< vector<int> > coordinates;
    FILE *file;
    
    if (images_to_extract.size() == 0) {
        printf("[ERROR] No images to be extracted.\n");
        exit(EXIT_FAILURE);
    }

    // First, read allboxes
    file = fopen(allboxes, "r");
    if (!file) {
        printf( "\n[ERROR] %s does not exist or cannot be read.\nExiting...\n", allboxes );
        exit(EXIT_FAILURE);
    }
    
    int num_particles = 0;
    int number_tilts = 1;
    int number_frames = 1;
    
    int curr_line = 1;
    int curr_idx_image_list = 0;
    

    while (fgets(line, sizeof(line), file) != NULL) {
        
        if (use_frames) {sscanf(line, "%d %d %d %d", &coord_x, &coord_y, &idx_tilt, &idx_frame);}
        else {sscanf(line, "%d %d %d", &coord_x, &coord_y, &idx_tilt);}
         
        // only extract images if they are in the list
        if (curr_line == images_to_extract[curr_idx_image_list]) {
            tmp.push_back(coord_x);
            tmp.push_back(coord_y);
            tmp.push_back(idx_tilt);
            tmp.push_back(idx_frame);
            
            coordinates.push_back(tmp);
            tmp.clear();
            
            num_particles ++;
            curr_idx_image_list ++;
           
            //if (idx_tilt+1 > number_tilts) {number_tilts = idx_tilt+1;}
            //if (idx_frame+1 > number_frames) {number_frames = idx_frame+1;} 
        } 
        if (idx_tilt+1 > number_tilts) {number_tilts = idx_tilt+1;}
        if (idx_frame+1 > number_frames) {number_frames = idx_frame+1;}
        curr_line ++;

        if ( curr_idx_image_list >= images_to_extract.size() ) {
            break;
        }
    }
    fclose(file);
    
    // convert to int array before passing it to extern C function 
    double coords[num_particles][6];

    for ( int i = 0; i < num_particles; i ++ ) {
        coords[i][0] = (double) coordinates[i][0];
        coords[i][1] = (double) coordinates[i][1];
        coords[i][2] = (double) coordinates[i][2];
        if (use_frames) {coords[i][3] = (double) coordinates[i][3];} 
        else {coords[i][3] = 0.0;}
        coords[i][4] = frame_shifts.at(i).first;
        coords[i][5] = frame_shifts.at(i).second;
        // printf("%f %f %f %f\n", coords[i][0], coords[i][1], coords[i][2], coords[i][3]);
    }
    // refine3dParams.numberFrames = number_frames;
    // Second, extract particles from tilt-series (extern C wrapper function)
    return ExtractParticlesFromMRCs( mrc_file, use_frames, stack_avg, coords, num_particles, number_tilts, number_frames, refine3dParams.particleBin, refine3dParams.boxsize, refine3dParams.outerMask, refine3dParams.scopePixel, write_to_disk, output_stack );
}   

void calcFrameWeights( vector< vector<double> > & frame_weights, int num_frames, double weight_width ) {
    
    double weight, sum;
    int counter;
    vector<double> curr_frame;
    
    for ( int i = 0; i < num_frames; i++ ) {
        counter = 0;
        sum = 0.0;
        curr_frame.clear();
        for ( int j = 0; j < num_frames; j++ ) {
            weight = exp( -pow( counter - i, 2 ) / weight_width ); 
            sum += weight;
            counter ++;
            curr_frame.push_back(weight);
        }
        for ( int j = 0; j < num_frames; j++ ) {
            curr_frame[j] = curr_frame[j] / ( sum / num_frames ) / num_frames;
        }
        frame_weights.push_back( curr_frame );
    }   
}

ImageProjectionComparison** PrepareComparisonObjects(LoopModeType mode,
                                                    int objectIndex,
                                                    int frameIndex,
                                                    bool frame_refine,
                                                    ReconstructedVolume* ref_3d,
                                                    CImage* image_set,
                                                    int & number_projections_per_refinement, 
                                                    map<int, float*> & priors_average, 
                                                    map<int, float*> & priors_variance
                                                    ) {
    
    float parameters[6];
    for (int i = 0; i < 6; i++) {
        parameters[i] = 0.0;
    }

    int ind, film, ptlind, scanor, mag, logp, ind_in_stack;
    EulerAngleType psi, theta, phi, tiltan;
    CoordinateType shx, shy;

    float df1, df2, angast, ppindex, occ, sigma, score, change, dosexx, cnfdnc, ptlccx;
    float avgScore_before = 0.0, avgScore_after = 0.0;
    float fx = 0.0, fy = 0.0, mx = 0.0, my = 0.0;
    CEulerAngles orientation;
    CPosition2D shift, micrographShift, frameShift;

    CEulerAngles eA;
    EulerAngleType newEulerAngles[3];
    CoordinateType vshift[2];
    EulerAngleType tiltAngle;
    EulerAngleType tiltAxisAngle;
    EulerAngleType normal[3];
    EulerAngleType matrix[16];
    int mindex;
    int pindex;
    int PPCounter;

    int prior_tilt = 0;
    
    int numberOfParticleProjections, numberOfLocalProjections;

    CTiltSeries* pTiltSeries;
    CParticleProjection* pPP;
    CMicrograph* pMicrograph;
    CParticle* pParticle;
    vector<CParticleProjection*> pParticleProjectionSet;
    
    // Get the tilt series
    pTiltSeries = CContainer::Instance()->mTiltSeriesPointer;

    vector<int> positionsList; 
    int projCounter;

    switch (mode) {
        case OnMicrographRotation: 
        case OnMicrographTranslation: 
        case OnMicrographDefocus: 
        case OnMicrograph: {
            prior_tilt = objectIndex;
            // FIXME
            if ( !frame_refine ) {
                // we just wanna refine micrograph i frame 0
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByMicrograph( objectIndex, frameIndex, positionsList );
            }
            else {    
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByMicrograph( objectIndex, positionsList );
            }
            break;
        }
        case OnParticleRotation:
        case OnParticleTranslation: 
        case OnParticle: {
            // FIXME
            if ( !frame_refine ) {
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByParticle( objectIndex, 0, positionsList );
            }
            else {
                numberOfParticleProjections = pTiltSeries->GetAllPositionsByParticle( objectIndex, positionsList );
            }
            break;
        }
        default:
        break;
    }

    if ( numberOfParticleProjections > 0 ) {
        
        float** p = new float*[numberOfParticleProjections];
        ImageProjectionComparison** comparison_objects = new ImageProjectionComparison*[numberOfParticleProjections];

        projCounter = 0;
        
        for ( int pos = 0; pos < positionsList.size(); pos++ ) {
            
            switch (mode) {
                
                case OnMicrographRotation: 
                case OnMicrographDefocus:
                case OnMicrographTranslation: 
                case OnMicrograph: {
                    pMicrograph = pTiltSeries->mMicrograph.at(positionsList[pos]);
                    pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
                    numberOfLocalProjections = pMicrograph->mParticleProjectionSet.size(); 
                    break;
                } 

                case OnParticleRotation: 
                case OnParticleTranslation: 
                case OnParticle: {
                    pParticle = pTiltSeries->mParticle.at(positionsList[pos]);
                    pParticleProjectionSet = pParticle->mParticleProjectionSet;
                    numberOfLocalProjections = pParticle->mParticleProjectionSet.size();
                    break;
                }

                default: {
                    printf( "[ERROR] ComputeAnglesInSPManifold: no correct mode selected (mode:%d)", mode );
                    exit( EXIT_FAILURE );
                }
                

            }        
            for ( PPCounter = 0; PPCounter < numberOfLocalProjections; PPCounter++ ) {
                
                switch (mode) {
                    case OnMicrographRotation: 
                    case OnMicrograph: {

                        // Get the particle projection
                        pPP = pParticleProjectionSet.at( PPCounter );

                        // Get the orientation of the particle.
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        // Get the normal and matrix refinement of the particle.
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );

                        // Get the angles of the micrograph.
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle() + parameters[2];
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle() + parameters[3];

                        break;
                    }
                    case OnMicrographDefocus: 
                    case OnMicrographTranslation: {
                        
                        pPP = pParticleProjectionSet.at( PPCounter );

                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        break;
                    }


                    case OnParticleRotation: {

                        pPP = pParticleProjectionSet.at( PPCounter );

                        // Get the angles of the micrograph.
                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();

                        // Get the orientation of the particle.
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA.SetTheta( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetTheta() + parameters[3] );
                        eA.SetPsi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPsi() + parameters[4] );
                        eA.SetPhi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPhi() + parameters[5] );

                        // Get the normal and matrix refinement of the particle.
                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );

                        break;
                    } 
                    case OnParticleTranslation: {
                        
                        pPP = pParticleProjectionSet.at( PPCounter );

                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        matrix[3] = matrix[3] + parameters[0]; // pX 
                        matrix[7] = matrix[7] + parameters[1]; // pY
                        matrix[11] = matrix[11] + parameters[2]; // pZ
                        break;
                    }
                    
                    case OnParticle: {
                        
                        pPP = pParticleProjectionSet.at( PPCounter );

                        mindex = pTiltSeries->GetMicrographPositionByIndex( pPP->GetMicrographIndex(), 0 );
                        tiltAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAngle();
                        tiltAxisAngle = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisAngle();
                        
                        pindex = pTiltSeries->GetParticlePositionByIndex( pPP->GetParticleIndex(), 0 );
                        eA = pTiltSeries->mParticle.at( pindex )->GetOrientation();

                        pTiltSeries->mParticle.at( pindex )->GetNormal( normal );
                        
                        pTiltSeries->mParticle.at( pindex )->GetMatrix( matrix );
                        
                        eA.SetTheta( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetTheta() + parameters[3] );
                        eA.SetPsi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPsi() + parameters[4] );
                        eA.SetPhi( pTiltSeries->mParticle.at( pindex )->GetOrientation().GetPhi() + parameters[5] );

                        matrix[3] = matrix[3] + parameters[0]; // pX 
                        matrix[7] = matrix[7] + parameters[1]; // pY
                        matrix[11] = matrix[11] + parameters[2]; // pZ
                        break;
                    }

                    

                    default: {
                        printf( "[ERROR] ComputeAnglesInSPManifold: no correct mode selected (mode:%d)", mode );
                        exit( EXIT_FAILURE );
                    }
                } 

                EulerAngleType ppsi, ptheta, pphi;
                ppsi = eA.GetPsi();
                ptheta = eA.GetTheta();
                pphi = eA.GetPhi();
                    
                // use local microgrpah geometry instead of global ones
                if ( frame_refine && mode == OnMicrographTranslation ) {
                    micrographShift = pPP->mMicrographTiltAxisShift;
                    tiltAngle = pPP->mMicrographTiltAngle;
                    tiltAxisAngle = pPP->mMicrographTiltAxisAngle;
                }   
                else {
                    micrographShift = pTiltSeries->mMicrograph.at( mindex )->GetTiltAxisShift();
                }

                // Map the angles for the particle projection
                SPAEulerAngles( eA.GetPsi(), eA.GetTheta(), eA.GetPhi(), tiltAngle, tiltAxisAngle, normal, matrix, newEulerAngles, vshift );

                // Save the angles into the data structure
                eA.SetAngles( newEulerAngles[0], newEulerAngles[1], newEulerAngles[2] );
                    
                // when refining defocus/astig, do NOT want to change current rot/shifts
                if ( mode != OnMicrographDefocus ) {
                    pPP->SetOrientation( eA );
                    pPP->SetShift( vshift[0], vshift[1] );
                }

                // Write the line
                ind = pPP->mIndexInStack; // ind
                orientation = pPP->mOrientation; // psi, theta, phi
                psi = orientation.GetPsi();
                theta = orientation.GetTheta();
                phi = orientation.GetPhi();
                shift = pPP->mShift; // shx, shy
                frameShift = pPP->mFrameShift; 

                // add micrograph shifts if OnMicrographShift mode
                if ( mode == OnMicrographTranslation || mode == OnMicrograph ) {
                    shx = shift.mCoordX + micrographShift.mCoordX + parameters[0]; //+ frameShift.mCoordX;
                    shy = shift.mCoordY + micrographShift.mCoordY + parameters[1]; // + frameShift.mCoordY;
                    if ( image_set->IsRaw() ) {
                        // frame shift should be already applied to (weighted) averages
                        // NOTE: frame shifts should be zero if no frames at all
                        shx += frameShift.mCoordX;
                        shy += frameShift.mCoordY;
                    }
                }
                else if ( mode == OnMicrographDefocus ) {
                    shx = shift.mCoordX;
                    shy = shift.mCoordY;
                }
                else {
                    shx = shift.mCoordX + micrographShift.mCoordX + frameShift.mCoordX;
                    shy = shift.mCoordY + micrographShift.mCoordY + frameShift.mCoordY;
                }
                
                mag = pPP->mMagnification; // mag
                film = pPP->mTiltSeriesIndex; // film
                    
                // add defocus offset and astig if OnMicrographDefocus mode
                if ( mode == OnMicrographDefocus ) {
                    df1 = pPP->mDefocus01 + parameters[3];
                    df2 = pPP->mDefocus02 + parameters[4];
                    angast = pPP->mAstigmatism + parameters[5];
                }
                else {
                    df1 = pPP->mDefocus01;
                    df2 = pPP->mDefocus02;
                    angast = pPP->mAstigmatism;
                }
        
                ppindex = pPP->mPpindex; // ppindex
                occ = pPP->mOcc; // occ
                logp = pPP->mLogp; // logp
                sigma = pPP->mSigma; // sigma
                score = pPP->mScore; // score
                change = pPP->mChange; // change
                ptlind = pPP->mParticleIndex; //ptlind
                tiltan = pPP->mMicrographTiltAxisAngle; // tiltan;
                dosexx = pPP->mDose; // dosexx
                scanor = pPP->mMicrographsIndex; // scanor
                cnfdnc = pPP->mFrameIndex; // cnfdnc
                ptlccx = pPP->mCCX; // ptlccx
                ind_in_stack = pPP->mIndexInLocalStack;

                p[projCounter] = new float[NUM_COL_REFINE3D];
                p[projCounter][0] = ind;
                p[projCounter][1] =  psi; p[projCounter][2] = theta; p[projCounter][3] = phi;  p[projCounter][4] = shx; p[projCounter][5] = shy;
                p[projCounter][6] = mag; p[projCounter][7] = film;  p[projCounter][8] = df1; p[projCounter][9] = df2; p[projCounter][10] = angast;  
                p[projCounter][11] = ppindex; p[projCounter][12] = occ; p[projCounter][13] = logp; p[projCounter][14] = sigma; p[projCounter][15] = score; p[projCounter][16] = change; 
                

                p[projCounter][17] = ind_in_stack;
                comparison_objects[projCounter] = CreateComparisonObject(ref_3d, 
                                                                        image_set->imageStack,
                                                                        p[projCounter], 
                                                                        refine3dParams.scopePixel,
                                                                        refine3dParams.boxsize,
                                                                        refine3dParams.beamEnergy,
                                                                        refine3dParams.sphericalAbberation,
                                                                        refine3dParams.amplitutdeContrast,
                                                                        refine3dParams.outerMask,
                                                                        refine3dParams.highResLimit,
                                                                        refine3dParams.lowResLimit,
                                                                        refine3dParams.molecularMass,
                                                                        refine3dParams.resLimitSignCC, 
                                                                        refine3dParams.useStatistics,
                                                                        refine3dParams.statisticFile, 
                                                                        refine3dParams.usePriors,
                                                                        priors_average.at(scanor), 
                                                                        priors_variance.at(scanor)
                                                                        );
                projCounter ++;
            }
        }
        number_projections_per_refinement = projCounter;
        for ( int i = 0; i < numberOfParticleProjections; i++ ) {
            delete [] p[i];
        }
        delete [] p;

        return comparison_objects;
    }
}

void DeleteComparisonObjects(ImageProjectionComparison** comparison_objects, int number_of_projections) {

    for (int i = 0; i < number_of_projections; i++) {
        DeleteComparisonObject(comparison_objects[i]);
    }
    delete [] comparison_objects;
}



/** Read a .par FREALIGN parameters file, for each particle, recover the
 * center (SHX and SHY) and the cross-correlation function (PRESA). */
int ReadParFile( LoopModeType mode, int objectIndex, char filename[200] ) {

    FILE *file;
    char line[2000];
    char CCFilename[200];
    int ind;
    int film;
    int counter;
    int mag;
    float df1;
    float df2;
    float angast;
    float ppindex;
    float occ;
    int logp;
    float sigma;
    float score;
    float change;
    // float presa;
    // float dpres;
    // float ccx;
    float toler = 1e-2;
    EulerAngleType psi;
    EulerAngleType theta;
    EulerAngleType phi;
    CoordinateType shx;
    CoordinateType shy;
    CEulerAngles orientation;
    CPosition2D shift;
    bool showInfo = false;

    CTiltSeries* pTiltSeries;
    CParticle* pParticle;
    CMicrograph* pMicrograph;
    CParticleProjection* pPP;
    vector<CParticleProjection*> pParticleProjectionSet;

    // Back up the file name
    strcpy( CCFilename, filename );

    file = fopen( filename, "r" );
    if (!file) {
        printf( "\n[ERROR] ReadParFile(): File %s does not exist.\nExiting...\n", filename );
        exit( EXIT_FAILURE );
    }
    PrintDebug( DebugData, "== Read file %s\n", 1, filename );

    

        // The rest of the variables must not change (check it).
        //      if (pPP->mIndexInStack != ind) {
        //        printf( "[ERROR] ReadParFile in ind: %d != %d\n", pPP->mIndexInStack, ind );
        //        showInfo = true;
        //      }
        if (( fabs( pPP->mOrientation.GetPsi() - psi ) > toler ) && //
            ( fabs( pPP->mOrientation.GetPsi() - psi - 360 ) > toler )) {
            printf( "f in psi: %lf != %lf\n", pPP->mOrientation.GetPsi(), psi );
            showInfo = true;
        }
        if (( fabs( pPP->mOrientation.GetTheta() - theta ) > toler ) && //
            ( fabs( pPP->mOrientation.GetTheta() - theta - 360 ) > toler )) {
            printf( "[ERROR] ReadParFile in theta: %lf != %lf\n", pPP->mOrientation.GetTheta(), theta );
            showInfo = true;
        }
        if (( fabs( pPP->mOrientation.GetPhi() - phi ) > toler ) && //
            ( fabs( pPP->mOrientation.GetPhi() - phi - 360 ) > toler )) {
            printf( "[ERROR] ReadParFile in phi: %lf != %lf\n", pPP->mOrientation.GetPhi(), phi );
            showInfo = true;
        }
        if (fabs( pPP->mMagnification - mag ) > toler) {
            printf( "[ERROR] ReadParFile in mag: %lf != %lf\n", pPP->mMagnification, mag );
            showInfo = true;
        }
        if (pPP->mTiltSeriesIndex != film) {
            // frealignx will change the film column to INCLUDE
            //printf( "[ERROR] ReadParFile in film: %d != %d\n", pPP->mTiltSeriesIndex, film );
            //showInfo = true;
        }
        if ( (fabs( pPP->mDefocus01 - df1 ) > toler) && (mode != OnMicrographDefocus) ) {
            printf( "[ERROR] ReadParFile in df1: %lf != %lf\n", pPP->mDefocus01, df1 );
            showInfo = true;
        }
        if ( (fabs( pPP->mDefocus02 - df2 ) > toler) && (mode != OnMicrographDefocus) ) {
            printf( "[ERROR] ReadParFile in df2: %lf != %lf\n", pPP->mDefocus02, df2 );
            showInfo = true;
        }
        if ((fabs( pPP->mAstigmatism - angast ) > toler) && (mode != OnMicrographDefocus)) {
            printf( "[ERROR] ReadParFile in angast: %lf != %lf\n", pPP->mAstigmatism, angast );
            showInfo = true;
        }

        if (( showInfo )) {
        printf( "\n[ERROR] ReadParFile() reading %s \n", filename );
            pPP->ShowInformation();
            printf( "Exiting...\n" );
            exit( EXIT_FAILURE );
        }

        // Count the number of particle projections read.
        counter++;
        
    
    fclose( file );

    return counter;
} // end ReadParFile()

/** Compute our cost function, using only the particles given by the
 * optimization <mode>, and the <index> (corresponding to the micrograph
 *  index or particle index depending on mode). First select the
 *  particles of interest and build the vector of values; then compute,
 *  in this case compute the summation.
 *  */
flrCostFunction::MeasureType ComputeCostFunction(LoopModeType mode, 
                                                int index, 
                                                int frameIndex, 
                                                bool frame_refine, 
                                                double phaseResidualThreshold, 
                                                double minScanOrderToUseForRefinement, 
                                                double maxScanOrderToUseForRefinement,
                                                int refineProjectionCutoff,  
                                                vector<double> weightList ) {

    flrCostFunction::MeasureType result;
    vector<float> CCX, tiltangles, scanorder, occ, confidence;
    float tValue;
    CTiltSeries* tiltSeries = CContainer::Instance()->mTiltSeriesPointer;
    CMicrograph* pM;
    CParticle* pP;
    int tInd = -1;
    unsigned int k; // for index
    unsigned int arraySize;
    flrCostFunction::MeasureType sum, norm, weight;
    vector<float>::iterator it;
    vector<CParticleProjection*>::iterator itPP;

    vector<int> positionsList; 
    int numberOfParticleProjections;

    switch (mode) {
        case OnMicrographRotation: 
        case OnMicrographTranslation: 
        case OnMicrographDefocus: 
        case OnMicrograph: {
            // FIXME
            if ( !frame_refine ) {
                numberOfParticleProjections = tiltSeries->GetAllPositionsByMicrograph( index, frameIndex, positionsList );
            }
            else {    
                numberOfParticleProjections = tiltSeries->GetAllPositionsByMicrograph( index, positionsList );
            }
            //pMicrograph = pTiltSeries->mMicrograph.at( pTiltSeries->GetMicrographPositionByIndex( objectIndex, 0 ) );
            //numberOfParticleProjections = pMicrograph->mParticleProjectionSet.size();
            //pParticleProjectionSet = pMicrograph->mParticleProjectionSet;
            break;
        }
        case OnParticleRotation:
        case OnParticleTranslation: 
        case OnParticle: {
            // FIXME
            if ( !frame_refine ) {
                numberOfParticleProjections = tiltSeries->GetAllPositionsByParticle( index, 0, positionsList );
            }
            else {
                numberOfParticleProjections = tiltSeries->GetAllPositionsByParticle( index, positionsList );
            }
            //pParticle = pTiltSeries->mParticle.at( pTiltSeries->GetParticlePositionByIndex( objectIndex, 0 ) );
            //numberOfParticleProjections = pParticle->mParticleProjectionSet.size();
            //pParticleProjectionSet = pParticle->mParticleProjectionSet;
            break;
        }
        default:
        break;
    }
    for ( int pos = 0; pos < positionsList.size(); pos++ ) {
        
        switch (mode) {
            case OnMicrographRotation:
            case OnMicrographTranslation: 
            case OnMicrographDefocus: 
            case OnMicrograph: {
                // FIXME
                //tInd = tiltSeries->GetMicrographPositionByIndex( index, 0 );
                pM = tiltSeries->mMicrograph.at(positionsList[pos]);
                arraySize = pM->mParticleProjectionSet.size();
                PrintDebug( DebugInfo, "== ComputeCostFunction() on M%d (%d/%d): ", 3, index, tInd + 1,
                            (int) tiltSeries->mListOfMicrographs.size() );
                for ( k = 0; k < arraySize; k++ ) {        
                    tValue = pM->mParticleProjectionSet.at( k )->mScore;
                    //printf("[%d]->%f, ", pM->GetIndex(), tValue);
                    CCX.push_back( tValue );
                    occ.push_back( pM->mParticleProjectionSet.at( k )->mOcc );
                    PrintDebug( DebugData, "P%d(%d/%d)=%f ", 4, pM->mParticleProjectionSet.at( k )->mParticleIndex, k + 1,
                                arraySize, tValue );
                } 
                break;
            } 
            case OnParticleRotation:
            case OnParticleTranslation: 
            case OnParticle: {
                // FIXME
                //tInd = tiltSeries->GetParticlePositionByIndex( index, 0 );
                pP = tiltSeries->mParticle.at(positionsList[pos]);
                arraySize = pP->mParticleProjectionSet.size();
                PrintDebug( DebugData, "== ComputeCostFunction() on P%d (%d/%d): ", 3, index, tInd + 1,
                            (int) tiltSeries->mListOfParticles.size() );
                for ( k = 0; k < arraySize; k++ ) {
                    
                    tValue = pP->mParticleProjectionSet.at( k )->mScore;
                    //printf("%f, ", tValue);
                    // retrieve tilt-angle of current particle projection
                    CCX.push_back( tValue );
                    tiltangles.push_back( tiltSeries->mMicrograph.at( tiltSeries->GetMicrographPositionByIndex( pP->mParticleProjectionSet.at( k )->GetMicrographIndex(), pP->mParticleProjectionSet.at( k )->GetFrameIndex() ) )->GetTiltAngle() );
                    scanorder.push_back( pP->mParticleProjectionSet.at( k )->mMicrographsIndex );
                    occ.push_back( pP->mParticleProjectionSet.at( k )->mOcc );
                    // confidence.push_back( pP->mParticleProjectionSet.at( k )->mConfidence );
                    PrintDebug( DebugData, "M%d(%d/%d)=%f ", 4, pP->mParticleProjectionSet.at( k )->mMicrographsIndex, k + 1,
                                arraySize, tValue );
                } 
                break;
            } 
            default: {
                printf( "[ERROR] On ComputeCostFunction: incorrect mode (%d).\n", mode );
                exit( EXIT_FAILURE );
            }
        }
    }
    //printf("\n");

    PrintDebug( DebugData, "CCX: [ " );

    // // Summation
    // sum = 0;
    // flrCostFunction::MeasureType constrainedSum = 0;
    // flrCostFunction::MeasureType constrainedCount = 0;
    // for ( it = CCX.begin(); it < CCX.end(); it++ ) {
        // sum = sum + *it;
        // PrintDebug( DebugData, "%f ", 1, *it );
        // if ( *it < phaseResidualThreshold ){
        // constrainedSum = constrainedSum + *it;
        // constrainedCount++;
        // }
    // }

    sum = 0.0;
    norm = 0.0;
    int valid = 0; 
    // Micrograph modes
    if ( tiltangles.size() == 0 ){
	    for ( int i = 0; i < CCX.size(); i++ ) {
            if ( occ[i] > 0.0 ) {
		        sum = sum + CCX[i];
                norm ++;
                valid++;
            }
	    }
    } 
    // Particle modes (with tilt-weighting)
    else {
	    // compute weights based on tilt-angle and scan-order
	    for ( int i = 0; i < tiltangles.size(); i++ ){
            weight = 0;
                  
            // use tilts from the starting point all the way to the end
            if ( maxScanOrderToUseForRefinement == -1 ) {
                if ( scanorder[i] >= minScanOrderToUseForRefinement ) {
                    valid++;
                    weight = 1; //weightList[ scanorder[i] ];
                }
                else {continue;}
            }
            // use tilts within user-specified range
            else {
                if ( (scanorder[i] >= minScanOrderToUseForRefinement) && (scanorder[i] <= maxScanOrderToUseForRefinement) ) {
                    valid++;
                    weight = 1; //weightList[ scanorder[i] ];
                }
                else {continue;}
            }
        
            /**     
            if (-15.0 <= tiltangles[i] && tiltangles[i] <= 15.0) {
                weight = 1;
                valid++;
            }
            */
            sum += weight * CCX[i];
		    norm += weight;
	    }
    }
    if ( norm > 0 ){
	    result = sum / norm;
    } 
    else {
	    /** 
	    result = -itk::NumericTraits<double>::infinity();
	    // use all available projections
	    sum = 0;
        for ( int i = 0; i < CCX.size(); i++ ) {
            if ( occ[i] > 0.0 ) {
                sum = sum + CCX[i];
            }
        }
	    norm = CCX.size();
	    result = sum / norm;
        */
        result = 0.0;
    }
    PrintDebug( DebugData, "] => %f\n", 1, result );
    PrintDebug( DebugData, " => cost: %f\n", 1, result );
    
    return (valid < refineProjectionCutoff ) ? -100 : result;
}

/** ReconstructDensityMap. This function is called from the main
 * function and calls the run_frealign_reconstruction.sh script. The
 * input is the filename of the .par file, ready for the
 * reconstruction. The script first refines the shifts for each
 * particle projection and then reconstruct the density map.
 * */
int ReconstructDensityMap( char pathfilename[200] ) {

  int returnValue;
  char filename[200];
  char command[200];
  char* p;

  // Parse the path and filename to obtain only the filename.
  if (( p = strrchr( pathfilename, '/' ) ) == NULL) {
    strcpy( filename, pathfilename );
  } else {
    strcpy( filename, p + 1 );
  }

  // Build the command to run.
  //sprintf( command, "./run_frealign_reconstruction.sh %s >> run_frealign_reconstruction.log", filename );
  sprintf( command, "${CSPDIR}/run_frealign_reconstruction.sh %s", filename );
  PrintDebug( DebugInfo, "Starts FREALIGN reconstruction (%s)\n", 1, filename );
  returnValue = system( command );
  PrintDebug( DebugInfo, "Ends FREALIGN reconstruction\n" );

  return returnValue;
}

/** Print the debug information, or not, depending in what kind of
 * information is.
 * */
int PrintDebug( DebugInfoType infoType, const char *format, int numarg, ... ) {
  int returnValue;
  va_list args;
  va_start(args, numarg);

  if (DebugOrNotDebug( infoType )) {
    returnValue = vprintf( format, args );
  } else {
    returnValue = -1;
  }

  return returnValue;
}

/** Determines if the debug type of information should be displayed or
 * not.
 * */
bool DebugOrNotDebug( DebugInfoType infoType ) {
  extern SParameters params;
  bool debugOrNot = false;

  switch (infoType) {
    case DebugInfo:
      debugOrNot = params.mDebugInfo;
      break;
    case DebugData:
      debugOrNot = params.mDebugData;
      break;
    case DebugBasic:
      debugOrNot = params.mDebugBasic;
      break;
    default:
      break;
  }

  if (params.mDebugFull) debugOrNot = true;
  if (params.mDebugNone) debugOrNot = false;

  return debugOrNot;
}
