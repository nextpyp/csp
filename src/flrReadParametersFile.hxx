/*=========================================================================

 Name:        flrReadParametersFile.hxx
 Category:    Header file
 Description: Headers for the parameters reading functions.
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrReadParametersFile.hxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/

#ifndef FLRREADPARAMETERSFILE_H_
#define FLRREADPARAMETERSFILE_H_

#include "flrClasses.hxx"
#include "flrMain.hxx"
#include <stdio.h>
#include <string.h>

// For the event observer
#include "itkCommand.h"

template <typename T>
std::string to_string(T value);

bool checkExistence(const char* filename);
std::string RemoveQuotations( char value[] );
double GetParamInListByIter( char value[], int iter, std::string delimiter=":" );

/** Structure with the parameters defined by the user. */
struct SParameters {
  /** Number of rings defined in the FT. */
  unsigned int mNumberOfRings;

  /** B-factor */
  float mBFactor;

  /** Folders */
  char mWorkingFolder[200];
  char mFrealignFolder[200];
  char mDataFolder[200];

  /** Par file from FREALIGN with the parameters for each image (particle projection) */
  char mParFile[80];
  char mBaseName[80];
  char mRootName[80];

  /** Number of images (particles projections) */
  unsigned int mTotalNumberOfParticlesProjections;
  float mInnerRadius, mOuterRadius;

  /** Optimizer parameters */
  unsigned int mOptimizerMaxIter; // Maximum number of iteration for the optimizer
  double mOptimizerValueTolerance; // Tolerance
  double mOptimizerStepTolerance; // Tolerance
  double mOptimizerStepLength; // Step length

  /** Image size */
  unsigned int mImageSizeX, mImageSizeY;

  int mNumberOfIterations;
  int mNumberOfRandomIterations;
  int mNumberOfThreads;

  /** Select variable to refine */
  bool mRefineMicrographTiltAngles;
  bool mRefineMicrographTiltAxisAngles;
  bool mRefineParticlesTheta;
  bool mRefineParticlesPsi;
  bool mRefineParticlesPhi;

  int mRefineProjectionCutoff; 

  double mToleranceMicrographTiltAngles;
  double mToleranceMicrographTiltAxisAngles;
  double mToleranceMicrographShifts;
  
  double mToleranceMicrographDefocus1;
  double mToleranceMicrographDefocus2;
  double mToleranceMicrographAstigmatism;
  
  double mToleranceParticlesTheta;
  double mToleranceParticlesPsi;
  double mToleranceParticlesPhi;
  double mToleranceParticlesShifts;

  int mUseImagesForRefinementMax;
  int mUseImagesForRefinementMin;

  int mUseImagesForReconstructionMax;
  int mUseImagesForReconstructionMin;


  /** Debug or not debug */
  bool mDebugInfo;
  bool mDebugData;
  bool mDebugFull;
  bool mDebugNone;
  bool mDebugBasic;

};

struct Refine3dParameters {
    /**
     * Parse parameters from .pyp_config (p) and frealign.config (f)
        Assume it only does"local" search without refining any geometry parameters (only evaluate scores)
     **/
    int currentIter = 0;
    
    bool isSPR = false;
    bool isTOMO = false;

    char* scratch = "/scratch"; 

    bool useStatistics = true; // f
    char statisticFile[200];
    
    bool usePriors = false;

    float particle_blur_sigma = -1.0;

    char symmetry[3]; // f 
    int firstParticleRefine; 
    int lastParticleRefine;
    double scopePixel; // p
    int dataBin; // p
    int particleBin; // p
    int boxsize;
    int numberFrames;
    double beamEnergy;  // p
    double sphericalAbberation; // p
    double amplitutdeContrast; // p
    double molecularMass; // p
    double outerMask; // p 
    double lowResLimit; // f
    double highResLimit; // f
    double resLimitSignCC = 30.0; // f 
    double maskRadiusGlobalSearch; // f 
    
    // The following parameters cannot be customized
    // as we only use refine3d to evaluate scores
    double ResLimitClassification = 8.0; // f 
    double AngularStep = 200; // f 
    double paddingFactor; // f
    bool globalSearch = 0; // f
    bool localSearch = 1; // f
    bool refinePsi = 0; // f
    bool refineTheta = 0; // f 
    bool refinePhi = 0; // f 
    bool refineShiftX = 0; // f 
    bool refineShiftY = 0; // f

};

extern Refine3dParameters refine3dParams;
extern int NUM_COL_REFINE3D;

void ProcessCommandLineOptions( int argc, char * argv[], SParameters *params );

/** Read the configuration (parameters) in <filename> and fill the
 * <params> structure. */
int LoadConfiguration( SParameters *params, char filename[80] );

/** Read the pyp and fyp configuration and fill the structure. */
int LoadParameters ( Refine3dParameters *refine_params, SParameters *csp_params, char pyp_config[80] );

/** Print the values of the parameters. */
void ShowConfiguration( SParameters *params );

/** Process the parameters, computing and reading the images, build the
 * data structure of CTiltSeries, etc. */
//int ProcessConfiguration( SParameters *params, CTiltSeries *tiltSerie );

/** Process the parameters, computing and reading the images, build the
 * data structure of CTiltSeries, etc. */
int ProcessConfiguration( SParameters *params, char parxfile[300], vector< CTiltSeries > & tiltSerie, int mode, int first_index, int last_index, int frame_index, vector<double> & weights, vector<int> & images_to_extract, vector< pair<double, double> > & frame_shifts, map<int, float*> & priors_average, map<int, float*> & priors_variance );

/** Write the extended .par file. */
//int WriteParFileOld();
//int WriteParFileOld( const flrCostFunction::ParametersType & parameters, LoopModeType mode );
//int WriteParFileOld( const flrCostFunction::ParametersType & parameters, LoopModeType mode, char filename[200] );
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
                  map<int, float*> priors_variance);

int WriteMultipleParFileNew( vector< CTiltSeries > &, char filename[200], char filename_new[200], LoopModeType mode, bool frame_refine, bool extract_particle_frames, int minScanOrderUsedForReconstruction, int maxScanOrderUsedForReconstruction );

void GetRefinedFrameShiftsAtTilt(vector< CTiltSeries > & dataset, int tilt_index, double ** refined_frame_shifts);
void GetRefinedFrameShifts(vector< CTiltSeries > & dataset, double ** refined_frame_shifts);

ImageProjectionComparison** PrepareComparisonObjects(LoopModeType mode,
                                                    int objectIndex,
                                                    int frameIndex,
                                                    bool frame_refine,
                                                    ReconstructedVolume* ref_3d,
                                                    CImage* image_set,
                                                    int & number_projections_per_refinement, 
                                                    map<int, float*> & priors_average, 
                                                    map<int, float*> & priors_variance
                                                    );
void DeleteComparisonObjects(ImageProjectionComparison** comparison_objects, int number_of_projections);

				  /** Call frealign. */
void RunFrealign( char pathfilename[200], char pathstackname[200], char pathreference[200] );

/** Read a .par FREALIGN parameters file, for each particle recover the
 * center (SHX and SHY) and the cross-correlation function (PRESA). */
int ReadParFile( LoopModeType mode, int objectIndex, char filename[200] );

Image** ExtractParticles( char mrc_file[1000], bool use_frames, bool stack_avg, char allboxes[1000], vector<int> images_to_extract, vector< pair<double, double> > frame_shifts, int write_to_disk, char output_stack[1000] );

void calcFrameWeights( vector< vector<double> > & frame_weights, int num_frames, double weight_width );


/** Compute our cost function, using only the particles given by the
 * optimization mode, micrograph index and particle index. */
flrCostFunction::MeasureType ComputeCostFunction( LoopModeType mode, int index, int frameIndex, bool frame_refine, double phaseResidualThreshold, double minScanOrderToUseForRefinement, double maxScanOrderToUseForRefinement, int refineProjectionCutoff, vector<double>  weightList );

/** Reconstruct the density map given the new variables. */
int ReconstructDensityMap( char pathfilename[200] );

/** Debug */
int PrintDebug( DebugInfoType infoType, const char *format, int numarg = 0, ... );
bool DebugOrNotDebug( DebugInfoType infoType );

/** Observer class */
class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self )

protected:
  CommandIterationUpdate() {
  }
  typedef const PowellOptimizerType *OptimizerPointer;
  //  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  //  typedef const OptimizerType *OptimizerPointer;

  void Execute( itk::Object *caller, const itk::EventObject & event ) {
    Execute( (const itk::Object *) caller, event );
  }

  void Execute( const itk::Object * object, const itk::EventObject & event ) {
    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer> ( object );
    if (!itk::IterationEvent().CheckEvent( &event )) {
      return;
    }
    std::cout << std::endl << optimizer->GetCurrentIteration() << " = ";
    std::cout << optimizer->GetValue() << " : ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};

#endif /* FLRREADPARAMETERSFILE_H_ */
