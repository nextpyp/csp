/*=========================================================================

Name:        flrMain.cxx
Category:    Implementation file
Description: This is intended to be the main() file.
Language:    C++
Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
Author:      $Author: fefo $
Revision:    $Rev: 492 $
Version:     $Id: flrMain.cxx 492 2010-08-21 12:43:09Z fefo $

===========================================================================*/
#define BZ_GENERATE_GLOBAL_INSTANCES

/**
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>
*/
#include "refine3d_cspt.h"

#include <dlfcn.h>

//using namespace blitz;

//#include <random/uniform.h>

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

#include "flrClasses.hxx"
#include "flrMain.hxx"
#include "flrCostFunctions.hxx"
#include "flrReadParametersFile.hxx"

#include <vtkTransform.h>

extern SParameters params;
extern int NUM_COL_REFINE3D = 18;

SParameters params = { 16, 7.6 };
Refine3dParameters refine3dParams;




int main( int argc, char * argv[] ) {

	// Usage
	if (argc < 8) {
		std::cerr << "Usage: " << argv[0] << " [parx_file(.parx)]  [mode]  [first_index]  [last_index]  [frame_refinement] [tiltseries(.mrc)]  [allboxes]  [stack_file]" << std::endl;
		return 1;
	}
    clock_t start; 
    clock_t begin = clock();
    double duration;

	struct tm *current;
	time_t now;
	time( &now );
	/**
    current = localtime( &now );
	printf( "\n %s\n Started on %02i-%02i-%02i at %02i:%02i:%02i.\n", argv[0], current->tm_year + 1900, current->tm_mon + 1, current->tm_mday,
			current->tm_hour, current->tm_min, current->tm_sec );
	cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    */
	

	ProcessCommandLineOptions( argc, argv, &params );

	vector< CTiltSeries > dataset;
	
	// File name string
	char filename[200];

	// Variables for the optimization loops.
	int loop;
	int loopM;
	int loopP;
	CEulerAngles eA;
    EulerAngleType matrix[16]; 

	// Read parameters from the configuration file.
	// LoadConfiguration( &params, "parameters.config" );
	// fill pyp, fyp, csp parameters in struct for refine3d
    LoadParameters( &refine3dParams, &params, ".pyp_config.toml" );

    // read the parx file (assume it is in relative folder of frealign/maps/) 
    char parxfile[300];
    sprintf( parxfile, "frealign/maps/%s", argv[1] );
    // mRootName is the same as parxfile but without .parx
    strncpy( params.mRootName, argv[1], strcspn( argv[1], "." ));
    

    // get the mode of operation
	// 0 - micrograph rotation mode (deprecated)
	// 1 - particles rotation mode (deprecated)
    // 2 - particles translation mode (deprecated)
    // 3 - frame refinement ( previous micrograph translation mode )
    // 4 - micrograph defocus refinement mode 
    // 5 - particle refinement 
    // 6 - micrograph refinement 

	int	mode = atoi( argv[2] );
	if ( ( mode <  -2 ) || ( mode > 6 ) ) {
		printf("CSP refinement only supports mode 0 - 4\nmode -1 : bypass refinement\n\mode -2 : only extract particle stacks\n");
		return 1;
	}
	// get what image in the parfile to be refined
	int first_index = atoi( argv[3] );
	int last_index = atoi( argv[4] );
    int frame_index = atoi( argv[5] );

	if ( first_index < 0 ) {
		printf( "First index should be an integer larger than zero.\n" );
		return 1;
	}
	if ( (last_index != -1) && (last_index < first_index) ) {
		printf( "Last index should be equal or greater than frist index, otherwise no images are refined.\n" );
		return 1;
	}

	// get weighting curve used for refinement by computing average score from previous parfile
    vector<double> weights;
    vector<int> images_to_extract;
    vector< pair<double, double> > frame_shifts;
    
    map <int, float*> priors_average; 
    map <int, float*> priors_variance;
    
    double** refined_frame_shifts;
    double** refined_frame_shifts_tilt;
    int number_of_images; 
    

	ProcessConfiguration(&params, 
                        parxfile, 
                        dataset, 
                        mode, 
                        first_index, 
                        last_index, 
                        frame_index, 
                        weights, 
                        images_to_extract, 
                        frame_shifts, 
                        priors_average, 
                        priors_variance);
    number_of_images = images_to_extract.size();
   
    // convert frame weights to array before passing it via extern c 
    double* frame_weights = new double[weights.size()];  
    for (int i = 0; i < weights.size(); i++) {
        frame_weights[i] = weights.at(i);
    }
    
    if (images_to_extract.size() == 0) {
        printf("No particle found using index %d - %d\n", first_index, last_index);
        return EXIT_SUCCESS;
    }

    /**
     * Particle extraction 
     */
    Image** particle_images; 
    Image** weighted_averages;
    CImage* raw_particles;
    CImage* weighted_particles;
    
    // werights for creating running averages for frame alignment
    // spr uses 15.0, tomo uses 0.05 (refine on almost raw frames)
    double sigma = (refine3dParams.isSPR) ? 15.0 : 0.05;

    bool frame_refine = atoi(argv[5]);
    bool write_to_disk = (mode == -2);
    char output_stack[300];
    sprintf( output_stack, "%s", argv[8] );

    bool extract_particle_frames = false;
    bool stack_avg = false;
    
    // assume if mode == 3, we use refine3d internal optimizor 
    bool cistem_optimizor = (mode == 3 || mode == 4);
    
    if ( true ) {    
        
        // check if using frames or not (based on argv[5]: .mrc -> normal ; .txt -> use frames)
        std::string extension, file(argv[6]), mrc("mrc"), txt("txt");
        std::string::size_type idx;
        idx = file.rfind('.');
        if ( idx != std::string::npos ) {
            extension = file.substr(idx+1);
        }
        else {
            printf( "[ERROR] File %s extension not found. Exit...\n" , argv[6]);
            return 1;
        }
        
        if ( extension.compare(mrc) == 0 ) {
            extract_particle_frames = false;
            frame_refine = false;
        }
        else if ( extension.compare(txt) == 0 ) {
            extract_particle_frames = true;
        }
        else {
            printf("File extension [.%s] not supported. Exit ...\n", extension.c_str());
            return EXIT_FAILURE;
        }

        if ( extract_particle_frames && ( !write_to_disk || !frame_refine ) ) {
            // perform weighted running average for refinement except mode -1 for generating particle stack
            stack_avg = true;
        }
        if ( mode != -1 ) {
            
            bool write_raw = write_to_disk && ( !extract_particle_frames || frame_refine );
            ifstream particle_stack(output_stack);
            // read particle from stack if stack exists
            // start = clock();
            if (mode != -2 && particle_stack.good()) {
                int images_to_extract_arr[number_of_images];
                for ( int i = 0; i < images_to_extract.size(); i++ ) { images_to_extract_arr[i] = images_to_extract.at(i); }
                particle_images = ReadParticlesFromStack( output_stack, images_to_extract_arr, number_of_images, refine3dParams.boxsize );
            }
            // otherwise extract particles from movies
            else {
                particle_images = ExtractParticles( argv[6], extract_particle_frames, stack_avg, argv[7], images_to_extract, frame_shifts, write_raw, output_stack );
            }
            // duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
            // printf("Time spent for particle extraction: %f s\n", duration);
            
            if (mode != 2 && mode != 3 && refine3dParams.particle_blur_sigma > 0) ApplyGaussianToImages(particle_images, number_of_images, refine3dParams.boxsize, refine3dParams.scopePixel, refine3dParams.particle_blur_sigma);

            raw_particles = new CImage(particle_images, number_of_images, true);
            if (stack_avg) {

                refined_frame_shifts = new double*[number_of_images];
                refined_frame_shifts_tilt = new double*[refine3dParams.numberFrames];
                for (int i = 0; i < number_of_images; i++) { refined_frame_shifts[i] = new double[3]; } 
                for (int i = 0; i < refine3dParams.numberFrames; i++) { refined_frame_shifts_tilt[i] = new double[3]; }

                bool write_average = ( write_to_disk != write_raw );
                GetRefinedFrameShifts(dataset, refined_frame_shifts);
                
                weighted_averages = WeightedAverages(particle_images, refined_frame_shifts, number_of_images, refine3dParams.numberFrames, refine3dParams.boxsize, refine3dParams.scopePixel, frame_refine, frame_weights, sigma, write_average, output_stack);
                
                weighted_particles = new CImage(weighted_averages, number_of_images, false);
            }
        }
        if ( write_to_disk ) { 
            printf( "Write out %s successfully. Exit ...\n", argv[8] );
            duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
            printf("Total time spent: %f s\n", duration);

            return EXIT_SUCCESS; 
        }
    
    }
    
    // read external reconstruction (reference) for refinement
    char reconstruction[200];  
    sprintf( reconstruction, "%s/%s_frames_CSP_01.mrc", refine3dParams.scratch, params.mBaseName );
    
    MRCFile* ref_mrc = ReadReconstruction( reconstruction, false, true );
    
    ReconstructedVolume* ref_3d = Init3DVolume( ref_mrc, refine3dParams.scopePixel, refine3dParams.symmetry, refine3dParams.molecularMass, refine3dParams.outerMask, refine3dParams.lowResLimit, refine3dParams.highResLimit );
    
    CloseMRCFile( ref_mrc );

    // Get the instance of the CContainer, run the constructor, and
	// associate the tiltseries with it.
	CContainer::Instance();
	// CContainer::Instance()->SetTiltSeries( &dataset );
	CContainer::Instance()->SetTiltSeries( &dataset[0] );

	// Optimizer definition.
	PowellOptimizerType::Pointer optimizer = PowellOptimizerType::New();

	// Cost function definition.
	flrCostFunction::Pointer costFunction = flrCostFunction::New();
	const unsigned int spaceDimension = costFunction->GetNumberOfParameters();
	costFunction->SetScanOrderRangeToUseForRefinement( params.mUseImagesForRefinementMin, params.mUseImagesForRefinementMax );
    
    costFunction->SetRefineProjectionCutoff(params.mRefineProjectionCutoff);
    // save the pointer to volume object of reference (used by refine3d)
    costFunction->SetReference( ref_3d );

    costFunction->SetPriors( priors_average, priors_variance );
    
    // set if refining particles on a frame basis or just tilt basis
    costFunction->SetFrameRefine( frame_refine );

    costFunction->UseCistemOptimizor( cistem_optimizor );
	
    // Points in the optimization space.
	typedef flrCostFunction::ParametersType ParametersType;
	ParametersType initialPosition( spaceDimension );
	ParametersType finalPosition( spaceDimension );
	ParametersType dummyPosition( spaceDimension );

	// Initial and final cost.
	flrCostFunction::MeasureType initialCost;
	flrCostFunction::MeasureType finalCost;

	LoopModeType loopMode;

	// It's a maximization.
	optimizer->SetMaximize( true );

	// Set the cost function to optimize.
	optimizer->SetCostFunction( costFunction.GetPointer() );

	// Set the configuration parameters. If they are negative the
	// default values are used.
	if (params.mOptimizerStepLength >= 0) optimizer->SetStepLength( params.mOptimizerStepLength );
	if (params.mOptimizerStepTolerance >= 0) optimizer->SetStepTolerance( params.mOptimizerStepTolerance );
	if (params.mOptimizerMaxIter >= 0) optimizer->SetMaximumIteration( params.mOptimizerMaxIter );
	if (params.mOptimizerValueTolerance >= 0) optimizer->SetValueTolerance( params.mOptimizerValueTolerance );
	optimizer->GlobalWarningDisplayOn();

    set<int>::iterator it; 
    int max_ind_frame = 0;
    int number_projections_per_refinement = 0;
    int mInd, pInd;
    
    // If micrograph modes 
	if ( mode == 0 || mode == 3 || mode == 4 || mode == 6 ){
        double time_refinement = 0.0;
		for ( int tiltSeriesIndex = 0; tiltSeriesIndex < dataset.size(); tiltSeriesIndex++ ) {
			dataset[tiltSeriesIndex].ComputeMicrographTiltAngleBounds(false);
			
			CContainer::Instance()->SetTiltSeries( &dataset[tiltSeriesIndex] );
            
            //start = clock();
            for ( it = dataset[tiltSeriesIndex].mListOfMicrographs.begin(); it != dataset[tiltSeriesIndex].mListOfMicrographs.end(); it++ ) {
                
                // no need to sperate frames since we need to refine all frames together
                max_ind_frame = 0;
                
                for ( int ind_frame = 0; ind_frame <= max_ind_frame; ind_frame++ ) {
                    
                    ImageProjectionComparison** comparison_objects;
                     
                    if (extract_particle_frames && !frame_refine) {
                        costFunction->SetImageSet( weighted_particles );
                    }
                    else {
                        costFunction->SetImageSet( raw_particles );
                    }

                    mInd = dataset[ tiltSeriesIndex ].GetMicrographPositionByIndex(*it, ind_frame);
                    if (mInd == -1) {continue;}

                    if ( mode == 0 ) {
            		    loopMode = OnMicrographRotation;
                    }
                    else if ( mode == 3 ) {
                        loopMode = OnMicrographTranslation;
                        // use weighted average instead 
                        costFunction->SetImageSet( weighted_particles );
                        // normal refine3d refinement (spr) 
                        if (!frame_refine) {refine3dParams.resLimitSignCC = 30.0;}
                    }
                    else if ( mode == 4 ) {
                        loopMode = OnMicrographDefocus;
                    }
                    else if ( mode == 6 ) {
                        loopMode = OnMicrograph;
                        comparison_objects = PrepareComparisonObjects(loopMode,
                                                                    *it,
                                                                    ind_frame,
                                                                    frame_refine,
                                                                    ref_3d,
                                                                    (extract_particle_frames && !frame_refine) ? weighted_particles : raw_particles,
                                                                    number_projections_per_refinement,
                                                                    priors_average, 
                                                                    priors_variance);
                        costFunction->SetComparisonObjects(comparison_objects);
                    }
                    
                    costFunction->SetMode( loopMode );
                    costFunction->SetCurrentOptimizationIndex( *it );
                    costFunction->SetCurrentOptimizationFrameIndex( ind_frame );

                    initialPosition[0] = 0.0; 
                    initialPosition[1] = 0.0; 
                    initialPosition[2] = 0.0; 
                    initialPosition[3] = 0.0; 
                    initialPosition[4] = 0.0;
                    initialPosition[5] = 0.0;
                    

                    optimizer->SetInitialPosition( initialPosition );
                    CContainer::Instance()->SetInitialPoint( initialPosition );
                    CContainer::Instance()->SetPreviousPoint( initialPosition );
                    
                    CContainer::Instance()->mNumberOfGetValueCalled = 1;
                    // if the micrograph is not in the refinement range, only evaluate once
                    if (!( ( *it >= params.mUseImagesForRefinementMin ) && \
                        ( ( params.mUseImagesForRefinementMax == -1 ) || ( *it <= params.mUseImagesForRefinementMax ) ) )) {
                        costFunction->EvaluateAllProjections(true);
                        costFunction->SetImageSet( raw_particles );
                        initialCost = costFunction->Evaluate( initialPosition );
                        CContainer::Instance()->SetInitialValue( initialCost );
                        CContainer::Instance()->SetPreviousValue( initialCost );
                        costFunction->EvaluateAllProjections(false);
                        continue;
                    }
                    
                    if ( mode == OnMicrographDefocus ) {

					    flrCostFunction::MeasureType bestCost = costFunction->Evaluate( initialPosition );
					    ParametersType bestInitialPosition( spaceDimension ), currentPosition( spaceDimension );
					    bestInitialPosition[0] = currentPosition[0] = initialPosition[0];
					    bestInitialPosition[1] = currentPosition[1] = initialPosition[1];
					    bestInitialPosition[2] = currentPosition[2] = initialPosition[2];
					    bestInitialPosition[3] = currentPosition[3] = initialPosition[3];
					    bestInitialPosition[4] = currentPosition[4] = initialPosition[4];
                        bestInitialPosition[5] = currentPosition[5] = initialPosition[5];

                        const float DEF1_TOL = params.mToleranceMicrographDefocus1; 
                        const float DEF2_TOL = params.mToleranceMicrographDefocus2;
                        const float STEP = 50.0;
                        for ( float i = -DEF1_TOL; i <= DEF1_TOL; i += STEP ) {
                            for ( float j = -DEF2_TOL; j <= DEF2_TOL; j += STEP ) {
                                currentPosition[3] = i;
                                currentPosition[4] = j;
                                flrCostFunction::MeasureType currentCost = costFunction->Evaluate( currentPosition );
                                // printf("%f <= ( df1 = %f, df2 = %f )\n", currentCost, currentPosition[3], currentPosition[4]);
                                if ( currentCost > bestCost ){
                                    bestCost = currentCost;
                                    bestInitialPosition[0] = currentPosition[0];
                                    bestInitialPosition[1] = currentPosition[1];
                                    bestInitialPosition[2] = currentPosition[2];
                                    bestInitialPosition[3] = currentPosition[3];
                                    bestInitialPosition[4] = currentPosition[4];
                                    bestInitialPosition[5] = currentPosition[5];
                                }
                            }
                        }
                        initialPosition[0] = bestInitialPosition[0];
                        initialPosition[1] = bestInitialPosition[1];
                        initialPosition[2] = bestInitialPosition[2];
                        initialPosition[3] = bestInitialPosition[3];
                        initialPosition[4] = bestInitialPosition[4];
                        initialPosition[5] = bestInitialPosition[5];

                    }

                    optimizer->SetInitialPosition( initialPosition );
                    CContainer::Instance()->mNumberOfGetValueCalled = 1;

				    CContainer::Instance()->SetInitialPoint( initialPosition );
				    CContainer::Instance()->SetPreviousPoint( initialPosition );
				    initialCost = costFunction->Evaluate( initialPosition );
			        CContainer::Instance()->SetInitialValue( initialCost );
				    CContainer::Instance()->SetPreviousValue( initialCost );

                    if (initialCost < -50) {
                        continue;
                    }


                    switch ( loopMode ) {
                        
                        case OnMicrographRotation: {
                            // compute the Powell step based on half the distance to the closest neighboring micrograph
                            float minb = fabs( dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() - dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngleLowerBound() );
                            float maxb = fabs( dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() - dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngleUpperBound() );
                            float step = min( minb, maxb ) / 2;
                            //optimizer->SetStepLength( step );
                        
                            if ( fabs( CContainer::Instance()->mTolerance[2] ) + fabs( CContainer::Instance()->mTolerance[3] ) > 0 ) {
                                // Run the optimization.
                                try {
                                    optimizer->StartOptimization();
                                } catch ( itk::ExceptionObject & e ) {
                                    printf( "[ERROR] Optimizing micrograph %d: %s\n", *it, e.GetDescription() );
                                    return EXIT_FAILURE;
                                }
                        
                                // Get the solution.
                                finalPosition = optimizer->GetCurrentPosition();
                                finalCost = optimizer->GetCurrentCost();
                                
                            }
                            else {
                                finalPosition = initialPosition;
                                finalCost = costFunction->Evaluate( initialPosition );
                            }
                            costFunction->Evaluate( finalPosition );

                            // print outputs in logs
                            PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
                            tiltSeriesIndex, *it );
                            PrintDebug( DebugBasic, "\nRotation refined from : f( Tilt = %8.4f, Axis = %8.4f ) = %8.4f =>", 3,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() + initialPosition[2], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisAngle() + initialPosition[3], initialCost );
                            PrintDebug( DebugBasic, "\nRotation refined to   : f( Tilt = %8.4f, Axis = %8.4f ) = %8.4f\n", 3,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() + finalPosition[2], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisAngle() + finalPosition[3], finalCost ); 
                            
                            // Update refined micrograph rotations
                            dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAngle( dataset[tiltSeriesIndex].mMicrograph.at( mInd )->GetTiltAngle() + finalPosition[2] );
                            dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAxisAngle( dataset[tiltSeriesIndex].mMicrograph.at( mInd )->GetTiltAxisAngle() + finalPosition[3] );
                                
                            break;
                        }
                        
                        case OnMicrographTranslation: {

                            // Update refined micrograph translations
                            CPosition2D imageShift;
                            if ( frame_refine ) {
                                int iteration = 5;
                                for ( int iter = 0; iter < iteration; iter++ ) {
                                    // update weighted average
                                    GetRefinedFrameShiftsAtTilt(dataset, *it, refined_frame_shifts_tilt);
                                    UpdateWeightedAverages(particle_images, weighted_averages, refined_frame_shifts_tilt, number_of_images, refine3dParams.numberFrames, refine3dParams.boxsize, refine3dParams.scopePixel, sigma);
                                    // run refinement again (it will update the shifts internally so no need to do it here)
                                    start = clock();
                                    initialCost = costFunction->Evaluate( initialPosition );
                                    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
                                    time_refinement += duration; 
                                }
                                costFunction->SetImageSet( raw_particles );
                                initialCost = costFunction->Evaluate( initialPosition );
                                
                                
                                /**
                                PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d frame %d ###\n", 3,
                                tiltSeriesIndex, *it, ind_frame );
                                PrintDebug( DebugBasic, "\nTranslation refined from : f( X = %8.4f, Y = %8.4f ) = %8.4f =>", 3,
                                dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordX + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordY + initialPosition[1], initialCost );
                                PrintDebug( DebugBasic, "\nTranslation refined to   : f( X = %8.4f, Y = %8.4f ) = %8.4f\n", 3,
                                dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordY + finalPosition[1], finalCost );
                                imageShift.SetPosition( dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetFrameShift().mCoordY + finalPosition[1] );
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetFrameShift( imageShift );
                                */
                            }
                            else {
                                /**
                                try {
                                    optimizer->StartOptimization();
                                } catch ( itk::ExceptionObject & e ) {
                                    printf( "[ERROR] Optimizing micrograph %d: %s\n", *it, e.GetDescription() );
                                    return EXIT_FAILURE;
                                }
                    
                                finalPosition = optimizer->GetCurrentPosition();
                                finalCost = optimizer->GetCurrentCost();
                                costFunction->Evaluate( finalPosition );


                                PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
                                tiltSeriesIndex, *it );
                                PrintDebug( DebugBasic, "\nTranslation refined from : f( X = %8.4f, Y = %8.4f ) = %8.4f =>", 3,
                                dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + initialPosition[1], initialCost );
                                PrintDebug( DebugBasic, "\nTranslation refined to   : f( X = %8.4f, Y = %8.4f ) = %8.4f\n", 3,
                                dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + finalPosition[1], finalCost );
                                imageShift.SetPosition( dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + finalPosition[1] );
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAxisShift( imageShift );
                                */
                            }
                            break;
                        }
                        
                        case OnMicrographDefocus: {
                            try {
                                optimizer->StartOptimization();
                            } catch ( itk::ExceptionObject & e ) {
                                printf( "[ERROR] Optimizing micrograph %d: %s\n", *it, e.GetDescription() );
                                return EXIT_FAILURE;
                            }

                            finalPosition = optimizer->GetCurrentPosition();
                            finalCost = optimizer->GetCurrentCost();
                            costFunction->Evaluate( finalPosition );
                            // print outputs in logs
                            PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
                            tiltSeriesIndex, *it );
                            PrintDebug( DebugBasic, "\nDefocus refined from : f( DF1 Offset = %8.4f, DF2 Offset = %8.4f, Astig = %8.4f ) = %8.4f =>", 3,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltDefocus1Offset() + initialPosition[3], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltDefocus2Offset() + initialPosition[4], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAstigmatism() + initialPosition[5], initialCost );
                            PrintDebug( DebugBasic, "\nDefocus refined to   : f( DF1 Offset = %8.4f, DF2 Offset = %8.4f, Astig = %8.4f ) = %8.4f\n", 3,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltDefocus1Offset() + finalPosition[3], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltDefocus2Offset() + finalPosition[4], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAstigmatism() + finalPosition[5], finalCost );

                            // Update micrograph defocus and astismatism
                            if (finalCost > initialCost) {
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltDefocusOffset( finalPosition[3], finalPosition[4] );
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAstigmatism( finalPosition[5] );
                            }
                            
                            break;
                        }       
                        
                        case OnMicrograph: {
                            try {
                                optimizer->StartOptimization();
                            } catch ( itk::ExceptionObject & e ) {
                                printf( "[ERROR] Optimizing micrograph %d: %s\n", *it, e.GetDescription() );
                                return EXIT_FAILURE;
                            }
                            CPosition2D imageShift;

                            finalPosition = optimizer->GetCurrentPosition();
                            finalCost = optimizer->GetCurrentCost();
                            costFunction->Evaluate( finalPosition );

                            PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3, tiltSeriesIndex, *it );
                            PrintDebug( DebugBasic, "\nTranslation refined from : f( X = %8.4f, Y = %8.4f, Tilt = %8.4f, Axis = %8.4f ) = %8.4f =>", 5,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + initialPosition[1], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() + initialPosition[2], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisAngle() + initialPosition[3], initialCost );
                            PrintDebug( DebugBasic, "\nTranslation refined to   : f( X = %8.4f, Y = %8.4f, Tilt = %8.4f, Axis = %8.4f ) = %8.4f\n", 5,
                            dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + finalPosition[1], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAngle() + finalPosition[2], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisAngle() + finalPosition[3], finalCost );
                           

                            if (finalCost > initialCost) {
                                 imageShift.SetPosition( dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( mInd )->GetTiltAxisShift().mCoordY + finalPosition[1] );
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAxisShift( imageShift );

                                // Update refined micrograph rotations
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAngle( dataset[tiltSeriesIndex].mMicrograph.at( mInd )->GetTiltAngle() + finalPosition[2] );
                                dataset[tiltSeriesIndex].mMicrograph.at( mInd )->SetTiltAxisAngle( dataset[tiltSeriesIndex].mMicrograph.at( mInd )->GetTiltAxisAngle() + finalPosition[3] );
                            }
                            DeleteComparisonObjects(comparison_objects, number_projections_per_refinement);

                            break;
                        }                  
                    }
                }
			}
		}
        //duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
        printf("Time spent for frame refinement: %f s\n", time_refinement);
        if (stack_avg) {
            //DeleteImageArray( weighted_averages, number_of_images );
            delete weighted_particles;
        }
        delete raw_particles;
	    //DeleteImageArray( particle_images, number_of_images );
    }

	// If particle mode
	else if ( (mode == 1) || (mode == 2) || (mode == 5) ) { 
		for ( int tiltSeriesIndex = 0; tiltSeriesIndex < dataset.size(); tiltSeriesIndex++ ) {
			
			CContainer::Instance()->SetTiltSeries( &dataset[tiltSeriesIndex] );
            for ( it = dataset[tiltSeriesIndex].mListOfParticles.begin(); it != dataset[tiltSeriesIndex].mListOfParticles.end(); it++ ) {
                
                ImageProjectionComparison** comparison_objects;
                pInd = dataset[ tiltSeriesIndex ].GetParticlePositionByIndex(*it, 0);
                
                // skip if occ is zero 
                if ( dataset[tiltSeriesIndex].mParticle.at( pInd )->GetOcc() <= 0.0 ) { continue; }

                // Set the mode of optimization, and the index of the current
				// particle being optimized.
                if ( mode == 1 ) {
				    loopMode = OnParticleRotation;
                }
                else if ( mode == 2 ) {
                    loopMode = OnParticleTranslation;  
                }
                else if ( mode == 5 ) {
                    loopMode = OnParticle;
                    comparison_objects = PrepareComparisonObjects(loopMode,
                                                                *it,
                                                                0,
                                                                frame_refine,
                                                                ref_3d,
                                                                (extract_particle_frames && !frame_refine) ? weighted_particles : raw_particles,
                                                                number_projections_per_refinement, 
                                                                priors_average, 
                                                                priors_variance);
                    
                    costFunction->SetComparisonObjects(comparison_objects);
                }
                
                costFunction->SetMode( loopMode );
				//costFunction->SetCurrentOptimizationIndex( dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
                costFunction->SetCurrentOptimizationIndex( *it );
                costFunction->SetCurrentOptimizationFrameIndex( 0 );
                costFunction->SetWeightList( weights );
                costFunction->EvaluateAllProjections(false);
                
                if (extract_particle_frames && !frame_refine) {
                    costFunction->SetImageSet( weighted_particles );
                }
                else {
                    costFunction->SetImageSet( raw_particles );
                }
                
                // Set the initial point in the optimization to the current
				// euler angles of the particle. The coordinates corresponding
				// to micrograph will not be affected in the optimization and
				// are different for each micrograph.
                
				initialPosition[0] = 0.0; 
				initialPosition[1] = 0.0; 
				eA = dataset[ tiltSeriesIndex ].mParticle.at( pInd )->GetOrientation();
			    dataset[ tiltSeriesIndex ].mParticle.at( pInd )->GetMatrix( matrix );	
				initialPosition[2] = 0.0; 
				initialPosition[3] = 0.0; 
				initialPosition[4] = 0.0; 
                initialPosition[5] = 0.0;
                
                if ( loopMode == OnParticle ){
					// random initial search
					flrCostFunction::MeasureType bestCost = costFunction->Evaluate( initialPosition );
					ParametersType bestInitialPosition( spaceDimension ), currentPosition( spaceDimension );
					bestInitialPosition[0] = currentPosition[0] = initialPosition[0];
					bestInitialPosition[1] = currentPosition[1] = initialPosition[1];
					bestInitialPosition[2] = currentPosition[2] = initialPosition[2];
					bestInitialPosition[3] = currentPosition[3] = initialPosition[3];
					bestInitialPosition[4] = currentPosition[4] = initialPosition[4];
                    bestInitialPosition[5] = currentPosition[5] = initialPosition[5];
					
                    // Redo global alignment if overall score is too low
                    if ( true ) {
                        
                        srand (time(NULL));

                        const int TRANS_TOL =  params.mToleranceParticlesShifts; // refine3dParams.scopePixel * refine3dParams.boxsize / 4.0; // use 1/4 of the boxsize
                        const int PSI_TOL = params.mToleranceParticlesPsi;
                        const int THETA_TOL = params.mToleranceParticlesTheta;
                        const int PHI_TOL = params.mToleranceParticlesPhi;

                        for ( int iter = 0; iter < params.mNumberOfRandomIterations; iter ++ ) {
                            currentPosition[0] = rand() % (TRANS_TOL*2) - TRANS_TOL;
                            currentPosition[1] = rand() % (TRANS_TOL*2) - TRANS_TOL;
                            currentPosition[2] = rand() % (TRANS_TOL*2) - TRANS_TOL;
                            currentPosition[3] = (params.mRefineParticlesTheta) ? rand() % (THETA_TOL*2) - THETA_TOL : 0.0;
                            currentPosition[4] = (params.mRefineParticlesPsi) ? rand() % (PSI_TOL*2) - PSI_TOL : 0.0;
                            currentPosition[5] = (params.mRefineParticlesPhi) ? rand() % (PHI_TOL*2) - PHI_TOL : 0.0;
                            flrCostFunction::MeasureType currentCost = costFunction->Evaluate( currentPosition );
                            //printf("%f <= ( x = %f, y = %f, z = %f, theta = %f, psi = %f, phi = %f )\n", currentCost, currentPosition[0], currentPosition[1], currentPosition[2], currentPosition[3], currentPosition[4], currentPosition[5]);
                            if ( currentCost > bestCost ){
                                bestCost = currentCost;
                                bestInitialPosition[0] = currentPosition[0];
                                bestInitialPosition[1] = currentPosition[1];
                                bestInitialPosition[2] = currentPosition[2];
                                bestInitialPosition[3] = currentPosition[3];
                                bestInitialPosition[4] = currentPosition[4];
                                bestInitialPosition[5] = currentPosition[5];
                            }
                        }
                        initialPosition[0] = bestInitialPosition[0];
                        initialPosition[1] = bestInitialPosition[1];
                        initialPosition[2] = bestInitialPosition[2];
                        initialPosition[3] = bestInitialPosition[3];
                        initialPosition[4] = bestInitialPosition[4];
                        initialPosition[5] = bestInitialPosition[5];
                    }
				}
				
                optimizer->SetInitialPosition( initialPosition );
				CContainer::Instance()->mNumberOfGetValueCalled = 1;

				CContainer::Instance()->SetInitialPoint( initialPosition );
				CContainer::Instance()->SetPreviousPoint( initialPosition );
				initialCost = costFunction->Evaluate( initialPosition );
				CContainer::Instance()->SetInitialValue( initialCost );
				CContainer::Instance()->SetPreviousValue( initialCost );


                if (initialCost < -1.0) {
                    continue;
                }

                // Run the optimization.
				try {
					optimizer->StartOptimization();
				} catch ( itk::ExceptionObject & e ) {
					printf( "[ERROR] Optimizing particle %d: %s\n", *it, e.GetDescription() );
					return EXIT_FAILURE;
				}
				
                // Get the solution.
				finalPosition = optimizer->GetCurrentPosition();
				finalCost = optimizer->GetCurrentCost();
                
                // update the scores using best position (need to evaluate all projections to obtain correct score distribution)
                costFunction->EvaluateAllProjections(true);
                costFunction->Evaluate( finalPosition ); 
                    
				// make final results modulo 360
				switch ( loopMode ) {
                    case OnParticleRotation: {
                        finalPosition[3] += eA.GetTheta();
				        finalPosition[4] += eA.GetPsi();
				        finalPosition[5] += eA.GetPhi();

				        if ( finalPosition[3] < 0 ){
				            finalPosition[3] = fmod(finalPosition[3],360.0)+360.0;
				        } else {
				            finalPosition[3] = fmod(finalPosition[3],360.0);
				        }
				        if ( finalPosition[4] < 0 ){
				            finalPosition[4] = fmod(finalPosition[4],360.0)+360.0;
				        } else {
				            finalPosition[4] = fmod(finalPosition[4],360.0);
				        }
				        if ( finalPosition[5] < 0 ){
				            finalPosition[5] = fmod(finalPosition[5],360.0)+360.0;
				        } else {
				            finalPosition[5] = fmod(finalPosition[5],360.0);
				        }
						
						// Update particle rotation
						dataset[tiltSeriesIndex].mParticle.at( pInd )->SetOrientation( finalPosition[3], finalPosition[5], finalPosition[4] );

						// print outputs in logs
                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, *it );
                        PrintDebug( DebugBasic, "\nRotation refined from: f( PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f =>", 4,
                        eA.GetPsi() + initialPosition[4], eA.GetTheta() + initialPosition[3], eA.GetPhi() + initialPosition[5], initialCost );
                        PrintDebug( DebugBasic, "\nRotation refined to  : f( PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f\n", 4,
                        finalPosition[4], finalPosition[3], finalPosition[5], finalCost );
                        break;
                    }

                    case OnParticleTranslation: {
                        
                        finalPosition[0] += matrix[3];
                        finalPosition[1] += matrix[7];
                        finalPosition[2] += matrix[11];

						// Update particle translation
						dataset[tiltSeriesIndex].mParticle.at( pInd )->SetMatrixTranslation( finalPosition[0], finalPosition[1], finalPosition[2] ); 
                        
						// print outputs in logs
                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, *it );
                        PrintDebug( DebugBasic, "\nTranslation refined from: f( X = %8.3f, Y = %8.3f, Z = %8.3f ) = %8.4f =>", 4,
                        matrix[3] + initialPosition[0], matrix[7] + initialPosition[1], matrix[11] + initialPosition[2], initialCost );
                        PrintDebug( DebugBasic, "\nTranslation refined to  : f( X = %8.3f, Y = %8.3f, Z = %8.3f ) = %8.4f\n", 4,
                        finalPosition[0], finalPosition[1], finalPosition[2], finalCost );
                        
                        break;
                    }

                    case OnParticle: {
                        
                        finalPosition[0] += matrix[3];
                        finalPosition[1] += matrix[7];
                        finalPosition[2] += matrix[11];
                        dataset[tiltSeriesIndex].mParticle.at( pInd )->SetMatrixTranslation( finalPosition[0], finalPosition[1], finalPosition[2] );
                        
                        finalPosition[3] += eA.GetTheta();
				        finalPosition[4] += eA.GetPsi();
				        finalPosition[5] += eA.GetPhi();

				        if ( finalPosition[3] < 0 ){
				            finalPosition[3] = fmod(finalPosition[3],360.0)+360.0;
				        } else {
				            finalPosition[3] = fmod(finalPosition[3],360.0);
				        }
				        if ( finalPosition[4] < 0 ){
				            finalPosition[4] = fmod(finalPosition[4],360.0)+360.0;
				        } else {
				            finalPosition[4] = fmod(finalPosition[4],360.0);
				        }
				        if ( finalPosition[5] < 0 ){
				            finalPosition[5] = fmod(finalPosition[5],360.0)+360.0;
				        } else {
				            finalPosition[5] = fmod(finalPosition[5],360.0);
				        }
						
						// Update particle rotation
                        if (finalCost > initialCost) { 
                            dataset[tiltSeriesIndex].mParticle.at( pInd )->SetOrientation( finalPosition[3], finalPosition[5], finalPosition[4] );
                        }
 
                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 2, tiltSeriesIndex, *it );
                        PrintDebug( DebugBasic, "\nAlignments refined from: f( X = %8.3f, Y = %8.3f, Z = %8.3f, PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f =>", 7,
                        matrix[3] + initialPosition[0], matrix[7] + initialPosition[1], matrix[11] + initialPosition[2], eA.GetPsi() + initialPosition[4], eA.GetTheta() + initialPosition[3], eA.GetPhi() + initialPosition[5], initialCost);
                        PrintDebug( DebugBasic, "\nAlignments refined to  : f( X = %8.3f, Y = %8.3f, Z = %8.3f, PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f\n", 7,
                        finalPosition[0], finalPosition[1], finalPosition[2], finalPosition[4], finalPosition[3], finalPosition[5], finalCost );

                        DeleteComparisonObjects(comparison_objects, number_projections_per_refinement);

                        break;
                    }
                }
			}
		}
        delete raw_particles;
	}
    map<int, float*>::iterator mit;
    for ( mit = priors_average.begin(); mit != priors_average.end(); mit++ ) {
        delete [] priors_average.at(mit->first);
        delete [] priors_variance.at(mit->first);
    }

    delete [] frame_weights; 
    
    if (mode != -1 && stack_avg) {
        for (int i = 0; i < number_of_images; i++) { 
            delete [] refined_frame_shifts[i];
        }
        for (int i = 0; i < refine3dParams.numberFrames; i++) {
            delete [] refined_frame_shifts_tilt[i];
        }
        delete [] refined_frame_shifts;
        delete [] refined_frame_shifts_tilt;
    }

    // delete reference
    Delete3DVolume( ref_3d );
    
	// Write parfile at the end
	sprintf( filename, "frealign/%s.series", params.mBaseName );
							
	char outfilename[200];
	sprintf( outfilename, "frealign/maps/%s_%06d_%06d.parx", params.mRootName, first_index, last_index );

    switch (mode) {
        case 0: 
            loopMode = OnMicrographRotation;
            break;
        case 1:
            loopMode = OnParticleRotation;
            break;
        case 2: 
            loopMode = OnParticleTranslation;
            break;
        case 3: 
            loopMode = OnMicrographTranslation;
            break;
        case 4:
            loopMode = OnMicrographDefocus;
            break;
        case 5:
            loopMode = OnParticle;
            break;
        case 6: 
            loopMode = OnMicrograph;
            break;
        default:
            loopMode = OnSkip;
            break; 
	}
    //start = clock();
    WriteMultipleParFileNew( dataset, filename, outfilename, loopMode, frame_refine, extract_particle_frames, params.mUseImagesForReconstructionMin, params.mUseImagesForReconstructionMax );
    //duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    //printf("Time spent for writing parfile: %f s\n", duration);
    
    duration = ( clock() - begin ) / (double) CLOCKS_PER_SEC;
    printf("Total time spent: %f s\n", duration);

	return EXIT_SUCCESS;
}
