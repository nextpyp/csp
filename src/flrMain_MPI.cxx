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

#include "mpi.h"

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <random/uniform.h>

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

// #include "local_include/itkFFTShiftImageFilter.h"
// #include "itkComplexToModulusImageFilter.h"
// #include "itkImageLinearIteratorWithIndex.h"
// #include "itkSphereSpatialFunction.h"
// #include "itkFloodFilledSpatialFunctionConditionalIterator.h"

#include <vtkTransform.h>

extern SParameters params;
SParameters params = { 16, 7.6 };
Refine3dParameters refine3dParams;
int main( int argc, char * argv[] ) {

	//  Initialize MPI.
	MPI_Init ( &argc, &argv );

	int my_id, num_procs;

	//  Get this processes's rank.
	MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

	//  Find out how many processes are available.
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	if ( my_id == 0 ){
		// Usage
		if (argc < 2) {
			std::cerr << "Usage: " << argv[0] << " parameters_file [iteration] [mode]" << std::endl;
			return 1;
		}

		struct tm *current;
		time_t now;
		time( &now );
		current = localtime( &now );
		printf( "\n %s\n Started on %02i-%02i-%02i at %02i:%02i:%02i.\n", argv[0], current->tm_year + 1900, current->tm_mon + 1, current->tm_mday,
				current->tm_hour, current->tm_min, current->tm_sec );
		cout << " Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
		cout << " Number of processes is " << num_procs << ".\n";
	}

	// cout << "Process " << my_id << " is active.\n";
	// Process the options given in the command line.
	ProcessCommandLineOptions( argc, argv, &params );

	// Create the tilt series where all the data belongs.
	// TODO: create an empty array of CTiltSeries and fill it at the same
	// time is filled the data.
	// CTiltSeries tiltSeries;

	vector< CTiltSeries > dataset;
	
	// File name string
	char filename[200];
	char MPIStateFilename[200];

	// For OMP
	int threadId = 0;

	// Variables for the optimization loops.
	int loop;
	int loopM;
	int loopP;
	CEulerAngles eA;
    EulerAngleType matrix[16]; 

	// Read parameters from the configuration file.
	LoadConfiguration( &params, argv[1] );
	
	// fill pyp, fyp, csp parameters in struct for refine3d
    LoadRefine3dParameters( &refine3dParams, "/scratch/.pyp_config", "/scratch/frealign.config", argv[1], "/scratch/.csp_current_iteration" );
    
    if ( my_id == 0 ) {
		// ShowConfiguration( &params );
	}

	// Number of iteration (if available)
	int starting_iteration = 0;
	if ( argc > 2 ){
		starting_iteration = atoi( argv[2] );
		params.mNumberOfIterations = starting_iteration + 1;
	}
	
	// get the mode of operation
	// 0 - micrograph mode
	// 1 - particles rotation mode
    // 2 - particles translation mode
    // 3 - micrograph translation mode
    // 4 - micrograph defocus refinement mode
	// 5 - both micrographs and particles
	int mode = 5;
	if ( argc > 3 ){
		mode = atoi( argv[3] );
	}
	
	// get what images to use for refinement
	double minScanOrderToUseForRefinement = 1;
	if ( argc > 4 ){
		minScanOrderToUseForRefinement = atof( argv[4] );
	}

	// get what images to use for refinement
	double maxScanOrderToUseForRefinement = 100;
	if ( argc > 5 ){
		maxScanOrderToUseForRefinement = atof( argv[5] );
	}

	double phaseResidualThreshold = 90.0;

	// Process the parameters, read .par file. Build tilt series and
	// micrograph. Read and attach the particle projections to the
	// micrograph. The created structure is accessed through the object
	// CTiltSeries tiltSeries.
    vector<double> weights;
	ProcessConfiguration( &params, dataset, mode, 0, -1, weights );

    // Get the instance of the CContainer, run the constructor, and
	// associate the tiltseries with it.
	CContainer::Instance();
	// CContainer::Instance()->SetTiltSeries( &dataset );
	CContainer::Instance()->SetTiltSeries( &dataset[0] );

	// // Set the number of available threads
	// omp_set_num_threads( params.mNumberOfThreads );

	// Optimizer definition.
	PowellOptimizerType::Pointer optimizer = PowellOptimizerType::New();

	// Cost function definition.
	flrCostFunction::Pointer costFunction = flrCostFunction::New();
	const unsigned int spaceDimension = costFunction->GetNumberOfParameters();
	costFunction->SetPhaseResidualThreshold( phaseResidualThreshold );
	costFunction->SetScanOrderRangeToUseForRefinement( minScanOrderToUseForRefinement, maxScanOrderToUseForRefinement );

	// Points in the optimization space.
	typedef flrCostFunction::ParametersType ParametersType;
	ParametersType initialPosition( spaceDimension );
	ParametersType finalPosition( spaceDimension );
	ParametersType dummyPosition( spaceDimension );

	// Initial and final cost.
	flrCostFunction::MeasureType initialCost;
	flrCostFunction::MeasureType finalCost;

	// loopMode defines in which mode is working the framework, it could
	// be OnMicrograph, OnParticleRotation, OnParticleTranslation or OnReconstruction.
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
	// if ( my_id == 0 ) {
		// optimizer->Print( std::cout );
	// }

	// Add an observer to the optimizer.
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	// optimizer->AddObserver( itk::AnyEvent(), observer );

	int parameterSize = 4;
	int * mpi_params = new int[ parameterSize ]; 	// mpi_params[0] - micrograph or particle
													// mpi_params[1] - tilt series index
													// mpi_params[2] - micrograph or particle index within tilt series
													// mpi_params[3] - cummulative job index
	int fparameterSize = 3;
	float * mpi_fparams = new float[ fparameterSize ];
	int tag_processing = 0;
	int tag_done = 1;
	int tag_update = 2;
	int source;
	
	// Create a dummy point to be used for writing .parx files
	for ( unsigned int dk = 0; dk < spaceDimension; dk++ )
		dummyPosition[dk] = 0;

	int totalNumberOfMicrographs = 0;
	int totalNumberOfParticles = 0;
	for ( int i = 0; i < dataset.size(); i++ ){
		totalNumberOfMicrographs += dataset[i].GetNumberOfMicrographs();
		totalNumberOfParticles += dataset[i].GetNumberOfParticles();
	}
	
	// // file for sharing tilt series parameters with all nodes
	// sprintf( MPIStateFilename, "%s/scratch/%s_MPI_state.bin", params.mFrealignFolder, params.mBaseName );

	// // load state from file
	// if ( mode < 2 ){
		// std :: ifstream src( MPIStateFilename, std :: ofstream :: binary );
		// if ( src.good() ){	// only if file open succesfully
			// stringstream inputString;
			// src >> inputString.rdbuf();
			// Array< double, 1 > A;
			// inputString >> A;
			
			// if ( my_id == 0 ){
				// cout << "\nLoading state from " << MPIStateFilename << endl;
				// // cout << "Loading A =\n" << A << endl;
			// }

			// if ( A.rows() != 2 * totalNumberOfMicrographs + 3 * totalNumberOfParticles ){
				// cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.rows() << " != " << 2 * totalNumberOfMicrographs + 3 * totalNumberOfParticles << endl;
			// }
			// int counter = 0;
			// for ( int t = 0; t < dataset.size(); t++ ){
				// for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
					// dataset[t].mMicrograph.at( i )->SetTiltAngle( A( counter ) );
					// dataset[t].mMicrograph.at( i )->SetTiltAxisAngle( A( counter + 1 ) );
					// counter+=2;
				// }
			// }
			// for ( int t = 0; t < dataset.size(); t++ ){
				// for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
					// dataset[t].mParticle.at( i )->SetOrientation( A( counter + 1 ), A( counter + 2 ), A( counter ) );
					// counter+=3;
				// }
			// }

			// // if ( A.cols() == 2 ){
				// // if ( A.rows() != totalNumberOfMicrographs ){
					// // cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.rows() << " != " << totalNumberOfMicrographs << endl;
				// // }
				// // int counter = 0;
				// // for ( int t = 0; t < dataset.size(); t++ ){
					// // for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
						// // dataset[t].mMicrograph.at( i )->SetTiltAngle( A( counter, 0 ) );
						// // dataset[t].mMicrograph.at( i )->SetTiltAxisAngle( A( counter, 1 ) );
						// // counter++;
					// // }
				// // }
			// // } else {
				// // if ( A.cols() != 3 ){
					// // cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.cols() << " != 3" << endl;
				// // }
				// // if ( A.rows() != totalNumberOfParticles ){
					// // cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.rows() << " != " << totalNumberOfParticles << endl;
				// // }
				// // int counter = 0;
				// // for ( int t = 0; t < dataset.size(); t++ ){
					// // for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
						// // dataset[t].mParticle.at( i )->SetOrientation( A( counter, 0 ), A( counter, 1 ), A( counter, 2 ) );
						// // counter++;
					// // }
				// // }
			// // }
		// } else {
			// if ( my_id == 0 ){
				// // cout << "File " << MPIStateFilename << " not found. Skip loading previous state." << endl;
			// }
		// }
		// src.close();
	// }
	
	// First update the lower and upper values for the tilt angle of each
	// micrographs given the separation between them.
	// tiltSeries.ComputeMicrographTiltAngleBounds();
	for ( int i = 0; i < dataset.size(); i++ ){
		dataset[i].ComputeMicrographTiltAngleBounds(true);
	}

	// Show the information from the tilt series.
	if ( my_id == 0 ) {
		for ( int i = 0; i < dataset.size(); i++ ){
			// dataset[i].ShowInformation();
		}
	}
	
	//// General Initialization
	char command[2000];
	//// sprintf( command, "clearscratch; cp microscope_parameters frealign/frealign_parameters_* parameters.config .parameters.config .csp_current_inner_iteration /scratch");
	//sprintf( command, "clearscratch; cp parameters.config .parameters.config frealign/frealign.config .csp_current_inner_iteration /scratch");
	//// sprintf( command, "clearscratch; cp microscope_parameters frealign/frealign_parameters_* %s .csp_current_inner_iteration /scratch", argv[1]);
	//system( command );
	//sprintf( command, "cp frealign/maps/%s_CSP_%02d.mrc /scratch; echo %d > /scratch/.csp_current_iteration", params.mBaseName, starting_iteration + 1, starting_iteration + 2 );
	//system( command );
		
	if ( my_id == 0 ){ // Master process
		
		// cout << "totalNumberOfMicrographs = " << totalNumberOfMicrographs << endl;
		// cout << "totalNumberOfParticles = " << totalNumberOfParticles << endl;

		cout << "\nRunning in mode " << mode << endl;		
	
		////////////////////////////////////////////////////////////////////
		/// Main loop, running through the micrograph and the particles
		/// optimizing their parameters, and then updating the density map.
		////////////////////////////////////////////////////////////////////
		for ( loop = starting_iteration; loop < params.mNumberOfIterations; loop++ ) {

			// store current iteration value in file
			char command[200];
			sprintf( command, "echo %d > .csp_current_iteration", loop + 2 );
			system( command );
			
			PrintDebug( DebugBasic, "\n[ Main loop iteration %d ]\n", 1, loop );

			// First update the lower and upper values for the tilt angle of each
			// micrographs given the separation between them.
			for ( int i = 0; i < dataset.size(); i++ ){
				dataset[i].ComputeMicrographTiltAngleBounds(false);
			}
			
			int jobNotSubmitted = 0;
			int jobSubmitted = 1;
			int jobDone = 2;

			Array< int, 1 > queueStatus;
			
			// do the micrographs loop first;
			// once all processes done,
			// do the particles loop
			for ( int type = 0; type < 5; type++ ){

				// override with mode selection
				if ( mode != 5 ){
					type = mode;
				}
			
				// initialize status queue and MPI parameters
				switch (type)
				{
				case 0: // micrographs rotation loop
					mpi_params[0] = 0;
					// queueStatus.resize( tiltSeries.GetNumberOfMicrographs() );
					queueStatus.resize( totalNumberOfMicrographs );
					break;
				case 1: // particles rotation loop
                     mpi_params[0] = 1;
                     queueStatus.resize( totalNumberOfParticles );
                     break;
                case 2:
					mpi_params[0] = 2; // particle translation loop
					// queueStatus.resize( tiltSeries.GetNumberOfParticles() );
					queueStatus.resize( totalNumberOfParticles );
					break;
                case 3: 
                    mpi_params[0] = 3; //micrograph translation loop
                    queueStatus.resize( totalNumberOfMicrographs );
                    break;
                case 4:
                    mpi_params[0] = 4; // micrograph defocus loop
                    queueStatus.resize( totalNumberOfMicrographs );
                    break;
				default:
					break;
				}
				queueStatus = jobNotSubmitted;

				// submit first 'num_procs' runs of optimization
				int nextJob = 0;
				int series = 0;
				int indexInSeries = 0;
				for ( int i = 1; i < num_procs; i++ ){
					if ( sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) > 0 ){

						mpi_params[1] = series;
						mpi_params[2] = indexInSeries;
						mpi_params[3] = nextJob;
						
						// cout << "MASTER: sending request for process = " << i << endl;
						MPI_Send ( mpi_params, parameterSize, MPI_INT, i, tag_processing, MPI_COMM_WORLD );
						
						// update queue status
						queueStatus( nextJob ) = jobSubmitted;
						nextJob++;
						
						indexInSeries++;
						
						switch (type)
						{
						case 0:
                        case 3:
                        case 4:
							if ( dataset[ series ].GetNumberOfMicrographs() == indexInSeries ){
								series++; indexInSeries = 0;
							}
							break;
						case 1:
                        case 2:
							if ( dataset[ series ].GetNumberOfParticles() == indexInSeries ){
								series++; indexInSeries = 0;
							}
							break;
						}
						
					} else {
						// all jobs have been submitted
						break;
					}
				}

				// submit remaining jobs
				while ( sum( where( queueStatus != jobDone, 1, 0 ) ) > 0 ){

					MPI_Status status;

					// cout << "MASTER: waiting for results from slaves ... " << endl;
					
					// receive notification from slaves
					MPI_Recv( mpi_params, parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

					int source = status.MPI_SOURCE;
					
					// retrieve object size and tag information
					MPI_Recv( mpi_fparams, fparameterSize, MPI_FLOAT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );		
					
					//cout << "MASTER: received results from slave " << source << endl;

					int tiltseriesIndex = mpi_params[1];
					if ( mpi_params[0] == 0 ){
						// Update micrograph's angles.
						int loopM = mpi_params[2];
						dataset[tiltseriesIndex].mMicrograph.at( loopM )->SetTiltAngle( dataset[tiltseriesIndex].mMicrograph.at( loopM )->GetTiltAngle() + mpi_fparams[0] );
						dataset[tiltseriesIndex].mMicrograph.at( loopM )->SetTiltAxisAngle( dataset[tiltseriesIndex].mMicrograph.at( loopM )->GetTiltAxisAngle() + mpi_fparams[1] );
						//printf("Received Tilt Series %d, micrograph %d, [ %f %f 0 0 0 ]", tiltseriesIndex, loopM, mpi_fparams[0], mpi_fparams[1] );
					} else if ( mpi_params[0] == 1 ) {
						// Update particles's orientation.
						int loopP = mpi_params[2];
						//dataset[tiltseriesIndex].mParticle.at( loopP )->SetOrientation( dataset[tiltseriesIndex].mParticle.at( loopP )->GetOrientation().GetTheta() + mpi_fparams[0], 
						//                                                                dataset[tiltseriesIndex].mParticle.at( loopP )->GetOrientation().GetPhi() + mpi_fparams[1],
						//								dataset[tiltseriesIndex].mParticle.at( loopP )->GetOrientation().GetPsi() + mpi_fparams[2] );
						dataset[tiltseriesIndex].mParticle.at( loopP )->SetOrientation( mpi_fparams[0], mpi_fparams[1], mpi_fparams[2] );
						//printf("Received Tilt Series %d, micrograph %d, [ 0 0 %f %f %f ]", tiltseriesIndex, loopP, mpi_fparams[0], mpi_fparams[1], mpi_fparams[2] );
					}
                    // Update particle's translation in 3D
                    else if ( mpi_params[0] == 2 ) {
                        int loopP = mpi_params[2];
                        dataset[tiltseriesIndex].mParticle.at( loopP )->SetMatrixTranslation( mpi_fparams[0], mpi_fparams[2], mpi_fparams[1] ); 
                    }
                    // Update micrograph translation
                    else if ( mpi_params[0] == 3 ) {
                        int loopM = mpi_params[2];
                        CPosition2D micrographTiltAxisShift;
                        micrographTiltAxisShift.SetPosition( dataset[ tiltseriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordX + mpi_fparams[0], dataset[ tiltseriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordY + mpi_fparams[1] );
                        dataset[tiltseriesIndex].mMicrograph.at( loopM )->SetTiltAxisShift( micrographTiltAxisShift );
                    }
                    // Update micrograph defocus offset and astigmatism
                    else if (  mpi_params[0] == 4 ) {
                        int loopM = mpi_params[2];
                        dataset[tiltseriesIndex].mMicrograph.at( loopM )->SetTiltDefocusOffset( mpi_fparams[0] );
                        dataset[tiltseriesIndex].mMicrograph.at( loopM )->SetTiltAstigmatism( mpi_fparams[1] );
                    }
					
					// update queue status
					queueStatus( mpi_params[3] ) = jobDone;

					// submit remaining jobs (if any)
					if ( nextJob < queueStatus.size() ){
						switch (type)
						{
							case 0:
								// refine micrograph
								mpi_params[0] = 0;
								break;
							case 1:
								// refine particle
								mpi_params[0] = 1;
								break;
                            case 2:
                                mpi_params[0] = 2;
                                break;
                            case 3: 
                                mpi_params[0] = 3;
                                break;
                            case 4:
                                mpi_params[0] = 4;
                                break;
						}

						mpi_params[1] = series;
						mpi_params[2] = indexInSeries;
						mpi_params[3] = nextJob;
						
						// cout << "MASTER: sending request for process = " << source << endl;

						// send parameters
						MPI_Send ( mpi_params, parameterSize, MPI_INT, source, tag_processing, MPI_COMM_WORLD );

						// update queue status
						queueStatus( nextJob ) = jobSubmitted;
						nextJob++;				

						indexInSeries++;
						
						switch (type)
						{
						case 0:
                        case 3:
                        case 4:
							if ( dataset[ series ].GetNumberOfMicrographs() == indexInSeries ){
								series++; indexInSeries = 0;
							}
							break;
						case 1:
                        case 2:
							if ( dataset[ series ].GetNumberOfParticles() == indexInSeries ){
								series++; indexInSeries = 0;
							}
							break;
						}
						
					}
				}
				
				// // diseminate new tilt series parameters to all nodes through file

				// Array< double, 1 > A;

				// A.resize( 2 * totalNumberOfMicrographs + 3 * totalNumberOfParticles );
				// int counter = 0;
				// for ( int t = 0; t < dataset.size(); t++ ){
					// for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
						// A( counter ) = dataset[t].mMicrograph.at( i )->GetTiltAngle();
						// A( counter + 1 ) = dataset[t].mMicrograph.at( i )->GetTiltAxisAngle();
						// counter+=2;
					// }
				// }
				// for ( int t = 0; t < dataset.size(); t++ ){
					// for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
						// A( counter ) = dataset[t].mParticle.at( i )->GetOrientation().GetPsi();
						// A( counter + 1 ) = dataset[t].mParticle.at( i )->GetOrientation().GetTheta();
						// A( counter + 2 ) = dataset[t].mParticle.at( i )->GetOrientation().GetPhi();
						// counter+=3;
					// }
				// }

				// switch ( type ){
					// case 0:
					// {
						// // update micrograph parameters
						// A.resize( totalNumberOfMicrographs, 2 );
						// int counter = 0;
						// for ( int t = 0; t < dataset.size(); t++ ){
							// for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
								// A( counter, 0 ) = dataset[t].mMicrograph.at( i )->GetTiltAngle();
								// A( counter, 1 ) = dataset[t].mMicrograph.at( i )->GetTiltAxisAngle();
								// counter++;
							// }
						// }
						// break;
					// }
					// case 1:
					// {
						// // update particle parameters
						// A.resize( totalNumberOfParticles, 3 );
						// int counter = 0;
						// for ( int t = 0; t < dataset.size(); t++ ){
							// for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
								// A( counter, 0 ) = dataset[t].mParticle.at( i )->GetOrientation().GetPsi();
								// A( counter, 1 ) = dataset[t].mParticle.at( i )->GetOrientation().GetTheta();
								// A( counter, 2 ) = dataset[t].mParticle.at( i )->GetOrientation().GetPhi();
								// counter++;
							// }
						// }
						// break;
					// }
				// }
				
				//cout << "Saving A =\n" << A << endl;

				// // dump data to file
				// std :: ofstream dst( MPIStateFilename, std :: ofstream :: binary );
				// if ( dst.good() ){
					// dst << A;
					// dst.close();
					// cout << "\nSaving state in file " << MPIStateFilename << endl;
				// } else {
					// cout << "ERROR - Cannot open file " << MPIStateFilename << endl;
					// break;
				// }

				// if in mode 0 or 1, we are done
				if ( mode != 5 ){
					// cout << "We're done with mode " << mode << ". Breaking loop.." << endl;
					break;
				}
				
				mpi_params[0] = type;
				
				// notify slaves to update
				for ( int i = 1; i < num_procs; i++ ){
					MPI_Send ( mpi_params, parameterSize, MPI_INT, i, tag_update, MPI_COMM_WORLD );
				}
			}
			
			if ( ( dataset.size() == 1 ) && ( mode == 5 ) ) { // do the reconstruction within same process
				
				//////////////////////////////////////////////////////////////////
				/// Reconstruct the density map with the new orientations
				/// determined in the previous loops. Call Frealign to perform
				/// this reconstruction (write par file mapping the angles and run
				/// Frealign).
				//////////////////////////////////////////////////////////////////
				PrintDebug( DebugBasic, "\n\n### Reconstruct density map %d ###\n", 1, loop );

				// Write the .parx file for the reconstruction step.
				sprintf( filename, "frealign/scratch/%s_CSP_%02d.parx", params.mBaseName, loop + 2 );

				loopMode = OnReconstruction;
				WriteParFile( dummyPosition, loopMode, 0, filename );
				// Reconstruct the density map with the new parameters.
				ReconstructDensityMap( filename );

				// retrieve new PR and PRES values
				sprintf( filename, "frealign/maps/%s_CSP_%02d.parx", params.mBaseName, loop + 2 );
			
				cout << "\nRetrieving PR and DPRES from " << filename << endl;
				ReadParFile( loopMode, 0, filename );
				
			} else {

				// if ( mode != 0 ){
					// Write the .parx file for the reconstruction step.
					//cout << "Write the .parx file for the reconstruction step." << endl;
					sprintf( filename, "frealign/%s.series", params.mBaseName );
					//WriteMultipleParFile( dataset, loop + 1, filename );

					// Write the .parx file for the reconstruction step.
					char outfilename[200];
					sprintf( outfilename, "frealign/maps/%s_CSP_%02d.parx", params.mBaseName, loop + 1 );
					

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
                            loopMode = OnMicrographDefocus;
                            break;
                        case 4:
                            loopMode = OnMicrographDefocus;
                            break;
                        default:
                            break; 
                    }
                    WriteMultipleParFileNew( dataset, loop + 1, filename, outfilename, loopMode );

				// }
			}

			// printf("\nMASTER: The main loop is ready yet since we are working with few\nimages, then the density map will not be good enough.\n");
		}

		// send termination message to all nodes
		for ( int i = 1; i < num_procs; i++ ){
			// cout << "MASTER - sending termination signal to proc " << i << endl;
			MPI_Send ( mpi_params, parameterSize, MPI_INT, i, tag_done, MPI_COMM_WORLD );
		}
		
		// if ( mode == 2 ){
			// // delete state file
			// char command[200];
			// sprintf( command, "rm -f %s .csp_current_iteration", MPIStateFilename);
			// system( command );
		// }

	} else {
		// slave processes
		
		int tag = tag_processing;
		int source = 0;

		while ( tag != tag_done ){

			// cout << "SLAVE " << my_id << " waiting communication from master." << endl;
			
			// retrieve object size and tag information
			MPI_Status status;
			MPI_Recv( mpi_params, parameterSize, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );		
			int tag = status.MPI_TAG;

			if ( tag == tag_done ){
				// cout << "SLAVE " << my_id << " received DONE signal." << endl;
				break;
			}
			
			// if ( tag == tag_update ){
			
				// std :: ifstream src( MPIStateFilename, std :: ofstream :: binary );
				// stringstream inputString;
				// src >> inputString.rdbuf();
				// src.close();
				// Array< double, 1 > A;
				// inputString >> A;
				
				// // if ( my_id == 1 ){
					// // cout << "Loading A =\n" << A << endl;
				// // }

				// if ( A.rows() != 2 * totalNumberOfMicrographs + 3 * totalNumberOfParticles ){
				// }
				// int counter = 0;
				// for ( int t = 0; t < dataset.size(); t++ ){
					// for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
						// dataset[t].mMicrograph.at( i )->SetTiltAngle( A( counter ) );
						// dataset[t].mMicrograph.at( i )->SetTiltAxisAngle( A( counter + 1 ) );
						// counter+=2;
					// }
				// }
				// for ( int t = 0; t < dataset.size(); t++ ){
					// for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
						// dataset[t].mParticle.at( i )->SetOrientation( A( counter + 1 ), A( counter + 2 ), A( counter ) );
						// counter+=3;
					// }
				// }
				
				// // if ( mpi_params[0] == 0 ){
					// // if ( A.rows() != totalNumberOfMicrographs ){
						// // cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.rows() << " != " << totalNumberOfMicrographs << endl;
					// // }
					// // int counter = 0;
					// // for ( int t = 0; t < dataset.size(); t++ ){
						// // for ( int i = 0; i < dataset[t].GetNumberOfMicrographs(); i++ ){
							// // dataset[t].mMicrograph.at( i )->SetTiltAngle( A( counter, 0 ) );
							// // dataset[t].mMicrograph.at( i )->SetTiltAxisAngle( A( counter, 1 ) );
							// // counter++;
						// // }
					// // }
				// // } else {
					// // if ( A.rows() != totalNumberOfParticles ){
						// // cout << "ERROR - Reading state file " << MPIStateFilename << ". Dimensions do not match: " << A.rows() << " != " << totalNumberOfParticles << endl;
					// // }
					// // int counter = 0;
					// // for ( int t = 0; t < dataset.size(); t++ ){
						// // for ( int i = 0; i < dataset[t].GetNumberOfParticles(); i++ ){
							// // dataset[t].mParticle.at( i )->SetOrientation( A( counter, 1 ), A( counter, 2 ), A( counter, 0 ) );
							// // counter++;
						// // }
					// // }
				// // }

				// // Update the lower and upper values for the tilt angle of each
				// // micrographs given the separation between them.
				// for ( int t = 0; t < dataset.size(); t++ ){
					// dataset[t].ComputeMicrographTiltAngleBounds();
				// }

				// continue;
			// }
			
			// cout << "SLAVE " << my_id << " received request from master for MODE = " << mpi_params[0] << ", LOOP = " << mpi_params[1] << endl;

			int tiltSeriesIndex = mpi_params[1];
				
			CContainer::Instance()->SetTiltSeries( &dataset[tiltSeriesIndex] );
			
			if ( mpi_params[0] == 0 || mpi_params[0] == 3 || mpi_params[0] == 4 ){
		
				// loop in micrographs
				int loopM = mpi_params[2];

				// Set the mode of optimization, and the index (number in the
				// parx file) of the current micrograph being optimized.
				if (  mpi_params[0] == 0 ) {
                    loopMode = OnMicrographRotation;
                }
                else if ( mpi_params[0] == 3 ) {
                    loopMode = OnMicrographTranslation;
                }
                else if ( mpi_params[0] == 4 ) {
                    loopMode = OnMicrographDefocus;
                }
                
				costFunction->SetMode( loopMode );
				costFunction->SetCurrentOptimizationIndex( dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ) );

				// Set the initial point in the optimization to the current tilt
				// angle and tilt axis angle. The coordinates corresponding to
				// particles will not be affected in the optimization and are
				// different for each particle.
				initialPosition[0] = dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngle(); // Micrograph tilt angle
				initialPosition[1] = dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisAngle(); // Micrograph tilt axis angle
				initialPosition[0] = 0; // Particle theta angle
				initialPosition[1] = 0; // Particle theta angle
				initialPosition[2] = 0; // Particle theta angle
				initialPosition[3] = 0; // Particle psi angle
				initialPosition[4] = 0; // Particle phi angle
				optimizer->SetInitialPosition( initialPosition );

				// Save the necessary information for the constrained
				// optimization. This will be used inside the cost function
				// (GetValue()) to make a candidate a point feasible or not.
				CContainer::Instance()->mNumberOfGetValueCalled = 1;

				CContainer::Instance()->SetInitialPoint( initialPosition );
				CContainer::Instance()->SetPreviousPoint( initialPosition );
				initialCost = costFunction->Evaluate( initialPosition );
				CContainer::Instance()->SetInitialValue( initialCost );
				CContainer::Instance()->SetPreviousValue( initialCost );
                
                switch ( loopMode ) {
                    case OnMicrographRotation: {
				        // compute the Powell step based on half the distance to the closest neighboring micrograph
				        float minb = fabs( dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngle() - dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngleLowerBound() );
				        float maxb = fabs( dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngle() - dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngleUpperBound() );
				        float step = min( minb, maxb ) / 2;
				        optimizer->SetStepLength( step );
				        // cout << "Setting step-lenght for micrograph " <<  loopM << " refinement to " << step << endl;
					
				        if ( fabs( CContainer::Instance()->mTolerance[0] ) + fabs( CContainer::Instance()->mTolerance[1] ) > 0 ) {
					        // Run the optimization.
					        try {
						        optimizer->StartOptimization();
					        } catch ( itk::ExceptionObject & e ) {
						        printf( "[ERROR] Optimizing micrograph %d: %s\n", dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ), e.GetDescription() );
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
				//#pragma omp ordered
				        
					    PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
					    tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ) );
					    PrintDebug( DebugBasic, "\nRotation refined from : f( Tilt = %8.4f, Axis = %8.4f ) = %8.4f =>", 3,
					    dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngle() + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisAngle() + initialPosition[1], initialCost );
					    PrintDebug( DebugBasic, "\nRotation refined to   : f( Tilt = %8.4f, Axis = %8.4f ) = %8.4f\n", 3,
					    dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAngle() + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisAngle() + finalPosition[1], finalCost );
				        
                        break;
                    }
                    case OnMicrographTranslation: {
                        try {
					        optimizer->StartOptimization();
				        } catch ( itk::ExceptionObject & e ) {
					        printf( "[ERROR] Optimizing micrograph %d: %s\n", dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ), e.GetDescription() );
					        return EXIT_FAILURE;
				        }
                  
                        finalPosition = optimizer->GetCurrentPosition();
                        finalCost = optimizer->GetCurrentCost();
                        
                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
                        tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ) );
                        PrintDebug( DebugBasic, "\nTranslation refined from : f( X = %8.4f, Y = %8.4f ) = %8.4f =>", 3,
                        dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordX + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordY + initialPosition[1], initialCost );
                        PrintDebug( DebugBasic, "\nTranslation refined to   : f( X = %8.4f, Y = %8.4f ) = %8.4f\n", 3,
                        dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordX + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAxisShift().mCoordY + finalPosition[1], finalCost );
                        break;
                    }
                    case OnMicrographDefocus: {
                        try {
                            optimizer->StartOptimization();
                        } catch ( itk::ExceptionObject & e ) {
                            printf( "[ERROR] Optimizing micrograph %d: %s\n", dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ), e.GetDescription() );
                            return EXIT_FAILURE;
                        }

                        finalPosition = optimizer->GetCurrentPosition();
                        finalCost = optimizer->GetCurrentCost();

                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, micrograph %d ###\n", 3,
                        tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfMicrographs.at( loopM ) );
                        PrintDebug( DebugBasic, "\nDefocus refined from : f( Offset = %8.4f, Astig = %8.4f ) = %8.4f =>", 3,
                        dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltDefocusOffset() + initialPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAstigmatism() + initialPosition[1], initialCost );
                        PrintDebug( DebugBasic, "\nDefocus refined to   : f( Offset = %8.4f, Astig = %8.4f ) = %8.4f\n", 3,
                        dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltDefocusOffset() + finalPosition[0], dataset[ tiltSeriesIndex ].mMicrograph.at( loopM )->GetTiltAstigmatism() + finalPosition[1], finalCost );
                        break;
                    }
                }
                // Update micrograph's angles or shifts.
                mpi_fparams[0] = finalPosition[0];
                mpi_fparams[1] = finalPosition[1];

			} else {
				
				// loop in particles
				int loopP = mpi_params[2];

				// Set the mode of optimization, and the index of the current
				// particle being optimized.
                if ( mpi_params[0] == 1 ) {
				    loopMode = OnParticleRotation;
                }
                else {
                    loopMode = OnParticleTranslation;
                }

				costFunction->SetMode( loopMode );
				costFunction->SetCurrentOptimizationIndex( dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
                costFunction->SetWeightList( weights );
				// Set the initial point in the optimization to the current
				// euler angles of the particle. The coordinates corresponding
				// to micrograph will not be affected in the optimization and
				// are different for each micrograph.
				initialPosition[0] = 0; // Micrograph tilt angle
				initialPosition[1] = 0; // Micrograph tilt axis angle
				eA = dataset[ tiltSeriesIndex ].mParticle.at( loopP )->GetOrientation();
			    dataset[ tiltSeriesIndex ].mParticle.at( loopP )->GetMatrix( matrix );	
                // initialPosition[2] = eA.GetTheta(); // Particle theta angle
				// initialPosition[3] = eA.GetPsi(); // Particle psi angle
				// initialPosition[4] = eA.GetPhi(); // Particle phi angle
				initialPosition[2] = 0; // Particle theta angle
				initialPosition[3] = 0; // Particle theta angle
				initialPosition[4] = 0; // Particle theta angle

				if ( fabs( CContainer::Instance()->mTolerance[2] ) > 0 ){

					// randomized initial position search
					flrCostFunction::MeasureType bestCost = costFunction->Evaluate( initialPosition );
					ParametersType bestInitialPosition( spaceDimension ), currentPosition( spaceDimension );
					bestInitialPosition[0] = initialPosition[0];
					bestInitialPosition[1] = initialPosition[1];
					bestInitialPosition[2] = initialPosition[2];
					bestInitialPosition[3] = initialPosition[3];
					bestInitialPosition[4] = initialPosition[4];
	
					// cout << "Starting from f(" << eA.GetTheta() << "," << eA.GetPsi() << "," << eA.GetPhi() << ") = " << bestCost << endl;
	
					currentPosition[0] = initialPosition[0];
					currentPosition[1] = initialPosition[1];
					ranlib::Uniform<flrCostFunction::MeasureType> rnd;
					rnd.seed((unsigned int)time(0));
					
					/**
					//int random_points = 50;
					//if ( starting_iteration == 0 ){
					//	random_points = 100;
					//}
					//for ( int i = 0; i < random_points; i++ ){
				    if ( loopMode == OnParticleRotation ) {	
					    for ( int i = 0; i < params.mNumberOfRandomIterations; i++ ){
						    //currentPosition[2] = fabs( CContainer::Instance()->mTolerance[2] ) * 2 * ( rnd.random() -.5 );
						    currentPosition[2] = 360 * ( rnd.random() - .5 );
						    //currentPosition[3] = fabs( CContainer::Instance()->mTolerance[3] ) * 2 * ( rnd.random() -.5 );
						    currentPosition[3] = 360 * ( rnd.random() - .5 );
						    //currentPosition[4] = fabs( CContainer::Instance()->mTolerance[4] ) * 2 * ( rnd.random() -.5 );
						    currentPosition[4] = 360 * ( rnd.random() - .5 );
						    //currentPosition[2] =  2 * ( rnd.random() -.5 );
						    //currentPosition[3] =  2 * ( rnd.random() -.5 );
						    //currentPosition[4] =  2 * ( rnd.random() -.5 );
						    //if ( currentPosition[2] > 0 ){ // Generate point in [-90,90] range
						    //	currentPosition[2] = eA.GetTheta() + currentPosition[2] * ( 90 - eA.GetTheta() );
						    //} else {
						    //	currentPosition[2] = eA.GetTheta() + currentPosition[2] * ( 90 + eA.GetTheta() );
						    //}
						    //if ( currentPosition[3] > 0 ){ // Generate point in [-180,180] range
						    //	currentPosition[3] = eA.GetPsi() + currentPosition[3] * ( 180 - eA.GetPsi() );
						    //} else {
						    //	currentPosition[3] = eA.GetPsi() + currentPosition[3] * ( 180 + eA.GetPsi() );
						    //}
						    //if ( currentPosition[4] > 0 ){ // Generate point in [-180,180] range
						    //	currentPosition[4] = eA.GetPhi() + currentPosition[4] * ( 180 - eA.GetPhi() );
						    //} else {
						    //	currentPosition[4] = eA.GetPhi() + currentPosition[4] * ( 180 + eA.GetPhi() );
						    //}
						    flrCostFunction::MeasureType currentCost = costFunction->Evaluate( currentPosition );
						    // flrCostFunction::MeasureType currentCost = costFunction->GetValue( currentPosition );
						    // cout << i << ": f(" << currentPosition[2] << "," << currentPosition[3] << "," << currentPosition[4] << ") = " << currentCost << endl;
						    if ( currentCost > bestCost ){
							    bestCost = currentCost;
							    bestInitialPosition[2] = currentPosition[2];
							    bestInitialPosition[3] = currentPosition[3];
							    bestInitialPosition[4] = currentPosition[4];
							    // cout << i << "\t: found better starting point at f(" << bestInitialPosition[2] << "," << bestInitialPosition[3] << "," << bestInitialPosition[4] << ") = " << bestCost << endl;
						    }
					    }

                    }
					initialPosition[2] = bestInitialPosition[2];
					initialPosition[3] = bestInitialPosition[3];
					initialPosition[4] = bestInitialPosition[4];
					*/

                    /** 
                    if ( loopMode == OnParticleRotation ) {
					    double STEP = 20;
					    for ( double i = 0; i < 360; i+= STEP ) {
						    for ( double j = 0; j < 360; j+= STEP ) {
							    for ( double k = 0; k < 360; k += STEP ) {
								    currentPosition[2] = i;
								    currentPosition[3] = j;
								    currentPosition[4] = k;

								    flrCostFunction::MeasureType currentCost = costFunction->Evaluate( currentPosition );
								    // printf("( %f, %f, %f ) => %f\n", currentPosition[2], currentPosition[3], currentPosition[4], currentCost);
								    if ( currentCost > bestCost ){
							    	    bestCost = currentCost;
							    	    bestInitialPosition[2] = currentPosition[2];
							    	    bestInitialPosition[3] = currentPosition[3];
							    	    bestInitialPosition[4] = currentPosition[4];
                                    }
                                }
							}
						}
					
					    initialPosition[2] = bestInitialPosition[2];
					    initialPosition[3] = bestInitialPosition[3];
					    initialPosition[4] = bestInitialPosition[4];
                    }
                    */
				}
				
				optimizer->SetInitialPosition( initialPosition );
				// Save the necessary information for the constrained
				// optimization. This will be used inside the cost function
				// (GetValue()) to make a candidate a point feasible or not.
				CContainer::Instance()->mNumberOfGetValueCalled = 1;

				CContainer::Instance()->SetInitialPoint( initialPosition );
				CContainer::Instance()->SetPreviousPoint( initialPosition );
				initialCost = costFunction->Evaluate( initialPosition );
				CContainer::Instance()->SetInitialValue( initialCost );
				CContainer::Instance()->SetPreviousValue( initialCost );

				// Run the optimization.
				try {
					optimizer->StartOptimization();
				} catch ( itk::ExceptionObject & e ) {
					printf( "[ERROR] Optimizing particle %d: %s\n", dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ), e.GetDescription() );
					return EXIT_FAILURE;
				}

				// optimizer->Print( std::cout );

				// Get the solution.
				finalPosition = optimizer->GetCurrentPosition();
				finalCost = optimizer->GetCurrentCost();

				// make final results modulo 360
                
				switch ( loopMode ) {
                    case OnParticleRotation: {
                        finalPosition[2] += eA.GetTheta();
				        finalPosition[3] += eA.GetPsi();
				        finalPosition[4] += eA.GetPhi();

				        if ( finalPosition[3] < 0 ){
				            finalPosition[3] = fmod(finalPosition[3],360.0)+360.0;
				        } else {
				            finalPosition[3] = fmod(finalPosition[3],360.0);
				        }
				        if ( finalPosition[2] < 0 ){
				            finalPosition[2] = fmod(finalPosition[2],360.0)+360.0;
				        } else {
				            finalPosition[2] = fmod(finalPosition[2],360.0);
				        }
				        if ( finalPosition[4] < 0 ){
				            finalPosition[4] = fmod(finalPosition[4],360.0)+360.0;
				        } else {
				            finalPosition[4] = fmod(finalPosition[4],360.0);
				        }

                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
                        PrintDebug( DebugBasic, "\nRotation refined from: f( PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f =>", 4,
                        eA.GetPsi() + initialPosition[3], eA.GetTheta() + initialPosition[2], eA.GetPhi() + initialPosition[4], initialCost );
                        PrintDebug( DebugBasic, "\nRotation refined to  : f( PSI = %8.3f, THETA = %8.3f, PHI = %8.3f ) = %8.4f\n", 4,
                        finalPosition[3], finalPosition[2], finalPosition[4], finalCost );
                        break;
                    }
                    case OnParticleTranslation: {
                        
                        finalPosition[2] += matrix[3];
                        finalPosition[3] += matrix[7];
                        finalPosition[4] += matrix[11];
                        
                        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
                        PrintDebug( DebugBasic, "\nTranslation refined from: f( X = %8.3f, Y = %8.3f, Z = %8.3f ) = %8.4f =>", 4,
                        matrix[3] + initialPosition[2], matrix[7] + initialPosition[3], matrix[11] + initialPosition[4], initialCost );
                        PrintDebug( DebugBasic, "\nTranslation refined to  : f( X = %8.3f, Y = %8.3f, Z = %8.3f ) = %8.4f\n", 4,
                        finalPosition[2], finalPosition[3], finalPosition[4], finalCost );
                        
                        break;
                    }
                }
				//#pragma omp ordered
				{
			//        PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
			//		PrintDebug( DebugBasic, "\nRefined from ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f =>", 8,
			//		dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ), threadId, initialPosition[0], initialPosition[1],
			//		eA.GetPsi() + initialPosition[3], eA.GetTheta() + initialPosition[2], eA.GetPhi() + initialPosition[4], initialCost );
			//		PrintDebug( DebugBasic, "\nRefined to   ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f\n", 8,
			//		dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ), threadId, finalPosition[0], finalPosition[1],
			//		eA.GetPsi() + finalPosition[3], eA.GetTheta() + finalPosition[2], eA.GetPhi() + finalPosition[4], finalCost );
			        //PrintDebug( DebugBasic, "\n\n### Refining tilt series %d, particle %d ###\n", 3, tiltSeriesIndex, dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ) );
					//PrintDebug( DebugBasic, "\nRefined from ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f =>", 8,
					//dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ), threadId, initialPosition[0], initialPosition[1],
					//eA.GetPsi() + initialPosition[3], eA.GetTheta() + initialPosition[2], eA.GetPhi() + initialPosition[4], initialCost );
					//PrintDebug( DebugBasic, "\nRefined to   ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f\n", 8,
					//dataset[ tiltSeriesIndex ].mListOfParticles.at( loopP ), threadId, finalPosition[0], finalPosition[1],
					//finalPosition[3], finalPosition[2], finalPosition[4], finalCost );
				}	
				
				// Update particles's orientation.
				mpi_fparams[0] = finalPosition[2];
				mpi_fparams[1] = finalPosition[4];
				mpi_fparams[2] = finalPosition[3];
			}

			// cout << "SLAVE " << my_id << " communicating result to master for MODE = " << mpi_params[0] << ", LOOP = " << mpi_params[1] << endl;

			MPI_Send ( mpi_params, parameterSize, MPI_INT, 0, tag_done, MPI_COMM_WORLD );
			MPI_Send ( mpi_fparams, fparameterSize, MPI_FLOAT, 0, tag_done, MPI_COMM_WORLD );
		}
	}
	
	delete [] mpi_params;
	delete [] mpi_fparams;

	// clean up scratch
	// sprintf( command, "clearscratch");
	// system( command );

	MPI_Finalize();

	if ( my_id == 0 ){
		cout << "\nNormal program termination." << endl;
	}

	/// Nos vamooooos...
	return EXIT_SUCCESS;
}
