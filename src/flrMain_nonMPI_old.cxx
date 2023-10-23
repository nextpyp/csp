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

#include "flrClasses.hxx"
#include "flrMain.hxx"
#include "flrCostFunctions.hxx"
#include "flrReadParametersFile.hxx"

#include "nbfTimer.h"

// #include "local_include/itkFFTShiftImageFilter.h"
// #include "itkComplexToModulusImageFilter.h"
// #include "itkImageLinearIteratorWithIndex.h"
// #include "itkSphereSpatialFunction.h"
// #include "itkFloodFilledSpatialFunctionConditionalIterator.h"

extern SParameters params;
SParameters params = { 16, 7.6 };

int main( int argc, char * argv[] ) {

  // Usage
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " parameters_file" << std::endl;
    return 1;
  }

  // Process the options given in the command line.
  ProcessCommandLineOptions( argc, argv, &params );

  // Create the tilt series where all the data belongs.
  // TODO: create an empty array of CTiltSeries and fill it at the same
  // time is filled the data.
  CTiltSeries tiltSeries;

  // File name string
  char filename[200];
  char MPIStateFilename[200];

  // For OMP
  int threadId;

  // Variables for the optimization loops.
  int loop;
  int loopM;
  int loopP;
  CEulerAngles eA;

  // Read parameters from the configuration file.
  LoadConfiguration( &params, argv[1] );
  ShowConfiguration( &params );

  // Process the parameters, read .par file. Build tilt series and
  // micrograph. Read and attach the particle projections to the
  // micrograph. The created structure is accessed through the object
  // CTiltSeries tiltSeries.
  ProcessConfiguration( &params, &tiltSeries );

  // Show the information from the tilt series.
  tiltSeries.ShowInformation();

  // Get the instance of the CContainer, run the constructor, and
  // associate the tiltseries with it.
  CContainer::Instance();
  CContainer::Instance()->SetTiltSeries( &tiltSeries );

  // // Set the number of available threads
  // omp_set_num_threads( params.mNumberOfThreads );

  // Optimizer definition.
  PowellOptimizerType::Pointer optimizer = PowellOptimizerType::New();

  // Cost function definition.
  flrCostFunction::Pointer costFunction = flrCostFunction::New();
  const unsigned int spaceDimension = costFunction->GetNumberOfParameters();

  // Points in the optimization space.
  typedef flrCostFunction::ParametersType ParametersType;
  ParametersType initialPosition( spaceDimension );
  ParametersType finalPosition( spaceDimension );
  ParametersType dummyPosition( spaceDimension );

  // Initial and final cost.
  flrCostFunction::MeasureType initialCost;
  flrCostFunction::MeasureType finalCost;

  // loopMode defines in which mode is working the framework, it could
  // be OnMicrograph, OnParticle or OnReconstruction.
  LoopModeType loopMode;

  // It's a minimization.
  optimizer->SetMaximize( false );

  // Set the cost function to optimize.
  optimizer->SetCostFunction( costFunction.GetPointer() );

  // Set the configuration parameters. If they are negative the
  // default values are used.
  if (params.mOptimizerStepLength >= 0) optimizer->SetStepLength( params.mOptimizerStepLength );
  if (params.mOptimizerStepTolerance >= 0) optimizer->SetStepTolerance( params.mOptimizerStepTolerance );
  if (params.mOptimizerMaxIter >= 0) optimizer->SetMaximumIteration( params.mOptimizerMaxIter );
  if (params.mOptimizerValueTolerance >= 0) optimizer->SetValueTolerance( params.mOptimizerValueTolerance );
  optimizer->GlobalWarningDisplayOn();
  //  optimizer->Print( std::cout );

  // Add an observer to the optimizer.
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::AnyEvent(), observer );

  ////////////////////////////////////////////////////////////////////
  /// Main loop, running through the micrograph and the particles
  /// optimizing their parameters, and then updating the density map.
  ////////////////////////////////////////////////////////////////////
  for ( loop = 0; loop < params.mNumberOfIterations; loop++ ) {

    PrintDebug( DebugBasic, "\n[ Main loop iteration %d ]\n", 1, loop );

    // First update the lower and upper values for the tilt angle of each
    // micrographs given the separation between them.
    tiltSeries.ComputeMicrographTiltAngleBounds();
    //////////////////////////////////////////////////////////////////
    /// Optimize over the parameters of each of the micrographs (tilt
    /// angle and tilt axis angle). For all the particles with
    /// projection in the micrograph refine the tilt angles.
    //////////////////////////////////////////////////////////////////   
    //#pragma omp parallel for private(threadId) ordered
    for ( loopM = 0; loopM < tiltSeries.GetNumberOfMicrographs(); loopM++ ) {

      //threadId = omp_get_thread_num();
      threadId = 0;
      PrintDebug( DebugBasic, "\n\n### Refining micrograph %d (th:%d) ###\n", 2,
                  tiltSeries.mListOfMicrographs.at( loopM ), threadId );

      // Set the mode of optimization, and the index (number in the
      // parx file) of the current micrograph being optimized.
      loopMode = OnMicrograph;
      costFunction->SetMode( loopMode );
      costFunction->SetCurrentOptimizationIndex( tiltSeries.mListOfMicrographs.at( loopM ) );

      // Set the initial point in the optimization to the current tilt
      // angle and tilt axis angle. The coordinates corresponding to
      // particles will not be affected in the optimization and are
      // different for each particle.
      initialPosition[0] = tiltSeries.mMicrograph.at( loopM )->GetTiltAngle(); // Micrograph tilt angle
      initialPosition[1] = tiltSeries.mMicrograph.at( loopM )->GetTiltAxisAngle(); // Micrograph tilt axis angle
      initialPosition[2] = 0; // Particle theta angle
      initialPosition[3] = 0; // Particle psi angle
      initialPosition[4] = 0; // Particle phi angle
      optimizer->SetInitialPosition( initialPosition );

      // Save the necessary information for the constrained
      // optimization. This will be used inside the cost function
      // (GetValue()) to make a candidate a point feasible or not.
      CContainer::Instance()->SetInitialPoint( initialPosition );
      CContainer::Instance()->SetPreviousPoint( initialPosition );
      initialCost = costFunction->Evaluate( initialPosition );
      CContainer::Instance()->SetInitialValue( initialCost );
      CContainer::Instance()->SetPreviousValue( initialCost );

      // Run the optimization.
      try {
        optimizer->StartOptimization();
      } catch ( itk::ExceptionObject & e ) {
        printf( "[ERROR] Optimizing micrograph %d: %s\n", tiltSeries.mListOfMicrographs.at( loopM ), e.GetDescription() );
        return EXIT_FAILURE;
      }

      // Get the solution.
      finalPosition = optimizer->GetCurrentPosition();
      finalCost = optimizer->GetCurrentCost();

      // Update micrograph's angles.
      tiltSeries.mMicrograph.at( loopM )->SetTiltAngle( finalPosition[0] );
      tiltSeries.mMicrograph.at( loopM )->SetTiltAxisAngle( finalPosition[1] );

      //#pragma omp ordered
      {
        PrintDebug( DebugBasic, "\nRefined from ( 0 %02d %02d ) : f( %f %f %f %f %f ) = %f =>", 8,
                    tiltSeries.mListOfMicrographs.at( loopM ), threadId, initialPosition[0], initialPosition[1],
                    initialPosition[3], initialPosition[2], initialPosition[4], initialCost );
        PrintDebug( DebugBasic, "\nRefined to   ( 0 %02d %02d ) : f( %f %f %f %f %f ) = %f\n", 8,
                    tiltSeries.mListOfMicrographs.at( loopM ), threadId, finalPosition[0], finalPosition[1],
                    finalPosition[3], finalPosition[2], finalPosition[4], finalCost );
      }
		
	  return 0;
		
    } /// end loop on micrographs

	// Create a dummy point to be used for writing .parx files
	for ( unsigned int dk = 0; dk < spaceDimension; dk++ )
	dummyPosition[dk] = 0;
	
	// file for sharing tilt series parameteres with all nodes
	sprintf( MPIStateFilename, "%s/%s/%s_MPI_state.parx", params.mFrealignFolder, params.mBaseName, params.mBaseName );
	// Switch mode to OnReconstruction
	loopMode = OnReconstruction;
	// Write the .parx file so that all slaves can update their tilt-series
	WriteParFile( dummyPosition, loopMode, 0, MPIStateFilename );


    //////////////////////////////////////////////////////////////////
    /// Optimize over the parameters of each particle
    /// (orientation). For each particle refine the orientation using
    /// the constraint.
    //////////////////////////////////////////////////////////////////
    //#pragma omp parallel for private(threadId) ordered
    for ( loopP = 0; loopP < tiltSeries.GetNumberOfParticles(); loopP++ ) {

      //threadId = omp_get_thread_num();
      threadId = 0;
      PrintDebug( DebugBasic, "\n\n### Refining particle %d (th:%d)###\n", 2, tiltSeries.mListOfParticles.at( loopP ),
                  threadId );

      // Set the mode of optimization, and the index of the current
      // particle being optimized.
      loopMode = OnParticle;
      costFunction->SetMode( loopMode );
      costFunction->SetCurrentOptimizationIndex( tiltSeries.mListOfParticles.at( loopP ) );

      // Set the initial point in the optimization to the current
      // euler angles of the particle. The coordinates corresponding
      // to micrograph will not be affected in the optimization and
      // are different for each micrograph.
      initialPosition[0] = 0; // Micrograph tilt angle
      initialPosition[1] = 0; // Micrograph tilt axis angle
      eA = tiltSeries.mParticle.at( loopP )->GetOrientation();
      initialPosition[2] = eA.GetTheta(); // Particle theta angle
      initialPosition[3] = eA.GetPsi(); // Particle psi angle
      initialPosition[4] = eA.GetPhi(); // Particle phi angle
      optimizer->SetInitialPosition( initialPosition );

      // Save the necessary information for the constrained
      // optimization. This will be used inside the cost function
      // (GetValue()) to make a candidate a point feasible or not.
      CContainer::Instance()->SetInitialPoint( initialPosition );
      CContainer::Instance()->SetPreviousPoint( initialPosition );
      initialCost = costFunction->Evaluate( initialPosition );
      CContainer::Instance()->SetInitialValue( initialCost );
      CContainer::Instance()->SetPreviousValue( initialCost );

      // Run the optimization.
      try {
        optimizer->StartOptimization();
      } catch ( itk::ExceptionObject & e ) {
        printf( "[ERROR] Optimizing particle %d: %s\n", tiltSeries.mListOfParticles.at( loopP ), e.GetDescription() );
        return EXIT_FAILURE;
      }

      // Get the solution.
      finalPosition = optimizer->GetCurrentPosition();
      finalCost = optimizer->GetCurrentCost();

      // Update particles's orientation.
      tiltSeries.mParticle.at( loopP )->SetOrientation( finalPosition[2], finalPosition[4], finalPosition[3] );

      //#pragma omp ordered
      {
        PrintDebug( DebugBasic, "\nRefined from ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f =>", 8,
                    tiltSeries.mListOfParticles.at( loopP ), threadId, initialPosition[0], initialPosition[1],
                    initialPosition[3], initialPosition[2], initialPosition[4], initialCost );
        PrintDebug( DebugBasic, "\nRefined to   ( 1 %02d %02d ) : f( %8.3f %8.3f %8.3f %8.3f %8.3f ) = %8.4f\n", 8,
                    tiltSeries.mListOfParticles.at( loopP ), threadId, finalPosition[0], finalPosition[1],
                    finalPosition[3], finalPosition[2], finalPosition[4], finalCost );
      }
    } /// end loop on particles

    //////////////////////////////////////////////////////////////////
    /// Reconstruct the density map with the new orientations
    /// determined in the previous loops. Call Frealign to perform
    /// this reconstruction (write par file mapping the angles and run
    /// Frealign).
    //////////////////////////////////////////////////////////////////
    PrintDebug( DebugBasic, "\n\n### Reconstruct density map %d ###\n", 1, loop );

    // Switch mode to OnReconstruction
    loopMode = OnReconstruction;

    // Create a dummy point with nil coordinates. This is used in the
    // OnReconstruction step after the refinement.
    for ( unsigned int dk = 0; dk < spaceDimension; dk++ )
      dummyPosition[dk] = 0;

    // Write the .parx file for the reconstruction step.
    sprintf( filename, "%s/%s/%s_rec_%07d.parx", params.mFrealignFolder, params.mBaseName, params.mBaseName, loop );
    WriteParFile( dummyPosition, loopMode, 0, filename );

    // Reconstruct the density map with the new parameters.
    ReconstructDensityMap( filename );

    printf("\nThe main loop is ready yet since we are working with few\nimages, then the density map will not be good enough.\n");

  }

  /// Nos vamooooos...
  return EXIT_SUCCESS;
}
