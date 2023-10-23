/*=========================================================================

 Name:        flrCostFunctions.cxx
 Category:    Implementation file
 Description:
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrCostFunctions.cxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "flrClasses.hxx"
#include "flrCostFunctions.hxx"
#include "flrReadParametersFile.hxx"

#include "nbfTimer.h"

extern SParameters params;

flrCostFunction::flrCostFunction() {
  /** B-factor */
  mRBFactor = params.mBFactor;

  // Set mode
  SetMode( OnMicrograph );

  /** Temporal, then should be defined with the constructor. */
  ImageType::SizeType inputSize;
  inputSize[0] = params.mImageSizeX;
  inputSize[1] = params.mImageSizeY;

  /** Fill the ringsRadius array. */
  for ( unsigned int rC = 0; rC < params.mNumberOfRings; rC++ ) {
    mRingsWeigths.push_back( 1 );
    mRingsRadius.push_back( rC * inputSize[0] / 2 / params.mNumberOfRings );
  }
  mRingsRadius.push_back( inputSize[0] ); // One more radius bound, the outermost.

  this->phaseResidualThreshold = 90.0;
}

flrCostFunction::flrCostFunction( WeightArrayType ringsWeigths ) {

  /** B-factor */
  mRBFactor = params.mBFactor;

  SetMode( OnMicrograph );

  /** Temporal, then should be defined with the constructor. */
  ImageType::SizeType inputSize;
  inputSize[0] = 128;
  inputSize[1] = 128;

  /** Fill the ringsRadius array. */
  for ( unsigned int rC = 0; rC < params.mNumberOfRings; rC++ ) {
    mRingsWeigths.push_back( ringsWeigths.at( rC ) );
    mRingsRadius.push_back( rC * inputSize[0] / 2 / params.mNumberOfRings );
  }
  mRingsRadius.push_back( inputSize[0] ); // One more radius bound, the outermost.

  this->phaseResidualThreshold = 90.0;
}

flrCostFunction::~flrCostFunction() {
}

void flrCostFunction::GetDerivative( const ParametersType &, DerivativeType & ) const {
}

flrCostFunction::MeasureType flrCostFunction::GetValue( const ParametersType & parameters ) const {

  MeasureType result;
  MeasureType notFeasibleResult = itk::NumericTraits<double>::infinity();
  ParametersType previousPoint = CContainer::Instance()->GetPreviousPoint();

  CContainer::Instance()->SetCurrentPoint( parameters );

  struct tm *current;
  time_t now;
  time( &now );
  current = localtime( &now );

  //  // Find the descent direction, given by the position where the previous
  //  // and current points differ. Remember that only one dimension is
  //  // optimized per iteration.
  //  descentDirection = CContainer::Instance()->GetDescentDirection();

   // printf( "\n[%02i-%02i-%02i@%02i:%02i:%02i] %7d: ", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday,
           // current->tm_hour, current->tm_min, current->tm_sec, CContainer::Instance()->mNumberOfGetValueCalled );
   // PrintParameters( CContainer::Instance()->GetCurrentPoint() );
   // cout << endl;

  // PrintDebug( DebugInfo, "\n= GetValue() call #%d", 1, CContainer::Instance()->mNumberOfGetValueCalled );
  // if (DebugOrNotDebug( DebugData )) {
    // PrintDebug( DebugData, "\n  Current point:  " );
    // PrintParameters( CContainer::Instance()->GetCurrentPoint() );
    // PrintDebug( DebugData, "\n  Previous point: " );
    // PrintParameters( CContainer::Instance()->GetPreviousPoint() );
    // if (CContainer::Instance()->GetPreviousValue() == notFeasibleResult) {
      // PrintDebug( DebugData, " valued %e\n", 1, CContainer::Instance()->GetPreviousValue() );
    // } else {
      // PrintDebug( DebugData, " valued %f\n", 1, CContainer::Instance()->GetPreviousValue() );
    // }
  // }

  // If the previous and the current points are the same, return the
  // previous value and do not recompute the cost function.
  if (CContainer::Instance()->IsSamePoint()) {
    //PrintDebug( DebugInfo, " == Same point. Return previous value = " );
    result = CContainer::Instance()->GetPreviousValue();
    //printf( " S %f", result );
  }
  // If is not a feasible point return \infty.
  //  else if (!this->( descentDirection )) {
  else if (!this->IsFeasiblePoint()) {
    //PrintDebug( DebugInfo, " == Point not feasible. Return minus infinite = " );
    result = notFeasibleResult;
    //printf( " X %f", result );
    // printf( " X " );
  }
  // Finally if it is feasible and is a new point, compute the result.
else {
	//PrintDebug( DebugInfo, "\n== Call Evaluate() " );
	++CContainer::Instance()->mNumberOfGetValueCalled;
	result = Evaluate( parameters );
	
	// //PrintDebug( DebugInfo, " = %f\n", result );
	// printf( "\n[%02i-%02i-%02i@%02i:%02i:%02i] %7d: ", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday,
	// current->tm_hour, current->tm_min, current->tm_sec, CContainer::Instance()->mNumberOfGetValueCalled );

#if 0
	printf( "\n[%02i-%02i-%02i@%02i:%02i:%02i] ", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday,
	current->tm_hour, current->tm_min, current->tm_sec );
// #endif

	// retrieve iteration number from file
	std :: ifstream src( ".csp_current_iteration" );
	stringstream inputString;
	src >> inputString.rdbuf();
	src.close();
	int loop_iteration;
	inputString >> loop_iteration;

	char object;  
	switch (mMode) {
	case OnMicrograph: {
			object = 'M';
			break;
		}
	case OnParticle: {
			object = 'P';
			break;
		}
	default: {
			object = 'X';
			break;
		}
	}
// #if 0
	printf( "T%02d\t%c%04d\tLI%04d\tPI%07d\t", CContainer::Instance()->GetTiltSeries()->GetIndex(), object, mCurrentOptimizationIndex, loop_iteration, CContainer::Instance()->mNumberOfGetValueCalled );
	PrintParameters( CContainer::Instance()->GetCurrentPoint() );
	printf( " = %f\n", result );
#endif
}

  CContainer::Instance()->SetPreviousPoint( CContainer::Instance()->GetCurrentPoint() );
  CContainer::Instance()->SetPreviousValue( result );

  // if (result == notFeasibleResult) {
    // PrintDebug( DebugInfo, "=> Returned value: %e\n", 1, result );
  // } else {
    // PrintDebug( DebugInfo, "=> Returned value: %f\n", 1, result );
  // }

  return result;

}

void flrCostFunction::SetMode( LoopModeType mode ) {
  mMode = mode;
  switch (mMode) {
    case OnMicrograph: {
      CContainer::Instance()->mTolerance[0] = params.mToleranceMicrographTiltAngles;
      CContainer::Instance()->mTolerance[1] = params.mToleranceMicrographTiltAxisAngles;
      CContainer::Instance()->mTolerance[2] = 0;
      CContainer::Instance()->mTolerance[3] = 0;
      CContainer::Instance()->mTolerance[4] = 0;
      break;
    }
    case OnParticle: {
      CContainer::Instance()->mTolerance[0] = 0;
      CContainer::Instance()->mTolerance[1] = 0;
      CContainer::Instance()->mTolerance[2] = params.mToleranceParticlesTheta;
      CContainer::Instance()->mTolerance[3] = params.mToleranceParticlesPsi;
      CContainer::Instance()->mTolerance[4] = params.mToleranceParticlesPhi;
      break;
    }
    default: {
      // This mode should never run
      CContainer::Instance()->mTolerance[0] = 0;
      CContainer::Instance()->mTolerance[1] = 0;
      CContainer::Instance()->mTolerance[2] = 0;
      CContainer::Instance()->mTolerance[3] = 0;
      CContainer::Instance()->mTolerance[4] = 0;
      break;
    }
  }
}

LoopModeType flrCostFunction::GetMode() {
  return mMode;
}

bool flrCostFunction::IsFeasiblePoint() const {
  bool condition;
  bool result;
  double step;
  double tolerance;
  double micrographTiltAngleUpperBound;
  double micrographTiltAngleLowerBound;
  int position;

  ParametersType initialPoint = CContainer::Instance()->GetInitialPoint();
  ParametersType currentPoint = CContainer::Instance()->GetCurrentPoint();
  CTiltSeries* pTiltSeries = CContainer::Instance()->GetTiltSeries();

  //// compute particle angle tolerance with respect to no perturbation
  //initialPoint[2] = 0;
  //initialPoint[3] = 0;
  //initialPoint[4] = 0;
  
  // Is the current point in the L1-feasible region?
  // TODO: The tolerance of the particle orientation should be
  // computed in 3D given the three euler angles.
  condition = true;
  for ( unsigned int k = 0; k < GetNumberOfParameters(); k++ ) {
    step = std::fabs( initialPoint[k] - currentPoint[k] );
    tolerance = CContainer::Instance()->mTolerance[k];
    ( step <= tolerance ) ? condition = condition : condition = false;	
    // For the micrograph's tilt angles the tolerance range is defined by
    //and upper and lower values.
	// printf("k = %d, step = %f, tolerance = %f, condition = %d\n", k, step, tolerance, condition );
    if ( ( k == 0 ) && ( tolerance > 0 ) && ( mMode == OnMicrograph ) ) {
      position = pTiltSeries->GetMicrographPositionByIndex( mCurrentOptimizationIndex );
      micrographTiltAngleUpperBound = pTiltSeries->mMicrograph.at( position )->GetTiltAngleUpperBound();
	  if ( pTiltSeries->mMicrograph.at( position )->GetTiltAngle() + currentPoint[k] >= micrographTiltAngleUpperBound) condition = false;

      micrographTiltAngleLowerBound = pTiltSeries->mMicrograph.at( position )->GetTiltAngleLowerBound();
      if ( pTiltSeries->mMicrograph.at( position )->GetTiltAngle() + currentPoint[k] <= micrographTiltAngleLowerBound) condition = false;
	  // printf("POS = %d; TOL = %f: %f, %f; cond = %d\n",position,currentPoint[k],micrographTiltAngleLowerBound,micrographTiltAngleUpperBound,condition);
    }
  }

  result = condition;
  //printf( "\tIsFeasible: %d\n", result );
  return result;
}

/**  */
void flrCostFunction::PrintParameters( const ParametersType & parameters ) const {
  printf( "( %10.6f", parameters[0] );
  for ( unsigned int k = 1; k < GetNumberOfParameters(); k++ ) {
    printf( ", %10.6f", parameters[k] );
  }
  printf( " )" );
}

/**  */
void flrCostFunction::SetCurrentOptimizationIndex( int index ) {
  mCurrentOptimizationIndex = index;
}

/**  */
int flrCostFunction::GetCurrentOptimizationIndex() {
  return mCurrentOptimizationIndex;
}

/**  */
flrCostFunction::MeasureType flrCostFunction::Evaluate( const ParametersType & parameters ) const {

  MeasureType result;
  unsigned int numberOfReadParticleProjections;
  unsigned int numberOfWrittenParticleProjections;
  extern SParameters params;
  char filename[200];
  int evaluationTime;
  int currentTiltSeries;
  char object;
  CTiltSeries* pTiltSeries;

  pTiltSeries = CContainer::Instance()->GetTiltSeries();
  currentTiltSeries = pTiltSeries->GetIndex();
  CContainer::Instance()->SetCurrentOptimizationIndex( mCurrentOptimizationIndex );

  switch (mMode) {
    case OnMicrograph: {
      object = 'M';
      break;
    }
    case OnParticle: {
      object = 'P';
      break;
    }
    default: {
      object = 'X';
      break;
    }
  }

  // Build the name of the .parx file to use with Frealign. Get the number
  // of evaluation of the cost function and concatenate strings.
  evaluationTime = CContainer::Instance()->mNumberOfGetValueCalled;
  sprintf( filename, "/scratch/%s_CSP_T%02d_%c%04d_I%07d.parx", params.mBaseName,
           currentTiltSeries, object, mCurrentOptimizationIndex, evaluationTime );
 
  // 1. For each particle map the point <parameters> belonging to the
  //    CSP manifold to the SP manifold. This will compute the
  //    orientation (3 euler angles) corresponding to THETA, PSI and
  //    PHI in the SP manifold, a triple for each particle.
  //ComputeAnglesInSPManifold( parameters, mMode );

  //CContainer::Instance()->mTiltSeriesPointer->ShowInformation();

  // 2. Run SP. We are using FREALIGN, then, first, we need to write
  // the .par file.

  // 2.1. Write the .par file
  //WriteParFile( parameters, mMode, filename );
  numberOfWrittenParticleProjections = WriteParFile( parameters, mMode, mCurrentOptimizationIndex, filename );

  // 2.2. Run FREALIGN
  RunFrealign( filename );

  // 3. Get the new values of the cross-correlation for each particle
  //    and the new computed center of the particle projection (SHX
  //    and SHY).

  // 3.1. Build the name of the .parx file to read, replacing _1_ for _2_
  //      in the filename
  // strncpy( strstr( filename, "_1_" ), "_2_", 3 );

  // 3.2. Read the .parx file
  numberOfReadParticleProjections = ReadParFile( mMode, mCurrentOptimizationIndex, filename );

  // 3.3. Check the number of written and read particle projections
  if (numberOfWrittenParticleProjections != numberOfReadParticleProjections) {
    printf( "[ERROR] ReadParFile(%s) reads only %d of %d particle projections.\nExiting...\n",
            filename,numberOfReadParticleProjections, numberOfWrittenParticleProjections );
    exit( EXIT_FAILURE );
  }

  // 3.2. Save the updated Position2D in each ParticleProjection

  // 3.3. Compute the cost function, adding the values of the
  //      cross-correlation for the particle of interest.
  result = ComputeCostFunction( mMode, mCurrentOptimizationIndex, this->phaseResidualThreshold, this->minScanOrderToUseForRefinement, this->maxScanOrderToUseForRefinement );
  
  // delete .parx file
  char command[200];
  sprintf( command, "rm -f %s", filename );
  system( command );

  //PrintParameters( CContainer::Instance()->GetCurrentPoint() );
  //printf( " = %f\n", result );

  //printf("result= %f",result);
  //result = 0;

  // 4. Return the cost function value.
  return result;
}

unsigned int flrCostFunction::GetNumberOfParameters( void ) const {
  return mSpaceDimension;
}