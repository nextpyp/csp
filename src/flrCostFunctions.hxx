/*=========================================================================

 Name:        flrCostFunctions.hxx
 Category:    Headers files
 Description: This is the header file with the definitions for the cost
 function.
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrCostFunctions.hxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/
#ifndef __flrcostfunctions_hxx
#define __flrcostfunctions_hxx

#include "flrClasses.hxx"
#include "flrMain.hxx"
#include "itkCostFunction.h"
#include <refine3d_cspt.h>
/**
 * Cost function computation. This is the cost function associated with
 * the optimizer. The parameters' (or variables) vector to minimize has
 * the following fields (in this order):
 *    0: Micrograph tilt angle
 *    1: Micrograph tilt axis angle
 *    2: Particle theta angle
 *    3: Particle psi angle
 *    4: Particle phi angle
 *
 * The optimization is performed on a micrograph (changing the first two
 * parameters) or on a particle (changing the last three parameters).
 * This is selected by the mMode member. However the optimization is
 * performed over a manifold with these five dimension, defining a
 * feasible region depending the mode.
 *
 */
class flrCostFunction : public itk::SingleValuedCostFunction {

private:
  // Dimension of the space of variables (points, parameters)
  static const unsigned int mSpaceDimension = VARIABLES_SPACE_DIMENSION;

  // Mode on which the optimizer is running. Mode has two possible
  // values given by the LoopModeType: OnMicrograph or OnParticle
  LoopModeType mMode;
  // Index (number inthe parx file) of the current object being
  // optimized. Should be a micrograph or a particle, that is
  // given by mMode.
  int mCurrentOptimizationIndex;
  int mCurrentOptimizationFrameIndex;

  double phaseResidualThreshold;
  double minScanOrderToUseForRefinement;
  double maxScanOrderToUseForRefinement;
  int refineProjectionCutoff; 
  vector<double> weightList;

  ReconstructedVolume* ref_3d; 
  MRCFile* stack;
  class CImage* image_set;
  ImageProjectionComparison** comparison_objects;

  map<int, float*> priors_average; 
  map<int, float*> priors_variance;

  bool mFrameRefine = false;
  bool mCistemOptimizor = false;
  bool mEvaluateAllProjections = false;
  //Particle** particle_set;

public:
	
  typedef flrCostFunction Self;
  typedef itk::SingleValuedCostFunction Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  itkNewMacro( Self )
  itkTypeMacro( flrCostFunction, SingleValuedCostFunction )
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::DerivativeType DerivativeType;
  typedef Superclass::MeasureType MeasureType;
  typedef vector<double> MyParametersType;

  // .............
  double mRBFactor;
  vector<float> mRingsRadius;
  WeightArrayType mRingsWeigths;
  // Tilt series with all the data
  class CTiltSeries* mTiltSeries;

  // Constructor
  flrCostFunction();
  flrCostFunction( WeightArrayType ringsWeigths );

  // Destructor
  ~flrCostFunction();

  // Return the dimension of the space of variables
  unsigned int GetNumberOfParameters( void ) const;

  // Get the derivatives: do nothing
  void GetDerivative( const ParametersType &, DerivativeType & ) const;

  // Get the value of the cost function at the point <parameters>
  MeasureType GetValue( const ParametersType & parameters ) const;

  // Set and get the mode of optimization.
  void SetMode( LoopModeType mode );
  LoopModeType GetMode();
  
  // set if refining on a frame basis
  void SetFrameRefine( bool );
    
  void UseCistemOptimizor( bool );
  void EvaluateAllProjections( bool );

  // Set and get the current index of the object being optimized.
  void SetCurrentOptimizationIndex( int );
  void SetCurrentOptimizationFrameIndex( int );
  int GetCurrentOptimizationIndex();
  int GetCurrentOptimizationFrameIndex();

  MeasureType Evaluate( const ParametersType & parameters ) const;

  // Check if the current point (CContainer::GetCurrentPoint()) is a
  // feasible point, checking the descent direction and the tolerance
  // for this paramenter
  bool IsFeasiblePoint() const;

  // Print the point <parameters>
  void PrintParameters( const ParametersType & parameters ) const;

  // Set Phase Residual Threshold
  void SetPhaseResidualThreshold( double d ){ this->phaseResidualThreshold = d; }

  // Set Scan Order range to use for refinement
  void SetScanOrderRangeToUseForRefinement( double min, double max ){ this->minScanOrderToUseForRefinement = min; this->maxScanOrderToUseForRefinement = max; }
  void SetRefineProjectionCutoff(int cutoff) {this->refineProjectionCutoff = cutoff;}
  void SetWeightList( vector<double> & weights ) { this->weightList = weights; }
  void SetReference( ReconstructedVolume* vol ) { this->ref_3d = vol; }
  void SetStack( MRCFile* particle_stack ) { this->stack = particle_stack; }
  void SetImageSet( CImage* images ) { this->image_set = images; }
  void SetPriors( map<int, float*> & priors_average, map<int, float*> & priors_variance ) {
      this->priors_average = priors_average;
      this->priors_variance = priors_variance;
  }
  void SetComparisonObjects(ImageProjectionComparison** comparison_objects) {
      this->comparison_objects = comparison_objects;
  }
  //void SetParticleSet( Particle** particles ) { this->particle_set = particles; }
};
#endif /* FLRCOSTFUNCTIONS_HXX_ */
