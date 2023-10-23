/*=========================================================================

 Name:        flrClasses.hxx
 Category:    Headers files
 Description: This header file contains the definition of the classes
 used in the <framework>
 Language:    C++
 Date:        $Date: 2010-08-30 20:03:51 -0300 (Mon, 30 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 507 $
 Version:     $Id: flrClasses.hxx 507 2010-08-30 23:03:51Z fefo $

 ===========================================================================*/

/** std::vector or std::list?
 * 
 * They store the data in different ways. std::vector stores the data
 * as an array internally, whereas std::list stores it as a doubly
 * linked list.
 *
 * In std::list, you can efficiently insert and remove elements at
 * arbitrary positions, while a vector only provides that for adding
 * or removing elements at the end. Vector provides better locality of
 * data (usually meaning better cache efficiency) and supports direct
 * access to the n-th element for any n between 0 and the size of the
 * vector, while a list must be traversed for that.xx 
 **/

// TODO: Revisit this classes and rewrite them in a correct way, with
// no redundant information and more compact. Also try too use the
// same names to the same kind of members.

#ifndef __flrclasses_hxx
#define __flrclasses_hxx

#include "flrMain.hxx"
#include "flrCostFunctions.hxx"

using namespace std;

/************************************************************************
 * Euler angles class
 ************************************************************************/
class CEulerAngles {

public:
  EulerAngleType mPsi, mTheta, mPhi;

  CEulerAngles();
  CEulerAngles( EulerAngleType psi, EulerAngleType theta, EulerAngleType phi );
  ~CEulerAngles();
  CEulerAngles& operator=( const CEulerAngles& eulerAngle );
  void SetPsi( EulerAngleType psi );
  void SetTheta( EulerAngleType theta );
  void SetPhi( EulerAngleType phi );
  void SetAngles( EulerAngleType psi, EulerAngleType theta, EulerAngleType phi );
  EulerAngleType GetTheta();
  EulerAngleType GetPhi();
  EulerAngleType GetPsi();
  void PrintAngles();

};

/************************************************************************
 * Position in 3D class
 ************************************************************************/
// TODO: templatizar CPosition en N dimensiones
class CPosition3D {

public:
  CoordinateType mCoordX, mCoordY, mCoordZ;
  CPosition3D();
  CPosition3D( CoordinateType x, CoordinateType y, CoordinateType z );
  ~CPosition3D();

};

/************************************************************************
 * Position in 2D class
 ************************************************************************/
class CPosition2D {

public:
  CoordinateType mCoordX, mCoordY;
  CPosition2D();
  CPosition2D( CoordinateType x, CoordinateType y );
  ~CPosition2D();
  void PrintCoordinates();
  void SetPosition( CoordinateType x, CoordinateType y );
  void GetPosition( CoordinateType *x, CoordinateType *y );
};

/************************************************************************
 * Particle projection class definition
 ************************************************************************/
class CParticleProjection {

public:
  char mFilename[120];
  int counter;

  // The following variables are the ones in the frealign .par file,
  // in this order.
  int mIndexInStack; // ind
  CEulerAngles mOrientation; // psi, theta, phi
  CPosition2D mShift; //shx, shy
  CPosition2D mFrameShift;
  int mMagnification; // mag
  int mTiltSeriesIndex; // film
  float mDefocus01; // df1
  float mDefocus02; // df2
  float mAstigmatism; // angast
  float mPpindex; // ppindex
  float mOcc; // occ
  int mLogp; // logp
  float mSigma; // sigma
  float mScore; // score
  float mPreviousScore; // score before csp refinement
  float mChange; // change
  //float mPhaseResidual; // presa
  //float mDeltaPhaseResidual; // dpres
  int mParticleIndex; //ptlind
  EulerAngleType mMicrographTiltAngle; // tiltan
  EulerAngleType mMicrographTiltAxisAngle; // tiltaixsan;
  CPosition2D mMicrographTiltAxisShift; // tilt shift
  float mDose; // dosexx
  int mMicrographsIndex; // scanor
  int mFrameIndex; // cnfdnc
  float mCCX; // ptlccx
  int mIndexInLocalStack; 

  CParticleProjection();
  CParticleProjection( int indexInStack, // ind
                       CEulerAngles orientation, // psi, theta, phi
                       CPosition2D shift, //shx, shy
                       CPosition2D frameShift,
                       int magnification, // mag
                       int tiltSeriesIndex, // film
                       float defocus01, // df1
                       float defocus02, // df2
                       float astigmatism, // angast
                       float ppindex, // ppindex
                       float occ, // occ
                       int logp, // logp
                       float sigma, // sigma
                       float score, // score
                       float change, // change
                       //float phaseResidual, // presa
                       //float deltaPhaseResidual, // dpres
                       int particleIndex, //ptlind
                       EulerAngleType micrographTiltAngle, // tiltan
                       EulerAngleType micrographTiltAxisAngle, // tiltaxisan;
                       float dose, // dosexx
                       int micrographsIndex, // scanor
                       int frameIndex, // cnfdnc
                       float particleCCX, // ptlccx
                       CPosition2D micrographTiltAxisShift,
                       int ind_local_stack
      );
  CParticleProjection( CEulerAngles orientation, CPosition2D position, int particleIndex, int micrographsIndex );
  ~CParticleProjection();

  void ShowInformation();
  int GetParticleIndex();
  int GetMicrographIndex();
  int GetFrameIndex();
  void SetIndex( int index );
  void SetOrientation( CEulerAngles orientation );
  void SetShift( CPosition2D shift );
  void SetShift( CoordinateType x, CoordinateType y );
  CEulerAngles GetOrientation();
  CPosition2D GetShift();

};

/************************************************************************
 * Particle class definition
 ************************************************************************/
class CParticle {

public:
  int mIndex, mFrameIndex;
  int mTiltSeriesIndex;
  float mOcc;
  CEulerAngles mOrientation;
  CPosition3D mPosition3D;
  vector< pair<int,int> > mListOfMicrographs;
  vector<CParticleProjection*> mParticleProjectionSet;
  EulerAngleType mNormal[3];
  EulerAngleType mMatrix[16];

  CParticle();
  CParticle( CEulerAngles eulerAngles,
             CPosition3D position3D,
             int mNumber,
             int frameIndex, 
             int mTiltSeriesIndex,
             EulerAngleType normal[3],
             EulerAngleType matrix[16],
             float occ);
  ~CParticle();

  int GetIndex();
  int GetFrameIndex();

  void SetNumber( int number );
  int GetTiltSeriesIndex();
  CEulerAngles GetOrientation();
  void GetNormal( EulerAngleType normal[3] );
  void GetMatrix( EulerAngleType matrix[16] );
  void SetMatrixTranslation( EulerAngleType x, EulerAngleType y, EulerAngleType z );
  CPosition3D GetPosition();

  float GetOcc();

  void SetOrientation( CEulerAngles* eulerAngles );
  void SetOrientation( EulerAngleType theta, EulerAngleType phi, EulerAngleType psi );
  void SetPosition( CPosition3D* position3D );
  void SetPosition( CoordinateType x, CoordinateType y, CoordinateType z );
  void AddProjection( CParticleProjection* p_Projection );
  bool HasProjectionInMicrograph( int micrographIndex, int frameIndex );
  void ShowInformation();
};

/************************************************************************
 * Micrograph class definition
 ************************************************************************/
class CMicrograph {

private:
  EulerAngleType mTiltAngle;
  EulerAngleType mTiltAxisAngle;
  CPosition2D mTiltAxisShift;
  CPosition2D mFrameShift;
  double mTiltSeriesDefocus1Offset;
  double mTiltSeriesDefocus2Offset;
  double mTiltSeriesAstigmatism;

public:

  // Index (number or Id) readed from the .par file
  int mMicrographIndex, mFrameIndex;

  // Index of the tilt series where the micrograph belongs.
  int mTiltSeriesIndex;

  // Tolerance in the tilt angle
  EulerAngleType mTiltAngleTolerance;

  // Maximum and minimum values that can be assigned to the tilt angle
  // given the separation of the micrograph in the tilt series.
  EulerAngleType mTiltAngleUpperBound;
  EulerAngleType mTiltAngleLowerBound;

  // Particles with projections into this micrograph, a vector of
  // indexes (ordered as they are readed) and the vector of pointers
  // to the CParticleProjection
  vector< pair<int,int> > mListOfParticleProjections;
  vector<CParticleProjection*> mParticleProjectionSet;

  CMicrograph();
  CMicrograph( EulerAngleType tiltAngle,
               EulerAngleType tiltAxisAngle,
               CPosition2D tiltAxisShift,
               CPosition2D frameShift,
               int numberID,
               int frameIndex,
               int tiltSeriesIndex,
               double tiltDefocus1Offset,
               double tiltDefocus2Offset,
               double tiltAstigmatism );
  ~CMicrograph();

  CMicrograph& operator=( const CMicrograph& cMicrograph );

  void PrintListOfParticleProjections();
  int GetNumberOfParticleProjections();
  int GetIndex();
  int GetFrameIndex();

  int AddParticleProjection( CParticleProjection* particleProjection );
  bool HasParticleProjection( int particleProjectionIndex, int frameIndex );

  EulerAngleType GetTiltAngle();
  void SetTiltAngle( EulerAngleType angle );

  EulerAngleType GetTiltAngleUpperBound();
  EulerAngleType GetTiltAngleLowerBound();

  EulerAngleType GetTiltAxisAngle();
  void SetTiltAxisAngle( EulerAngleType angle );

  CPosition2D GetTiltAxisShift();
  void SetTiltAxisShift( CPosition2D position );

  CPosition2D GetFrameShift();
  void SetFrameShift( CPosition2D position );

  EulerAngleType GetTiltAngleTolerance();
  void SetTiltAngleTolerance( EulerAngleType tolerance );

  void SetTiltDefocusOffset( double tiltDefocus1Offset, double tiltDefocus2Offset );
  double GetTiltDefocus1Offset();
  double GetTiltDefocus2Offset();

  void SetTiltAstigmatism( double tiltAstigmatism );
  double GetTiltAstigmatism();
};

/************************************************************************
 * TiltSeries class definition
 ************************************************************************/
class CTiltSeries {

public:
  // Number or ID of the tilt series in the data set.
  int mTiltSeriesIndex;

  // Max index of frames in the data set.
  int mMaxFrameIndex = 0;


  // List of particles in the 3D volume scanned in this tilt series.
  set<int> mListOfParticles;
  //vector<int> mListOfParticlesFrames;

  // List of micrographs in this tilt series.
  set<int> mListOfMicrographs;
  //set<int> mListOfMicrographsFrames;
  // List of particle projections in the tilt series.
  vector<int> mListOfParticlesProjections;

  // The following are containers for micrographs, particles and
  // particles projections. These are vectors load with push-back,
  // then the order of the objects in each one is as they come from
  // the parx file. In order to recover the N-th object, we need to
  // access to the place where it was pushed; this number is given by
  // the function GetMicrographPositionByIndex() and
  // GetParticlePositionByIndex() where the argument is the number
  // (index) of the object (micrograph or particle) in the par file
  // and returns the position in the lists or sets where to find the
  // object.
  vector<CMicrograph*> mMicrograph;
  vector<CParticle*> mParticle;
  vector<CParticleProjection*> mParticleProjectionSet;

  // Index of particles. These maps are used by the functions
  // GetMicrographPositionByIndex() and
  // GetParticlePositionByIndex(). Each one stores for the index-th
  // object the position in the lists or sets where to find the
  // object. For example, if the micrograph index=23 is the third to
  // appear in the parx file mMicrographByIndex[23]=2=3-1 Then,
  // GetMicrographPositionByIndex(23) = 2, mMicrograph.at(2) is the
  // micrograph number 23, and mListOfMicrographs.at(2) = 23.
  IntMapType mParticleByIndex;
  IntMapType mMicrographByIndex;

public:
  CTiltSeries();
  CTiltSeries( int mTiltSeriesIndex );
  ~CTiltSeries();

  int GetIndex();
  void SetIndex( int index );
  void SetMaxFrameIndex( int frameIndex );

  int AddParticleProjection( CParticleProjection* particleProjection );

  int GetNumberOfMicrographs();
  void PrintListOfMicrographs();
  int AddMicrograph( CMicrograph* p_Micrograph );
  bool HasMicrograph( int micrographIndex, int frameIndex );
  void ComputeMicrographTiltAngleBounds(bool=false);

  int GetNumberOfParticles();
  void PrintListOfParticles();
  int AddParticle( CParticle* p_Particle );
  bool HasParticle( int particleIndex, int frameIndex );

  void ShowInformation();

  // Returns the position in the vector mMicrograph of the micrograph
  // given by a particular mMicrographIndex. The same for the
  // Particles.
  int GetMicrographPositionByIndex( int micrographIndex, int frameIndex );
  int GetParticlePositionByIndex( int particleIndex, int frameIndex );

  int GetAllPositionsByMicrograph( int mInd, vector<int> & positions );
  int GetAllPositionsByMicrograph( int mInd, int frameIndex, vector<int> & positions );
  int GetAllPositionsByParticle( int pInd, vector<int> & positions );
  int GetAllPositionsByParticle( int pInd, int frameIndex, vector<int> & positions );
};

/************************************************************************
 * Container class definition
 *
 * This class is a singleton and is used as a container where the
 * information is saved and can be accessed everywhere. Is generally
 * used to recover the tilt series pointer.
 *
 ************************************************************************/
class CContainer {
public:

  int mNumberOfGetValueCalled;
  int mCurrentOptimizationIndex;

  // Set and get the current index of the object being optimized.
  void SetCurrentOptimizationIndex( int );
  int GetCurrentOptimizationIndex();

  // The mSelectionMask selects which of parameters will be refined.
  // The order of the parameters in the mask is the same as the
  // parameters of the optimizer (see the definition of the
  // flrCostFunction.
  //
  //  IS NOT BEEN USED NOW, THE mMode CHANGES THE VALUES OF mTolerance
  //  SETTING IN ZERO THE UNDESIRED DIRECTIONS
  //
  flrCostFunction::ParametersType mSelectionMask;

  // Tolerance (or feasible region) for each dimension, same order as
  // the selection mask.
  double mTolerance[6];

  // Previous and current point in the optimization
  flrCostFunction::ParametersType mPreviousPoint;
  flrCostFunction::ParametersType mCurrentPoint;
  flrCostFunction::ParametersType mInitialPoint;
  flrCostFunction::ParametersType mBestPoint;

  // Previous and current vale in the previous and current point in the optimization
  flrCostFunction::MeasureType mPreviousValue;
  flrCostFunction::MeasureType mCurrentValue;
  flrCostFunction::MeasureType mInitialValue;
  flrCostFunction::MeasureType mBestValue;

  CTiltSeries* mTiltSeriesPointer;
  //  CDensityMap* mDensityMapPointer;

  static CContainer* Instance();

  void SetTiltSeries( CTiltSeries* tiltSeries );
  CTiltSeries* GetTiltSeries();

  //  void SetDensityMap( CDensityMap* densityMap );
  //  CDensityMap* GetDensityMap();

  void SetPreviousPoint( flrCostFunction::ParametersType point );
  flrCostFunction::ParametersType GetPreviousPoint();

  void SetCurrentPoint( flrCostFunction::ParametersType point );
  flrCostFunction::ParametersType GetCurrentPoint();

  void SetInitialPoint( flrCostFunction::ParametersType point );
  flrCostFunction::ParametersType GetInitialPoint();

  void SetBestPoint( flrCostFunction::ParametersType point );
  flrCostFunction::ParametersType GetBestPoint();

  void SetPreviousValue( flrCostFunction::MeasureType value );
  flrCostFunction::MeasureType GetPreviousValue();

  void SetInitialValue( flrCostFunction::MeasureType value );
  flrCostFunction::MeasureType GetInitialValue();

  void SetCurrentValue( flrCostFunction::MeasureType value );
  flrCostFunction::MeasureType GetCurrentValue();

  int GetDescentDirection();

  bool IsSamePoint();

private:

  // Constructor is private so that it can not be called.
  CContainer();

  // Copy constructor is private too.
  CContainer( CContainer const& ) {
  }

  // Assignment operator is private too.
  // CContainer& operator=( CContainer const& ){}

  // The only one instance of the singleton, private too.
  static CContainer* mpInstance;

};


class CImage {

public:
    Image ** imageStack;
    int numImages = 0;
    bool isRaw; 

    CImage();
    CImage(Image** image, int number_of_images, bool is_raw);
    ~CImage();

    void SetImage(Image** image);
    bool CImage::IsRaw();
};


#endif
