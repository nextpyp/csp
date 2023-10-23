/*=========================================================================

 Name:        flrClasses.cxx

 Category:    Implementation files
 Description: Implementation of the functions defined in flrClasses.cxx
 Language:    C++
 Date:        $Date: 2010-08-30 23:59:03 -0300 (Mon, 30 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 508 $
 Version:     $Id: flrClasses.cxx 508 2010-08-31 02:59:03Z fefo $

 ===========================================================================*/
#include "flrClasses.hxx"
#include "flrReadParametersFile.hxx"
#include "vtkMath.h"

/************************************************************************
 * Euler angles class
 ************************************************************************/
/** Default constructor */
CEulerAngles::CEulerAngles() {
    mPsi = 0;
    mTheta = 0;
    mPhi = 0;
}

/** Alternate constructor  */
CEulerAngles::CEulerAngles( EulerAngleType psi, EulerAngleType theta, EulerAngleType phi ) {
    mPsi = psi;
    mTheta = theta;
    mPhi = phi;
}

/** Destructor */
CEulerAngles::~CEulerAngles() {
}

/** Operator copy */
CEulerAngles& CEulerAngles::operator=( const CEulerAngles& in ) {
    mPsi = in.mPsi;
    mTheta = in.mTheta;
    mPhi = in.mPhi;
    return *this;
}

/**  */
void CEulerAngles::SetAngles( EulerAngleType psi, EulerAngleType theta, EulerAngleType phi ) {
    this->SetPsi( psi );
    this->SetTheta( theta );
    this->SetPhi( phi );
}

/**  */
void CEulerAngles::SetTheta( EulerAngleType theta ) {
    mTheta = theta;
}

/**  */
void CEulerAngles::SetPhi( EulerAngleType phi ) {
    mPhi = phi;
}

/**  */
void CEulerAngles::SetPsi( EulerAngleType psi ) {
    mPsi = psi;
}

/**  */
EulerAngleType CEulerAngles::GetTheta() {
    return mTheta;
}

/**  */
EulerAngleType CEulerAngles::GetPhi() {
    return mPhi;
}

/**  */
EulerAngleType CEulerAngles::GetPsi() {
    return mPsi;
}

/**  */
void CEulerAngles::PrintAngles() {
    printf( "(%f,%f,%f)", mTheta, mPhi, mPsi );
}

/************************************************************************
 * Position in 3D class
 ************************************************************************/

/** Default constructor */
CPosition3D::CPosition3D() {
    mCoordX = 0;
    mCoordY = 0;
    mCoordZ = 0;
}

/** Alternate constructor  */
CPosition3D::CPosition3D( CoordinateType x, CoordinateType y, CoordinateType z ) {
    mCoordX = x;
    mCoordY = y;
    mCoordZ = z;
}

/** Destructor */
CPosition3D::~CPosition3D() {
}

/************************************************************************
 * Position in 2D class
 ************************************************************************/

/** Default constructor */
CPosition2D::CPosition2D() {
    mCoordX = 0;
    mCoordY = 0;
}

/** Alternate constructor  */
CPosition2D::CPosition2D( CoordinateType x, CoordinateType y ) {
    mCoordX = x;
    mCoordY = y;
}

/** Destructor */
CPosition2D::~CPosition2D() {
}

void CPosition2D::PrintCoordinates() {
    printf( "(%f,%f)", mCoordX, mCoordY );
}

void CPosition2D::SetPosition( CoordinateType x, CoordinateType y ) {
    mCoordX = x;
    mCoordY = y;
}

void CPosition2D::GetPosition( CoordinateType *x, CoordinateType *y ) {
    *x = mCoordX;
    *y = mCoordY;
}

/************************************************************************
 * Particle class definition
 ************************************************************************/

/** Default constructor */
CParticle::CParticle() {
    CEulerAngles mOrientation();
    CPosition3D mPosition3D();
    mIndex = 0;
    mTiltSeriesIndex = 0;
}

/** Alternate constructor  */
CParticle::CParticle( CEulerAngles eulerAngles,
                      CPosition3D position3D,
                      int number,
                      int frameIndex, 
                      int tiltSeriesIndex,
                      EulerAngleType normal[3],
                      EulerAngleType matrix[16],
                      float occ ) {
    this->mOrientation = eulerAngles;
    this->mPosition3D = position3D;
    this->mIndex = number;
    this->mTiltSeriesIndex = tiltSeriesIndex;
    this->mNormal[0] = normal[0];
    this->mNormal[1] = normal[1];
    this->mNormal[2] = normal[2];
    this->mMatrix[0] = matrix[0];
    this->mMatrix[1] = matrix[1];
    this->mMatrix[2] = matrix[2];
    this->mMatrix[3] = matrix[3];
    this->mMatrix[4] = matrix[4];
    this->mMatrix[5] = matrix[5];
    this->mMatrix[6] = matrix[6];
    this->mMatrix[7] = matrix[7];
    this->mMatrix[8] = matrix[8];
    this->mMatrix[9] = matrix[9];
    this->mMatrix[10] = matrix[10];
    this->mMatrix[11] = matrix[11];
    this->mMatrix[12] = matrix[12];
    this->mMatrix[13] = matrix[13];
    this->mMatrix[14] = matrix[14];
    this->mMatrix[15] = matrix[15];
    this->mOcc = occ;
    this->mFrameIndex = frameIndex; 

    // The transformation matrix may not be orthonormal because of the truncations 
    // done when transfering the 3DAVG matrices to the CSP framework.
    // We need to make sure transformation matrix is norm=1
    EulerAngleType A[3][3], U[3][3], w[3], VT[3][3];
    A[0][0] = matrix[0]; A[0][1] = matrix[1]; A[0][2] = matrix[2];
    A[1][0] = matrix[4]; A[1][1] = matrix[5]; A[1][2] = matrix[6];
    A[2][0] = matrix[8]; A[2][1] = matrix[9]; A[2][2] = matrix[10];
    vtkMath::SingularValueDecomposition3x3(A,U,w,VT);
    EulerAngleType wmax = max( w[0], max( w[1], w[2] ) );

    this->mMatrix[0] /= wmax;
    this->mMatrix[1] /= wmax;
    this->mMatrix[2] /= wmax;
    this->mMatrix[4] /= wmax;
    this->mMatrix[5] /= wmax;
    this->mMatrix[6] /= wmax;
    this->mMatrix[8] /= wmax;
    this->mMatrix[9] /= wmax;
    this->mMatrix[10] /= wmax;
}

/** Destructor */
CParticle::~CParticle() {
}

/**  */
int CParticle::GetIndex() {
    return mIndex;
}

/** */
int CParticle::GetFrameIndex() {
    return mFrameIndex;
}

/**  */
void CParticle::SetNumber( int number ) {
    mIndex = number;
}

/**  */
int CParticle::GetTiltSeriesIndex() {
    return mTiltSeriesIndex;
}

/** Get the euler angles  */
CEulerAngles CParticle::GetOrientation() {
    return mOrientation;
}

/** Set the euler angles */
void CParticle::SetOrientation( CEulerAngles* eulerAngles ) {
    mOrientation.mTheta = eulerAngles->mTheta;
    mOrientation.mPhi = eulerAngles->mPhi;
    mOrientation.mPsi = eulerAngles->mPsi;
}

/** Set the euler angles, alternative */
void CParticle::SetOrientation( EulerAngleType theta, EulerAngleType phi, EulerAngleType psi ) {
    mOrientation.mTheta = theta;
    mOrientation.mPhi = phi;
    mOrientation.mPsi = psi;
}

void CParticle::GetNormal( EulerAngleType normal[3] ) {
    normal[0] = mNormal[0];
    normal[1] = mNormal[1];
    normal[2] = mNormal[2];
}

void CParticle::GetMatrix( EulerAngleType matrix[16] ) {
    matrix[0] = mMatrix[0];
    matrix[1] = mMatrix[1];
    matrix[2] = mMatrix[2];
    matrix[3] = mMatrix[3];
    matrix[4] = mMatrix[4];
    matrix[5] = mMatrix[5];
    matrix[6] = mMatrix[6];
    matrix[7] = mMatrix[7];
    matrix[8] = mMatrix[8];
    matrix[9] = mMatrix[9];
    matrix[10] = mMatrix[10];
    matrix[11] = mMatrix[11];
    matrix[12] = mMatrix[12];
    matrix[13] = mMatrix[13];
    matrix[14] = mMatrix[14];
    matrix[15] = mMatrix[15];
}
void CParticle::SetMatrixTranslation( EulerAngleType x, EulerAngleType y, EulerAngleType z ){
    mMatrix[3] = x;
    mMatrix[7] = y;
    mMatrix[11] = z;
}

float CParticle::GetOcc() {
    return this->mOcc;
}

/** Get the 3D position */
CPosition3D CParticle::GetPosition() {
    return mPosition3D;
}

/** Set the 3D position */
void CParticle::SetPosition( CPosition3D* position3D ) {
    mPosition3D.mCoordX = position3D->mCoordX;
    mPosition3D.mCoordY = position3D->mCoordY;
    mPosition3D.mCoordZ = position3D->mCoordZ;
}

void CParticle::SetPosition( CoordinateType x, CoordinateType y, CoordinateType z ) {
    mPosition3D.mCoordX = x;
    mPosition3D.mCoordY = y;
    mPosition3D.mCoordZ = z;
}

void CParticle::AddProjection( CParticleProjection* p_Projection ) {
    mListOfMicrographs.push_back( {p_Projection->GetMicrographIndex(), p_Projection->GetFrameIndex() } );
    mParticleProjectionSet.push_back( p_Projection );
}

bool CParticle::HasProjectionInMicrograph( int index, int frameIndex ) {
    bool isPresent = false;
    vector< pair<int,int> >::iterator it;
    for ( it = mListOfMicrographs.begin(); it != mListOfMicrographs.end(); it++ ) {
        if ( (it->first == index) && (it->second == frameIndex) ) {
            isPresent = true;
            break;
        }
    }
    return isPresent;
}

void CParticle::ShowInformation() {
    printf( "[Particle %d]\n", mIndex );
    printf( "Belongs to tilt series: %d. ", mTiltSeriesIndex );
    printf( "Orientation: " );
    mOrientation.PrintAngles();
    printf( ". " );
    printf( "Has projections into micrographs: " );
    for ( unsigned int ind = 0; ind < mListOfMicrographs.size(); ind++ ) {
        printf( "%d ", mListOfMicrographs.at( ind ) );
    }
    printf( "\n" );
    //  printf( "Normal: (%lf, %lf, %lf)\n", mNormal[0], mNormal[1], mNormal[2] );
    //  printf( "Refinement: \t%lf\t%lf\t%lf\t%lf\n", mMatrix[0], mMatrix[1], mMatrix[2], mMatrix[3] );
    //  printf( "\t\t%lf\t%lf\t%lf\t%lf\n", mMatrix[4], mMatrix[5], mMatrix[6], mMatrix[7] );
    //  printf( "\t\t%lf\t%lf\t%lf\t%lf\n", mMatrix[8], mMatrix[9], mMatrix[10], mMatrix[11] );
    //  printf( "\t\t%lf\t%lf\t%lf\t%lf\n", mMatrix[12], mMatrix[13], mMatrix[14], mMatrix[15] );
    printf( "Normal: (%lf, %lf, %lf)\n", mNormal[0], mNormal[1], mNormal[2] );
    printf( "Refinement matrix: [ %.3lf %.3lf %.3lf %.3lf;", mMatrix[0], mMatrix[1], mMatrix[2], mMatrix[3] );
    printf( " %.3lf %.3lf %.3lf %.3lf;", mMatrix[4], mMatrix[5], mMatrix[6], mMatrix[7] );
    printf( " %.3lf %.3lf %.3lf %.3lf;", mMatrix[8], mMatrix[9], mMatrix[10], mMatrix[11] );
    printf( " %.3lf %.3lf %.3lf %.3lf ]\n", mMatrix[12], mMatrix[13], mMatrix[14], mMatrix[15] );
}

/************************************************************************
 * Particle projection class definition
 ************************************************************************/

/** Default constructor */
CParticleProjection::CParticleProjection() {
    CEulerAngles orientation( 0, 0, 0 );
    CPosition2D position( 0, 0 );
    mOrientation = orientation;
    mShift = position;
    mParticleIndex = 0;
    mMicrographsIndex = 0;
    //  mData = ImageType::New();
}

/** Destructor */
CParticleProjection::~CParticleProjection() {
}

CParticleProjection::CParticleProjection( int indexInStack, // ind
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
                                          EulerAngleType micrographTiltAngle, // tiltan;
                                          EulerAngleType micrographTiltAxisAngle,
                                          float dose, // dosexx
                                          int micrographsIndex, // scanor
                                          int frameIndex, // cnfdnc
                                          float CCX, // ptlccx
                                          CPosition2D micrographTiltAxisShift, 
                                          int ind_local_stack
) {
    mIndexInStack = indexInStack; // ind
    mOrientation = orientation; // psi, theta, phi
    mShift = shift; //shx, shy
    mFrameShift = frameShift;
    mMagnification = magnification; // mag
    mTiltSeriesIndex = tiltSeriesIndex; // film
    mDefocus01 = defocus01; // df1
    mDefocus02 = defocus02; // df2
    mAstigmatism = astigmatism; // angast
    mPpindex = ppindex; // occ
    mOcc = occ; // occ
    mLogp = logp; // logp
    mSigma = sigma; // sigma
    mScore = score; // score
    mPreviousScore = score;
    mChange = change; // change
    //mPhaseResidual = phaseResidual; // presa
    //mDeltaPhaseResidual = deltaPhaseResidual; // dpres
    mParticleIndex = particleIndex; //ptlind
    mMicrographTiltAngle = micrographTiltAngle;
    mMicrographTiltAxisAngle = micrographTiltAxisAngle; // tiltaxisan;
    mDose = dose; // dosexx
    mMicrographsIndex = micrographsIndex; // scanor
    mFrameIndex = frameIndex; // cnfdnc
    mCCX = CCX; // ptlccx
    mMicrographTiltAxisShift = micrographTiltAxisShift,
    mIndexInLocalStack = ind_local_stack;
}

/** Alternate constructor */
CParticleProjection::CParticleProjection( CEulerAngles orientation,
                                          CPosition2D position,
                                          int particleIndex,
                                          int micrographsIndex ) {

    mOrientation = orientation;
    mShift = position;
    mParticleIndex = particleIndex;
    mMicrographsIndex = micrographsIndex;

}

void CParticleProjection::ShowInformation() {
    printf( "[Particle Projection %d] ", mParticleIndex );
    printf( "Belongs to micrograph: %d. ", mMicrographsIndex );
    printf( "Orientation: " );
    mOrientation.PrintAngles();
    printf( ". " );
    printf( "2D position: " );
    mShift.PrintCoordinates();
    printf( ". " );
    //  ImageType::SizeType size = mData->GetLargestPossibleRegion().GetSize();
    //  printf( "Image size: (%lu,%lu). ", size[0], size[1] );
    printf( "\n" );
}

int CParticleProjection::GetParticleIndex() {
    return mParticleIndex;
}

int CParticleProjection::GetMicrographIndex() {
    return mMicrographsIndex;
}

int CParticleProjection::GetFrameIndex() {
    return mFrameIndex;
}

void CParticleProjection::SetIndex( int index ) {
    mParticleIndex = index;
}

/**  */
void CParticleProjection::SetOrientation( CEulerAngles orientation ) {
    mOrientation = orientation;
}

/**  */
void CParticleProjection::SetShift( CPosition2D shift ) {
    mShift = shift;
}

/**  */
void CParticleProjection::SetShift( CoordinateType x, CoordinateType y ) {
    mShift.mCoordX = x;
    mShift.mCoordY = y;
}

/**  */
CEulerAngles CParticleProjection::GetOrientation() {
    return mOrientation;
}

/**  */
CPosition2D CParticleProjection::GetShift() {
    return mShift;
}

///**  */
//ImageType::SizeType CParticleProjection::GetImageSize() {
//  ImageType::SizeType imageSize;
//  imageSize = this->mData->GetLargestPossibleRegion().GetSize();
//  return imageSize;
//}

/************************************************************************
 * Micrograph class definition
 ************************************************************************/

/** Default constructor */
CMicrograph::CMicrograph() {
}

/** Destructor */
CMicrograph::~CMicrograph() {
}

/** Alternate constructor */
CMicrograph::CMicrograph( EulerAngleType tiltAngle,
                          EulerAngleType tiltAxisAngle,
                          CPosition2D tiltAxisShift,
                          CPosition2D frameShift,
                          int numberID,
                          int frameIndex,
                          int tiltSeriesIndex,
                          double tiltDefocus1Offset,
                          double tiltDefocus2Offset,
                          double tiltAstigmatism ) {

    mTiltAngle = tiltAngle;
    mTiltAxisAngle = tiltAxisAngle;
    mTiltAxisShift = tiltAxisShift;
    mFrameShift = frameShift;
    mMicrographIndex = numberID;
    mFrameIndex = frameIndex;
    mTiltSeriesIndex = tiltSeriesIndex;
    mTiltSeriesDefocus1Offset = tiltDefocus1Offset;
    mTiltSeriesDefocus2Offset = tiltDefocus2Offset;
    mTiltSeriesAstigmatism = tiltAstigmatism;
    mTiltAngleTolerance = itk::NumericTraits<double>::infinity();
}

/**  */
int CMicrograph::GetNumberOfParticleProjections() {
    return mListOfParticleProjections.size();
}

/**  */
int CMicrograph::GetIndex() {
    return mMicrographIndex;
}

/** */
int CMicrograph::GetFrameIndex() {
    return mFrameIndex;
}

/**  */
void CMicrograph::PrintListOfParticleProjections() {
    vector< pair<int,int> >::iterator it;
    for ( it = mListOfParticleProjections.begin(); it != mListOfParticleProjections.end(); it++ ) {
        printf( "%d ", *it );
    }
}

/**  */
int CMicrograph::AddParticleProjection( CParticleProjection* particleProjection ) {
    mListOfParticleProjections.push_back( {particleProjection->GetParticleIndex(), particleProjection->GetFrameIndex() } );
    mParticleProjectionSet.push_back( particleProjection );

    return mListOfParticleProjections.size();
}

/**  */
CMicrograph& CMicrograph::operator=( const CMicrograph& in ) {
    this->mTiltAxisAngle = in.mTiltAxisAngle;
    this->mTiltAngle = in.mTiltAngle;
    this->mTiltAxisAngle = in.mTiltAxisAngle;
    this->mTiltAxisShift = in.mTiltAxisShift;
    this->mMicrographIndex = in.mMicrographIndex;
    this->mTiltSeriesIndex = in.mTiltSeriesIndex;
    this->mTiltSeriesDefocus1Offset = in.mTiltSeriesDefocus1Offset;
    this->mTiltSeriesDefocus2Offset = in.mTiltSeriesDefocus2Offset;
    this->mTiltSeriesAstigmatism = in.mTiltSeriesAstigmatism;
    this->mListOfParticleProjections = in.mListOfParticleProjections;
    this->mParticleProjectionSet = in.mParticleProjectionSet;
    return *this;
}

/**  */
bool CMicrograph::HasParticleProjection( int index, int frameIndex ) {
    bool isPresent = false;
    vector< pair<int,int> >::iterator it;
    for ( it = mListOfParticleProjections.begin(); it != mListOfParticleProjections.end(); it++ ) {
        if ( (it->first == index) && (it->second == frameIndex) ) {
            isPresent = true;
            break;
        }
    }
    return isPresent;
}

/**  */
EulerAngleType CMicrograph::GetTiltAngle() {
    return mTiltAngle;
}

/** */
void CMicrograph::SetTiltAngle( EulerAngleType angle ) {
    mTiltAngle = angle;
}

/**  */
EulerAngleType CMicrograph::GetTiltAngleUpperBound() {
    return mTiltAngleUpperBound;
}

/**  */
EulerAngleType CMicrograph::GetTiltAngleLowerBound() {
    return mTiltAngleLowerBound;
}

/**  */
EulerAngleType CMicrograph::GetTiltAxisAngle() {
    return mTiltAxisAngle;
}

/**  */
void CMicrograph::SetTiltAxisAngle( EulerAngleType angle ) {
    mTiltAxisAngle = angle;
}

/**  */
CPosition2D CMicrograph::GetTiltAxisShift() {
    return mTiltAxisShift;
}

/**  */
void CMicrograph::SetTiltAxisShift( CPosition2D position ) {
    mTiltAxisShift.SetPosition( position.mCoordX, position.mCoordY );
}

/** */
CPosition2D CMicrograph::GetFrameShift() {
    return mFrameShift;
}

/** */
void CMicrograph::SetFrameShift( CPosition2D position ) {
    mFrameShift.SetPosition( position.mCoordX, position.mCoordY );
}

void CMicrograph::SetTiltAngleTolerance( EulerAngleType tolerance ) {
    mTiltAngleTolerance = tolerance;
}

EulerAngleType CMicrograph::GetTiltAngleTolerance() {
    return mTiltAngleTolerance;
}
void CMicrograph::SetTiltDefocusOffset( double tiltDefocus1Offset, double tiltDefocus2Offset ) {
    mTiltSeriesDefocus1Offset = tiltDefocus1Offset;
    mTiltSeriesDefocus2Offset = tiltDefocus2Offset;
}
double CMicrograph::GetTiltDefocus1Offset() {
    return mTiltSeriesDefocus1Offset;
}
double CMicrograph::GetTiltDefocus2Offset() {
    return mTiltSeriesDefocus2Offset;
}
void CMicrograph::SetTiltAstigmatism( double tiltAstigmatism ) {
    mTiltSeriesAstigmatism = tiltAstigmatism;
}
double CMicrograph::GetTiltAstigmatism() {
    return mTiltSeriesAstigmatism;
}
/************************************************************************
 * TiltSeries class definition
 ************************************************************************/

/** Default constructor */
CTiltSeries::CTiltSeries() {
  mTiltSeriesIndex = -1;
}

/** Destructor */
CTiltSeries::~CTiltSeries() {
}

/** Alternate constructor */
CTiltSeries::CTiltSeries( int index ) {
    mTiltSeriesIndex = index;
}

/**  */
int CTiltSeries::GetNumberOfMicrographs() {
    return mListOfMicrographs.size();
}

/**  */
int CTiltSeries::GetNumberOfParticles() {
    return mListOfParticles.size();
}

/**  */
int CTiltSeries::GetIndex() {
    return mTiltSeriesIndex;
}

/**  */
void CTiltSeries::SetIndex( int index ) {
    mTiltSeriesIndex = index;
}

/** */
void CTiltSeries::SetMaxFrameIndex( int frameIndex ) {
    if ( frameIndex > mMaxFrameIndex ) {  
        mMaxFrameIndex = frameIndex;
    }
}

/**  */
void CTiltSeries::PrintListOfParticles() {
    
    set<int>::iterator it;
    
    for ( it = mListOfParticles.begin(); it != mListOfParticles.end(); it++ ) {
        printf( "%d ", *it );
    }
}

/**  */
void CTiltSeries::PrintListOfMicrographs() {
    
    set<int>::iterator it;
    
    for ( it = mListOfMicrographs.begin(); it != mListOfMicrographs.end(); it++ ) {
        printf( "%d ", *it );
    }

}

/**  */
int CTiltSeries::AddMicrograph( CMicrograph* p_Micrograph ) {
    mMicrograph.push_back( p_Micrograph );
    mListOfMicrographs.insert( p_Micrograph->GetIndex() );
    //mListOfMicrographsFrames.push_back( p_Micrograph->GetFrameIndex() );
    mMicrographByIndex[ {p_Micrograph->GetIndex(), p_Micrograph->GetFrameIndex()} ] = mMicrograph.size() - 1;
    return mListOfMicrographs.size();
}

/**  */
int CTiltSeries::AddParticle( CParticle* p_Particle ) {
    
    int index = p_Particle->GetIndex();
    int frameIndex = p_Particle->GetFrameIndex();
    if (!this->HasParticle( index, frameIndex )) {
        mParticle.push_back( p_Particle );
        mListOfParticles.insert( index );
        //mListOfParticlesFrames.push_back( frameIndex );
        mParticleByIndex[ {index, frameIndex} ] = mParticle.size() - 1;
    }
    return mListOfParticles.size();
}


int CTiltSeries::GetAllPositionsByMicrograph( int mInd, vector<int> & positions ) {
    int numberOfProjections = 0;
    map<pair<int,int>, int>::iterator it;
    for ( it = mMicrographByIndex.begin(); it != mMicrographByIndex.end(); it++ ) {
        if ( it->first.first == mInd ) {
            positions.push_back( mMicrographByIndex[it->first] );
            numberOfProjections += mMicrograph.at( mMicrographByIndex[it->first] )->mParticleProjectionSet.size();
        }
    }
    return numberOfProjections;
}

int CTiltSeries::GetAllPositionsByMicrograph( int mInd, int frameIndex, vector<int> & positions ) {
    int numberOfProjections = 0;
    map<pair<int,int>, int>::iterator it;
    for ( it = mMicrographByIndex.begin(); it != mMicrographByIndex.end(); it++ ) {
        if ( (it->first.first == mInd) && (it->first.second == frameIndex) ) {
            positions.push_back( mMicrographByIndex[it->first] );
            numberOfProjections += mMicrograph.at( mMicrographByIndex[it->first] )->mParticleProjectionSet.size();
        }
    }
    return numberOfProjections;
}

int CTiltSeries::GetAllPositionsByParticle( int pInd, vector<int> & positions ) {
    int numberOfProjections = 0;
    map<pair<int,int>, int>::iterator it;
    for ( it = mParticleByIndex.begin(); it != mParticleByIndex.end(); it++ ) {
        if ( it->first.first == pInd ) {
            positions.push_back( mParticleByIndex[it->first] );
            numberOfProjections += mParticle.at( mParticleByIndex[it->first] )->mParticleProjectionSet.size();
        }
    }
    return numberOfProjections;
}

int CTiltSeries::GetAllPositionsByParticle( int pInd, int frameIndex, vector<int> & positions ) {
    int numberOfProjections = 0;
    map<pair<int,int>, int>::iterator it;
    for ( it = mParticleByIndex.begin(); it != mParticleByIndex.end(); it++ ) {
        if ( (it->first.first == pInd) && (it->first.second == frameIndex) ) {
            positions.push_back( mParticleByIndex[it->first] );
            numberOfProjections += mParticle.at( mParticleByIndex[it->first] )->mParticleProjectionSet.size();
        }
    }
    return numberOfProjections;
}


/**
 * Add a particle projection, given by <particleProjection>, into the
 * micrograph <micrographIndex> of this tilt series. If the micrographs is not
 * included in the tilt series, it is added and then the particle projection
 * is added to it. If the micrograph exists and already has the particle
 * projection given by particleProjection.mParticleIndex then send a warning
 * message. In other case, add the particleProjection to the selected
 * micrograph. [REVISAR --fefo]
 * */
int CTiltSeries::AddParticleProjection( CParticleProjection* p_ParticleProjection ) {

    CMicrograph* p_Micrograph;
    CParticle* p_Particle;
    int pInd = p_ParticleProjection->GetParticleIndex();
    int mInd = p_ParticleProjection->GetMicrographIndex();
    int frameInd = p_ParticleProjection->GetFrameIndex();
    
    //Check if the tilt series has the micrograph <micrographIndex>
    if (!this->HasMicrograph( mInd, frameInd )) {
        return -1;
    } 

    //Check if the tilt series has the particle.
    if (!this->HasParticle( pInd, frameInd )) {
        return -1;
    }

    mListOfParticlesProjections.push_back( p_ParticleProjection->mIndexInStack );
    mParticleProjectionSet.push_back( p_ParticleProjection );

    // Get the index micrograph and the particle.
    p_Micrograph = this->mMicrograph.at( this->GetMicrographPositionByIndex( mInd, frameInd ) );
    p_Particle = this->mParticle.at( this->GetParticlePositionByIndex( pInd, frameInd) );
    //printf("(%d, %d) is in %d\n", pInd, frameInd, this->GetParticlePositionByIndex( pInd, frameInd));
    // If the micrograph does not have this particle projection add it. If does
    // send a warning message.
    if (!p_Micrograph->HasParticleProjection( pInd, frameInd )) {
        p_Micrograph->AddParticleProjection( p_ParticleProjection );
    } 
    else {
        printf( "[WARNING] Particle %d already has a projection into micrograph %d."
        " Problems with the .par file?\n", pInd, mInd );
    }

    if (!p_Particle->HasProjectionInMicrograph( mInd, frameInd )) {
        p_Particle->AddProjection( p_ParticleProjection );
    } 
    else {
        printf( "[WARNING] Particle %d already has a projection into micrograph %d.\n", pInd, mInd );
    }

    return pInd;
}

/**  */
void CTiltSeries::ShowInformation() {
    printf( "[ Tilt series %03d ]\n", mTiltSeriesIndex );

    printf( "Particles present (%lu): ", mListOfParticles.size() );
    this->PrintListOfParticles();
    printf( "\n" );

    printf( "Micrographs present (%lu): ", mListOfMicrographs.size() );
    this->PrintListOfMicrographs();
    printf( "\n" );

    for ( unsigned int indM = 0; indM < this->mMicrograph.size(); indM++ ) {
        //this->PrintListOfParticles();
        printf( "\n# [ Micrograph %03d ]\n", this->mMicrograph.at( indM )->mMicrographIndex );
        printf( "Tilt Angle: %lf\n", this->mMicrograph.at( indM )->GetTiltAngle() );
        printf( "Tilt Axis Angle: %lf\n", this->mMicrograph.at( indM )->GetTiltAxisAngle() );
        // printf( "Tilt Angle Tolerance: %lf\n", this->mMicrograph.at( indM )->GetTiltAngleTolerance() );
        printf( "Tilt Angle Tolerance: [%lf , %1f]\n", this->mMicrograph.at( indM )->GetTiltAngleLowerBound(), this->mMicrograph.at( indM )->GetTiltAngleUpperBound() );
        printf( "Particles present in this micrograph (%lu): ",
                this->mMicrograph.at( indM )->mListOfParticleProjections.size() );
        this->mMicrograph.at( indM )->PrintListOfParticleProjections();
        printf( "\n" );

        for ( int indPP = 0; indPP < this->mMicrograph.at( indM )->GetNumberOfParticleProjections(); indPP++ ) {
        printf( "## " );
        this->mMicrograph.at( indM )->mParticleProjectionSet.at( indPP )->ShowInformation();
        }
    }

    //printf( "\n\n# [ Particles information ]\n" );
    for ( unsigned int indP = 0; indP < this->mParticle.size(); indP++ ) {
        printf( "\n# " );
        this->mParticle.at( indP )->ShowInformation();
    }
    printf( "\n\n" );
}

/**  */
bool CTiltSeries::HasMicrograph( int index, int frameIndex ) {
    bool isPresent = false;
    map<pair<int,int>, int>::iterator it;
    /**
    for ( it = mListOfMicrographs.begin(); it < mListOfMicrographs.end(); it++ ) {
        if (index == *it) isPresent = true;
    }*/
    for ( it = mMicrographByIndex.begin(); it != mMicrographByIndex.end(); it++ ) {
        if ( (it->first.first == index) && (it->first.second == frameIndex) ) {
            isPresent = true;
            break;
        }
    }
    return isPresent;
}

/**
 * Process micrograph's tilt angles setting the lower and upper
 * bounds for their variations given the angles between them.
 * */
void CTiltSeries::ComputeMicrographTiltAngleBounds( bool display ) {
    // TODO(5): This implementation assumes that the micrograph at 0 in
    // pTiltSeries->mMicrograph has the smallest tilt angle and they
    // are ordered increasing the tilt angle.
    unsigned int index;
    EulerAngleType alpha;
    EulerAngleType beta;
    float fraction = 0.45;
    CMicrograph* pMicrograph;
    
    if ( display ) PrintDebug( DebugInfo, "Tilt angle search range for Micrographs\n" );

    // cout << "Micrograph size = " << this->mMicrograph.size() << endl;

    index = 0;
    //pMicrograph = this->mMicrograph.at( index );
    pMicrograph = this->mMicrograph.at( mMicrographByIndex[{index, 0}] );

    if ( this->mMicrograph.size() == 1 ) {
       pMicrograph->mTiltAngleUpperBound = pMicrograph->GetTiltAngle() + fraction * 2;
        pMicrograph->mTiltAngleLowerBound = 1.1 * pMicrograph->GetTiltAngle();
        return;
    }

    //alpha = abs( pMicrograph->GetTiltAngle() - this->mMicrograph.at( index + 1 )->GetTiltAngle() );
    alpha = abs( pMicrograph->GetTiltAngle() - this->mMicrograph.at( mMicrographByIndex[{index+1, 0}] )->GetTiltAngle() );
    pMicrograph->mTiltAngleUpperBound = pMicrograph->GetTiltAngle() + fraction * alpha;
    pMicrograph->mTiltAngleLowerBound = 1.1 * pMicrograph->GetTiltAngle();
    
    // cout << "Micrograph index = " << pMicrograph->mMicrographIndex << endl;
    
    // PrintDebug( DebugInfo, "M = %02f\n", pMicrograph->mMicrographIndex );
    
    if ( display ) PrintDebug( DebugInfo, "M = %02d: ( %7.3f, %7.3f ) starting from %7.3f\n", pMicrograph->mMicrographIndex,
            pMicrograph->GetTiltAngleLowerBound(), pMicrograph->GetTiltAngleUpperBound(), pMicrograph->GetTiltAngle() );
    
    for ( int unsigned cnt = 1; cnt < this->mListOfMicrographs.size()-1; cnt++ ) {
    //for ( int unsigned cnt = 1; cnt < ( this->mMicrograph.size() - 1 ); cnt++ ) {
        index = cnt; //pTiltSeries->GetMicrographPositionByIndex( cnt );
        pMicrograph = this->mMicrograph.at( mMicrographByIndex[{index, 0}] );

        alpha = abs( pMicrograph->GetTiltAngle() - this->mMicrograph.at( mMicrographByIndex[{index+1, 0}] )->GetTiltAngle() );
        pMicrograph->mTiltAngleUpperBound = pMicrograph->GetTiltAngle() + fraction * alpha;

        beta = abs( pMicrograph->GetTiltAngle() - this->mMicrograph.at( mMicrographByIndex[{index-1, 0}] )->GetTiltAngle() );
        pMicrograph->mTiltAngleLowerBound = pMicrograph->GetTiltAngle() - fraction * beta;
        //printf( "%d = ( %f -> %f  )\n", pMicrograph->mMicrographIndex, pMicrograph->mTiltAngleUpperBound , pMicrograph->mTiltAngleLowerBound );
        if ( display ) PrintDebug( DebugInfo, "M = %02d: ( %7.3f, %7.3f ) starting from %7.3f\n", pMicrograph->mMicrographIndex,
                pMicrograph->GetTiltAngleLowerBound(), pMicrograph->GetTiltAngleUpperBound(), pMicrograph->GetTiltAngle() );
    }

    //index = this->mMicrograph.size() - 1;
    //pMicrograph = this->mMicrograph.at( index );
    index = this->mListOfMicrographs.size() - 1;
    pMicrograph = this->mMicrograph.at( mMicrographByIndex[{index, 0}] );
    
    beta = abs( pMicrograph->GetTiltAngle() - this->mMicrograph.at( mMicrographByIndex[{index-1, 0}] )->GetTiltAngle() );
    pMicrograph->mTiltAngleUpperBound = 1.1 * pMicrograph->GetTiltAngle();
    pMicrograph->mTiltAngleLowerBound = pMicrograph->GetTiltAngle() - fraction * beta;

    if ( display ) PrintDebug( DebugInfo, "M = %02d: ( %7.3f, %7.3f ) starting from %7.3f\n", pMicrograph->mMicrographIndex,
            pMicrograph->GetTiltAngleLowerBound(), pMicrograph->GetTiltAngleUpperBound(), pMicrograph->GetTiltAngle() );

}

/** Returns the index in the vector mMicrograph of the micrograph given
 * by a particular mMicrographIndex
 * */
int CTiltSeries::GetMicrographPositionByIndex( int index, int frameIndex ) {
    //  if (( index < 1 ) || ( index > mListOfMicrographs.size() )) {
    //    printf( "[ERROR] GetMicrographPositionByIndex: Index (%d) is out of bound (1-%d).\n", index, mListOfMicrographs.size() );
    //  }
    //
    if ( this->HasMicrograph( index, frameIndex ) ) {
        return mMicrographByIndex[{index, frameIndex}];
    }
    else {
        return -1;
    }
}

/** Returns true if the particle <index> is present in the tilt series. */
bool CTiltSeries::HasParticle( int index, int frameIndex ) {
    bool isPresent = false;
    map<pair<int, int>, int>::iterator it;
    /**
    for ( it = mListOfParticles.begin(); it < mListOfParticles.end(); it++ ) {
        if (index == *it) isPresent = true;
    }
    */
    for ( it = mParticleByIndex.begin(); it != mParticleByIndex.end(); it++ ) {
        if ( (it->first.first == index) && (it->first.second == frameIndex) ) {
            isPresent = true;
            break;
        }
    }
    return isPresent;
}

/** Returns the index in the vector mParticle of the particles given
 * by a particular particle projection.
 * */
int CTiltSeries::GetParticlePositionByIndex( int index, int frameIndex ) {
    //  if (( index < 1 ) || ( index > mListOfParticles.size() )) {
    //    printf( "[ERROR] GetParticlePositionByIndex: Index (%d) is out of bound (1-%d).\n", index, mListOfParticles.size() );
    //  }
    if ( this->HasParticle( index, frameIndex ) ) {
        return mParticleByIndex[{index, frameIndex}];
    }
    else {
        return -1;
    }
}

/************************************************************************
 * CContainer class definition
 ************************************************************************/
// Global static pointer used to ensure a single instance of the class.
CContainer* CContainer::mpInstance = NULL;
/** This function is called to create an instance of the class.
 Calling the constructor publicly is not allowed. The constructor
 is private and is only called by this Instance function.
 */
CContainer* CContainer::Instance() {
    if (!mpInstance) // Only allow one instance of class to be generated.
    mpInstance = new CContainer();
    return mpInstance;
}

CContainer::CContainer() {
    extern SParameters params;
    // Count the number of times GetValue() is called.
    mNumberOfGetValueCalled = 0;

    //
    mCurrentOptimizationIndex = -1;

    // Initialize the parameter selection mask
    mSelectionMask = flrCostFunction::ParametersType( VARIABLES_SPACE_DIMENSION );
    if (params.mRefineMicrographTiltAngles) {
        mSelectionMask[0] = 1;
    } else {
        mSelectionMask[0] = 0;
    }
    if (params.mRefineMicrographTiltAxisAngles) {
        mSelectionMask[1] = 1;
    } else {
        mSelectionMask[1] = 0;
    }
    if (params.mRefineParticlesTheta) {
        mSelectionMask[2] = 1;
    } else {
        mSelectionMask[2] = 0;
    }
    if (params.mRefineParticlesPsi) {
        mSelectionMask[3] = 1;
    } else {
        mSelectionMask[3] = 0;
    }
    if (params.mRefineParticlesPhi) {
        mSelectionMask[4] = 1;
    } else {
        mSelectionMask[4] = 0;
    }

    // Initialize the tolerances, same order as the selection mask
    mTolerance[0] = params.mToleranceMicrographShifts;
    mTolerance[1] = params.mToleranceMicrographShifts;
    mTolerance[2] = params.mToleranceMicrographShifts;
    mTolerance[3] = params.mToleranceParticlesTheta;
    mTolerance[4] = params.mToleranceParticlesPsi;
    mTolerance[5] = params.mToleranceParticlesPhi;

    // Initialize the previous and current values
    mPreviousValue = itk::NumericTraits<double>::max();
    mCurrentValue = itk::NumericTraits<double>::max();

    // Initialize the previous and current points
    mPreviousPoint = flrCostFunction::ParametersType( VARIABLES_SPACE_DIMENSION );
    mCurrentPoint = flrCostFunction::ParametersType( VARIABLES_SPACE_DIMENSION );
    mInitialPoint = flrCostFunction::ParametersType( VARIABLES_SPACE_DIMENSION );
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        mPreviousPoint[k] = itk::NumericTraits<double>::max();
        mCurrentPoint[k] = itk::NumericTraits<double>::max();
        mInitialPoint[k] = itk::NumericTraits<double>::max();
    }

}

// Set and get the current index of the object being optimized.
void SetCurrentOptimizationIndex( int );
int GetCurrentOptimizationIndex();

void CContainer::SetTiltSeries( CTiltSeries* tiltSeries ) {
    mTiltSeriesPointer = tiltSeries;
}

CTiltSeries* CContainer::GetTiltSeries() {
    return mTiltSeriesPointer;
}

void CContainer::SetPreviousPoint( flrCostFunction::ParametersType point ) {
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        mPreviousPoint[k] = point[k];
    }
}

void CContainer::SetCurrentPoint( flrCostFunction::ParametersType point ) {
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        mCurrentPoint[k] = point[k];
    }
}

flrCostFunction::ParametersType CContainer::GetPreviousPoint() {
    return mPreviousPoint;
}

flrCostFunction::ParametersType CContainer::GetCurrentPoint() {
    return mCurrentPoint;
}

void CContainer::SetInitialPoint( flrCostFunction::ParametersType point ) {
  for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
    mInitialPoint[k] = point[k];
  }
}

flrCostFunction::ParametersType CContainer::GetInitialPoint() {
  return mInitialPoint;
}

void CContainer::SetBestPoint( flrCostFunction::ParametersType point ) {
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        mBestPoint[k] = point[k];
    }
}

flrCostFunction::ParametersType CContainer::GetBestPoint() {
    return mBestPoint;
}

void CContainer::SetPreviousValue( flrCostFunction::MeasureType value ) {
    mPreviousValue = value;
}

flrCostFunction::MeasureType CContainer::GetPreviousValue() {
    return mPreviousValue;
}

void CContainer::SetInitialValue( flrCostFunction::MeasureType value ) {
  mInitialValue = value;
}

flrCostFunction::MeasureType CContainer::GetInitialValue() {
    return mInitialValue;
}

void CContainer::SetCurrentValue( flrCostFunction::MeasureType value ) {
    mCurrentValue = value;
}

flrCostFunction::MeasureType CContainer::GetCurrentValue() {
    return mCurrentValue;
}

/** Gets in which direction is it moving (varies) respect to the initial point */
int CContainer::GetDescentDirection() {
    int direction = 0;
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        if (mCurrentPoint[k] != mInitialPoint[k]) {
        direction = k + 1;
        }
    }
    return direction;
}

/** True if the previous and current point are the same. */
bool CContainer::IsSamePoint() {
    bool same = true;
    for ( unsigned int k = 0; k < VARIABLES_SPACE_DIMENSION; k++ ) {
        if (mCurrentPoint[k] != mPreviousPoint[k]) {
        same = false;
        }
    }
    return same;
}

void CContainer::SetCurrentOptimizationIndex( int index ) {
    mCurrentOptimizationIndex = index;
}

int CContainer::GetCurrentOptimizationIndex() {
    return mCurrentOptimizationIndex;
}



/************************************************************************
 * Image class
 ************************************************************************/

/** Default constructor */
CImage::CImage() {
}

/** Alternate constructor  */
CImage::CImage(Image** image, int number_of_images, bool is_raw) {
    imageStack = image; 
    numImages = number_of_images;
    isRaw = is_raw;
}

/** Destructor */
CImage::~CImage() {
    DeleteImageArray( imageStack, numImages );
}

void CImage::SetImage(Image** image) {
    imageStack = image; 
}

bool CImage::IsRaw() {
    return isRaw; 
}

