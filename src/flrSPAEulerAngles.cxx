/*=========================================================================

 Name:        flrSPAEulerAngles.cxx
 Category:    Implementation file
 Description: This is the file with the function that computes the mapping
 between the CSP and SP manifolds.
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrSPAEulerAngles.cxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/

#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <limits>
#include <ctime>
#include <iostream>


#include "flrMain.hxx"
#include "flrReadParametersFile.hxx"

vtkMatrix4x4* eulerZXZtoZYZ( vtkMatrix4x4* matrixZXZ ) {
    
    double z1, x, z2, y;
    double pi = atan2(0, -1);

    if ( matrixZXZ->GetElement( 2, 2 ) < 1 ) {
        if ( matrixZXZ->GetElement( 2, 2 ) > -1 ) {
            x = acos( matrixZXZ->GetElement( 2, 2 ) );
            z1 = atan2( (matrixZXZ->GetElement( 0, 2 )/sin(x)), (matrixZXZ->GetElement( 1, 2 )/sin(x)) );
            z2 = atan2( (matrixZXZ->GetElement( 2, 0 )/sin(x)), - (matrixZXZ->GetElement( 2, 1 )/sin(x)) );
        }
        else {
            x = pi;
            z2 = atan2( matrixZXZ->GetElement( 0, 1 ), matrixZXZ->GetElement( 0, 0 ) );
            z1 = 0;
        }
    }
    else {
        x = 0;
        z1 = 0;
        z2 = atan2(  matrixZXZ->GetElement( 0, 1 ), matrixZXZ->GetElement( 0, 0 ) );
    }
    
    double z1_degree = vtkMath::DegreesFromRadians( z1 );
    double z2_degree = vtkMath::DegreesFromRadians( z2 );
    double x_degree = vtkMath::DegreesFromRadians( x );

    vtkTransform * t = vtkTransform::New();
    t->PostMultiply();

    t->RotateZ( z2_degree );
    t->RotateX( x_degree );
    t->RotateZ( z1_degree );

    vtkMatrix4x4 * m = t->GetMatrix();
    
    if ( m->GetElement( 2, 2 ) < 1 ) {
        if ( m->GetElement( 2, 2 ) > -1 ) {
            y = acos( m->GetElement( 2, 2 ) );
            z1 = atan2( (m->GetElement( 1, 2 )/sin(y)), (m->GetElement( 0, 2 )/sin(y)) );
            z2 = atan2( (m->GetElement( 2, 1 )/sin(y)), - (m->GetElement( 2, 0 )/sin(y)) );
        }
        else {
            y = pi;
            z1 = - atan2( m->GetElement( 1, 0 ), m->GetElement( 1, 1 ) );
            z2 = 0;
        }
    }
    else {
        y = 0;
        z2 = 0;
        z1 = atan2(  m->GetElement( 1, 0 ), matrixZXZ->GetElement( 1, 1 ) );
    }
    
    z1_degree = vtkMath::DegreesFromRadians( z1 );
    z2_degree = vtkMath::DegreesFromRadians( z2 );
    double y_degree = vtkMath::DegreesFromRadians( y );


    vtkTransform * t2 = vtkTransform::New();
    t2->PostMultiply();
    
    t2->RotateZ( - z2_degree );
    t2->RotateY( - y_degree );
    t2->RotateZ( - z1_degree );
    
    //vtkMatrix4x4 * ret = t2->GetMatrix();
    
    matrixZXZ->DeepCopy( t2->GetMatrix() );
    t->Delete();
    t2->Delete();

    //return ret;
}

vtkMatrix4x4* eulerTwoZYZtoOneZYZ( vtkMatrix4x4* matrixZYZZYZ ) {

    double y, z1, z2; 
    double pi = atan2( 0, -1 );

    if ( matrixZYZZYZ->GetElement( 2, 2 ) < 1 ) {
        if ( matrixZYZZYZ->GetElement( 2, 2 ) > -1 ) {
            y = acos( matrixZYZZYZ->GetElement( 2, 2 ) );
            z2 = atan2( (matrixZYZZYZ->GetElement( 2, 1 )/sin(y)), (matrixZYZZYZ->GetElement( 2, 0 )/sin(y)) );
            z1 = atan2( (matrixZYZZYZ->GetElement( 1, 2 )/sin(y)), - (matrixZYZZYZ->GetElement( 0, 2 )/sin(y)) );
        }
        else {
            y = pi;
            z1 = atan2( - matrixZYZZYZ->GetElement( 0, 1 ), - matrixZYZZYZ->GetElement( 0, 0 ) );
            z2 = 0;
        }
    }
    else {
        y = 0;
        z2 = 0;
        z1 = atan2(  matrixZYZZYZ->GetElement( 0, 1 ), matrixZYZZYZ->GetElement( 0, 0 ) );
    }

    double z1_degree = vtkMath::DegreesFromRadians( z1 );
    double z2_degree = vtkMath::DegreesFromRadians( z2 );
    double y_degree = vtkMath::DegreesFromRadians( y );

    vtkTransform * t = vtkTransform::New();
    t->PostMultiply();

    t->RotateZ( - z2_degree );
    t->RotateY( - y_degree );
    t->RotateZ( - z1_degree );
    
    //vtkMatrix4x4 * ret = t->GetMatrix();
    matrixZYZZYZ->DeepCopy( t->GetMatrix() );
    t->Delete();

    //return ret;


}

void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] ) {
    /**
    std::clock_t start;
    double duration;
    start = std::clock();
    */

    EulerAngleType refinematrix[] = { matrix[0], //
                                      matrix[1], //
                                      matrix[2], //
                                      matrix[3], //
                                      matrix[4], //
                                      matrix[5], //
                                      matrix[6], //
                                      matrix[7], //
                                      matrix[8], //
                                      matrix[9], //
                                      matrix[10], //
                                      matrix[11], //
                                      0, //
                                      0, //
                                      0, //
                                      1  };
    
    
    EulerAngleType rotationmatrix_reverse[] = { matrix[0], //
                                      -matrix[4], //
                                      matrix[8], //
                                      0, //
                                      -matrix[1], //
                                      matrix[5], //
                                      -matrix[9], //
                                      0, //
                                      matrix[2], //
                                      -matrix[6], //
                                      matrix[10], //
                                      0, //
                                      0, //
                                      0, //
                                      0, //
                                      1 };
    /**
    EulerAngleType rotationmatrix_reverse[] = { matrix[0], //
                                      -matrix[4], //
                                      -matrix[8], //
                                      0, //
                                      -matrix[1], //
                                      matrix[5], //
                                      matrix[9], //
                                      0, //
                                      -matrix[2], //
                                      matrix[6], //
                                      matrix[10], //
                                      0, //
                                      0, //
                                      0, //
                                      0, //
                                      1 };
    */
    EulerAngleType rotationmatrix[] = { matrix[0], //
                                      matrix[1], //
                                      matrix[2], //
                                      0, //
                                      matrix[4], //
                                      matrix[5], //
                                      matrix[6], //
                                      0, //
                                      matrix[8], //
                                      matrix[9], //
                                      matrix[10], //
                                      0, //
                                      0, //
                                      0, //
                                      0, //
                                      1 };
    
    // 3DAVG refinement transformation (rotation component only)
    vtkMatrix4x4 * refinement = vtkMatrix4x4::New();
    refinement->DeepCopy( refinematrix );
    
    vtkMatrix4x4 * refinementRotation = vtkMatrix4x4::New();
    refinementRotation->DeepCopy( rotationmatrix );
    
    vtkMatrix4x4 * refinementRotation_reverse = vtkMatrix4x4::New();
    refinementRotation_reverse->DeepCopy( rotationmatrix_reverse );
    
    // Rotation
    vtkTransform * imod = vtkTransform::New();
    imod->PostMultiply();

    imod->RotateZ( -tilt_axis_angle );
    imod->RotateY( tilt_angle );
    
    //
    vtkMatrix4x4 * imod_matrix = vtkMatrix4x4::New();
    imod_matrix->DeepCopy( imod->GetMatrix() );


    vtkTransform * n = vtkTransform::New();
    n->PostMultiply();
    
    n->RotateZ( normal[2] );
    n->RotateX( normal[0] );
    n->RotateZ( normal[1] );
    
    //
    vtkMatrix4x4 * norm_matrix = vtkMatrix4x4::New();
    norm_matrix->DeepCopy( n->GetMatrix() );

    vtkMatrix4x4 * norm_and_rot = vtkMatrix4x4::New();
    vtkMatrix4x4 * rot_invert =  vtkMatrix4x4::New();
    vtkMatrix4x4 * rot_invert_trans =  vtkMatrix4x4::New();
    
    vtkMatrix4x4::Invert( refinementRotation, rot_invert );  
    rot_invert_trans->DeepCopy( rot_invert );


    eulerZXZtoZYZ( rot_invert );
    eulerZXZtoZYZ( norm_matrix );

    //vtkMatrix4x4::Multiply4x4( eulerZXZtoZYZ( rot_invert ), eulerZXZtoZYZ( n->GetMatrix() ), norm_and_rot );
    vtkMatrix4x4::Multiply4x4( rot_invert , norm_matrix, norm_and_rot );
    

    vtkMatrix4x4 * tmp = vtkMatrix4x4::New();
    vtkTransform * local = vtkTransform::New();
    local->PostMultiply();

    
    local->RotateZ( psi );
    local->RotateY( theta );
    local->RotateZ( phi );
    
    //
    vtkMatrix4x4 * local_matrix = vtkMatrix4x4::New();
    local_matrix->DeepCopy( local->GetMatrix() );
    // eulerZXZtoZYZ( local_matrix );
    
    //vtkMatrix4x4::Multiply4x4( eulerZXZtoZYZ( local->GetMatrix() ), norm_and_rot, tmp );
    vtkMatrix4x4::Multiply4x4( local_matrix, norm_and_rot, tmp );

    vtkMatrix4x4 * r = vtkMatrix4x4::New();
    eulerTwoZYZtoOneZYZ( tmp );
    //vtkMatrix4x4::Multiply4x4( eulerTwoZYZtoOneZYZ( tmp ), imod->GetMatrix(), r );
    vtkMatrix4x4::Multiply4x4( tmp, imod_matrix, r );
    

    // Translation 
    vtkMatrix4x4 * refinementTranslation = vtkMatrix4x4::New();
    vtkMatrix4x4::Multiply4x4( rot_invert_trans, refinement, refinementTranslation );
    refinementTranslation->SetElement( 0, 3, -refinementTranslation->GetElement( 0, 3 ) ); 
    refinementTranslation->SetElement( 2, 3, -refinementTranslation->GetElement( 2, 3 ) );
    
    vtkMatrix4x4 * trans_invert = vtkMatrix4x4::New();
    vtkMatrix4x4::Invert( refinementTranslation, trans_invert );
    vtkMatrix4x4 * refinement_rotation_reverse_invert = vtkMatrix4x4::New();
    vtkMatrix4x4::Invert( refinementRotation_reverse, refinement_rotation_reverse_invert );
    
    vtkMatrix4x4 * trans_tmp =  vtkMatrix4x4::New();
    vtkMatrix4x4 * trans_tmp1 =  vtkMatrix4x4::New();
    vtkTransform * trans_local = vtkTransform::New();
    trans_local->PostMultiply();
    trans_local->RotateZ( 0 );
    trans_local->RotateY( 0 );
    trans_local->RotateZ( 0 );
    vtkMatrix4x4::Multiply4x4( trans_local->GetMatrix(), trans_invert, trans_tmp );
    vtkMatrix4x4::Multiply4x4( refinement_rotation_reverse_invert, trans_tmp, trans_tmp1 );
    
    vtkTransform * nt = vtkTransform::New();
    nt->PostMultiply();
    nt->RotateZ( normal[1] );
    nt->RotateX( normal[0] );
    nt->RotateZ( normal[2] );

    vtkMatrix4x4 * trans_tmp2 =  vtkMatrix4x4::New();
    vtkMatrix4x4::Multiply4x4( nt->GetMatrix(), trans_tmp1, trans_tmp2);
    
    vtkMatrix4x4 * t =  vtkMatrix4x4::New();
    vtkTransform * trans_imod = vtkTransform::New();
    trans_imod->PostMultiply();
    trans_imod->RotateY( tilt_angle );
    trans_imod->RotateZ( - tilt_axis_angle );
    vtkMatrix4x4::Multiply4x4( trans_imod->GetMatrix(), trans_tmp2, t);

    // decompose rotation matrix into three Euler angles
    double pi = atan2( 0, -1 );

    if ( r->GetElement( 2, 2 ) < 1 ) {
        if ( r->GetElement( 2, 2 ) > -1 ) {
            theta = acos( r->GetElement( 2, 2 ) );
            psi = atan2( (r->GetElement( 2, 1 )/sin(theta)), (r->GetElement( 2, 0 )/sin(theta)) );
            phi = atan2( (r->GetElement( 1, 2 )/sin(theta)), -(r->GetElement( 0, 2 )/sin(theta)) );
        }
        else {
            theta = pi;
            phi = atan2( - r->GetElement( 0, 1 ), - r->GetElement( 0, 0 ) );
            psi = 0;
        }
    }
    else {
        theta = 0;
        phi = atan2( r->GetElement( 0, 1 ), r->GetElement( 0, 0 ) );
        psi = 0;
    }
    
    frealign[0] = vtkMath::DegreesFromRadians( psi );
    frealign[1] = vtkMath::DegreesFromRadians( theta );
    frealign[2] = vtkMath::DegreesFromRadians( phi );
    


    // frealign does not use negative angles, so we add 360 to each negative angle
    frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
    frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
    frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];
    frealign[0] == 360.0 ? frealign[0] = 0.0 : frealign[0];
    frealign[1] == 360.0 ? frealign[1] = 0.0 : frealign[1];
    frealign[2] == 360.0 ? frealign[2] = 0.0 : frealign[2];

    shift[0] = t->GetElement( 0, 3 );
    shift[1] = t->GetElement( 1, 3 );
    
    
    // delete intermediate matrix
    refinement->Delete();
    refinementRotation->Delete();
    refinementRotation_reverse->Delete();
    imod->Delete();
    n->Delete();
    norm_and_rot->Delete();
    rot_invert->Delete();
    tmp->Delete();
    local->Delete();
    
    refinementTranslation->Delete();
    trans_invert->Delete();
    refinement_rotation_reverse_invert->Delete();
    trans_tmp->Delete();
    trans_tmp1->Delete();
    trans_local->Delete();
    nt->Delete();
    trans_tmp2->Delete();
    trans_imod->Delete();

    norm_matrix->Delete();
    local_matrix->Delete();
    rot_invert_trans->Delete();
    imod_matrix->Delete();

    // delete product rotation and translation matrix
    t->Delete();
    r->Delete();
    

}
/**
void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] ) {

  EulerAngleType rotationmatrix[] = { matrix[0], //
                                      -matrix[4], //
                                      matrix[8], //
                                      0, //
                                      -matrix[1], //
                                      matrix[5], //
                                      -matrix[9], //
                                      0, //
                                      matrix[2], //
                                      -matrix[6], //
                                      matrix[10], //
                                      0, //
                                      matrix[12], //
                                      matrix[13], //
                                      matrix[14], //
                                      matrix[15] };
     
    //EulerAngleType rotationmatrix[] = { matrix[0], //
    //                                  matrix[1], //
    //                                  matrix[2], //
     //                                 0, //
    //                                  matrix[4], //
    //                                  matrix[5], //
    //                                  matrix[6], //
    //                                  0, //
    //                                  matrix[8], //
    //                                  matrix[9], //
    //                                  matrix[10], //
    //                                  0, //
    //                                  matrix[12], //
    //                                  matrix[13], //
    //                                  matrix[14], //
    //                                  matrix[15] };
  // 3DAVG refinement transformation (rotation component only)
  vtkMatrix4x4 * refinementRotation = vtkMatrix4x4::New();
  refinementRotation->DeepCopy( rotationmatrix );

  // transformation matrix
  vtkTransform * t = vtkTransform::New();
  t->PostMultiply();

  // tilt axis angle rotation
  t->RotateZ( tilt_axis_angle );

  // apply tilt angle rotation
  t->RotateY( -tilt_angle );

  // apply spike euler angles (pure rotation R1)
  t->RotateZ( -normal[2] + phi );
  t->RotateX( -normal[0] + theta );
  t->RotateZ( -normal[1] + psi );
  
  // apply 3DAVG refinement transformation (rotation only, R2)
  vtkMatrix4x4 * local = vtkMatrix4x4::New();
  vtkMatrix4x4::Multiply4x4( refinementRotation, t->GetMatrix(), local );
  t->SetMatrix( local );

  vtkMatrix4x4 * final = t->GetMatrix();

  if ( fabs( final->GetElement( 2, 2 ) ) < 1 ) {
    theta = acos( final->GetElement( 2, 2 ) );
    double num = final->GetElement( 2, 1 ) / sin( theta );
    double den = -final->GetElement( 2, 0 ) / sin( theta );
    psi = atan2( num, den );

    num = final->GetElement( 1, 2 ) / sin( theta );
    den = final->GetElement( 0, 2 ) / sin( theta );
    phi = atan2( num, den );
  } else {
	// check if m00 is valid:
	double m00 = final->GetElement( 0, 0 );
    if ( m00 < -1 ){
		m00 = -1;
		cerr << "WARNING : Non-orthogonal matrix, m00 =" << final->GetElement( 0, 0 ) << *refinementRotation << endl;
	}
	if ( m00 >  1 ){
		m00 =  1;
		cerr << "WARNING : Non-orthogonal matrix, m00 =" << final->GetElement( 0, 0 ) << *refinementRotation << endl;
	}
	
	phi = 0.0;
	theta = 0.0;
	psi = acos( m00 );
	if ( final->GetElement( 2, 2 ) < fabs( final->GetElement( 2, 2 ) ) ){
		theta = vtkMath::Pi();
		psi = acos( -m00 );
	}
  }

  frealign[0] = vtkMath::DegreesFromRadians( psi );
  frealign[1] = vtkMath::DegreesFromRadians( theta );
  frealign[2] = vtkMath::DegreesFromRadians( phi );

  // frealign does not use negative angles, so we add 360 to each negative angle
  frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
  frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
  frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];
  frealign[0] == 360.0 ? frealign[0] = 0.0 : frealign[0];
  frealign[1] == 360.0 ? frealign[1] = 0.0 : frealign[1];
  frealign[2] == 360.0 ? frealign[2] = 0.0 : frealign[2];

  shift[0] = 0;
  shift[1] = 0;
}
**/
