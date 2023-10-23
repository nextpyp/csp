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

#include "flrMain.hxx"
#include "flrReadParametersFile.hxx"

//using namespace std;

void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] ) {

  PrintDebug( DebugData, "\t[%10.6f %10.6f %10.6f %10.6f %10.6f ] => ", 5, tilt_angle, tilt_axis_angle, theta, psi, phi );

  /********************************************************************************************
   * ******************************************************************************************
   * ******************************************************************************************/
  //    // spike normal
  //    double normal[] = { 0, 0, 0 };
  //
  //    // refinement matrix
  //    double matrix[] = { 1, 0, 0, 0, //
  //                        0, 1, 0, 0, //
  //                        0, 0, 1, 0, //
  //                        0, 0, 0, 1 }; //
  //
  //    // Refinement matrix
  //    vtkMatrix4x4 * refinement = vtkMatrix4x4::New();
  //    refinement->DeepCopy( matrix );
  //
  //    // Transformation matrix
  //    vtkTransform * t = vtkTransform::New();
  //
  //    t->PostMultiply();
  //
  //    // Apply tilt axis angle rotation
  //    //t->RotateZ( tilt_axis_angle );
  //    t->RotateZ( tilt_axis_angle );
  //
  //    // Apply tilt angle rotation
  //    //t->RotateY( -tilt_angle );
  //    t->RotateY( -tilt_angle );
  //
  //    // Apply spike euler angles
  //    t->RotateZ( -normal[2] );
  //    t->RotateX( -normal[0] );
  //    t->RotateZ( -normal[1] );
  //
  //    // Apply refinement matrix
  //    //  t->Concatenate( refinement );
  //    t->RotateZ( phi );
  //    t->RotateY( theta );
  //    t->RotateZ( psi );
  //
  //    // Get the transformation matrix
  //    vtkMatrix4x4 * final = t->GetMatrix();
  //
  //    // Recover euler angles from transformation
  //    if (fabs( final->GetElement( 2, 2 ) ) < 1 - numeric_limits<float>::min()) {
  //      theta = acos( final->GetElement( 2, 2 ) );
  //      double num = final->GetElement( 2, 1 ) / sin( theta );
  //      double den = -final->GetElement( 2, 0 ) / sin( theta );
  //      psi = atan2( num, den );
  //      num = final->GetElement( 1, 2 ) / sin( theta );
  //      den = final->GetElement( 0, 2 ) / sin( theta );
  //      phi = atan2( num, den );
  //    } else {
  //      double sign = final->GetElement( 2, 2 ) / fabs( final->GetElement( 2, 2 ) );
  //      theta = vtkMath::Pi() * ( 1 - sign ) / 2.0;
  //      phi = 0.0;
  //      float arg = final->GetElement( 0, 0 ) / cos( theta );
  //      psi = acos( arg );
  //    }
  //
  //    // Convert to degrees
  //    frealign[0] = psi * vtkMath::RadiansToDegrees();
  //    frealign[1] = theta * vtkMath::RadiansToDegrees();
  //    frealign[2] = phi * vtkMath::RadiansToDegrees();
  //
  //    // Frealign does not use negative angles, so we add 360 to each
  //    // negative angle
  //    frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
  //    frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
  //    frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];
  //    frealign[0] == 360.0 ? frealign[0] = 0.0 : frealign[0];
  //    frealign[1] == 360.0 ? frealign[1] = 0.0 : frealign[1];
  //    frealign[2] == 360.0 ? frealign[2] = 0.0 : frealign[2];
  //
  //    // Retrieve shifts from final refinement
  //    double shifts[] = { final->GetElement( 0, 3 ), final->GetElement( 1, 3 ), final->GetElement( 2, 3 ) };
  //
  //    // Set translation component to zero
  //    final->SetElement( 0, 3, 0 );
  //    final->SetElement( 1, 3, 0 );
  //    final->SetElement( 2, 3, 0 );
  //
  //    double xaxis[] = { 1, 0, 0, 1 };
  //    double nxaxis[4];
  //    final->MultiplyPoint( xaxis, nxaxis );
  //
  //    double yaxis[] = { 0, 1, 0, 1 };
  //    double nyaxis[4];
  //    final->MultiplyPoint( yaxis, nyaxis );
  //
  //    double zaxis[] = { 0, 0, 1, 1 };
  //    double nzaxis[4];
  //    final->MultiplyPoint( zaxis, nzaxis );
  //
  //    final->Delete();
  //
  //    double x3axis[] = { nxaxis[0], nxaxis[1], nxaxis[2] };
  //    double y3axis[] = { nyaxis[0], nyaxis[1], nyaxis[2] };
  //
  //    double sx = vtkMath::Dot( shifts, x3axis );
  //    double sy = vtkMath::Dot( shifts, y3axis );
  //    shift[0] = sx;
  //    shift[1] = sy;

  /********************************************************************************************
   * ******************************************************************************************
   * ******************************************************************************************/
  // tilt angle
  //double tilt_angle = atof( argv[1] );

  // tilt axis rotation
  //double tilt_axis_angle = atof( argv[2] );

  //  cout << "Tilt-angle, tilt axis=" << tilt_angle << "," << tilt_axis_angle << endl;

  //double axis_correction = 0; //atof( argv[3] );
  double tilt_X_correction = 0; //atof( argv[4] );
  double tilt_Y_correction = 0; //atof( argv[5] );
  //  cout << "Shift corrections =[" << axis_correction << "," << tilt_X_correction << "," << tilt_Y_correction << "]" << endl;

  //// spike normal
  //double normal[] = { atof( argv[6] ), atof( argv[7] ), atof( argv[8] ) };
  //cout << "Normal = " << normal[0] << "," << normal[1] << "," << normal[2] << endl;

  // 3DAVG refinement
  // This transformation matrix is different from the matrix[0-15] in the 3DAVG refinement file.
  // The ordering of the rotations has to be changed to match the new geometry.
  //  - The ordering in 3DAVG is: RotZ1 * RotX * RotZ2.
  //  - The ordering here is: RotZ2 * RotX * RotZ1
  // The new transformation matrix has the following expression:

  //  double matrix[] = { atof( argv[9] ), //
  //                      -atof( argv[13] ), //
  //                      atof( argv[17] ), //
  //                      atof( argv[12] ), //
  //                      -atof( argv[10] ), //
  //                      atof( argv[14] ), //
  //                      -atof( argv[18] ), //
  //                      atof( argv[16] ), //
  //                      atof( argv[11] ), //
  //                      -atof( argv[15] ), //
  //                      atof( argv[19] ), //
  //                      atof( argv[20] ), //
  //                      atof( argv[21] ), //
  //                      atof( argv[22] ), //
  //                      atof( argv[23] ), //
  //                      atof( argv[24] ) };

  //  double rotationmatrix[] = { atof( argv[9] ), //
  //                              -atof( argv[13] ), //
  //                              atof( argv[17] ), //
  //                              0, //
  //                              -atof( argv[10] ), //
  //                              atof( argv[14] ), //
  //                              -atof( argv[18] ), //
  //                              0, //
  //                              atof( argv[11] ), //
  //                              -atof( argv[15] ), //
  //                              atof( argv[19] ), //
  //                              0, //
  //                              atof( argv[21] ), //
  //                              atof( argv[22] ), //
  //                              atof( argv[23] ), //
  //                              atof( argv[24] ) };

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
                                      matrix[12], //
                                      matrix[13], //
                                      matrix[14], //
                                      matrix[15] };

  // 3DAVG refinement transformation
  vtkMatrix4x4 * refinement = vtkMatrix4x4::New();
  refinement->DeepCopy( matrix );

  // 3DAVG refinement transformation (rotation component only)
  vtkMatrix4x4 * refinementRotation = vtkMatrix4x4::New();
  refinementRotation->DeepCopy( rotationmatrix );

  // transformation matrix
  vtkTransform * t = vtkTransform::New();
  t->PostMultiply();

  // make auxiliary 2D rotation transformation
  vtkTransform * r2D = vtkTransform::New();
  r2D->RotateZ( tilt_axis_angle );

  // correction in the direction perpendicular to the tilt axis
  double correction[] = { .5, 0, 0, 1 };
  double tcorrection[4];
  r2D->MultiplyPoint( correction, tcorrection );

  // apply .box coordinate discretization error
  t->Translate( tilt_X_correction - tcorrection[0], tilt_Y_correction - tcorrection[1], 0 );

  // tilt axis angle rotation
  t->RotateZ( tilt_axis_angle );

  // Difference vector between 2D rotation origins (IMOD vs. FREALING)
  double diff2D[] = { .5, .5, 0, 1 };

  // Compute: t = Rot(C1) wrt C2 - C1 = Rot * ( C1 - C2 ) - ( C1 - C2 )
  double tdiff2D[4];
  r2D->MultiplyPoint( diff2D, tdiff2D );
  tdiff2D[0] -= diff2D[0];
  tdiff2D[1] -= diff2D[1];
  tdiff2D[2] -= diff2D[2];
  r2D->Delete();

  // printf("tdiff2D = [ %.5f, %.5f, %.5f ]\n",tdiff2D[0],tdiff2D[1],tdiff2D[2]);

  // apply rotation center and refinement translations
  t->Translate( -tdiff2D[0], -tdiff2D[1], -tdiff2D[2] );

  // correct for center of tilt axis
  t->Translate( -.5, 0, 0 );

  // apply tilt angle rotation
  t->RotateY( -tilt_angle );

  // Convert to IMOD's tilt axis location
  t->Translate( .5, 0, .5 );

  // The remaining transformations are the spike normal and the 3DAVG refinement transformation.
  // Spike normals are a pure rotation R1. 3DAVG refinement is a full rotation and translation matrix F = R2 * T2
  // If the two transformations are composed, the net rotation is then: R = R1 * R2 and the translation component is T2.
  // The origin of the net rotation R is different from the image rotation, so to account for this we express
  // as R * T3, where T3 corrects for the difference in the origin of the rotation.

  // apply spike euler angles (pure rotation R1)
  t->RotateZ( -normal[2] );
  t->RotateX( -normal[0] );
  t->RotateZ( -normal[1] );

  // apply 3DAVG refinement transformation (rotation only, R2)
  vtkMatrix4x4 * local = vtkMatrix4x4::New();
  vtkMatrix4x4::Multiply4x4( refinementRotation, t->GetMatrix(), local );
  t->SetMatrix( local );

  // Apply rotations defined by the CSP euler angles
  t->RotateZ( phi );
  t->RotateY( theta );
  t->RotateZ( psi );

  // compute translation due to change in rotation origin for R1 * R2
  vtkTransform * r = vtkTransform::New();
  r->PostMultiply();
  r->RotateZ( -normal[2] );
  r->RotateX( -normal[0] );
  r->RotateZ( -normal[1] );
  vtkMatrix4x4 * local1 = vtkMatrix4x4::New();
  vtkMatrix4x4::Multiply4x4( r->GetMatrix(), refinementRotation, local1 );
  r->SetMatrix( local1 );

  // Difference vector between rotation origins C1=[51,50,50] and C2=[51,51,50]
  double diff[] = { 0, 1, 0, 1 };

  // Compute: t = Rot(C1) wrt C2 - C1 = Rot * ( C1 - C2 ) - ( C1 - C2 )
  double tdiff[4];
  r->MultiplyPoint( diff, tdiff );
  tdiff[0] -= diff[0];
  tdiff[1] -= diff[1];
  tdiff[2] -= diff[2];

  // printf("tdiff = [ %.5f, %.5f, %.5f ]\n",tdiff[0],tdiff[1],tdiff[2]);

  // compute post-multiplying translation component from refinement matrix, T2
  refinementRotation->Invert();
  vtkMatrix4x4 * f = vtkMatrix4x4::New();
  vtkMatrix4x4::Multiply4x4( refinementRotation, refinement, f );

  // printf("f translations = [ %.5f, %.5f, %.5f ]\n",f->GetElement(0,3),f->GetElement(1,3),f->GetElement(2,3));

  // apply rotation center and refinement translations
  t->Translate( -tdiff[0] + f->GetElement( 0, 3 ), -tdiff[1] - f->GetElement( 1, 3 ), -tdiff[2] + f->GetElement( 2, 3 ) );

  vtkMatrix4x4 * final = t->GetMatrix();

  //  // my calculations
  //  double psi, theta, phi;

  if (fabs( final->GetElement( 2, 2 ) ) < 1 - numeric_limits<float>::min()) {
    theta = acos( final->GetElement( 2, 2 ) );
    double num = final->GetElement( 2, 1 ) / sin( theta );
    double den = -final->GetElement( 2, 0 ) / sin( theta );
    psi = atan2( num, den );

    num = final->GetElement( 1, 2 ) / sin( theta );
    den = final->GetElement( 0, 2 ) / sin( theta );
    phi = atan2( num, den );
  } else {
    cout << "Using alternate" << endl;
    double sign = final->GetElement( 2, 2 ) / fabs( final->GetElement( 2, 2 ) );
    theta = vtkMath::Pi() * ( 1 - sign ) / 2.0;
    phi = 0.0;
    // THIS ORIGINAL CALCULATION WAS NOT CORRECT
    // psi = -  sign * atan2( final->GetElement(0,1) , final->GetElement(0,0) );
    float arg = final->GetElement( 0, 0 ) / cos( theta );
    psi = acos( arg );
  }

  //double frealign[3];

  frealign[0] = psi * vtkMath::RadiansToDegrees();
  frealign[1] = theta * vtkMath::RadiansToDegrees();
  frealign[2] = phi * vtkMath::RadiansToDegrees();

  // frealign does not use negative angles, so we add 360 to each negative angle
  frealign[0] < 0.0 ? frealign[0] = 360.0 + frealign[0] : frealign[0];
  frealign[1] < 0.0 ? frealign[1] = 360.0 + frealign[1] : frealign[1];
  frealign[2] < 0.0 ? frealign[2] = 360.0 + frealign[2] : frealign[2];
  frealign[0] == 360.0 ? frealign[0] = 0.0 : frealign[0];
  frealign[1] == 360.0 ? frealign[1] = 0.0 : frealign[1];
  frealign[2] == 360.0 ? frealign[2] = 0.0 : frealign[2];

  // Now project xyz shifts onto view plane for use in FREALIGN

  vtkTransform * rt = vtkTransform::New();
  rt->RotateZ( phi * vtkMath::RadiansToDegrees() );
  rt->RotateY( theta * vtkMath::RadiansToDegrees() );
  rt->RotateZ( psi * vtkMath::RadiansToDegrees() );

  vtkMatrix4x4 * l = vtkMatrix4x4::New();
  vtkMatrix4x4 * nt = vtkMatrix4x4::New();
  rt->GetInverse( l );
  vtkMatrix4x4::Multiply4x4( l, final, nt );

  //  double sx = 0;//nt->GetElement( 0, 3 );
  //  double sy = 0;//nt->GetElement( 1, 3 );
  double sx = nt->GetElement( 0, 3 );
  double sy = nt->GetElement( 1, 3 );
  //double sz = nt->GetElement( 2, 3 );

  shift[0] = sx;
  shift[1] = sy;

  PrintDebug( DebugData, "[ %10.6f %10.6f %10.6f %10.6f %10.6f ]\n", 5, frealign[0], frealign[1], frealign[2], sx, sy );

}
