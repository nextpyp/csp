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

vtkMatrix4x4* eulerZXZtoZYZ( vtkMatrix4x4* matrixZXZ );
vtkMatrix4x4* eulerTwoZYZtoOneZYZ( vtkMatrix4x4* matrixZYZZYZ);
    

void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] );
    

