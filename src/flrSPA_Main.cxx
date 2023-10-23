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
#include "flrSPAEulerAngles.hxx"
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
    double phi = 270;
    double theta = -120;
    double psi = -60;
    double tilt_angle = 0.0;
    double tilt_axis = 0.0;
    double norm[3] = {0,0,0};
    double matrix[16] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    double frealign[3] = {0,0,0};
    double shifts[2] = {0,0};

    SPAEulerAngles( psi, theta, phi, tilt_angle, tilt_axis, norm, matrix, frealign, shifts  );
    printf("PSI = %f, THETA = %f, PHI = %f \n", frealign[0], frealign[1], frealign[2]);
    printf( "x = %f, y = %f\n", shifts[0], shifts[1] );
}
