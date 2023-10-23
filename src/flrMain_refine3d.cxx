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

/**
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>
*/
#include "refine3d.h"

#include <dlfcn.h>

//using namespace blitz;

//#include <random/uniform.h>

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

#include <vtkTransform.h>

extern SParameters params;
SParameters params = { 16, 7.6 };
Refine3dParameters refine3dParams;




int main( int argc, char * argv[] ) {

    RunNormalRefine3d();
    return EXIT_SUCCESS;
}
