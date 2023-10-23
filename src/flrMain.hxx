/*=========================================================================

 Name:        flrMain.hxx
 Category:    Headers files
 Description: This is the main header file, contains most of the
 includes, typedef, constant definitions, etc. of the
 <framework>.
 Language:    C++
 Date:        $Date: 2010-08-21 09:43:09 -0300 (Sat, 21 Aug 2010) $
 Author:      $Author: fefo $
 Revision:    $Rev: 492 $
 Version:     $Id: flrMain.hxx 492 2010-08-21 12:43:09Z fefo $

 ===========================================================================*/
#ifndef __flrmain_hxx
#define __flrmain_hxx

#define VARIABLES_SPACE_DIMENSION 6
//#define _DEBUG 3

#include <list>
#include <vector>
//#include "flrClasses.hxx"

#include "itkImage.h"
// #include "itkBioRadImageIO.h"
#include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"
// #include "itkCastImageFilter.h"
// #include "itkRescaleIntensityImageFilter.h"
// #include "itkImageRegionConstIterator.h"
// #include "itkImageRegionIterator.h"
// #include "itkVnlFFTRealToComplexConjugateImageFilter.h"
// //#include "local_include/itkFFTShiftImageFilter.h"
// #include "itkComplexToModulusImageFilter.h"

// #include "itkSphereSpatialFunction.h"
// #include "itkFloodFilledSpatialFunctionConditionalIterator.h"
// #include "itkFloodFilledSpatialFunctionConditionalConstIterator.h"

#include "itkPowellOptimizer.h"
#include "itkLBFGSOptimizer.h"

using namespace std;

/** Image dimension. */
const unsigned int Dimension = 2;

/** Input and output pixel type. */
typedef float PixelType;
typedef unsigned short WritePixelType;

/** Input and output image types. Filter for casting before writing. */
typedef itk::Image<PixelType, Dimension> ImageType;
// typedef itk::Image<WritePixelType, Dimension> WriteImageType;
// typedef itk::CastImageFilter<ImageType, WriteImageType> CastToWriteImageFilterType;

// /** Density map type, a 3D image. */
// typedef itk::Image<PixelType, 3> DensityMapType;
// typedef itk::ImageFileReader<DensityMapType> DensityMapReaderType;
// typedef itk::ImageFileWriter<DensityMapType> DensityMapWriterType;

// /** Image intensity re-scaler */
// typedef itk::RescaleIntensityImageFilter<WriteImageType, WriteImageType> RescaleFilterType;

// /** Filter for computing the FFT */
// typedef itk::FFTRealToComplexConjugateImageFilter<PixelType, Dimension> FFTFilterType;
// typedef FFTFilterType::OutputImageType FFTImageType;
// //typedef itk::FFTShiftImageFilter<FFTFilterType::OutputImageType, FFTFilterType::OutputImageType>  FFTShiftFilterType;

// /** Iterators for the defined image types. */
// typedef itk::ImageRegionConstIterator<ImageType> imageConstIteratorType;
// typedef itk::ImageRegionIterator<ImageType> imageIteratorType;
// typedef itk::ImageRegionConstIterator<WriteImageType> writeImageConstIteratorType;
// typedef itk::ImageRegionIterator<WriteImageType> writeImageIteratorType;
// typedef itk::ImageRegionConstIterator<FFTImageType> FFTImageConstIteratorType;
// typedef itk::ImageRegionIterator<FFTImageType> FFTImageIteratorType;

// /** Image reader and writer. */
// typedef itk::ImageFileReader<ImageType> ReaderType;
// typedef itk::ImageFileWriter<WriteImageType> WriterType;

// /** Filter to compute the modulus of a complex */
// typedef itk::ComplexToModulusImageFilter<FFTFilterType::OutputImageType, WriteImageType>
//     ModulusFilterType;

// /** Conditional iterators */
// typedef itk::SphereSpatialFunction<Dimension> FunctionType;
// typedef FunctionType::InputType FunctionPositionType;
// typedef itk::FloodFilledSpatialFunctionConditionalIterator<ImageType, FunctionType>
//     imageConditionalIteratorType;
// typedef itk::FloodFilledSpatialFunctionConditionalConstIterator<ImageType, FunctionType>
//     imageConditionalConstIteratorType;
// typedef itk::FloodFilledSpatialFunctionConditionalIterator<WriteImageType, FunctionType>
//     writeImageConditionalIteratorType;
// typedef itk::FloodFilledSpatialFunctionConditionalConstIterator<WriteImageType, FunctionType>
//     WriteImageConditionalConstIteratorType;
// typedef itk::FloodFilledSpatialFunctionConditionalIterator<FFTFilterType::OutputImageType,
//     FunctionType> FFTConditionalIteratorType;
// typedef itk::FloodFilledSpatialFunctionConditionalConstIterator<FFTFilterType::OutputImageType,
//     FunctionType> FFTConditionalConstIteratorType;

/** Array of weights for the Fourier frequency rings. */
typedef vector<float> WeightArrayType;
typedef double EulerAngleType;
typedef double CoordinateType;
typedef map<pair<int,int>, int> IntMapType;

/** Optimizer type definition */
typedef itk::PowellOptimizer PowellOptimizerType;
typedef itk::LBFGSOptimizer LBFGSOptimizerType;

enum LoopModeType {
    OnMicrographRotation, OnParticleRotation, OnParticleTranslation, OnMicrographTranslation, OnMicrographDefocus, OnParticle, OnMicrograph, OnSkip, 
};

enum DebugInfoType {
    DebugFull, DebugData, DebugInfo, DebugBasic
};

#endif 
