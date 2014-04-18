/*=============================================================================
//  --- Airway Segmenter ---+
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Authors: Marc Niethammer, Yi Hong, Johan Andruejol, Cory Quammen
=============================================================================*/

#ifndef __itkPhysicalSpaceBinaryDilateImageFilter_hxx
#define __itkPhysicalSpaceBinaryDilateImageFilter_hxx

#include "itkPhysicalSpaceBinaryDilateImageFilter.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryContourImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkFastMarchingImageFilter.h>

#include <itkImageFileWriter.h>

namespace itk {

template< typename TInputImage, typename TOutputImage >
PhysicalSpaceBinaryDilateImageFilter< TInputImage, TOutputImage >
::PhysicalSpaceBinaryDilateImageFilter()
{
  m_DilationDistance = 0.0;
}

template< typename TInputImage, typename TOutputImage >
PhysicalSpaceBinaryDilateImageFilter< TInputImage, TOutputImage >
::~PhysicalSpaceBinaryDilateImageFilter()
{
}

template< typename TInputImage, typename TOutputImage >
void
PhysicalSpaceBinaryDilateImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  os << indent << "DilationDistance: " << m_DilationDistance << "\n";
}

template< typename TInputImage, typename TOutputImage >
void
PhysicalSpaceBinaryDilateImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  const InputImageType * input = this->GetInput();
  if ( !input )
    {
    return;
    }

  // Check whether we are dealing with isotropic voxels. If not, we must use a
  // fast marching-based approach.
  bool useFastMarching = false;
  for ( typename InputSizeType::SizeValueType i = 0; i < InputImageDimension; ++i )
    {
    if ( std::fabs( input->GetSpacing()[ 0 ] - input->GetSpacing()[ i ] ) > 1e-5 )
      {
      useFastMarching = true;
      break;
      }
    }

  if ( useFastMarching )
    {

    std::cout << "Using fast marching\n";

    // Compute the outer contour of the binary input image. The voxels on this contour
    // serve as trial seed points for the fast marching algorithm.
    typedef BinaryContourImageFilter< InputImageType, InputImageType > ContourFilterType;
    typename ContourFilterType::Pointer contourFilter = ContourFilterType::New();
    contourFilter->SetForegroundValue( 1 );
    contourFilter->SetBackgroundValue( 0 );
    contourFilter->SetInput( input );
    try
      {
      contourFilter->UpdateLargestPossibleRegion();
      }
    catch ( ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not compute binary contour of input\n" );
      return;
      }

    typedef itk::Image< float, InputImageDimension >                       FloatImageType;
    typedef itk::FastMarchingImageFilter< FloatImageType, FloatImageType > FastMarchingFilterType;
    typedef typename FastMarchingFilterType::NodeContainer                 NodeContainer;
    typedef typename FastMarchingFilterType::NodeType                      NodeType;
    typedef itk::ImageRegionConstIterator< InputImageType >                InputIteratorType;

    typename NodeContainer::Pointer aliveSeeds = NodeContainer::New();
    typename NodeContainer::Pointer trialSeeds = NodeContainer::New();
    aliveSeeds->Initialize();
    trialSeeds->Initialize();

    NodeType node;
    node.SetValue( 0.0 ); // Seed value is 0 for all of them, because these are all starting nodes

    InputIteratorType it( input, input->GetLargestPossibleRegion() );
    InputIteratorType contourIt( contourFilter->GetOutput(),
                                 contourFilter->GetOutput()->GetLargestPossibleRegion() );

    for ( it.GoToBegin(), contourIt.GoToBegin(); !it.IsAtEnd(); ++it, ++contourIt )
      {
      node.SetIndex( it.GetIndex() );
      if ( (it.Get() > 0) && (contourIt.Get() <= 0) )
        {
        // Completely inside object, so mark as an alive seed.
        aliveSeeds->InsertElement( aliveSeeds->Size(), node );
        }
      else if ( contourIt.Get() > 0 )
        {
        // On contour, so a trial seed.
        trialSeeds->InsertElement( trialSeeds->Size(), node );
        }
      }

    typename FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
    fastMarching->SetInput( NULL );
    fastMarching->SetAlivePoints( aliveSeeds );
    fastMarching->SetTrialPoints( trialSeeds );
    fastMarching->SetSpeedConstant( 1.0 ); // To solve a simple Eikonal equation
    fastMarching->SetOutputRegion( input->GetBufferedRegion() );
    fastMarching->SetOutputSpacing( input->GetSpacing() );
    fastMarching->SetOutputOrigin( input->GetOrigin() );

    // Dilate a little past the desired distance
    double maxSpacing = input->GetSpacing()[ 0 ];
    for ( typename InputSizeType::SizeValueType i = 1; i < InputImageDimension; ++i )
      {
      maxSpacing = std::max( maxSpacing, input->GetSpacing()[ i ] );
      }

    std::cout << "Dilation distance: " << m_DilationDistance << std::endl;
    fastMarching->SetStoppingValue( m_DilationDistance + maxSpacing );

    try
      {
      fastMarching->UpdateLargestPossibleRegion();
      }
    catch ( itk::ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not update the fast marching filter" );
      }

    typedef itk::ImageFileWriter< FloatImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( fastMarching->GetOutput() );
    writer->SetFileName( "fastMarchingDilate-filter.nrrd" );
    try {
      writer->Update();
    } catch ( itk::ExceptionObject & except ) {
      std::cout << except << std::endl;
    }

    typedef itk::BinaryThresholdImageFilter< FloatImageType, OutputImageType >
      ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetLowerThreshold( -1.0 );
    thresholdFilter->SetUpperThreshold( m_DilationDistance );
    thresholdFilter->SetOutsideValue( 0 );
    thresholdFilter->SetInsideValue( 1 );
    thresholdFilter->SetInput( fastMarching->GetOutput() );

    try
      {
      thresholdFilter->GraftOutput( this->GetOutput() );
      thresholdFilter->UpdateLargestPossibleRegion();
      this->GraftOutput( thresholdFilter->GetOutput() );
      }
    catch ( itk::ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not update threshold filter" );
      }
    }
  else
    {
    std::cout << "Using standard binary dilate filter\n";

    // Use the standard dilation filter for images with equal spacing
    // Compute radius in terms of pixels based on erosion distance and voxel size
    unsigned int radius = std::floor( m_DilationDistance / input->GetSpacing()[0] );
    std::cout << "Radius: " << radius << std::endl;

    typedef itk::BinaryBallStructuringElement< InputPixelType,
                                               InputImageDimension > StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius( radius );
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter<InputImageType, OutputImageType, StructuringElementType>
      BinaryDilateImageFilterType;

    typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetForegroundValue( 1 );
    dilateFilter->SetBackgroundValue( 0 );
    dilateFilter->SetKernel( structuringElement );
    dilateFilter->SetInput( input );

    try
      {
      dilateFilter->GraftOutput( this->GetOutput() );
      dilateFilter->UpdateLargestPossibleRegion();
      this->GraftOutput( dilateFilter->GetOutput() );
      }
    catch ( itk::ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not update binary dilation filter" );
      }
    }
}


} // end namespace itk

#endif // __itkPhysicalSpaceBinaryDilateImageFilter_hxx
