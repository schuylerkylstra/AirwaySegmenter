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

#ifndef __itkPhysicalSpaceBinaryErodeImageFilter_hxx
#define __itkPhysicalSpaceBinaryErodeImageFilter_hxx

#include "itkPhysicalSpaceBinaryErodeImageFilter.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryContourImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkFastMarchingImageFilter.h>

namespace itk {

template< typename TInputImage, typename TOutputImage >
PhysicalSpaceBinaryErodeImageFilter< TInputImage, TOutputImage >
::PhysicalSpaceBinaryErodeImageFilter()
{
  m_ErosionDistance = 0.0;
}

template< typename TInputImage, typename TOutputImage >
PhysicalSpaceBinaryErodeImageFilter< TInputImage, TOutputImage >
::~PhysicalSpaceBinaryErodeImageFilter()
{
}

template< typename TInputImage, typename TOutputImage >
void
PhysicalSpaceBinaryErodeImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  os << indent << "ErosionDistance: " << m_ErosionDistance << "\n";
}

template< typename TInputImage, typename TOutputImage >
void
PhysicalSpaceBinaryErodeImageFilter< TInputImage, TOutputImage >
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
      if ( it.Get() <= 0 )
        {
        // Completely outside object, so mark as an alive seed.
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

    // Erode a little past the desired distance
    double maxSpacing = input->GetSpacing()[ 0 ];
    for ( typename InputSizeType::SizeValueType i = 1; i < InputImageDimension; ++i )
      {
      maxSpacing = std::max( maxSpacing, input->GetSpacing()[ i ] );
      }

    fastMarching->SetStoppingValue( m_ErosionDistance + maxSpacing );

    try
      {
      fastMarching->UpdateLargestPossibleRegion();
      }
    catch ( itk::ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not update the fast marching filter" );
      }

    typedef itk::BinaryThresholdImageFilter< FloatImageType, OutputImageType >
      ThresholdFilterType;
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetLowerThreshold( 0.0 );
    thresholdFilter->SetUpperThreshold( m_ErosionDistance );
    thresholdFilter->SetOutsideValue( 1 );
    thresholdFilter->SetInsideValue( 0 );
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
    std::cout << "Using standard binary erode filter\n";

    // Use the standard erosion filter for images with equal spacing
    // Compute radius in terms of pixels based on erosion distance and voxel size
    unsigned int radius = std::ceil( m_ErosionDistance / input->GetSpacing()[0] );
    std::cout << "Radius: " << radius << std::endl;

    typedef itk::BinaryBallStructuringElement< InputPixelType,
                                               InputImageDimension > StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius( radius );
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryErodeImageFilter<InputImageType, OutputImageType, StructuringElementType>
      BinaryErodeImageFilterType;

    typename BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
    erodeFilter->SetForegroundValue( 1 );
    erodeFilter->SetBackgroundValue( 0 );
    erodeFilter->SetBoundaryToForeground( true );
    erodeFilter->SetKernel( structuringElement );
    erodeFilter->SetInput( input );

    try
      {
      erodeFilter->GraftOutput( this->GetOutput() );
      erodeFilter->UpdateLargestPossibleRegion();
      this->GraftOutput( erodeFilter->GetOutput() );
      }
    catch ( itk::ExceptionObject & itkNotUsed(exception) )
      {
      itkExceptionMacro( << "Could not update binary erosion filter" );
      }
    }
}


} // end namespace itk

#endif // __itkPhysicalSpaceBinaryErodeImageFilter_hxx
