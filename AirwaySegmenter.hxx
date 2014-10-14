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
#ifndef AirwaySegmenter_hxx_included
#define AirwaySegmenter_hxx_included

#include <cfloat>

/* ITK includes */
#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkFastMarchingImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkIdentityTransform.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSmartPointer.h>
#include <itkSpatialOrientationAdapter.h>

/* Local ITK includes. */
#include "itkAirwaySurfaceWriter.h"
#include "itkMaskedOtsuThresholdImageFilter.h"
#include "itkPhysicalSpaceBinaryDilateImageFilter.h"
#include "itkPhysicalSpaceBinaryErodeImageFilter.h"

/* VTK includes */
#include "vtkPolyData.h"
#include "vtkXMLPolyDataWriter.h"

#include "DebugMacros.h"
#include "FastMarchIt.hxx"
#include "LabelIt.hxx"
#include "ProgramArguments.h"

namespace AirwaySegmenter {

  /*******************************************************************/
  /** Compute Otsu threshold of input image.
  /*******************************************************************/
  template< typename TImage >
  typename TImage::Pointer OtsuThreshold( TImage* image, double insideValue, double outsideValue, double & threshold )
  {
    typedef itk::OtsuThresholdImageFilter< TImage, TImage > OtsuThresholdFilterType;
    typename OtsuThresholdFilterType::Pointer otsuThresholdFilter = OtsuThresholdFilterType::New();
    otsuThresholdFilter->SetInsideValue( insideValue );
    otsuThresholdFilter->SetOutsideValue( outsideValue );
    otsuThresholdFilter->SetInput( image );
    TRY_UPDATE( otsuThresholdFilter );

    threshold = otsuThresholdFilter->GetThreshold();

    typename TImage::Pointer output = otsuThresholdFilter->GetOutput();
    output->DisconnectPipeline();

    return output;
  }

  /*******************************************************************/
  /** Binary threshold. */
  /*******************************************************************/
  template< typename TInputImage, typename TOutputImage >
  typename TOutputImage::Pointer BinaryThreshold( TInputImage* input,
                                                  double lowerThreshold, double upperThreshold,
                                                  double outsideValue, double insideValue )
  {
    typedef itk::BinaryThresholdImageFilter< TInputImage, TOutputImage > ThresholdingFilterType;
    typename ThresholdingFilterType::Pointer thresholdFilter = ThresholdingFilterType::New();
    thresholdFilter->SetLowerThreshold( lowerThreshold );
    thresholdFilter->SetUpperThreshold( upperThreshold );
    thresholdFilter->SetOutsideValue( outsideValue );
    thresholdFilter->SetInsideValue( insideValue );
    thresholdFilter->SetInput( input );
    TRY_UPDATE( thresholdFilter );

    typename TOutputImage::Pointer output = thresholdFilter->GetOutput();
    output->DisconnectPipeline();

    return output;
  }


  /*******************************************************************/
  /** Run the algorithm on an input image and write it to the output
      image. */
  /*******************************************************************/
  template< class TInput, class TOutput >
  int Execute( const ProgramArguments & args,
               TInput * originalImage,
               TOutput & output,
               itk::SmartPointer< TInput > & resampledInput,
               typename TInput::PixelType & airwayThreshold )
  {
    /* Typedefs */
    typedef float                      TFloatType;
    typedef typename TInput::PixelType T;
    typedef T                          TPixelType;
    typedef T                          TLabelPixelType;

    const unsigned char DIMENSION = 3;

    typedef itk::Image<TPixelType, DIMENSION>      InputImageType;
    typedef itk::Image<TPixelType, DIMENSION>      OutputImageType;
    typedef itk::Image<TLabelPixelType, DIMENSION> LabelImageType;
    typedef itk::Image<TFloatType, DIMENSION>      FloatImageType;
    typedef itk::Image<unsigned char, DIMENSION>   UCharImageType;

    typedef itk::ImageFileReader<InputImageType>  ReaderType;
    typedef itk::ImageFileReader<LabelImageType>  ReaderLabelType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    typedef itk::ImageFileWriter<LabelImageType>  WriterLabelType;

    typedef typename LabelImageType::SizeType    TSize;
    typedef typename LabelImageType::SpacingType TSpacing;
    typedef typename LabelImageType::PointType   TOrigin;
    typedef typename LabelImageType::IndexType   TIndex;

    std::string sDebugFolder( args.sDebugFolder );

    if ( args.bDebug && (sDebugFolder.compare("None") == 0 ||
         sDebugFolder.compare("") == 0) ) {
      sDebugFolder = "."; //Outputing debug result to current folder if not precised otherwise
    }

    /*  Automatic Resampling to RAI */
    typename InputImageType::DirectionType originalImageDirection = originalImage->GetDirection();

    itk::SpatialOrientationAdapter adapter;
    typename InputImageType::DirectionType RAIDirection = adapter.ToDirectionCosines(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);

    bool shouldConvert = false;

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        if (abs(originalImageDirection[i][j] - RAIDirection[i][j]) > 1e-6) {
          shouldConvert = true;
          break;
        }
      }
    }

    typedef itk::ResampleImageFilter< InputImageType, InputImageType > ResampleImageFilterType;
    typename ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();

    if ( shouldConvert ) {
      typedef itk::IdentityTransform< double, DIMENSION > IdentityTransformType;

      // Figure out bounding box of rotated image
      double boundingBox[6] = { DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX };
      typedef typename InputImageType::IndexType  IndexType;
      typedef typename InputImageType::RegionType RegionType;
      typedef typename InputImageType::PointType  PointType;

      RegionType region = originalImage->GetLargestPossibleRegion();
      IndexType lowerUpper[2];
      lowerUpper[0] = region.GetIndex();
      lowerUpper[1] = region.GetUpperIndex();

      for ( unsigned int i = 0; i < 8; ++i ) {
        IndexType cornerIndex;
        cornerIndex[0] = lowerUpper[ (i & 1u) >> 0 ][0];
        cornerIndex[1] = lowerUpper[ (i & 2u) >> 1 ][1];
        cornerIndex[2] = lowerUpper[ (i & 4u) >> 2 ][2];
        std::cout << "cornerIndex: " << cornerIndex << std::endl;

        PointType point;
        originalImage->TransformIndexToPhysicalPoint( cornerIndex, point );
        boundingBox[0] = std::min( point[0], boundingBox[0] );
        boundingBox[1] = std::max( point[0], boundingBox[1] );
        boundingBox[2] = std::min( point[1], boundingBox[2] );
        boundingBox[3] = std::max( point[1], boundingBox[3] );
        boundingBox[4] = std::min( point[2], boundingBox[4] );
        boundingBox[5] = std::max( point[2], boundingBox[5] );
      }

      // Now transform the bounding box from physical space to index space
      PointType lowerPoint;
      lowerPoint[0] = boundingBox[0];
      lowerPoint[1] = boundingBox[2];
      lowerPoint[2] = boundingBox[4];

      PointType upperPoint;
      upperPoint[0] = boundingBox[1];
      upperPoint[1] = boundingBox[3];
      upperPoint[2] = boundingBox[5];

      typename InputImageType::Pointer dummyImage = InputImageType::New();
      dummyImage->SetOrigin( lowerPoint );
      dummyImage->SetSpacing( originalImage->GetSpacing() );
      dummyImage->SetLargestPossibleRegion( RegionType() );

      IndexType newLower, newUpper;
      dummyImage->TransformPhysicalPointToIndex( lowerPoint, newLower );
      dummyImage->TransformPhysicalPointToIndex( upperPoint, newUpper );

      RegionType outputRegion;
      outputRegion.SetIndex( newLower );
      outputRegion.SetUpperIndex( newUpper );

      // Find the minimum pixel value in the image. This will be used as the default value
      // in the resample filter.
      typedef itk::MinimumMaximumImageCalculator< InputImageType > MinMaxType;
      typename MinMaxType::Pointer minMaxCalculator = MinMaxType::New();
      minMaxCalculator->SetImage( originalImage );
      minMaxCalculator->Compute();

      resampleFilter->SetTransform( IdentityTransformType::New() );
      resampleFilter->SetInput( originalImage );
      resampleFilter->SetSize( outputRegion.GetSize() );
      resampleFilter->SetOutputOrigin( lowerPoint );
      resampleFilter->SetOutputSpacing( originalImage->GetSpacing() );
      resampleFilter->SetDefaultPixelValue( minMaxCalculator->GetMinimum() );
      TRY_UPDATE( resampleFilter );
      DEBUG_WRITE_LABEL_IMAGE( resampleFilter );

      originalImage = resampleFilter->GetOutput();
    }

    resampledInput = originalImage;

    /* Write RAI Image if asked to */
    if ( args.bRAIImage ) {
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( originalImage);
      writer->SetFileName( args.sRAIImagePath.c_str() );
      TRY_UPDATE( writer );
    }

    /* Initial Otsu threshold first to separate patient from background. */
    double initialThreshold = 0.0;
    typename InputImageType::Pointer otsuThreshold = OtsuThreshold( originalImage, 0, 1, initialThreshold );
    DEBUG_WRITE_LABEL_IMAGE_ONLY( otsuThreshold );

    std::cout << "Initial Otsu threshold: " << initialThreshold << std::endl;

    /* Dilation */
    typedef itk::PhysicalSpaceBinaryDilateImageFilter< LabelImageType, LabelImageType >
      DilateFilterType;
    typedef itk::PhysicalSpaceBinaryErodeImageFilter< LabelImageType, LabelImageType >
      ErodeFilterType;

    typename LabelImageType::Pointer initialBinaryImage;

    // Optionally remove thin objects such as breathing masks
    if ( args.bRemoveBreathingMask ) {
      // Initial erosion that removes narrow objects such as breathing masks
      typename ErodeFilterType::Pointer thinErosion = ErodeFilterType::New();
      thinErosion->SetErosionDistance( args.dBreathingMaskThickness );
      thinErosion->SetInput( otsuThreshold );
      TRY_UPDATE( thinErosion );
      DEBUG_WRITE_LABEL_IMAGE( thinErosion );

      typename DilateFilterType::Pointer airwayDilation = DilateFilterType::New();
      airwayDilation->SetDilationDistance( args.dBreathingMaskThickness );
      airwayDilation->SetInput( thinErosion->GetOutput() );
      TRY_UPDATE( airwayDilation );
      DEBUG_WRITE_LABEL_IMAGE( airwayDilation );

      initialBinaryImage = airwayDilation->GetOutput();
    } else {
      initialBinaryImage = otsuThreshold;
    }

    /* Custom Fast marching */
    typedef itk::FastMarchingImageFilter< FloatImageType, FloatImageType > FastMarchingFilterType;
    typedef typename FastMarchingFilterType::NodeContainer                 NodeContainer;
    typedef typename FastMarchingFilterType::NodeType                      NodeType;
    typedef itk::ImageRegionIterator<InputImageType>                       InputIteratorType;
    typedef itk::ImageRegionIterator<LabelImageType>                       IteratorType;
    typedef itk::ImageRegionConstIterator<LabelImageType>                  ConstIteratorType;

    /* Instantiations */
    typename FastMarchingFilterType::Pointer fastMarchingDilate = FastMarchingFilterType::New();
    typename NodeContainer::Pointer trialSeeds = NodeContainer::New();
    typename NodeContainer::Pointer aliveSeeds = NodeContainer::New();
    trialSeeds->Initialize();
    aliveSeeds->Initialize();

    /* Nodes are created as stack variables and initialized with
     * a value and an itk::Index position. NodeType node; */
    NodeType node;
    node.SetValue( 0.0 );

    /* Loop through the output image and set all voxels to 0 seed voxels */
    ConstIteratorType binaryImageIterator( initialBinaryImage,
                                           initialBinaryImage->GetLargestPossibleRegion() );
    ConstIteratorType imageIterator( originalImage,
                                     originalImage->GetLargestPossibleRegion() );

    unsigned int uiNumberOfTrialSeeds = 0;
    unsigned int uiNumberOfAliveSeeds = 0;
    imageIterator.GoToBegin();

    for ( binaryImageIterator.GoToBegin(); !binaryImageIterator.IsAtEnd(); ++binaryImageIterator ) {
      if ( binaryImageIterator.Get() > 0 ) {
        node.SetIndex( binaryImageIterator.GetIndex() );

        if (imageIterator.Get() > 60) {
          aliveSeeds->InsertElement( uiNumberOfAliveSeeds++, node ); // Alive seed
        }
        else {
          trialSeeds->InsertElement( uiNumberOfTrialSeeds++, node ); // Trial seed
        }
      }

      ++imageIterator;
    }

    if ( args.bDebug ) {
      std::cout << std::endl << std::endl;
      std::cout << "FastMarching Dilation: - Number of Alive Seeds: " << aliveSeeds->Size() << std::endl;
      std::cout << "             - Number of Trial Seeds: " << trialSeeds->Size() << std::endl;
    }

    fastMarchingDilate->SetTrialPoints( trialSeeds ); //The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()
    fastMarchingDilate->SetAlivePoints( aliveSeeds );
    fastMarchingDilate->SetInput( NULL );
    fastMarchingDilate->SetSpeedConstant( 1.0 );  // to solve a simple Eikonal equation
     // The FastMarchingImageFilter requires the user to specify the size of the image to be produced as
     // output. This is done using the SetOutputSize().
    fastMarchingDilate->SetOutputSize( otsuThreshold->GetBufferedRegion().GetSize() );
    fastMarchingDilate->SetOutputRegion( otsuThreshold->GetBufferedRegion() );
    fastMarchingDilate->SetOutputSpacing( otsuThreshold->GetSpacing() );
    fastMarchingDilate->SetOutputOrigin( otsuThreshold->GetOrigin() );
    fastMarchingDilate->SetStoppingValue( args.dErodeDistance + args.dMaxAirwayRadius+ 1 );
    TRY_UPDATE( fastMarchingDilate );
    DEBUG_WRITE_IMAGE( fastMarchingDilate );

    typename LabelImageType::Pointer thresholdDilation =
      BinaryThreshold< FloatImageType, LabelImageType >( fastMarchingDilate->GetOutput(),
                                                         0.0, args.dMaxAirwayRadius, 0, 1 );
    DEBUG_WRITE_LABEL_IMAGE_ONLY( thresholdDilation );

    /* Erosion (Thus creating a closing) */

    /* Instantiations */
    typename FastMarchingFilterType::Pointer fastMarchingClose = FastMarchingFilterType::New();
    trialSeeds->Initialize();
    aliveSeeds->Initialize();

    // loop through the output image
    // and set all voxels to 0 seed voxels
    typedef itk::ImageRegionConstIterator< FloatImageType >  ConstFloatIteratorType;
    ConstFloatIteratorType floatDilatedImageIterator( fastMarchingDilate->GetOutput(), fastMarchingDilate->GetOutput()->GetLargestPossibleRegion() );
    ConstIteratorType binaryDilatedImageIterator( thresholdDilation, thresholdDilation->GetLargestPossibleRegion() );
    uiNumberOfTrialSeeds = 0;
    uiNumberOfAliveSeeds = 0;
    floatDilatedImageIterator.GoToBegin();

    for ( binaryDilatedImageIterator.GoToBegin(); !binaryDilatedImageIterator.IsAtEnd(); ++binaryDilatedImageIterator ) {
      if ( binaryDilatedImageIterator.Get() == 0 ) {
        node.SetIndex( binaryDilatedImageIterator.GetIndex() );

        if (floatDilatedImageIterator.Get() > args.dMaxAirwayRadius + args.dErodeDistance) {
          aliveSeeds->InsertElement( uiNumberOfAliveSeeds++, node ); // Alive seed
        }
        else {
          trialSeeds->InsertElement( uiNumberOfTrialSeeds++, node ); // Trial seed
        }
      }

      ++floatDilatedImageIterator;
    }

    if (args.bDebug) {
      std::cout<<std::endl<<std::endl;
      std::cout<<"FastMarching Close:   - Number of Alive Seeds: "<<aliveSeeds->Size()<<" "<<uiNumberOfAliveSeeds << std::endl;
      std::cout<<"                      - Number of Trial Seeds: "<<trialSeeds->Size()<<" "<<uiNumberOfTrialSeeds << std::endl;
    }

    fastMarchingClose->SetTrialPoints( trialSeeds ); // The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()
    fastMarchingClose->SetAlivePoints( aliveSeeds );
    fastMarchingClose->SetInput( NULL );
    fastMarchingClose->SetSpeedConstant( 1.0 );  // To solve a simple Eikonal equation
    fastMarchingClose->SetOutputSize( thresholdDilation->GetBufferedRegion().GetSize() ); // The FastMarchingImageFilter requires the user to specify the size of the image to be produced as output. This is done using the SetOutputSize()
    fastMarchingClose->SetOutputRegion( thresholdDilation->GetBufferedRegion() );
    fastMarchingClose->SetOutputSpacing( thresholdDilation->GetSpacing() );
    fastMarchingClose->SetOutputOrigin( thresholdDilation->GetOrigin() );

    fastMarchingClose->SetStoppingValue( args.dErodeDistance + args.dMaxAirwayRadius+ 1 );
    TRY_UPDATE( fastMarchingClose );
    DEBUG_WRITE_IMAGE( fastMarchingClose );

    // Done with thresholdDilation
    thresholdDilation = NULL;

    typename LabelImageType::Pointer thresholdClosing =
      BinaryThreshold< FloatImageType, LabelImageType >( fastMarchingClose->GetOutput(),
                                                         0.0, args.dMaxAirwayRadius, 1, 0 );
    DEBUG_WRITE_LABEL_IMAGE_ONLY( thresholdClosing );

    typedef itk::MaskedOtsuThresholdImageFilter< InputImageType,
                                                 LabelImageType,
                                                 LabelImageType > MaskedOtsuThresholdFilterType;
    typedef itk::ConnectedComponentImageFilter< LabelImageType,
                                                LabelImageType >  ConnectedComponentType;
    typedef itk::RelabelComponentImageFilter< LabelImageType,
                                              LabelImageType >    RelabelComponentType;
    typedef itk::BinaryThresholdImageFilter< LabelImageType,
                                             LabelImageType >     FinalThresholdingFilterType;

    /* Difference between closed image and ostu-threshold of the original one */
    typedef itk::AbsoluteValueDifferenceImageFilter<LabelImageType, LabelImageType, LabelImageType > TAbsoluteValueDifferenceFilter;
    typename TAbsoluteValueDifferenceFilter::Pointer absoluteValueDifferenceFilter = TAbsoluteValueDifferenceFilter::New();
    absoluteValueDifferenceFilter->SetInput1( otsuThreshold );
    absoluteValueDifferenceFilter->SetInput2( thresholdClosing );
    TRY_UPDATE( absoluteValueDifferenceFilter );
    DEBUG_WRITE_LABEL_IMAGE( absoluteValueDifferenceFilter );

    // Done with thresholdClosing
    thresholdClosing = NULL;

    /* Create a slightly eroded version of the closed image
     * This is to prevent any weird effects at the outside of the face */
    typename LabelImageType::Pointer thresholdDifference =
      BinaryThreshold< FloatImageType, LabelImageType >( fastMarchingClose->GetOutput(), 0.0, args.dMaxAirwayRadius + args.dErodeDistance, 1, 0 );
    DEBUG_WRITE_LABEL_IMAGE_ONLY( thresholdDifference );

    /* The masking */
    typedef itk::MaskImageFilter<LabelImageType, LabelImageType, LabelImageType > TMaskImageFilter;
    typename TMaskImageFilter::Pointer absoluteValueDifferenceFilterMasked = TMaskImageFilter::New();
    absoluteValueDifferenceFilterMasked->SetInput1( absoluteValueDifferenceFilter->GetOutput() );
    absoluteValueDifferenceFilterMasked->SetInput2( thresholdDifference ); // Second input is the mask
    TRY_UPDATE( absoluteValueDifferenceFilterMasked );
    DEBUG_WRITE_LABEL_IMAGE( absoluteValueDifferenceFilterMasked );

    /* Extract largest component of the difference */
    if (args.bDebug) {
      std::cout << "Extracting largest connected component ... ";
    }

    typename ConnectedComponentType::Pointer connected = ConnectedComponentType::New();
    typename RelabelComponentType::Pointer relabel = RelabelComponentType::New();
    typename FinalThresholdingFilterType::Pointer largestComponentThreshold = FinalThresholdingFilterType::New();

    connected->SetInput ( absoluteValueDifferenceFilterMasked->GetOutput());
    TRY_UPDATE( connected );
    DEBUG_WRITE_LABEL_IMAGE( connected );

    // Label the components in the image and relabel them so that object numbers
    // increase as the size of the objects decrease.
    relabel->SetInput( connected->GetOutput() );
    relabel->SetNumberOfObjectsToPrint( 5 );
    TRY_UPDATE( relabel );
    DEBUG_WRITE_LABEL_IMAGE( relabel );

    int componentNumber = 0;

    if (args.iComponent <= 0) {
      componentNumber = LabelIt< T >( relabel->GetOutput(),
                                      args.upperSeed,
                                      args.upperSeedRadius,
                                      args.bDebug );
      if (componentNumber <= 0 ) {
        componentNumber = LabelIt< T >( relabel->GetOutput(),
                                        args.lowerSeed,
                                        args.lowerSeedRadius,
                                        args.bDebug );
      }
      std::cout << "Label found = " << componentNumber << std::endl;
    } else {
      componentNumber = args.iComponent;
    }

    largestComponentThreshold->SetInput( relabel->GetOutput() );
    largestComponentThreshold->SetLowerThreshold( componentNumber ); // object #1
    largestComponentThreshold->SetUpperThreshold( componentNumber ); // object #1
    largestComponentThreshold->SetInsideValue(1);
    largestComponentThreshold->SetOutsideValue(0);
    TRY_UPDATE( largestComponentThreshold );
    DEBUG_WRITE_LABEL_IMAGE( largestComponentThreshold );

    typedef itk::ConnectedThresholdImageFilter< LabelImageType, LabelImageType >
      ConnectedThresholdFilterType;
    typename ConnectedThresholdFilterType::Pointer firstAirwayFragmentFilter =
      ConnectedThresholdFilterType::New();
    firstAirwayFragmentFilter->SetInput( relabel->GetOutput() );
    firstAirwayFragmentFilter->SetLower( 1 );

    // Add in fragmented components of the airway marked with seeds
    LabelImageType * relabelImage = relabel->GetOutput();
    if ( args.bAddAirwayFragments ) {
      for ( size_t i = 0; i < args.airwayFragmentSeeds.size(); ++i ) {
        std::vector< float > fragmentSeed = args.airwayFragmentSeeds[i];

        typename ConnectedThresholdFilterType::InputImageType::PointType point;
        point[0] = -fragmentSeed[0];
        point[1] = -fragmentSeed[1];
        point[2] =  fragmentSeed[2];
        std::cout << "Fragment seed " << i << ": " << point << std::endl;
        typename ConnectedThresholdFilterType::IndexType index;
        relabelImage->TransformPhysicalPointToIndex( point, index );
        firstAirwayFragmentFilter->AddSeed( index );
      }
    }
    TRY_UPDATE( firstAirwayFragmentFilter );
    DEBUG_WRITE_LABEL_IMAGE( firstAirwayFragmentFilter );

    // Now combine the fragments with the largest connected component.
    typedef typename itk::AddImageFilter< LabelImageType > AddLabelImageFilterType;
    typename AddLabelImageFilterType::Pointer firstFragmentCombineFilter =
      AddLabelImageFilterType::New();
    firstFragmentCombineFilter->SetInput1( largestComponentThreshold->GetOutput() );
    firstFragmentCombineFilter->SetInput2( firstAirwayFragmentFilter->GetOutput() );
    TRY_UPDATE( firstFragmentCombineFilter );
    DEBUG_WRITE_LABEL_IMAGE( firstFragmentCombineFilter );

    typename FinalThresholdingFilterType::Pointer firstCombineThresholdFilter =
      FinalThresholdingFilterType::New();
    firstCombineThresholdFilter->SetLowerThreshold( 1 );
    firstCombineThresholdFilter->SetInsideValue( 1 );
    firstCombineThresholdFilter->SetOutsideValue( 0 );
    firstCombineThresholdFilter->SetInput( firstFragmentCombineFilter->GetOutput() );
    TRY_UPDATE( firstCombineThresholdFilter );
    DEBUG_WRITE_LABEL_IMAGE( firstCombineThresholdFilter );

    /* Now do another Otsu thresholding but just around the current segmentation. */
    /* For this, we need to make it first a little bit bigger */
    FloatImageType::Pointer fastMarchFirstCombineThresholdFilter =
      FastMarchIt< T >( firstCombineThresholdFilter->GetOutput(), "Out", args.dErodeDistance, args.dMaxAirwayRadius );

    // To make sure we get roughly twice the volume if the object would have a circular cross section
    double upperThreshold = (sqrt(2.0)-1)*args.dMaxAirwayRadius;
    typename LabelImageType::Pointer extendSegmentation =
      BinaryThreshold< FloatImageType, LabelImageType >( fastMarchFirstCombineThresholdFilter, 0.0, upperThreshold, 0, 1 );
    DEBUG_WRITE_LABEL_IMAGE_ONLY( extendSegmentation );

    /* Now do another Otsu thresholding but restrict the statistics to the currently obtained area  (custom ostu-threshold filter) */
    typename MaskedOtsuThresholdFilterType::Pointer maskedOtsuThresholdFilter = MaskedOtsuThresholdFilterType::New();
    maskedOtsuThresholdFilter->SetInsideValue( 1 );
    maskedOtsuThresholdFilter->SetOutsideValue( 0 );
    maskedOtsuThresholdFilter->SetMaskImage( extendSegmentation.GetPointer() );
    maskedOtsuThresholdFilter->SetInput( originalImage );
    TRY_UPDATE( maskedOtsuThresholdFilter );
    DEBUG_WRITE_LABEL_IMAGE( maskedOtsuThresholdFilter );
    T dThreshold = maskedOtsuThresholdFilter->GetThreshold();

    // Masked Otsu filter above does a whole-image thresholding, so mask it here
    // by the bounds of the patient.
    typename TMaskImageFilter::Pointer maskedOtsu = TMaskImageFilter::New();
    maskedOtsu->SetInput1( maskedOtsuThresholdFilter->GetOutput() );
    maskedOtsu->SetInput2( thresholdDifference ); // Second input is the  mask
    TRY_UPDATE( maskedOtsu );
    DEBUG_WRITE_LABEL_IMAGE( maskedOtsu );

    airwayThreshold = dThreshold;

    /***************************************************************/
    /*
     * Second part of the code : getting rid of lung automatically
     */
    /***************************************************************/
    TOrigin imageOrigin = originalImage->GetOrigin();
    TSpacing imageSpacing = originalImage->GetSpacing();
    TSize imageSize = originalImage->GetBufferedRegion().GetSize();

    /* Need to convert from RAS (Slicer) to ITK's coordinate system: LPS */

    float ballX, ballY, ballZ;

    ballX = -args.lowerSeed[0];
    ballY = -args.lowerSeed[1];
    ballZ = args.lowerSeed[2];

    if (args.bDebug) {
      std::cout << "(x, y, z): " << ballX  << " " << ballY  << " " << ballZ << ", Radius: " << args.lowerSeedRadius << std::endl;
    }

    /* Compute the cubic region around the ball */
    int ballRegion[6];

    ballRegion[0] = int(floor( ( ballX - args.lowerSeedRadius - imageOrigin[0] ) / imageSpacing[0] ));
    ballRegion[1] = int(floor( ( ballY - args.lowerSeedRadius - imageOrigin[1] ) / imageSpacing[1] ));
    ballRegion[2] = int(floor( ( ballZ - args.lowerSeedRadius - imageOrigin[2] ) / imageSpacing[2] ));
    ballRegion[3] = int(ceil( ( ballX + args.lowerSeedRadius - imageOrigin[0] ) / imageSpacing[0] ));
    ballRegion[4] = int(ceil( ( ballY + args.lowerSeedRadius - imageOrigin[1] ) / imageSpacing[1] ));
    ballRegion[5] = int(ceil( ( ballZ + args.lowerSeedRadius - imageOrigin[2] ) / imageSpacing[2] ));

    ballRegion[0] = ballRegion[0] > 0 ? ballRegion[0] : 0;
    ballRegion[1] = ballRegion[1] > 0 ? ballRegion[1] : 0;
    ballRegion[2] = ballRegion[2] > 0 ? ballRegion[2] : 0;

    ballRegion[3] = ballRegion[3] < (imageSize[0]-1) ? ballRegion[3] : (imageSize[0]-1);
    ballRegion[4] = ballRegion[4] < (imageSize[1]-1) ? ballRegion[4] : (imageSize[1]-1);
    ballRegion[5] = ballRegion[5] < (imageSize[2]-1) ? ballRegion[5] : (imageSize[2]-1);

    if (args.bDebug) {
      std::cout << "Origin: "  << imageOrigin[0] << " " << imageOrigin[1] << " "  << imageOrigin[2] << std::endl;
      std::cout << "Spacing: " << imageSpacing[0] << " " << imageSpacing[1] << " " << imageSpacing[2] << std::endl;
      std::cout << "size: "    << imageSize[0] << " " << imageSize[1] << " " << imageSize[2] << std::endl;
      std::cout << "Ball Region: "  << ballRegion[0] << " " << ballRegion[1] << " " << ballRegion[2] << " " << ballRegion[3] << " " << ballRegion[4] << " " << ballRegion[5] << std::endl;
    }

    typename InputImageType::Pointer imageBranch = InputImageType::New();
    typename InputImageType::SizeType sizeBranch;

    sizeBranch[0] = ballRegion[3]-ballRegion[0]+1;
    sizeBranch[1] = ballRegion[4]-ballRegion[1]+1;
    sizeBranch[2] = ballRegion[5]-ballRegion[2]+1;

    typename InputImageType::IndexType startBranch;

    startBranch.Fill(0);

    typename InputImageType::RegionType regionBranch;

    regionBranch.SetSize( sizeBranch );
    regionBranch.SetIndex( startBranch );
    imageBranch->SetRegions( regionBranch );

    try {
      imageBranch->Allocate();
    } catch (itk::ExceptionObject & excep ) {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      std::cerr << "Please verify your parameters, this is often caused by a misplaced trachea carina"<<std::endl;
    }

    for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ ) {
      for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ ) {
        for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ ) {
          double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
          double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
          double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;
          TIndex pixelIndexBranch;

          pixelIndexBranch[0] = iI - ballRegion[0];
          pixelIndexBranch[1] = iJ - ballRegion[1];
          pixelIndexBranch[2] = iK - ballRegion[2];
          imageBranch->SetPixel( pixelIndexBranch, 0 );

          if( iX * iX + iY * iY + iZ * iZ <= args.lowerSeedRadius * args.lowerSeedRadius ) // If within the radius of the ball
          {
            TIndex pixelIndex;
            pixelIndex[0] = iI;
            pixelIndex[1] = iJ;
            pixelIndex[2] = iK;

            if( maskedOtsu->GetOutput()->GetPixel(pixelIndex) )
            {
              maskedOtsu->GetOutput()->SetPixel(pixelIndex, 0);
              imageBranch->SetPixel( pixelIndexBranch, 1 );
            }
          }
        }
      }
    }

    if ( args.bDebug ) {
      try {
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( imageBranch );
        std::string fileName( args.sDebugFolder );
        fileName.append( "/imageBranch-binary.nrrd" );
        writer->SetFileName( fileName.c_str() );
        writer->Update();
      } catch ( itk::ExceptionObject & except ) {
        std::cerr << "Exception occurred when writing imageBranch.nrrd" << std::endl;
        std::cerr << except << std::endl;
      }
    }

    /* Clean up the ball region: First pass */

    /* First with simple threshold */
    typename ConnectedComponentType::Pointer connectedBranch = ConnectedComponentType::New();
    typename RelabelComponentType::Pointer relabelBranch = RelabelComponentType::New();

    connectedBranch->SetInput( imageBranch );
    TRY_UPDATE( connectedBranch );
    DEBUG_WRITE_LABEL_IMAGE( connectedBranch );

    relabelBranch->SetInput( connectedBranch->GetOutput() );
    relabelBranch->SetNumberOfObjectsToPrint( 5 );
    TRY_UPDATE( relabelBranch );
    DEBUG_WRITE_LABEL_IMAGE( relabelBranch );

    /* Get geometry statistics */
    typedef itk::LabelGeometryImageFilter<LabelImageType> LabelGeometryImageFilterType;
    typename LabelGeometryImageFilterType::Pointer labelBranchGeometry = LabelGeometryImageFilterType::New();
    labelBranchGeometry->SetInput( relabelBranch->GetOutput() );
    labelBranchGeometry->CalculateOrientedBoundingBoxOn();

    /* Just in case */
    labelBranchGeometry->CalculateOrientedLabelRegionsOff();
    labelBranchGeometry->CalculatePixelIndicesOff();
    labelBranchGeometry->CalculateOrientedIntensityRegionsOff();
    TRY_UPDATE( labelBranchGeometry );

    int nBranchParts = relabelBranch->GetNumberOfObjects();
    int nBranchId = 1;

    if( nBranchParts > 1 ) {
      if (args.bDebug) {
        std::cout << "Number of parts in branch: " << nBranchParts << std::endl;
      }

      double minDist2Ball;
      int minLabel;
      double dBallIndexX = ( ballRegion[3] - ballRegion[0] ) / 2.0;
      double dBallIndexY = ( ballRegion[4] - ballRegion[1] ) / 2.0;
      double dBallIndexZ = ( ballRegion[5] - ballRegion[2] ) / 2.0;

      for( int nNumParts=1; nNumParts<=nBranchParts; nNumParts++ ) {
        typename LabelGeometryImageFilterType::BoundingBoxType boundingBox = labelBranchGeometry->GetBoundingBox( nNumParts );
        double xTmp = ( boundingBox[0] + boundingBox[1] ) / 2.0 - dBallIndexX;
        double yTmp = ( boundingBox[2] + boundingBox[3] ) / 2.0 - dBallIndexY;
        double zTmp = ( boundingBox[4] + boundingBox[5] ) / 2.0 - dBallIndexZ;
        double distTmp = sqrt( xTmp * xTmp + yTmp * yTmp + zTmp * zTmp );

        if (args.bDebug) {
          std::cout << "( " << xTmp << ", " << yTmp << ", " << zTmp << "), " << distTmp << std::endl;
        }

        if( nNumParts == 1 || minDist2Ball > distTmp) {
          minDist2Ball = distTmp;
          minLabel = nNumParts;
        }

        if (args.bDebug) {
          std::cout << boundingBox << std::endl;
        }
      }

      nBranchId = minLabel;
    }

    /* Get the biggest element (i.e. lung + airway) */

    typename FinalThresholdingFilterType::Pointer branchThreshold = FinalThresholdingFilterType::New();
    branchThreshold->SetInput( relabelBranch->GetOutput() );
    branchThreshold->SetLowerThreshold( nBranchId );
    branchThreshold->SetUpperThreshold( nBranchId );
    branchThreshold->SetInsideValue(1);
    branchThreshold->SetOutsideValue(0);
    TRY_UPDATE( branchThreshold );
    DEBUG_WRITE_LABEL_IMAGE( branchThreshold );

    typename ConnectedComponentType::Pointer connectedFinalWithoutLung = ConnectedComponentType::New();
    typename RelabelComponentType::Pointer relabelFinalWithoutLung = RelabelComponentType::New();
    connectedFinalWithoutLung->SetInput( maskedOtsu->GetOutput() );
    TRY_UPDATE( connectedFinalWithoutLung );
    DEBUG_WRITE_LABEL_IMAGE( connectedFinalWithoutLung );

    relabelFinalWithoutLung->SetInput( connectedFinalWithoutLung->GetOutput() );
    relabelFinalWithoutLung->SetNumberOfObjectsToPrint( 5 );
    TRY_UPDATE( relabelFinalWithoutLung );
    DEBUG_WRITE_LABEL_IMAGE( relabelFinalWithoutLung );

    /* Clean up the ball region: Second pass */
    if (args.bDebug) {
      std::cout << "Get rid of residual lungs in the ball region ... " << std::endl;
    }

    /* First get the original data in the lung+airway regions of the ball */
    for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ ) {
      for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ ) {
        for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ ) {
          double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
          double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
          double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;

          TIndex pixelIndexBranch;
          pixelIndexBranch[0] = iI - ballRegion[0];
          pixelIndexBranch[1] = iJ - ballRegion[1];
          pixelIndexBranch[2] = iK - ballRegion[2];
          imageBranch->SetPixel( pixelIndexBranch, -1024 );

          if( iX * iX + iY * iY + iZ * iZ <= args.lowerSeedRadius * args.lowerSeedRadius ) {
            TIndex pixelIndex;
            pixelIndex[0] = iI;
            pixelIndex[1] = iJ;
            pixelIndex[2] = iK;

            if( branchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) ) {
              imageBranch->SetPixel( pixelIndexBranch, originalImage->GetPixel( pixelIndex ) );
            }
          }
        }
      }
    }

    if ( args.bDebug ) {
      try {
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( imageBranch );
        std::string fileName( args.sDebugFolder );
        fileName.append( "/imageBranch-originalImage.nrrd" );
        writer->SetFileName( fileName.c_str() );
        writer->Update();
      } catch ( itk::ExceptionObject & except ) {
        std::cerr << "Exception occurred when writing imageBranch.nrrd" << std::endl;
        std::cerr << except << std::endl;
      }
    }

    // Now apply an Otsu threshold on it.
    double branchThresholdValue = 0.0;
    typename InputImageType::Pointer otsuThresholdBranch = OtsuThreshold( imageBranch.GetPointer(), 1, 0, branchThresholdValue );

    if (args.bDebug) {
      std::cout << "Getting rid of the small lungs parts ... " << std::endl;
    }

    for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ ) {
      for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ ) {
        for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ ) {
          double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
          double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
          double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;

          TIndex pixelIndexBranch;
          pixelIndexBranch[0] = iI - ballRegion[0];
          pixelIndexBranch[1] = iJ - ballRegion[1];
          pixelIndexBranch[2] = iK - ballRegion[2];
          imageBranch->SetPixel( pixelIndexBranch, 0 );

          if( iX * iX + iY * iY + iZ * iZ <= args.lowerSeedRadius * args.lowerSeedRadius ) {
            if( branchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) &&
                otsuThresholdBranch->GetPixel( pixelIndexBranch ) ) {
              imageBranch->SetPixel( pixelIndexBranch, 1 );
            }
          }
        }
      }
    }

    typename ConnectedComponentType::Pointer connectedCleanedBranch = ConnectedComponentType::New();
    typename RelabelComponentType::Pointer relabelCleanedBranch = RelabelComponentType::New();

    connectedCleanedBranch->SetInput( imageBranch );
    TRY_UPDATE( connectedCleanedBranch )
    DEBUG_WRITE_LABEL_IMAGE( connectedCleanedBranch );

    relabelCleanedBranch->SetInput( connectedCleanedBranch->GetOutput() );
    relabelCleanedBranch->SetNumberOfObjectsToPrint( 5 );
    TRY_UPDATE( relabelCleanedBranch );
    DEBUG_WRITE_LABEL_IMAGE( relabelCleanedBranch );

    typename FinalThresholdingFilterType::Pointer cleanedBranchThreshold = FinalThresholdingFilterType::New();
    cleanedBranchThreshold->SetInput( relabelCleanedBranch->GetOutput() );
    cleanedBranchThreshold->SetLowerThreshold( 1 );
    cleanedBranchThreshold->SetUpperThreshold( 1 );
    cleanedBranchThreshold->SetInsideValue( 1 );
    cleanedBranchThreshold->SetOutsideValue( 0 );
    TRY_UPDATE( cleanedBranchThreshold );
    DEBUG_WRITE_LABEL_IMAGE( cleanedBranchThreshold );

    if (args.bDebug) {
      std::cout << "Final airway label ... " << std::endl;
    }

    /* Find the airway label using the pyriform aperture */
    int nNumAirway = 0;

    if (args.iComponent <= 0) {
      nNumAirway = LabelIt<T>(relabelFinalWithoutLung->GetOutput(), args.upperSeed, args.upperSeedRadius, args.bDebug);
      if (nNumAirway <= 0 ) {
        nNumAirway = LabelIt<T>(relabelFinalWithoutLung->GetOutput(), args.lowerSeed, args.lowerSeedRadius, args.bDebug);
      }
      std::cout << "Label found = " << nNumAirway << std::endl;
    } else {
      nNumAirway = args.iComponent;
    }

    // Check if the maximum label found is 0, meaning that no label was found in the
    // nose region => Nasal cavity probably not segmented!
    if (nNumAirway == 0) {
      std::cerr << "WARNING !" << std::endl;
      std::cerr << "The maximum label found in the spherical region around the pyriform aperture was zero !" << std::endl;
      std::cerr << "This probably means that nasal cavity is not segmented (or the point is misplaced)." << std::endl;
      std::cerr << " Advice: use --debug to ouput and check all the labels found and/or increase the upperSeedRadius to cover more space" << std::endl;

      if (args.bNoWarning) {
        return EXIT_FAILURE;
      }
    }

    if (args.bDebug) {
      std::cout << "The label " << nNumAirway << " is picked as the airway." << std::endl;
    }

    typename FinalThresholdingFilterType::Pointer finalAirwayThreshold = FinalThresholdingFilterType::New();
    finalAirwayThreshold->SetInput( relabelFinalWithoutLung->GetOutput() );
    finalAirwayThreshold->SetLowerThreshold( nNumAirway );
    finalAirwayThreshold->SetUpperThreshold( nNumAirway );
    finalAirwayThreshold->SetInsideValue(1);
    finalAirwayThreshold->SetOutsideValue(0);
    TRY_UPDATE( finalAirwayThreshold );
    DEBUG_WRITE_LABEL_IMAGE( finalAirwayThreshold );

    /* Finally paste the ball back */
    if (args.bDebug) {
      std::cout << "Putting the branches back ... " << std::endl;
    }

    for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ ) {
      for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ ) {
        for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ ) {
          double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
          double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
          double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;

          if( iX * iX + iY * iY + iZ * iZ <= args.lowerSeedRadius * args.lowerSeedRadius ) {
            TIndex pixelIndex;

            pixelIndex[0] = iI;
            pixelIndex[1] = iJ;
            pixelIndex[2] = iK;

            TIndex pixelIndexBranch;

            pixelIndexBranch[0] = iI - ballRegion[0];
            pixelIndexBranch[1] = iJ - ballRegion[1];
            pixelIndexBranch[2] = iK - ballRegion[2];

            if( cleanedBranchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) ) {
              finalAirwayThreshold->GetOutput()->SetPixel(pixelIndex, 1);
            }
          }
        }
      }
    }

    // Add in fragments of the airway
    typename ConnectedThresholdFilterType::Pointer finalFragmentFilter =
      ConnectedThresholdFilterType::New();
    finalFragmentFilter->SetInput( relabelFinalWithoutLung->GetOutput() );
    finalFragmentFilter->SetLower( 1 );

    // Add in fragmented components of the airway marked with seeds
    relabelImage = relabelFinalWithoutLung->GetOutput();
    if ( args.bAddAirwayFragments ) {
      for ( size_t i = 0; i < args.airwayFragmentSeeds.size(); ++i ) {
        std::vector< float > fragmentSeed = args.airwayFragmentSeeds[i];

        typename ConnectedThresholdFilterType::InputImageType::PointType point;
        point[0] = -fragmentSeed[0];
        point[1] = -fragmentSeed[1];
        point[2] =  fragmentSeed[2];
        std::cout << "Fragment seed " << i << ": " << point << std::endl;
        typename ConnectedThresholdFilterType::IndexType index;
        relabelImage->TransformPhysicalPointToIndex( point, index );
        finalFragmentFilter->AddSeed( index );
      }
    }
    TRY_UPDATE( finalFragmentFilter );
    DEBUG_WRITE_LABEL_IMAGE( finalFragmentFilter );

    typename AddLabelImageFilterType::Pointer finalFragmentCombineFilter =
      AddLabelImageFilterType::New();
    finalFragmentCombineFilter->SetInput1( finalAirwayThreshold->GetOutput() );
    finalFragmentCombineFilter->SetInput2( finalFragmentFilter->GetOutput() );
    TRY_UPDATE( finalFragmentCombineFilter );
    DEBUG_WRITE_LABEL_IMAGE( finalFragmentCombineFilter );
    DEBUG_WRITE_LABEL_IMAGE( finalFragmentCombineFilter );

    typename FinalThresholdingFilterType::Pointer finalCombineThresholdFilter =
      FinalThresholdingFilterType::New();
    finalCombineThresholdFilter->SetLowerThreshold( 1 );
    finalCombineThresholdFilter->SetInsideValue( 1 );
    finalCombineThresholdFilter->SetOutsideValue( 0 );
    finalCombineThresholdFilter->SetInput( finalFragmentCombineFilter->GetOutput() );
    TRY_UPDATE( finalCombineThresholdFilter );
    DEBUG_WRITE_LABEL_IMAGE( finalCombineThresholdFilter );

    typename LabelImageType::Pointer finalSegmentation = finalCombineThresholdFilter->GetOutput();

    // Optionally remove the maxillary sinus(es)
    // NOTE!!! This has not been updated to support discontinuous airways
    // (those with fragments).
    if (args.bRemoveMaxillarySinuses) {
      std::cout << "maxillarySinusesSeeds "<<args.maxillarySinusesSeeds.size() << std::endl;

      /* First thing, Erode to severe the small connection
       * between the sinuses and the airway
       * Note that the erosion is very small
       * Using custom fast marching function
       */
      typename LabelImageType::Pointer thresholdSlightErosion =
        BinaryThreshold< FloatImageType, LabelImageType >( FastMarchIt< T >( finalSegmentation, "In", args.dErodeDistance, args.dMaxAirwayRadius),
                         args.dMaxAirwayRadius*args.erosionPercentage, args.dMaxAirwayRadius, 0, 1 );
      DEBUG_WRITE_LABEL_IMAGE_ONLY( thresholdSlightErosion );

      /* Create the image that we will substract from the segmentation
       * It should be all the maxillary sinuses
       */

      typename ConnectedComponentType::Pointer connectedSinuses = ConnectedComponentType::New();
      typename RelabelComponentType::Pointer relabelSinuses = RelabelComponentType::New();

      connectedSinuses->SetInput( thresholdSlightErosion );
      relabelSinuses->SetInput( connectedSinuses->GetOutput() );
      relabelSinuses->SetNumberOfObjectsToPrint( 5 );
      try {
        relabelSinuses->Update(); // First label the different part
                                  // that have been separated by the
                                  // erosion
      } catch ( itk::ExceptionObject & excep ) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
      }

      // Done with thresholdSlightErosion
      thresholdSlightErosion = NULL;

      typedef typename itk::AddImageFilter< LabelImageType > AddLabelImageFilterType;
      typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType > LabelThresholdFilterType;
      typename AddLabelImageFilterType::Pointer addFilter = AddLabelImageFilterType::New();

      typedef itk::ImageDuplicator< LabelImageType > DuplicatorType;
      typename DuplicatorType::Pointer duplicator = DuplicatorType::New(); // The blank image is declared with a duplicator>. Should probably be done otherwise

      duplicator->SetInputImage(relabelSinuses->GetOutput());
      TRY_UPDATE( duplicator );
      DEBUG_WRITE_LABEL_IMAGE( duplicator );

      typename LabelImageType::Pointer sinusesImage = duplicator->GetOutput();

      sinusesImage->FillBuffer(0);

      int airwayLabel = LabelIt<T>(relabelSinuses->GetOutput(), args.upperSeed, args.upperSeedRadius, args.bDebug); //Get the airway label

      for (int i = 0; i < args.maxillarySinusesSeeds.size(); ++i) { //For each seed
        int seedLabel = LabelIt<T>(relabelSinuses->GetOutput(), args.maxillarySinusesSeeds[i], args.maxillarySinusesSeedsRadius, args.bDebug); //Get the label of the maxillary sinus
        //std::cout<<"Seed Label: "<<seedLabel<<std::endl;

        if (airwayLabel == seedLabel) { //The airway label MUST be different than the seed label, Otherwise the airway would be masked out
          std::cerr <<"WARNING !"<<std::endl;
          std::cerr << "The airway label found is equal to the label found with seed #" << i << " (seed =" << args.maxillarySinusesSeeds[i][0] << ", " << args.maxillarySinusesSeeds[i][1]  << ", " << args.maxillarySinusesSeeds[i][2] << ")" << std::endl;
          std::cerr << "Review the seed position and/or the percentage used." << std::endl;

          if (args.bNoWarning) return EXIT_FAILURE; else continue;
        }

        typename LabelThresholdFilterType::Pointer thresholdOut =
          LabelThresholdFilterType::New();
        TRY_UPDATE( thresholdOut );
        DEBUG_WRITE_LABEL_IMAGE( thresholdOut );

        addFilter->SetInput1( sinusesImage );
        addFilter->SetInput2( thresholdOut->GetOutput() );
        TRY_UPDATE( addFilter );
        DEBUG_WRITE_LABEL_IMAGE( addFilter );

        sinusesImage = addFilter->GetOutput();
      }

      /* Dilate the undesired part
       * (So they approximately are their original size)
       * Note that the erosion is very small
       * Using custom fast marching function
       */
      typename LabelImageType::Pointer thresholdSlightDilation =
        BinaryThreshold< FloatImageType, LabelImageType >( FastMarchIt<T>( addFilter->GetOutput(), "Out", args.dErodeDistance, args.dMaxAirwayRadius),
                                                           0, args.dMaxAirwayRadius * args.erosionPercentage, 1, 0 );
      DEBUG_WRITE_LABEL_IMAGE_ONLY( thresholdSlightDilation );

      /* Mask all the undesired part from the segmentation */
      typename TMaskImageFilter::Pointer substractSinusesMask = TMaskImageFilter::New();
      substractSinusesMask->SetMaskImage( thresholdSlightDilation );
      substractSinusesMask->SetInput( finalSegmentation );
      substractSinusesMask->Update();
      TRY_UPDATE( substractSinusesMask );
      DEBUG_WRITE_LABEL_IMAGE( substractSinusesMask );

      // Done with thresholdSlightDilation
      thresholdSlightDilation = NULL;

      /* Clean up */

      typename ConnectedComponentType::Pointer connectedCleanUp = ConnectedComponentType::New();
      typename RelabelComponentType::Pointer relabelCleanUp = RelabelComponentType::New();

      connectedCleanUp->SetInput( substractSinusesMask->GetOutput() );
      relabelCleanUp->SetInput( connectedCleanUp->GetOutput() );
      relabelCleanUp->SetNumberOfObjectsToPrint( 5 );
      TRY_UPDATE( relabelCleanUp );
      DEBUG_WRITE_LABEL_IMAGE( relabelCleanUp );

      nNumAirway = LabelIt<T>(relabelCleanUp->GetOutput(), args.upperSeed, args.upperSeedRadius, args.bDebug);

      typename FinalThresholdingFilterType::Pointer thresholdCleanUp = FinalThresholdingFilterType::New();

      thresholdCleanUp->SetLowerThreshold( nNumAirway );
      thresholdCleanUp->SetUpperThreshold( nNumAirway );
      thresholdCleanUp->SetOutsideValue( 0 );
      thresholdCleanUp->SetInsideValue( 1 );
      thresholdCleanUp->SetInput( relabelCleanUp->GetOutput() );
      TRY_UPDATE( thresholdCleanUp );
      DEBUG_WRITE_LABEL_IMAGE( thresholdCleanUp );

      finalSegmentation = thresholdCleanUp->GetOutput();
    }

    if (args.bDebug) {
      std::cout << "Writing the final image ... " << std::endl;
    }

    // Optionally remove tracheal tube
    if ( args.bRemoveTrachealTube && args.trachealTubeSeed.size() >= 3 ) {
      // Do another Otsu threshold on the patient masked region only.
      // This should separate high-intensity objects such as bones and
      // tracheal tubes from the rest of the body.
      typename MaskedOtsuThresholdFilterType::Pointer trachealTubeFilter =
        MaskedOtsuThresholdFilterType::New();
      trachealTubeFilter->SetInsideValue( 0 );
      trachealTubeFilter->SetOutsideValue( 1 );
      trachealTubeFilter->SetMaskImage( otsuThreshold.GetPointer() );
      trachealTubeFilter->SetInput( originalImage );
      TRY_UPDATE( trachealTubeFilter );
      DEBUG_WRITE_LABEL_IMAGE( trachealTubeFilter );

      typename ConnectedComponentType::Pointer trachealConnected = ConnectedComponentType::New();
      trachealConnected->SetInput( trachealTubeFilter->GetOutput() );
      TRY_UPDATE( trachealConnected );
      DEBUG_WRITE_LABEL_IMAGE( trachealConnected );

      typename RelabelComponentType::Pointer trachealRelabel = RelabelComponentType::New();
      trachealRelabel->SetNumberOfObjectsToPrint( 5 );
      trachealRelabel->SetInput( trachealConnected->GetOutput() );
      TRY_UPDATE( trachealRelabel );
      DEBUG_WRITE_LABEL_IMAGE( trachealRelabel );

      int componentNumber = LabelIt< T >( trachealRelabel->GetOutput(),
                                          args.trachealTubeSeed,
                                          args.trachealTubeSeedRadius,
                                          args.bDebug );

      typename FinalThresholdingFilterType::Pointer trachealLargestComponentThreshold = FinalThresholdingFilterType::New();
      trachealLargestComponentThreshold->SetInput( trachealRelabel->GetOutput() );
      trachealLargestComponentThreshold->SetLowerThreshold( componentNumber );
      trachealLargestComponentThreshold->SetUpperThreshold( componentNumber );
      trachealLargestComponentThreshold->SetInsideValue( 1 );
      trachealLargestComponentThreshold->SetOutsideValue( 0 );
      TRY_UPDATE( trachealLargestComponentThreshold );
      DEBUG_WRITE_LABEL_IMAGE( trachealLargestComponentThreshold );

      // Morphological closing to fill in the center of the tracheal tube.
      double closingDistance = 5.0; // mm
      typename DilateFilterType::Pointer trachealDilate = DilateFilterType::New();
      trachealDilate->SetDilationDistance( closingDistance );
      trachealDilate->SetInput( trachealLargestComponentThreshold->GetOutput() );
      TRY_UPDATE( trachealDilate );
      DEBUG_WRITE_LABEL_IMAGE( trachealDilate );

      typename ErodeFilterType::Pointer trachealErode = ErodeFilterType::New();
      trachealErode->SetErosionDistance( closingDistance - args.trachealTubeDilationDistance );
      trachealErode->SetInput( trachealDilate->GetOutput() );
      TRY_UPDATE( trachealErode );
      DEBUG_WRITE_LABEL_IMAGE( trachealErode );

      // Iterate over the closed tracheal tube mask and make it part of the airway.
      ConstIteratorType trachealIterator( trachealErode->GetOutput(),
                                          trachealErode->GetOutput()->GetLargestPossibleRegion() );
      IteratorType finalIterator( finalSegmentation,
                                  finalSegmentation->GetLargestPossibleRegion() );
      InputIteratorType inputIterator( originalImage, originalImage->GetLargestPossibleRegion() );
      for ( trachealIterator.GoToBegin(),
            finalIterator.GoToBegin(),
            inputIterator.GoToBegin();
            !trachealIterator.IsAtEnd();
            ++trachealIterator,
            ++finalIterator,
            ++inputIterator ) {
        if ( trachealIterator.Get() > 0 ) {
          finalIterator.Set( 1 );
          inputIterator.Set( -1024 );
        }
      }
    }

    output = finalSegmentation;

    DEBUG_WRITE_LABEL_IMAGE_ONLY( finalSegmentation );

    if (args.bDebug) {
      std::cout << "done." << std::endl;
    }

    // Write the threshold used in the Otsu-thresholding
    std::cout << "Threshold computed: " << dThreshold << std::endl;
    if ( args.returnParameterFile.size() > 0 ) {
      std::ofstream parametersFile( args.returnParameterFile.c_str() );
      if ( parametersFile.is_open() ) {
        parametersFile << "computedThreshold = " << dThreshold << std::endl;
      }
    }

    return EXIT_SUCCESS;
  }

  /*******************************************************************/
  /** Execute the algorithm on an image read from a file. */
  /*******************************************************************/
  template <class T>
  int ExecuteFromFile( const ProgramArguments & args, T)
  {
    /* Typedefs */
    typedef float TFloatType;
    typedef T TPixelType;
    typedef T TLabelPixelType;

    const unsigned char DIMENSION = 3;

    typedef itk::Image<TPixelType, DIMENSION> InputImageType;
    typedef itk::Image<TPixelType, DIMENSION> OutputImageType;
    typedef itk::Image<TLabelPixelType, DIMENSION> LabelImageType;
    typedef itk::Image<TFloatType, DIMENSION> FloatImageType;

    typedef itk::Image<unsigned char, DIMENSION> UCharImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typedef itk::ImageFileReader<LabelImageType> ReaderLabelType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    typedef itk::ImageFileWriter<LabelImageType> WriterLabelType;

    // Read the input file
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( args.inputImage ); // rRead the input image
    TRY_UPDATE( reader );

    // Run the algorithm
    typename OutputImageType::Pointer algorithmOutput;
    typename InputImageType::Pointer resampledInput;
    typename InputImageType::PixelType airwayThreshold;
    int result = Execute( args, reader->GetOutput(),
                          algorithmOutput, resampledInput,
                          airwayThreshold );
    if ( result != EXIT_SUCCESS ) {
      return result;
    }

    // Write the result.
    typename WriterLabelType::Pointer writer = WriterLabelType::New();
    writer->SetInput( algorithmOutput );
    writer->SetFileName( args.outputImage.c_str() );
    TRY_UPDATE( writer );

    /* Write Surface */
    if ( args.outputGeometry != "" && args.createGeometry ) {
      typedef itk::AirwaySurfaceWriter<InputImageType, LabelImageType>
        SurfaceWriterType;
      typename SurfaceWriterType::Pointer surfaceWriter =
        SurfaceWriterType::New();

      surfaceWriter->SetFileName( args.outputGeometry.c_str() );
      surfaceWriter->SetUseFastMarching( true );
      surfaceWriter->SetMaskImage( algorithmOutput );
      surfaceWriter->SetInput( resampledInput );
      surfaceWriter->SetThreshold( airwayThreshold );
      TRY_UPDATE( surfaceWriter );
    } else {
      // Create empty polydata and save it.
      vtkPolyData * empty = vtkPolyData::New();
      vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
      writer->SetFileName( args.outputGeometry.c_str() );
#if VTK_MAJOR_VERSION <= 5
      writer->SetInput( empty );
#else
      writer->SetInputData( empty );
#endif
      writer->Update();

      empty->Delete();
      writer->Delete();
    }

    return EXIT_SUCCESS;
  }

} // end namespace AirwaySegmenter

#endif // AirwaySegmenter_hxx_included
