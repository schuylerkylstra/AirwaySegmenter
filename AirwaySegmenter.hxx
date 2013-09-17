#ifndef AirwaySegmenter_hxx_included
#define AirwaySegmenter_hxx_included

/* ITK includes */

#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkFastMarchingImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkIdentityTransform.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMatrix.h>
#include <itkMultiplyImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSmartPointer.h>
#include <itkSpatialOrientationAdapter.h>
#include <itkSubtractImageFilter.h>

#include "itkAirwaySurfaceWriter.h"
#include "itkMaskedOtsuThresholdImageFilter.h"


namespace AirwaySegmenter {

	/* This Fast marching is factorised for preventing the memory to 
	 * be overwhelmed by unused nodes. Now that 
	 * this function is for internal use, it assumes all 
	 * its parameters are given correctly
	 */
	
	template <class T> itk::Image<float, 3>::Pointer FastMarchIt(typename itk::Image<T, 3>::Pointer image,std::string type, double erodedDistance, double airwayRadius)
	{
		/* Necessaries typedefs */
		
		typedef itk::Image<float, 3> FloatImageType;
		typedef itk::Image<T, 3> LabelImageType;
		typedef itk::FastMarchingImageFilter<FloatImageType, FloatImageType> FastMarchingFilterType;
		typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
		typedef typename FastMarchingFilterType::NodeType NodeType;
		typedef itk::ImageRegionConstIterator<LabelImageType> ConstIteratorType;
	
		/* Instantiations */
		
		typename FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
		typename NodeContainer::Pointer seeds = NodeContainer::New();
		
		seeds->Initialize();
		NodeType node; //Nodes are created as stack variables and  initialized with a value and an itk::Index position. NodeType node;
		node.SetValue( 0.0 ); // Seed value is 0 for all of them, because these are all starting nodes
		ConstIteratorType it( image, image->GetLargestPossibleRegion() ); // Loop through the output image and set all voxels to 0 seed voxels
	
		unsigned int uiNumberOfSeeds = 0;
		
		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			if (type.compare("Out") == 0) // Dilatation
			{
				if ( it.Get() > 0 )
				{
					node.SetIndex( it.GetIndex() );
					seeds->InsertElement( uiNumberOfSeeds++, node );
				}
			}
			else if (type.compare("In") == 0)// Erosion
			{
				if ( it.Get() == 0 )
				{
					node.SetIndex( it.GetIndex() );
					seeds->InsertElement( uiNumberOfSeeds++, node );
				}
			}
		}
	
		fastMarching->SetTrialPoints( seeds ); // The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()
	
		/*
		 * The FastMarchingImageFilter requires the user to specify
		 * the size of the image to be produced as output
		 * This is done using the SetOutputSize()
		 */
	
		fastMarching->SetInput( NULL );
		fastMarching->SetSpeedConstant( 1.0 );  // To solve a simple Eikonal equation
	
		fastMarching->SetOutputSize( image->GetBufferedRegion().GetSize() );
		fastMarching->SetOutputRegion( image->GetBufferedRegion() );
		fastMarching->SetOutputSpacing( image->GetSpacing() );
		fastMarching->SetOutputOrigin( image->GetOrigin() );
	
		fastMarching->SetStoppingValue( airwayRadius + erodedDistance + 1 );
	
		try
		{
			fastMarching->Update();
		}
		catch(itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		return fastMarching->GetOutput();
	}
	
	/* Look for the labels in the sphercal region within radius of
	 * the ball center and return the most represented label
	 */
	
	template <class T> int LabelIt(typename itk::Image<T, 3>::Pointer image, std::vector<float> ballCenter, double radius, bool printLabels)
	{
		typedef itk::Image<T, 3> TImage;
		typedef typename itk::Image<T, 3>::SizeType TSize;
		typedef typename itk::Image<T, 3>::SpacingType TSpacing;
		typedef typename itk::Image<T, 3>::PointType TOrigin;
		typedef typename itk::Image<T, 3>::IndexType TIndex;
	
		TOrigin imageOrigin = image->GetOrigin();
		TSpacing imageSpacing = image->GetSpacing();
		TSize imageSize = image->GetBufferedRegion().GetSize();
	
		/* Convert to LPS system */
		
		float x, y, z;
		
		x = -ballCenter[0];
		y = -ballCenter[1];
		z =  ballCenter[2];
	
		/* Get the bounding box around the piryna aperture */
		
		int region[6];
		
		region[0] = int( floor( ( x - radius - imageOrigin[0] ) / imageSpacing[0] ) );
		region[1] = int( floor( ( y - radius - imageOrigin[1] ) / imageSpacing[1] ) );
		region[2] = int( floor( ( z - radius - imageOrigin[2] ) / imageSpacing[2] ) );
		region[3] = int(  ceil( ( x + radius - imageOrigin[0] ) / imageSpacing[0] ) );
		region[4] = int(  ceil( ( y + radius - imageOrigin[1] ) / imageSpacing[1] ) );
		region[5] = int(  ceil( ( z + radius - imageOrigin[2] ) / imageSpacing[2] ) );
	
		region[0] = region[0] > 0 ? region[0] : 0;
		region[1] = region[1] > 0 ? region[1] : 0;
		region[2] = region[2] > 0 ? region[2] : 0;
	
		region[3] = region[3] < (imageSize[0]-1) ? region[3] : (imageSize[0]-1);
		region[4] = region[4] < (imageSize[1]-1) ? region[4] : (imageSize[1]-1);
		region[5] = region[5] < (imageSize[2]-1) ? region[5] : (imageSize[2]-1);
	
		if (printLabels) std::cout << "Region: " << region[0] << ", " << region[1] << ", " << region[2] << ", " << region[3] << ", " << region[4] << ", "  << region[5] << std::endl;
	
		std::map<int, int> labels;
	
		for( int iI=region[0]; iI<=region[3]; iI++ )
		{
			for( int iJ=region[1]; iJ<=region[4]; iJ++ )
			{
				for( int iK=region[2]; iK<=region[5]; iK++ )
				{
					/* Get real space position */
					
					double iX = iI * imageSpacing[0] + imageOrigin[0] - x;
					double iY = iJ * imageSpacing[1] + imageOrigin[1] - y;
					double iZ = iK * imageSpacing[2] + imageOrigin[2] - z;
	
					/* If within the nose ball region */
					
					if( iX * iX + iY * iY + iZ * iZ <= radius * radius )
					{
						TIndex pixelIndex;
						
						pixelIndex[0] = iI;
						pixelIndex[1] = iJ;
						pixelIndex[2] = iK;
	
						int label = image->GetPixel( pixelIndex );
	
						if( label != 0 ) labels[label] += 1; // Ignore label 0, it's always the background
					}
				}
			}
		}
	
		int finalLabel = 0;
		int labelCount = 0;
		
		for(std::map<int, int>::const_iterator it = labels.begin(); it != labels.end(); ++it)
		{
			if (printLabels) std::cout<<"Labels "<<it->first<<" :   "<<it->second<<std::endl;
	
			if (it->second > labelCount)
			{
				labelCount = it->second;
				finalLabel = it->first;
			}
		}
	
		return finalLabel;
	}
	
	template <class T> int DoIt(int argc, char* argv[], T)
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
	
		typedef typename LabelImageType::SizeType TSize;
		typedef typename LabelImageType::SpacingType TSpacing;
		typedef typename LabelImageType::PointType TOrigin;
		typedef typename LabelImageType::IndexType TIndex;
	
		/* Parsing the input */
		
		PARSE_ARGS; // Parse the input arguments
	
		if(bDebug && (sDebugFolder.compare("None") == 0 || sDebugFolder.compare("") == 0)) sDebugFolder = "."; //Outputing debug result to current folder if not precised otherwise
	
		typename ReaderType::Pointer reader = ReaderType::New();
		
		reader->SetFileName( inputImage ); // rRead the input image
	
		try
		{
			reader->Update();  
		}
		catch ( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		/*  Automatic Resampling to RAI */
		
		typename InputImageType::Pointer originalImage = reader->GetOutput();
		typename InputImageType::DirectionType originalImageDirection = originalImage->GetDirection();
	
		itk::SpatialOrientationAdapter adapter;
		typename InputImageType::DirectionType RAIDirection = adapter.ToDirectionCosines(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
	
		bool shouldConvert = false;
		
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (abs(originalImageDirection[i][j] - RAIDirection[i][j]) > 1e-6)
				{
					shouldConvert = true;
					break;
				}
			}
		}
	
		typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleImageFilterType;
		typename ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
	
		if (shouldConvert)
		{
			typedef itk::IdentityTransform<double, DIMENSION> IdentityTransformType;
	
			resampleFilter->SetTransform(IdentityTransformType::New());
			resampleFilter->SetInput(originalImage);
			resampleFilter->SetSize(originalImage->GetLargestPossibleRegion().GetSize());
			resampleFilter->SetOutputOrigin(originalImage->GetOrigin());
			resampleFilter->SetOutputSpacing(originalImage->GetSpacing());
			resampleFilter->Update();
	
			originalImage = resampleFilter->GetOutput();
		}
	
		/* Write RAI Image if asked to */
		
		if (bRAIImage)
		{
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput( originalImage);
			writer->SetFileName( sRAIImagePath.c_str() );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
		
		/* Otsu thresholding first */
	
		typedef itk::OtsuThresholdImageFilter<InputImageType, LabelImageType > OtsuThresholdFilterType;
		typename OtsuThresholdFilterType::Pointer otsuThresholdFilter = OtsuThresholdFilterType::New();
	
		otsuThresholdFilter->SetInsideValue( 0 );
		otsuThresholdFilter->SetOutsideValue( 1 );
		otsuThresholdFilter->SetInput( originalImage );
		otsuThresholdFilter->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( otsuThresholdFilter->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/otsu.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl; std::cerr << excep << std::endl; 
			}
		}
	
		/* Dilatation */
	
		typedef itk::BinaryThresholdImageFilter<FloatImageType, LabelImageType > ThresholdingFilterType;
		typename ThresholdingFilterType::Pointer thresholdDilation = ThresholdingFilterType::New();
	
		thresholdDilation->SetLowerThreshold( 0.0 );
		thresholdDilation->SetUpperThreshold( dMaxAirwayRadius );
		thresholdDilation->SetOutsideValue( 0 );
		thresholdDilation->SetInsideValue( 1 );
	
		/* Custom Fast marching */
		
		typedef itk::FastMarchingImageFilter<FloatImageType, FloatImageType>  FastMarchingFilterType;
		typedef typename FastMarchingFilterType::NodeContainer  NodeContainer;
		typedef typename FastMarchingFilterType::NodeType NodeType;
		typedef itk::ImageRegionConstIterator<LabelImageType>  ConstIteratorType;
	
		/* Instantiations */
		
		typename FastMarchingFilterType::Pointer fastMarchingDilate = FastMarchingFilterType::New();
		typename NodeContainer::Pointer trialSeeds = NodeContainer::New();
		typename NodeContainer::Pointer aliveSeeds = NodeContainer::New();
		trialSeeds->Initialize();
		aliveSeeds->Initialize();
	
		/* Nodes are created as stack variables and initialized with
		 * a value and an itk::Index position. NodeType node;
		 */ 
		
		NodeType node;
		node.SetValue( 0.0 );
	
		/* Loop through the output image and set all voxels to 0 seed voxels */
		
		ConstIteratorType binaryImageIterator( otsuThresholdFilter->GetOutput(),otsuThresholdFilter->GetOutput()->GetLargestPossibleRegion() );
		ConstIteratorType imageIterator( originalImage,originalImage->GetLargestPossibleRegion() );
	
		unsigned int uiNumberOfTrialSeeds = 0;
		unsigned int uiNumberOfAliveSeeds = 0;
		imageIterator.GoToBegin();
		
		for ( binaryImageIterator.GoToBegin(); !binaryImageIterator.IsAtEnd(); ++binaryImageIterator )
		{
			if ( binaryImageIterator.Get() > 0 )
			{
				node.SetIndex( binaryImageIterator.GetIndex() );
	
				if (imageIterator.Get() > 60) aliveSeeds->InsertElement( uiNumberOfAliveSeeds++, node ); // Alive seed 
				else trialSeeds->InsertElement( uiNumberOfTrialSeeds++, node ); // Trial seed
			}
			
			++imageIterator;
		}
	
		if (bDebug)
		{
			std::cout<<std::endl<<std::endl;
			std::cout<< "FastMarching Dilatation:	- Number of Alive Seeds: "<< aliveSeeds->Size() <<std::endl;
			std::cout<< "							- Number of Trial Seeds: "<< trialSeeds->Size() <<std::endl;
		}
	
		fastMarchingDilate->SetTrialPoints( trialSeeds ); //The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()
		fastMarchingDilate->SetAlivePoints( aliveSeeds );
		fastMarchingDilate->SetInput( NULL );
		fastMarchingDilate->SetSpeedConstant( 1.0 );  // to solve a simple Eikonal equation
		fastMarchingDilate->SetOutputSize( otsuThresholdFilter->GetOutput()->GetBufferedRegion().GetSize() ); // The FastMarchingImageFilter requires the user to specify the size of the image to be produced as output. This is done using the SetOutputSize()
		fastMarchingDilate->SetOutputRegion( otsuThresholdFilter->GetOutput()->GetBufferedRegion() );
		fastMarchingDilate->SetOutputSpacing( otsuThresholdFilter->GetOutput()->GetSpacing() );
		fastMarchingDilate->SetOutputOrigin( otsuThresholdFilter->GetOutput()->GetOrigin() );
		fastMarchingDilate->SetStoppingValue( dErodeDistance + dMaxAirwayRadius+ 1 );
	
		try
		{
			fastMarchingDilate->Update();
		}
		catch(itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		if (bDebug)
		{
			typedef itk::ImageFileWriter<FloatImageType> WriterFloatType;
			typename WriterFloatType::Pointer writer = WriterFloatType::New();
			
			writer->SetInput( fastMarchingDilate->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/fmt-out.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		thresholdDilation->SetInput( fastMarchingDilate->GetOutput() );
		thresholdDilation->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			writer->SetInput( thresholdDilation->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/fmt-out.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl; std::cerr << excep << std::endl; 
			}
		}
	
		/* Erosion (Thus creating a closing) */
	
		typename ThresholdingFilterType::Pointer thresholdClosing = ThresholdingFilterType::New();
	
		thresholdClosing->SetLowerThreshold( 0.0 ); 
		thresholdClosing->SetUpperThreshold( dMaxAirwayRadius );
		thresholdClosing->SetOutsideValue( 1 ); 
		thresholdClosing->SetInsideValue( 0 );
	
		/* Instantiations */
		typename FastMarchingFilterType::Pointer fastMarchingClose = FastMarchingFilterType::New();
		trialSeeds->Initialize();
		aliveSeeds->Initialize();
	
		// loop through the output image
		// and set all voxels to 0 seed voxels
		typedef itk::ImageRegionConstIterator<FloatImageType>  ConstFloatIteratorType;
		ConstFloatIteratorType floatDilatedImageIterator( fastMarchingDilate->GetOutput(), fastMarchingDilate->GetOutput()->GetLargestPossibleRegion() );
		ConstIteratorType binaryDilatedImageIterator( thresholdDilation->GetOutput(), thresholdDilation->GetOutput()->GetLargestPossibleRegion() );
		uiNumberOfTrialSeeds = 0;
		uiNumberOfAliveSeeds = 0;
		floatDilatedImageIterator.GoToBegin();
		
		for ( binaryDilatedImageIterator.GoToBegin(); !binaryDilatedImageIterator.IsAtEnd(); ++binaryDilatedImageIterator )
		{
			if ( binaryDilatedImageIterator.Get() == 0 )
			{
				node.SetIndex( binaryDilatedImageIterator.GetIndex() );
	
				if (floatDilatedImageIterator.Get() > dMaxAirwayRadius + dErodeDistance) aliveSeeds->InsertElement( uiNumberOfAliveSeeds++, node ); // Alive seed
				else trialSeeds->InsertElement( uiNumberOfTrialSeeds++, node ); // Trial seed
			}
			
			++floatDilatedImageIterator;
		}
	
		if (bDebug)
		{
			std::cout<<std::endl<<std::endl;
			std::cout<<"FastMarching Close:		- Number of Alive Seeds: "<<aliveSeeds->Size()<<" "<<uiNumberOfAliveSeeds << std::endl;
			std::cout<<"                    	- Number of Trial Seeds: "<<trialSeeds->Size()<<" "<<uiNumberOfTrialSeeds << std::endl;
		}
	
		fastMarchingClose->SetTrialPoints( trialSeeds ); // The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()
		fastMarchingClose->SetAlivePoints( aliveSeeds );
	
		fastMarchingClose->SetInput( NULL );
		fastMarchingClose->SetSpeedConstant( 1.0 );  // To solve a simple Eikonal equation
	
		fastMarchingClose->SetOutputSize( thresholdDilation->GetOutput()->GetBufferedRegion().GetSize() ); // The FastMarchingImageFilter requires the user to specify the size of the image to be produced as output. This is done using the SetOutputSize()
		fastMarchingClose->SetOutputRegion( thresholdDilation->GetOutput()->GetBufferedRegion() );
		fastMarchingClose->SetOutputSpacing( thresholdDilation->GetOutput()->GetSpacing() );
		fastMarchingClose->SetOutputOrigin( thresholdDilation->GetOutput()->GetOrigin() );
	
		fastMarchingClose->SetStoppingValue( dErodeDistance + dMaxAirwayRadius+ 1 );
	
		try
		{
			fastMarchingClose->Update();
		}
		catch(itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		if (bDebug)
		{
			typedef itk::ImageFileWriter<FloatImageType> WriterFloatType;
			typename WriterFloatType::Pointer writer = WriterFloatType::New();
			writer->SetInput( fastMarchingClose->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/fmt-out.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		thresholdClosing->SetInput( fastMarchingClose->GetOutput() );
		thresholdClosing->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( thresholdClosing->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/fmtIn-out.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		/* Difference between closed image and ostu-threshold of the original one */
	
		typedef itk::AbsoluteValueDifferenceImageFilter<LabelImageType, LabelImageType, LabelImageType > TAbsoluteValueDifferenceFilter;
		typename TAbsoluteValueDifferenceFilter::Pointer absoluteValueDifferenceFilter = TAbsoluteValueDifferenceFilter::New();
	
		absoluteValueDifferenceFilter->SetInput1( otsuThresholdFilter->GetOutput() );
		absoluteValueDifferenceFilter->SetInput2( thresholdClosing->GetOutput() );
		absoluteValueDifferenceFilter->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			writer->SetInput( absoluteValueDifferenceFilter->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/avd.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		/* Create a sligthly eroded version of the closed image 
		 * This is to prevent any weird effects at the outside of the face
		 */
	
	
		typename ThresholdingFilterType::Pointer thresholdDifference = ThresholdingFilterType::New();
	
		thresholdDifference->SetLowerThreshold( 0.0 ); // We can get this simply by taking a slightly different threshold for the marching in case
		thresholdDifference->SetUpperThreshold( dMaxAirwayRadius+dErodeDistance );
		thresholdDifference->SetOutsideValue( 1 );
		thresholdDifference->SetInsideValue( 0 );
		
		thresholdDifference->SetInput( thresholdClosing->GetInput() ); // The closed image was the input of the closing threshold (we don't want to re-run a fast marching)
	
		/* The masking */
		
		typedef itk::MaskImageFilter<LabelImageType, LabelImageType, LabelImageType > TMaskImageFilter;
		typename TMaskImageFilter::Pointer absoluteValueDifferenceFilterMasked = TMaskImageFilter::New();
	
		absoluteValueDifferenceFilterMasked->SetInput1( absoluteValueDifferenceFilter->GetOutput() );
		absoluteValueDifferenceFilterMasked->SetInput2( thresholdDifference->GetOutput() ); // Second input is the mask
	
		/* Extract largest component of the difference */
	
		if (bDebug) std::cout << "Extracting largest connected component ... ";
	
		typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType > ConnectedComponentType;
		typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType > RelabelComponentType;
		typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType > FinalThresholdingFilterType;
	
		typename ConnectedComponentType::Pointer connected = ConnectedComponentType::New();
		typename RelabelComponentType::Pointer relabel = RelabelComponentType::New();
		typename FinalThresholdingFilterType::Pointer largestComponentThreshold = FinalThresholdingFilterType::New();
		
		//connected->SetFullyConnected( true );
		connected->SetInput ( absoluteValueDifferenceFilterMasked->GetOutput());
		relabel->SetInput( connected->GetOutput() );
		relabel->SetNumberOfObjectsToPrint( 5 );
		relabel->Update(); // Label the components in the image and relabel them so that object numbers increase as the size of the objects decrease
	
		int componentNumber = 0;
		
		if (iComponent <= 0)
		{
			componentNumber = LabelIt<T>(relabel->GetOutput(), upperSeed, upperSeedRadius, bDebug);
			if (componentNumber <= 0 ) componentNumber = LabelIt<T>(relabel->GetOutput(), lowerSeed, lowerSeedRadius, bDebug);
			std::cout<<"Label found = "<<componentNumber<<std::endl;
		}
		else
		{
			componentNumber = iComponent;
		}
	
		
		largestComponentThreshold->SetInput( relabel->GetOutput() );
		largestComponentThreshold->SetLowerThreshold( componentNumber ); // object #1
		largestComponentThreshold->SetUpperThreshold( componentNumber ); // object #1
		largestComponentThreshold->SetInsideValue(1);
		largestComponentThreshold->SetOutsideValue(0);
		largestComponentThreshold->Update(); // Pull out the largest object
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer lccWriter = WriterLabelType::New();
			lccWriter->SetInput( largestComponentThreshold->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/lcc.nhdr";
			lccWriter->SetFileName( filename );
	
			try
			{
				lccWriter->Update();
			}
			catch( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
	
			std::cout << "done." << std::endl;
		}
	
		/* Now do another Otsu thresholding but just around the current segmentation */
	
		/*  For this, we need to make it first a little bit bigger */
	
		typename ThresholdingFilterType::Pointer thresholdExtendedSegmentation =ThresholdingFilterType::New();
	
		thresholdExtendedSegmentation->SetInput( FastMarchIt<T>(largestComponentThreshold->GetOutput(), "Out", dErodeDistance, dMaxAirwayRadius) ); // There are enough few seeds that they can all be considered as part of the trial seeds
		thresholdExtendedSegmentation->SetLowerThreshold( 0.0 );
	
		thresholdExtendedSegmentation->SetUpperThreshold( (sqrt(2.0)-1)*dMaxAirwayRadius ); // To make sure we get roughly twice the volume if the object would have a circular cross section
	
		thresholdExtendedSegmentation->SetOutsideValue( 0 ); 
		thresholdExtendedSegmentation->SetInsideValue( 1 );
		
		thresholdExtendedSegmentation->Update(); // Force it so we have the mask available
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			writer->SetInput( thresholdExtendedSegmentation->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/seg-extended.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl; 
			}
		}
	
		/* Now do another Otsu thresholding but restrict the statistics to the currently obtained area  (custom ostu-threshold filter) */
	
		typedef itk::MaskedOtsuThresholdImageFilter<InputImageType, LabelImageType, LabelImageType >MaskedOtsuThresholdFilterType;
		typename MaskedOtsuThresholdFilterType::Pointer maskedOtsuThresholdFilter = MaskedOtsuThresholdFilterType::New();
	
		// TODO: not sure about these inside/outside settings, check!!
		
		maskedOtsuThresholdFilter->SetInsideValue( 1 );
		maskedOtsuThresholdFilter->SetOutsideValue( 0 );
		maskedOtsuThresholdFilter->SetMaskImage( thresholdExtendedSegmentation->GetOutput() );
		maskedOtsuThresholdFilter->SetInput( originalImage );
		maskedOtsuThresholdFilter->Update();
		
		T dThreshold = maskedOtsuThresholdFilter->GetThreshold(); // Get the threshold used in the otsu-thresholding
		std::cout << "Threshold computed: " << dThreshold << std::endl;
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer labelWriterMaskedOtsu = WriterLabelType::New();
			
			labelWriterMaskedOtsu->SetInput( maskedOtsuThresholdFilter->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/otst-out-masked.nhdr";
			labelWriterMaskedOtsu->SetFileName( filename );
	
			try
			{
				labelWriterMaskedOtsu->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		/* Now mask it again and extract the largest component */
		
		if (bDebug) std::cout << " mask it again and extract the largest component ... " << std::endl;
	
		typename TMaskImageFilter::Pointer maskedOtsu = TMaskImageFilter::New();
	
		maskedOtsu->SetInput1( maskedOtsuThresholdFilter->GetOutput() );
		maskedOtsu->SetInput2( thresholdDifference->GetOutput() ); // Second input is the mask
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			writer->SetInput( maskedOtsu->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/maskedOtsu-out_second.nhdr";
			writer->SetFileName( filename );
	
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
			
			std::cout << " Update ... " << std::endl;
		}
	
		try
		{
			maskedOtsu->Update();
		}
		catch ( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return EXIT_FAILURE;
		}
	
		/* Extract largest the airway with the lungs */
	
		if (bDebug) std::cout << "Extracting final largest connected component ... ";
	
		typename ConnectedComponentType::Pointer connectedFinal = ConnectedComponentType::New();
		typename RelabelComponentType::Pointer relabelFinal = RelabelComponentType::New();
		typename FinalThresholdingFilterType::Pointer finalThreshold = FinalThresholdingFilterType::New();
		
		connectedFinal->SetInput ( maskedOtsu->GetOutput());
		relabelFinal->SetInput( connectedFinal->GetOutput() );
		relabelFinal->SetNumberOfObjectsToPrint( 5 );
		relabelFinal->Update(); // Label the components in the image and relabel them so that object numbers increase as the size of the objects decrease
	
		if (iComponent <= 0)
		{
			componentNumber = LabelIt<T>(relabelFinal->GetOutput(), upperSeed, upperSeedRadius, bDebug);
			if ( componentNumber <= 0 ) componentNumber = LabelIt<T>(relabelFinal->GetOutput(), lowerSeed, lowerSeedRadius, bDebug);
			//std::cout<<"Label found = "<<componentNumber<<std::endl;
		}
		else
		{
			componentNumber = iComponent;
		}
	
		finalThreshold->SetInput( relabelFinal->GetOutput() );
		finalThreshold->SetLowerThreshold( componentNumber ); // object #1
		finalThreshold->SetUpperThreshold( componentNumber ); // object #1
		finalThreshold->SetInsideValue(1);
		finalThreshold->SetOutsideValue(0);
		finalThreshold->Update(); // Pull out the largest object
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			writer->SetInput( finalThreshold->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/airway-with-lung.nrrd";
			writer->SetFileName( filename );
			
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
	
			std::cout<<"..done"<<std::endl;
		}
	
		/* Second part of the code : getting rid of lung automatically */
	
		TOrigin imageOrigin = originalImage->GetOrigin();
		TSpacing imageSpacing = originalImage->GetSpacing();
		TSize imageSize = originalImage->GetBufferedRegion().GetSize();
	
		/* Need to convert from RAS (Slicer) to ITK's coordinate system: LPS */
	
		float ballX, ballY, ballZ;
		
		ballX = -lowerSeed[0];
		ballY = -lowerSeed[1];
		ballZ = lowerSeed[2];
	
		if (bDebug) std::cout << "(x, y, z): " << ballX  << " " << ballY  << " " << ballZ << ", Radius: " << lowerSeedRadius << std::endl;
	
		/* Compute the cubic region around the ball */
	
		int ballRegion[6];
		
		ballRegion[0] = int(floor( ( ballX - lowerSeedRadius - imageOrigin[0] ) / imageSpacing[0] ));
		ballRegion[1] = int(floor( ( ballY - lowerSeedRadius - imageOrigin[1] ) / imageSpacing[1] ));
		ballRegion[2] = int(floor( ( ballZ - lowerSeedRadius - imageOrigin[2] ) / imageSpacing[2] ));
		ballRegion[3] = int(ceil( ( ballX + lowerSeedRadius - imageOrigin[0] ) / imageSpacing[0] ));
		ballRegion[4] = int(ceil( ( ballY + lowerSeedRadius - imageOrigin[1] ) / imageSpacing[1] ));
		ballRegion[5] = int(ceil( ( ballZ + lowerSeedRadius - imageOrigin[2] ) / imageSpacing[2] ));
	
		ballRegion[0] = ballRegion[0] > 0 ? ballRegion[0] : 0;
		ballRegion[1] = ballRegion[1] > 0 ? ballRegion[1] : 0;
		ballRegion[2] = ballRegion[2] > 0 ? ballRegion[2] : 0;
	
		ballRegion[3] = ballRegion[3] < (imageSize[0]-1) ? ballRegion[3] : (imageSize[0]-1);
		ballRegion[4] = ballRegion[4] < (imageSize[1]-1) ? ballRegion[4] : (imageSize[1]-1);
		ballRegion[5] = ballRegion[5] < (imageSize[2]-1) ? ballRegion[5] : (imageSize[2]-1); 
	
		if (bDebug)
		{
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
	
		try
		{
			imageBranch->Allocate();
		}
		catch (itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			std::cerr << "Please verify your parameters, this is often caused by a misplaced trachea carina"<<std::endl;
		}
	
		for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ )
		{
			for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ )
			{
				for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ )
				{
					double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
					double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
					double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;
					TIndex pixelIndexBranch;
					
					pixelIndexBranch[0] = iI - ballRegion[0];
					pixelIndexBranch[1] = iJ - ballRegion[1];
					pixelIndexBranch[2] = iK - ballRegion[2];
					imageBranch->SetPixel( pixelIndexBranch, 0 );
	
					if( iX * iX + iY * iY + iZ * iZ <= lowerSeedRadius * lowerSeedRadius ) // If within the radius of the ball
					{
						TIndex pixelIndex;
						pixelIndex[0] = iI;
						pixelIndex[1] = iJ;
						pixelIndex[2] = iK;
	
						if( finalThreshold->GetOutput()->GetPixel(pixelIndex) )
						{
							finalThreshold->GetOutput()->SetPixel(pixelIndex, 0);
							imageBranch->SetPixel( pixelIndexBranch, 1 );
						}
					}
				}
			}
		}
	
		/* Clean up the ball region: First pass */
	
		/* First with simple threshold */
	
		typename ConnectedComponentType::Pointer connectedBranch = ConnectedComponentType::New();
		typename RelabelComponentType::Pointer relabelBranch = RelabelComponentType::New();
	
		connectedBranch->SetInput( imageBranch );
		connectedBranch->Update();
		relabelBranch->SetInput( connectedBranch->GetOutput() );
		relabelBranch->SetNumberOfObjectsToPrint( 5 );
		relabelBranch->Update();
	
		/* Get geometry statistics */
	
		typedef itk::LabelGeometryImageFilter<LabelImageType> LabelGeometryImageFilterType;
		typename LabelGeometryImageFilterType::Pointer labelBranchGeometry = LabelGeometryImageFilterType::New();
		
		labelBranchGeometry->SetInput( relabelBranch->GetOutput() );
		labelBranchGeometry->CalculateOrientedBoundingBoxOn();
		
		/* Just in case */
		
		labelBranchGeometry->CalculateOrientedLabelRegionsOff();
		labelBranchGeometry->CalculatePixelIndicesOff();
		labelBranchGeometry->CalculateOrientedIntensityRegionsOff();
		labelBranchGeometry->Update();
	
		int nBranchParts = relabelBranch->GetNumberOfObjects();
		int nBranchId = 1;
		
		if( nBranchParts > 1 )
		{
			if (bDebug) std::cout << "Number of parts in branch: " << nBranchParts << std::endl;
	
			double minDist2Ball;
			int minLabel;
			double dBallIndexX = ( ballRegion[3] - ballRegion[0] ) / 2.0;
			double dBallIndexY = ( ballRegion[4] - ballRegion[1] ) / 2.0;
			double dBallIndexZ = ( ballRegion[5] - ballRegion[2] ) / 2.0;
	
			for( int nNumParts=1; nNumParts<=nBranchParts; nNumParts++ )
			{
				typename LabelGeometryImageFilterType::BoundingBoxType boundingBox = labelBranchGeometry->GetBoundingBox( nNumParts );
				double xTmp = ( boundingBox[0] + boundingBox[1] ) / 2.0 - dBallIndexX;
				double yTmp = ( boundingBox[2] + boundingBox[3] ) / 2.0 - dBallIndexY;
				double zTmp = ( boundingBox[4] + boundingBox[5] ) / 2.0 - dBallIndexZ;
				double distTmp = sqrt( xTmp * xTmp + yTmp * yTmp + zTmp * zTmp );
	
				if (bDebug) std::cout << "( " << xTmp << ", " << yTmp << ", " << zTmp << "), " << distTmp << std::endl;
	
				if( nNumParts == 1 || minDist2Ball > distTmp)
				{
					minDist2Ball = distTmp;
					minLabel = nNumParts;
				}
	
				if (bDebug) std::cout << boundingBox << std::endl;
	
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
		branchThreshold->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( branchThreshold->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/FinalThresholdNoBranch.nrrd";
			writer->SetFileName( filename );
			
			try 
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		typename ConnectedComponentType::Pointer connectedFinalWithoutLung = ConnectedComponentType::New();
		typename RelabelComponentType::Pointer relabelFinalWithoutLung = RelabelComponentType::New();
	
		connectedFinalWithoutLung->SetInput( finalThreshold->GetOutput() );
		relabelFinalWithoutLung->SetInput( connectedFinalWithoutLung->GetOutput() );
		relabelFinalWithoutLung->SetNumberOfObjectsToPrint( 5 );
		relabelFinalWithoutLung->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( relabelFinalWithoutLung->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/relabelFinalWithLung.nrrd";
			writer->SetFileName( filename );
			
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
			
			std::cout << "Relabeled" << std::endl;
		}
	
		/* Clean up the ball region: Second pass */
	
		if (bDebug) std::cout << "Get rid of residual lungs in the ball region ... " << std::endl;
	
		/* First get the original data in the lung+airway regions of the ball */
		
		for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ )
		{
			for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ )
			{
				for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ )
				{
					double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
					double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
					double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;
	
					TIndex pixelIndexBranch;
					pixelIndexBranch[0] = iI - ballRegion[0];
					pixelIndexBranch[1] = iJ - ballRegion[1];
					pixelIndexBranch[2] = iK - ballRegion[2];
					imageBranch->SetPixel( pixelIndexBranch, -1024 );
	
					if( iX * iX + iY * iY + iZ * iZ <= lowerSeedRadius * lowerSeedRadius )
					{
						TIndex pixelIndex;
						pixelIndex[0] = iI;
						pixelIndex[1] = iJ;
						pixelIndex[2] = iK;
	
						if( branchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) ) imageBranch->SetPixel( pixelIndexBranch, originalImage->GetPixel( pixelIndex ) );
					}
				}
			}
		}
	
		//Now apply an otsu threshold on it
		typename OtsuThresholdFilterType::Pointer otsuThresholdBranchFilter = OtsuThresholdFilterType::New();
		
		otsuThresholdBranchFilter->SetInsideValue(1);
		otsuThresholdBranchFilter->SetOutsideValue(0);
		otsuThresholdBranchFilter->SetInput( imageBranch );
		otsuThresholdBranchFilter->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( relabelFinalWithoutLung->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "otsuBranchCleaning.nrrd";
			writer->SetFileName( filename );
			
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep ) 
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
			
		if (bDebug) std::cout << "Getting rid of the small lungs parts ... " << std::endl;
	
		for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ )
		{
			for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ )
			{
				for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ )
				{
					double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
					double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
					double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;
	
					TIndex pixelIndexBranch;
					pixelIndexBranch[0] = iI - ballRegion[0];
					pixelIndexBranch[1] = iJ - ballRegion[1];
					pixelIndexBranch[2] = iK - ballRegion[2];
					imageBranch->SetPixel( pixelIndexBranch, 0 );
	
					if( iX * iX + iY * iY + iZ * iZ <= lowerSeedRadius * lowerSeedRadius )
						if( branchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) && otsuThresholdBranchFilter->GetOutput()->GetPixel( pixelIndexBranch ) ) imageBranch->SetPixel( pixelIndexBranch, 1 );
				}
			}
		}
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( relabelFinalWithoutLung->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "CleanedBallRegion.nrrd";
			writer->SetFileName( filename );
			
			try 
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
		
		typename ConnectedComponentType::Pointer connectedCleanedBranch = ConnectedComponentType::New();
		typename RelabelComponentType::Pointer relabelCleanedBranch = RelabelComponentType::New();
	
		connectedCleanedBranch->SetInput( imageBranch );
		connectedCleanedBranch->Update();
		relabelCleanedBranch->SetInput( connectedCleanedBranch->GetOutput() );
		relabelCleanedBranch->SetNumberOfObjectsToPrint( 5 );
		relabelCleanedBranch->Update(); // Of course, get rid of any small residual effects
	
		typename FinalThresholdingFilterType::Pointer cleanedBranchThreshold = FinalThresholdingFilterType::New();
		
		cleanedBranchThreshold->SetInput( relabelCleanedBranch->GetOutput() );
		cleanedBranchThreshold->SetLowerThreshold( 1 );
		cleanedBranchThreshold->SetUpperThreshold( 1 );
		cleanedBranchThreshold->SetInsideValue( 1 );
		cleanedBranchThreshold->SetOutsideValue( 0 );
		cleanedBranchThreshold->Update();
	
		if (bDebug) std::cout << "Final airway label ... " << std::endl;
	
		/* Find the airway label using the pyryna apreture position */
	
		int nNumAirway = 0;
		
		if (iComponent <= 0)
		{
			nNumAirway = LabelIt<T>(relabelFinalWithoutLung->GetOutput(), upperSeed, upperSeedRadius, bDebug);
			if (nNumAirway <= 0 ) nNumAirway = LabelIt<T>(relabelFinalWithoutLung->GetOutput(), lowerSeed, lowerSeedRadius, bDebug);
			std::cout<<"Label found = "<<nNumAirway<<std::endl;
		}
		else
		{
			nNumAirway = iComponent;
		}
	
		/* Check if the maximum label found is 0, 
		 * meaning that no label was found in the nose region
		 * => Nasal cavity probably not segmented !
		*/
		
		if (nNumAirway == 0) 
		{
			std::cerr<<"WARNING !"<<std::endl;
			std::cerr<<"The maximum label found in the spherical region around the pyrina aperture was zero !"<<std::endl;
			std::cerr<<"This probably means that nasal cavity is not segmented (or the point is misplaced)."<<std::endl;
			std::cerr<<" Advice: use --debug to ouput and check all the labels found and/or increase the upperSeedRadius to cover more space"<<std::endl;
			
			if (bNoWarning) return EXIT_FAILURE;
		}
	
		if (bDebug) std::cout << "The label " << nNumAirway << " is picked as the airway." << std::endl;
	
		typename FinalThresholdingFilterType::Pointer finalAirwayThreshold = FinalThresholdingFilterType::New();
		
		finalAirwayThreshold->SetInput( relabelFinalWithoutLung->GetOutput() );
		finalAirwayThreshold->SetLowerThreshold( nNumAirway ); 
		finalAirwayThreshold->SetUpperThreshold( nNumAirway ); 
		finalAirwayThreshold->SetInsideValue(1);
		finalAirwayThreshold->SetOutsideValue(0); 
		finalAirwayThreshold->Update();
	
		if (bDebug)
		{
			typename WriterLabelType::Pointer writer = WriterLabelType::New();
			
			writer->SetInput( finalAirwayThreshold->GetOutput() );
			std::string filename = sDebugFolder;
			filename += "/final_threshold.nrrd";
			writer->SetFileName( filename );
			
			try
			{
				writer->Update();
			}
			catch ( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught !" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
	
		/* Finaly paste the ball back */
	
		if (bDebug) std::cout << "Putting the branches back ... " << std::endl;
	
		for( int iI=ballRegion[0]; iI<=ballRegion[3]; iI++ )
		{
			for( int iJ=ballRegion[1]; iJ<=ballRegion[4]; iJ++ )
			{
				for( int iK=ballRegion[2]; iK<=ballRegion[5]; iK++ )
				{
					double iX = iI * imageSpacing[0] + imageOrigin[0] - ballX;
					double iY = iJ * imageSpacing[1] + imageOrigin[1] - ballY;
					double iZ = iK * imageSpacing[2] + imageOrigin[2] - ballZ;
	
					if( iX * iX + iY * iY + iZ * iZ <= lowerSeedRadius * lowerSeedRadius )
					{
						TIndex pixelIndex;
						
						pixelIndex[0] = iI;
						pixelIndex[1] = iJ;
						pixelIndex[2] = iK;
	
						TIndex pixelIndexBranch;
						
						pixelIndexBranch[0] = iI - ballRegion[0];
						pixelIndexBranch[1] = iJ - ballRegion[1];
						pixelIndexBranch[2] = iK - ballRegion[2];
	
						if( cleanedBranchThreshold->GetOutput()->GetPixel( pixelIndexBranch ) ) finalAirwayThreshold->GetOutput()->SetPixel(pixelIndex, 1);
					}
				}
			}
		}
	
	
		typename LabelImageType::Pointer FinalSegmentation = finalAirwayThreshold->GetOutput();
	
		/* Optionnaly remove the maxillary sinus(es) */
	
		if (bRemoveMaxillarySinuses)
		{
			std::cout << "maxillarySinusesSeeds "<<maxillarySinusesSeeds.size() << std::endl;
			
			/* First thing, Erode to severe the small connection
			 * between the sinuses and the airway
			 * Note that the erosion is very small
			 * Using custom fast marching function
			 */
	
			typename ThresholdingFilterType::Pointer thresholdSlightErosion = ThresholdingFilterType::New();
			
			thresholdSlightErosion->SetLowerThreshold( dMaxAirwayRadius*erosionPercentage ); //Should be set ?
			thresholdSlightErosion->SetUpperThreshold( dMaxAirwayRadius );
			thresholdSlightErosion->SetOutsideValue( 0 );
			thresholdSlightErosion->SetInsideValue( 1 );
	
			thresholdSlightErosion->SetInput( FastMarchIt<T>( FinalSegmentation, "In", dErodeDistance, dMaxAirwayRadius));
			thresholdSlightErosion->Update();
	
			if (bDebug)
			{
				typename WriterLabelType::Pointer writer = WriterLabelType::New();
				
				writer->SetInput( thresholdSlightErosion->GetOutput() );
				std::string filename = sDebugFolder;
				filename += "SlightlyEroderSegmentation.nrrd";
				writer->SetFileName( filename );
	
				try
				{
					writer->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}
			}
	
			/* Create the image that we will substract from the segmentation
			 * It should be all the maxillary sinuses
			 */
	
			typename ConnectedComponentType::Pointer connectedSinuses = ConnectedComponentType::New();
			typename RelabelComponentType::Pointer relabelSinuses = RelabelComponentType::New();
	
			connectedSinuses->SetInput( thresholdSlightErosion->GetOutput() );
			relabelSinuses->SetInput( connectedSinuses->GetOutput() );
			relabelSinuses->SetNumberOfObjectsToPrint( 5 );
			relabelSinuses->Update(); // First label the different part that have been separated by the erosion
	
			typedef typename itk::AddImageFilter< LabelImageType > AddLabelImageFilterType;
			typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType > LabelThresholdFilterType;
			typename AddLabelImageFilterType::Pointer addFilter = AddLabelImageFilterType::New();
			
			typedef itk::ImageDuplicator< LabelImageType > DuplicatorType;
			typename DuplicatorType::Pointer duplicator = DuplicatorType::New(); // The blank image is declared with a duplicator>. Should probably be done otherwise
			
			duplicator->SetInputImage(relabelSinuses->GetOutput());
			duplicator->Update();
			
			typename LabelImageType::Pointer sinusesImage = duplicator->GetOutput();
			
			sinusesImage->FillBuffer(0);
	
			int airwayLabel = LabelIt<T>(relabelSinuses->GetOutput(), upperSeed, upperSeedRadius, bDebug); //Get the airway label
	
			for (int i = 0; i < maxillarySinusesSeeds.size(); ++i) //For each seed
			{
				int seedLabel = LabelIt<T>(relabelSinuses->GetOutput(), maxillarySinusesSeeds[i], maxillarySinusesSeedsRadius, bDebug); //Get the label of the maxillary sinus
				//std::cout<<"Seed Label: "<<seedLabel<<std::endl;
	
				if (airwayLabel == seedLabel) //The airway label MUST be different thant the seed label, Otherwise the airway would be maked out
				{
					std::cerr <<"WARNING !"<<std::endl;
					std::cerr << "The airway label found is equal to the label found with seed #" << i << " (seed =" << maxillarySinusesSeeds[i][0] << ", " << maxillarySinusesSeeds[i][1]  << ", " << maxillarySinusesSeeds[i][2] << ")" << std::endl;
					std::cerr << "Review the seed position and/or the percentage used." << std::endl;
					
					if (bNoWarning) return EXIT_FAILURE; else continue;
				}
	
				typename LabelThresholdFilterType::Pointer thresholdOut = LabelThresholdFilterType::New();
				
				thresholdOut->SetLowerThreshold( seedLabel );
				thresholdOut->SetUpperThreshold( seedLabel );
				thresholdOut->SetOutsideValue( 0 );
				thresholdOut->SetInsideValue( 1 );
				thresholdOut->SetInput( relabelSinuses->GetOutput() );
				thresholdOut->Update(); //Threshold out everything but the region given by the seed
	
				addFilter->SetInput1( sinusesImage );
				addFilter->SetInput2( thresholdOut->GetOutput() );
				addFilter->Update(); //Add it to the blank image
	
				sinusesImage = addFilter->GetOutput();
			}
	
			if (bDebug)
			{
				typename WriterLabelType::Pointer writer = WriterLabelType::New();
				
				writer->SetInput( addFilter->GetOutput() );
				std::string filename = sDebugFolder;
				filename += "UndesiredParts.nrrd";
				writer->SetFileName( filename );
	
				try
				{
					writer->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}
			}
	
			/* Dilate the undesired part
			 * (So they approxemately are their original size)
			 * Note that the erosion is very small
			 * Using custom fast marching function
			 */
	
			typename ThresholdingFilterType::Pointer thresholdSlighDilatation = ThresholdingFilterType::New();
			
			thresholdSlighDilatation->SetLowerThreshold( 0 );
			thresholdSlighDilatation->SetUpperThreshold( dMaxAirwayRadius*erosionPercentage ); //Should be set ?
			thresholdSlighDilatation->SetOutsideValue( 1 );
			thresholdSlighDilatation->SetInsideValue( 0 );
			thresholdSlighDilatation->SetInput( FastMarchIt<T>( addFilter->GetOutput(), "Out", dErodeDistance, dMaxAirwayRadius));
			thresholdSlighDilatation->Update();
	
			if (bDebug)
			{
				typename WriterLabelType::Pointer writer = WriterLabelType::New();
				
				writer->SetInput( thresholdSlighDilatation->GetOutput() );
				std::string filename = sDebugFolder;
				filename += "UndesiredParts_NormalSize.nrrd";
				writer->SetFileName( filename );
	
				try
				{
					writer->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}
			}
	
			/* Mask all the undesired part from the segmentation */
	
			typename TMaskImageFilter::Pointer substractSinusesMask = TMaskImageFilter::New();
			
			substractSinusesMask->SetMaskImage( thresholdSlighDilatation->GetOutput() );
			substractSinusesMask->SetInput( FinalSegmentation );
			substractSinusesMask->Update();
	
			if (bDebug)
			{
				typename WriterLabelType::Pointer writer = WriterLabelType::New();
				
				writer->SetInput( substractSinusesMask->GetOutput() );
				std::string filename = sDebugFolder;
				filename += "TrimmedAirway_Dirty.nrrd";
				writer->SetFileName( filename );
	
				try
				{
					writer->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}
			}
	
			/* Clean up */
			
			typename ConnectedComponentType::Pointer connectedCleanUp = ConnectedComponentType::New();
			typename RelabelComponentType::Pointer relabelCleanUp = RelabelComponentType::New();
	
			connectedCleanUp->SetInput( substractSinusesMask->GetOutput() );
			relabelCleanUp->SetInput( connectedCleanUp->GetOutput() );
			relabelCleanUp->SetNumberOfObjectsToPrint( 5 );
			relabelCleanUp->Update();
	
			nNumAirway = LabelIt<T>(relabelCleanUp->GetOutput(), upperSeed, upperSeedRadius, bDebug);
	
			typename FinalThresholdingFilterType::Pointer thresholdCleanUp = FinalThresholdingFilterType::New();
			
			thresholdCleanUp->SetLowerThreshold( nNumAirway );
			thresholdCleanUp->SetUpperThreshold( nNumAirway );
			thresholdCleanUp->SetOutsideValue( 0 );
			thresholdCleanUp->SetInsideValue( 1 );
			thresholdCleanUp->SetInput( relabelCleanUp->GetOutput() );
			thresholdCleanUp->Update();
	
			FinalSegmentation = thresholdCleanUp->GetOutput();
		}
	
		if (bDebug) std::cout << "Writing the final image ... " << std::endl;
	
		/* Write final image */
		
		typename WriterLabelType::Pointer lccWriterFinal = WriterLabelType::New();
		lccWriterFinal->SetInput( FinalSegmentation );
		lccWriterFinal->SetFileName( outputImage.c_str() );
	
		try
		{
			lccWriterFinal->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		/* Write Surface */
	
		typename itk::AirwaySurfaceWriter<InputImageType, LabelImageType>::Pointer surfaceWriter= itk::AirwaySurfaceWriter<InputImageType, LabelImageType>::New();
		
		surfaceWriter->SetFileName( outputGeometry.c_str() );
		surfaceWriter->SetUseFastMarching(true);
		surfaceWriter->SetMaskImage( FinalSegmentation );
		surfaceWriter->SetInput( originalImage );
		surfaceWriter->SetThreshold( dThreshold );
	
		try
		{
			surfaceWriter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	
		if (bDebug) std::cout << "done." << std::endl;
	
		return EXIT_SUCCESS;
	}

} // end namespace AirwaySegmenter

#endif // AirwaySegmenter_hxx_included
