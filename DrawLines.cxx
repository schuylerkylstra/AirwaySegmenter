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
#include <cstdlib>
#include <vector>

#include "DrawLinesCLP.h"
#include "AirwaySegmenterConfig.h"
#include "AirwaySegmenterCore.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

/*******************************************************************
 ** Test wether a point is in a cylinder.
 ** Derived from http://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml
 *******************************************************************/
bool InCylinder( const std::vector< float > & testPt,
                 const std::vector< float > & pt1,
                 const std::vector< float > & pt2,
                 double radius )
{
  float dx = pt2[0] - pt1[0];
  float dy = pt2[1] - pt1[1];
  float dz = pt2[2] - pt1[2];
  
  float pdx = testPt[0] - pt1[0];
  float pdy = testPt[1] - pt1[1];
  float pdz = testPt[2] - pt1[2];
  
  float dot = pdx*dx + pdy*dy + pdz*dz;
  
  float lengthSquared = dx*dx + dy*dy + dz*dz;
  
  if ( dot < 0.0 || dot > lengthSquared ) {
    // Point lies outside the ends of the cylinder
    return false;
  } else {
    // Check if distance from point to central axis is <= radius
    float dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot / lengthSquared;
    return ( dsq <= radius*radius );
  }
  
  return false;
}

/*******************************************************************
 ** Draw a line segment in an image.
 *******************************************************************/
template< typename TImage >
void DrawLine( TImage * image,
               const std::vector< float > & pt1,
               const std::vector< float > & pt2,
               double radius, typename TImage::PixelType value )
{
  std::cout << "Drawing line from (" << pt1[0] << ", " << pt1[1] << ", " << pt1[2]
    << ") to (" << pt2[0] << ", " << pt2[1] << "," << pt2[2] << ")" << std::endl;

  typedef typename TImage::PixelType PixelType;
  typedef TImage                     ImageType;
  typedef typename TImage::IndexType IndexType;
  typedef typename TImage::PointType PointType;

  // Simply iterate over all pixels in the image and perform point-in-cylinder test
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  IteratorType iter( image, image->GetBufferedRegion() );
  for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter ) {
    IndexType index = iter.GetIndex();
    PointType physicalPt;
    image->TransformIndexToPhysicalPoint( index, physicalPt );
    
    std::vector< float > testPt( 3, 0.0f );
    testPt[0] = physicalPt[0];
    testPt[1] = physicalPt[1];
    testPt[2] = physicalPt[2];
    
    if ( InCylinder( testPt, pt1, pt2, radius ) ) {
      iter.Set( value );
    }
  }
}

/*******************************************************************
 ** The primary function, templated by pixel type.
 *******************************************************************/
template< typename TPixel >
int Execute( const std::string & inputImage,
             const std::string & outputImage,
             const std::vector< std::vector< float > > & endpoints,
             double radius, TPixel value )
{
  // First, check that there are an even number of endpoints.
  // If not, give an error.
  if ( endpoints.size() % 2 != 0 ) {
    std::cerr << "Draw Lines requires an even number of endpoints.\n";
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef TPixel                             PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage.c_str() );
  try {
    reader->Update();
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << "Exception caught when reading input file.\n";
    std::cerr << except << std::endl;
    return EXIT_FAILURE;
  }

  ImageType * input = reader->GetOutput();

  // Draw lines...
  for ( size_t i = 0; i < endpoints.size() / 2; ++i ) {
    // Transform point from RAS to LPS
    std::vector< float > pt1 = endpoints[2*i];
    pt1[0] = -pt1[0];
    pt1[1] = -pt1[1];
    
    std::vector< float > pt2 = endpoints[2*i + 1];
    pt2[0] = -pt2[0];
    pt2[1] = -pt2[1];

    DrawLine( input, pt1, pt2, radius, value );
  }

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( reader->GetOutput() );
  writer->SetFileName( outputImage.c_str() );
  try {
    writer->Update();
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << "Exception caught when writing output file.\n";
    std::cerr << except << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     inputPixelType;
  itk::ImageIOBase::IOComponentType inputComponentType;

  int ret = EXIT_FAILURE;

  try {
    AirwaySegmenter::GetImageType( inputImage,
                                   inputPixelType,
                                   inputComponentType );

    switch ( inputComponentType ) {
#if defined(SUPPORT_SHORT_PIXEL)
      case itk::ImageIOBase::SHORT:
        ret = Execute( inputImage, outputImage, endpoints, radius, static_cast<short>(value) );
        break;
#endif
#if defined(SUPPORT_INT_PIXEL)
      case itk::ImageIOBase::INT:
        ret = Execute( inputImage, outputImage, endpoints, radius, static_cast<int>(value) );
        break;
#endif
    }
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << except << std::endl;

    return EXIT_FAILURE;
  }

  return ret;
}
