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

#include "AirwayGeometryCreatorCLP.h"
#include "AirwaySegmenterConfig.h"

#include "GetImageType.h"

#include <itkImageFileReader.h>

#include "itkAirwaySurfaceWriter.hxx"

template< typename TPixelType >
int Execute( const std::string & inputImage, const std::string & segmentationImage,
             const std::string & outputGeometry, double threshold, TPixelType )
{
  typedef TPixelType TLabelPixelType;

  const unsigned char DIMENSION = 3;

  typedef itk::Image<TPixelType, DIMENSION>      InputImageType;
  typedef itk::Image<TPixelType, DIMENSION>      OutputImageType;
  typedef itk::Image<TLabelPixelType, DIMENSION> LabelImageType;

  typedef itk::ImageFileReader<InputImageType>   ReaderType;
  typedef itk::ImageFileReader<LabelImageType>   ReaderLabelType;

  typename ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( inputImage );
  try {
    inputReader->Update();
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << "Exception caught when reading input" << std::endl;
    std::cerr << except << std::endl;
  }

  typename ReaderLabelType::Pointer segmentationReader = ReaderLabelType::New();
  segmentationReader->SetFileName( segmentationImage );
  try {
    segmentationReader->Update();
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << "Exception caught when reading segmentation" << std::endl;
    std::cerr << except << std::endl;
  }

  typedef itk::AirwaySurfaceWriter<InputImageType, LabelImageType>
    SurfaceWriterType;
  typename SurfaceWriterType::Pointer surfaceWriter =
    SurfaceWriterType::New();

  surfaceWriter->SetFileName( outputGeometry.c_str() );
  surfaceWriter->SetUseFastMarching( true );
  surfaceWriter->SetMaskImage( segmentationReader->GetOutput() );
  surfaceWriter->SetInput( inputReader->GetOutput() );
  surfaceWriter->SetThreshold( threshold );
  try {
    surfaceWriter->Update();
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << "Exception caught when writing airway surface" << std::endl;
    std::cerr << except << std::endl;
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
    AirwaySegmenter::GetImageType( inputImage, inputPixelType, inputComponentType );

    switch ( inputComponentType ) {
#if defined(SUPPORT_UCHAR_PIXEL)
      case itk::ImageIOBase::UCHAR:
        std::cout << "Unsigned char images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_CHAR_PIXEL)
      case itk::ImageIOBase::CHAR:
        std::cout << "Char images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_USHORT_PIXEL)
      case itk::ImageIOBase::USHORT:
        std::cout << "Unsigned short images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_SHORT_PIXEL)
      case itk::ImageIOBase::SHORT:
        ret = Execute( inputImage, segmentationImage, outputGeometry,
                       threshold, static_cast<short>(0) );
        break;
#endif
#if defined(SUPPORT_UINT_PIXEL)
      case itk::ImageIOBase::UINT:
        std::cout << "Unsigned int images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_INT_PIXEL)
      case itk::ImageIOBase::INT:
        ret = Execute( inputImage, segmentationImage, outputGeometry,
                       threshold, static_cast<int>(0) );
        break;
#endif
#if defined(SUPPORT_ULONG_PIXEL)
      case itk::ImageIOBase::ULONG:
        std::cout << "Unsigned long images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_LONG_PIXEL)
      case itk::ImageIOBase::LONG:
        ret = Execute( inputImage, segmentationImage, outputGeometry,
                       threshold, static_cast<long>(0) );
        break;
#endif
#if defined(SUPPORT_FLOAT_PIXEL)
      case itk::ImageIOBase::FLOAT:
        std::cout << "Float images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_DOUBLE_PIXEL)
      case itk::ImageIOBase::DOUBLE:
        std::cout << "Double images not supported" << std::endl;
        break;
#endif
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
    }
  } catch( itk::ExceptionObject & excep ) {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;

    return EXIT_FAILURE;
  }

  return 0;
}