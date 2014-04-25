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

/* STL includes */
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

/* Local includes */
#include "AirwaySegmenterCLP.h"
#include "AirwaySegmenterConfig.h"
#include "AirwaySegmenter.hxx"
#include "ProgramArguments.h"

#include "itkAirwaySurfaceWriter.h"
#include "itkMaskedOtsuThresholdImageFilter.h"

namespace
{
  /*******************************************************************/
  /** Print all the settings. */
  /*******************************************************************/
  int OutputAllSettings(int argc, char* argv[])
  {
    PARSE_ARGS;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Parameter settings:" << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "input image                  = " << inputImage << std::endl;
    std::cout << "output image                 = " << outputImage << std::endl;
    std::cout << "output geometry              = " << outputGeometry << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "lowerSeed                    = " << lowerSeed[0] << ", "
                                                   << lowerSeed[1] << ", "
                                                   << lowerSeed[2] << std::endl;
    std::cout << "lowerSeedRadius              = " << lowerSeedRadius << std::endl;
    std::cout << "upperSeed                    = " << upperSeed[0] << ", "
                                                   << upperSeed[1] << ", "
                                                   << upperSeed[2]  << std::endl;
    std::cout << "upperSeedRadius              = " << upperSeedRadius << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    for (size_t i = 0; i < airwayFragmentSeeds.size(); ++i) {
      std::cout << "airwayFragmentSeed " << i << " = "
        << airwayFragmentSeeds[i][0] << ", "
        << airwayFragmentSeeds[i][1] << ", "
        << airwayFragmentSeeds[i][2] << std::endl;
    }
    if (airwayFragmentSeeds.size()) {
      std::cout << "----------------------------------------------------------------------------------" << std::endl;
    }

    if (trachealTubeSeed.size() >= 3) {
      std::cout << "trachealTubeSeed             = " << trachealTubeSeed[0] << ", "
                                                     << trachealTubeSeed[1] << ", "
                                                     << trachealTubeSeed[2] << std::endl;
    }
    std::cout << "trachealTubeSeedRadius       = " << trachealTubeSeedRadius << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "bRemoveMaxillarySinuses      = " << bRemoveMaxillarySinuses << std::endl;
    for (size_t i = 0; i < maxillarySinusesSeeds.size(); ++i) {
      std::cout << "maxillarySinusesSeed " << i << " = "
        << maxillarySinusesSeeds[i][0] << ", "
        << maxillarySinusesSeeds[i][1] << ", "
        << maxillarySinusesSeeds[i][2] << std::endl;
    }
    std::cout << "maxillarySinusesSeedsRadius  = " << maxillarySinusesSeedsRadius << std::endl;
    std::cout << "erosionPercentage            = " << erosionPercentage << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "dMaxAirwayRadius             = " << dMaxAirwayRadius << std::endl;
    std::cout << "dErodeDistance               = " << dErodeDistance << std::endl;
    std::cout << "iComponent                   = " << iComponent << std::endl;
    std::cout << "bRemoveBreathingMask         = " << bRemoveBreathingMask << std::endl;
    std::cout << "dBreathingMaskThickness      = " << dBreathingMaskThickness << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "bNoWarning                   = " << bNoWarning <<std::endl;
    std::cout << "bDebug                       = " << bDebug << std::endl;
    std::cout << "sDebugFolder                 = " << sDebugFolder << std::endl;

    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "bRAIImage                    = " << bRAIImage << std::endl;
    std::cout << "sRAIImagePath                = " << sRAIImagePath << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;

    return 0;
  }

  /*******************************************************************/
  /** Query the image type. */
  /*******************************************************************/
  void GetImageType (std::string fileName,itk::ImageIOBase::IOPixelType &pixelType,itk::ImageIOBase::IOComponentType &componentType)
  {
    typedef itk::Image<unsigned char, 3> ImageType;
    itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();

    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();
  }

} // End namespace

/*******************************************************************/
int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  AirwaySegmenter::ProgramArguments args;
  args.inputImage           = inputImage;
  args.outputImage          = outputImage;
  args.outputGeometry       = outputGeometry;
  args.returnParameterFile = returnParameterFile;

  args.lowerSeed       = lowerSeed;
  args.lowerSeedRadius = lowerSeedRadius;
  args.upperSeed       = upperSeed;
  args.upperSeedRadius = upperSeedRadius;

  args.bAddAirwayFragments = bAddAirwayFragments;
  args.airwayFragmentSeeds = airwayFragmentSeeds;

  args.bRemoveBreathingMask    = bRemoveBreathingMask;
  args.dBreathingMaskThickness = dBreathingMaskThickness;

  args.bRemoveTrachealTube          = bRemoveTrachealTube;
  args.trachealTubeSeed             = trachealTubeSeed;
  args.trachealTubeSeedRadius       = trachealTubeSeedRadius;
  args.trachealTubeDilationDistance = trachealTubeDilationDistance;

  args.bRemoveMaxillarySinuses     = bRemoveMaxillarySinuses;
  args.maxillarySinusesSeeds       = maxillarySinusesSeeds;
  args.maxillarySinusesSeedsRadius = maxillarySinusesSeedsRadius;
  args.erosionPercentage           = erosionPercentage;

  args.dMaxAirwayRadius        = dMaxAirwayRadius;
  args.dErodeDistance          = dErodeDistance;
  args.iComponent              = iComponent;

  args.bNoWarning   = bNoWarning;
  args.bDebug       = bDebug;
  args.sDebugFolder = sDebugFolder;

  args.bRAIImage     = bRAIImage;
  args.sRAIImagePath = sRAIImagePath;

  if (bDebug) {
    // Output the arguments
    OutputAllSettings( argc, argv );
  }

  itk::ImageIOBase::IOPixelType     inputPixelType;
  itk::ImageIOBase::IOComponentType inputComponentType;

  int ret = EXIT_FAILURE;

  try {
    GetImageType(inputImage, inputPixelType, inputComponentType);

    switch( inputComponentType ) {
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
        ret = AirwaySegmenter::ExecuteFromFile( args, static_cast<short>(0) );
        break;
#endif
#if defined(SUPPORT_UINT_PIXEL)
      case itk::ImageIOBase::UINT:
        std::cout << "Unsigned int images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_INT_PIXEL)
      case itk::ImageIOBase::INT:
        ret = AirwaySegmenter::ExecuteFromFile( args, static_cast<int>(0) );
        break;
#endif
#if defined(SUPPORT_ULONG_PIXEL)
      case itk::ImageIOBase::ULONG:
        std::cout << "Unsigned long images not supported" << std::endl;
        break;
#endif
#if defined(SUPPORT_LONG_PIXEL)
      case itk::ImageIOBase::LONG:
        ret = AirwaySegmenter::ExecuteFromFile( args, static_cast<long>(0) );
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

  return ret;
}
