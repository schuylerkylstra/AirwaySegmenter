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
#ifndef AirwaySegmenterArgs_hxx_included
#define AirwaySegmenterArgs_hxx_included

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
#include "ProgramArgumentsArgs.h"

namespace AirwaySegmenterArgs 
{

  int Execute( const AirwaySegmenterArgs::ProgramArguments & args )
  {
    std::string argsFileName = args.argsFile;
    std::ofstream argsFile;

    argsFile.open(argsFileName.c_str()); 

    
    argsFile << " --createOutputGeometry";

    argsFile << " --lowerSeed " << args.lowerSeed[0] << ","
                                << args.lowerSeed[1] << ","
                                << args.lowerSeed[2];
    argsFile << " --lowerSeedRadius " << args.lowerSeedRadius;

    argsFile << " --upperSeed " << args.upperSeed[0] << ","
                                << args.upperSeed[1] << ","
                                << args.upperSeed[2] ;
    argsFile << " --upperSeedRadius " << args.upperSeedRadius;

    if (args.airwayFragmentSeeds.size()){
      argsFile << " --addAirwayFragments";
    }
    for (size_t i = 0; i < args.airwayFragmentSeeds.size(); ++i) {
      argsFile << " --airwayFragmentSeed "
        << args.airwayFragmentSeeds[i][0] << ","
        << args.airwayFragmentSeeds[i][1] << ","
        << args.airwayFragmentSeeds[i][2];
    }

    if (args.trachealTubeSeed.size() >= 3) {
      argsFile << " --trachealTubeSeed " << args.trachealTubeSeed[0] << ","
                                         << args.trachealTubeSeed[1] << ","
                                         << args.trachealTubeSeed[2];
    }
    argsFile << " --trachealTubeSeedRadius " << args.trachealTubeSeedRadius;

    if (args.bRemoveMaxillarySinuses) {
      argsFile << " --removeMaxillarySinuses";
    }
    for (size_t i = 0; i < args.maxillarySinusesSeeds.size(); ++i) {
      argsFile << " --maxillarySinusesSeed "
        << args.maxillarySinusesSeeds[i][0] << ","
        << args.maxillarySinusesSeeds[i][1] << ","
        << args.maxillarySinusesSeeds[i][2];
    }
    argsFile << " --maxillarySinusesSeedsRadius " << args.maxillarySinusesSeedsRadius;

    argsFile << " --erosionPercentage " << args.erosionPercentage;
    argsFile << " --maxAirwayRadius " << args.dMaxAirwayRadius;
    argsFile << " --erodeDistance " << args.dErodeDistance;
    argsFile << " --component " << args.iComponent;
    
    if(args.bRemoveBreathingMask){
      argsFile << " --breathingMaskThickness " << args.dBreathingMaskThickness;
    }
    
    argsFile << " --noWarning";
    if(args.bRAIImage){
      argsFile << " --RAIImage";
      argsFile << " --RAIImagePath " << args.sRAIImagePath << "nrrd";
    }
    argsFile << args.inputImage << " ";
    argsFile << args.outputImage<< " ";
    argsFile << args.outputGeometry;
    
    argsFile.close();
    return EXIT_SUCCESS;
  }

  /*******************************************************************/
  /** Execute the algorithm on an image read from a file. */
  /*******************************************************************/
  template <class T>
  int ExecuteFromFile( const AirwaySegmenterArgs::ProgramArguments & args, T)
  {
    return Execute( args );
  }

} // end namespace AirwaySegmenterArgs

#endif // AirwaySegmenterArgs_hxx_included
