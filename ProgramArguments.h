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
//  Authors: Marc Niethammer, Yi Hong, Johan Andruejol
=============================================================================*/

#ifndef AirwaySegmenter_ProgramArguments_h_included
#define AirwaySegmenter_ProgramArguments_h_included

namespace AirwaySegmenter {

class ProgramArguments {
public:
  std::string inputImage;
  std::string outputImage;
  std::string outputGeometry;

  std::vector< float > lowerSeed;
  double               lowerSeedRadius;
  std::vector< float > upperSeed;
  double               upperSeedRadius;

  double dMaxAirwayRadius;
  double dErodeDistance;
  int    iMaximumNumberOfCVIterations;
  double dCVLambda;
  int    iComponent;

  std::vector< std::vector< float > > maxillarySinusesSeeds;
  double                              maxillarySinusesSeedsRadius;
  double                              erosionPercentage;
  bool                                bRemoveMaxillarySinuses;

  bool        bNoWarning;
  bool        bDebug;
  std::string sDebugFolder;

  bool        bRAIImage;
  std::string sRAIImagePath;
};

} // end namespace AirwaySegmenter


#endif // AirwaySegmenter_ProgramArguments_h_included
