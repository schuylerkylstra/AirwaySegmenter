#include "AirwayWriteArgs.h"

namespace AirwaySegmenter {

  void StoreArgs( const ProgramArguments & args )
  {
    std::ofstream argsFile;

    argsFile.open(args.argsFile.c_str()); 

    
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
      argsFile << " --trachealTubeSeedRadius " << args.trachealTubeSeedRadius;
    }
    

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
      argsFile << " --RAIImagePath " << args.sRAIImagePath;
    }

    argsFile.close();
  }

}


