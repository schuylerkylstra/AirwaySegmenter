#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "ProgramArguments.h"

#ifndef AirwaySegmenter_AirwayWriteArgs_h_included
#define AirwaySegmenter_AirwayWriteArgs_h_included

namespace AirwaySegmenter {

  void WriteArgsToFile( const ProgramArguments & args );

}

#endif