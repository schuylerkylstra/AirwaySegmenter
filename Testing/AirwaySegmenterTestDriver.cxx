#include <cstdlib>
#include <iostream>
#include <string>

#include <itkImageFileReader.h>
#include <itkTestingComparisonImageFilter.h>

#include "itksys/Process.h"

#if 0
#define DEBUG(txt) \
  { \
  std::cout << txt << std::endl; \
  }
#else
#define DEBUG(txt)
#endif

// Finds the argument index that separates arguments
// meant for the actual test and arguments meant
// for comparing test results
int SeparateArguments(int argc, char* const argv[])
{
  DEBUG("In SeparateArguments");

  for (int i = 0; i < argc; ++i) {
    if ( strcmp( argv[i], "--compare" ) == 0 ) {
      return i;
    }
  }

  DEBUG("Leaving SeparateArguments");

  return argc;
}

// Invoke the test process
int TestDriverInvokeProcess( int argc, char* argv[] )
{
  DEBUG("In TestDriverInvokeProcess");

  // Copy only the argv we want to send to this process.
  // Skip the first argument, which is the test driver
  // process itself.
  DEBUG("Processing args");
  char **processArgv = new char*[argc];
  for ( int i = 1; i < argc; ++i ) {
    processArgv[i-1] = argv[i];
  }
  processArgv[argc-1] = NULL;
  DEBUG("Done processing args");

  itksysProcess *process = itksysProcess_New();
  itksysProcess_SetCommand(process, processArgv);
  itksysProcess_SetPipeShared(process, itksysProcess_Pipe_STDOUT, true);
  itksysProcess_SetPipeShared(process, itksysProcess_Pipe_STDERR, true);
  itksysProcess_Execute(process);
  itksysProcess_WaitForExit(process, NULL);
  delete[] processArgv;

  DEBUG("Process launched");

  int state = itksysProcess_GetState(process);
  switch ( state ) {
    case itksysProcess_State_Error:
    {
      std::cerr << "AirwaySegmenterTestDriver: Process error: "
        << itksysProcess_GetErrorString(process) << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }

    case itksysProcess_State_Exception:
    {
      std::cerr << "AirwaySegmenterTestDriver: Process exception: "
        << itksysProcess_GetExceptionString(process) << std::endl;
      itksysProcess_Delete(process);
    }

    case itksysProcess_State_Executing:
    {
      // this is not a possible state after itksysProcess_WaitForExit
      std::cerr << "AirwaySegmenterTestDriver: Internal error: process can't be in Executing State." << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }
    case itksysProcess_State_Exited:
    {
      // this is the normal case - it is treated later
      break;
    }
    case itksysProcess_State_Expired:
    {
      // this is not a possible state after itksysProcess_WaitForExit
      std::cerr << "AirwaySegmenterTestDriver: Internal error: process can't be in Expired State." << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }
    case itksysProcess_State_Killed:
    {
      std::cerr << "AirwaySegmenterTestDriver: The process has been killed." << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }
    case itksysProcess_State_Disowned:
    {
      std::cerr << "AirwaySegmenterTestDriver: Process disowned." << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }
    default:
    {
      // this is not a possible state after itksysProcess_WaitForExit
      std::cerr << "AirwaySegmenterTestDriver: Internal error: unknown State." << std::endl;
      itksysProcess_Delete(process);
      return 1;
    }
  }

  int retCode = itksysProcess_GetExitValue(process);
  if ( retCode != 0 ) {
    std::cerr << "itkTestDriver: Process exited with return value: " << retCode << std::endl;
    }
  itksysProcess_Delete(process);

  DEBUG("Leaving TestDriverInvokeProcess");

  return retCode;
}

int main(int argc, char* argv[])
{
  int testArgc = SeparateArguments(argc, argv);

  // Process --comparison options
  std::string testFile;
  std::string baselineFile;
  for ( int i = testArgc; i < argc; ++i ) {
    if ( strcmp( argv[i], "--compare" ) == 0 && i + 2 < argc ) {
      testFile = std::string( argv[i+1] );
      baselineFile = std::string( argv[i+2] );
      // Allow only two files to be compared at the moment
      break;
    }
  }

  // Launch the test command as a subprocess
  if ( TestDriverInvokeProcess( testArgc, argv ) != EXIT_SUCCESS ) {
    std::cerr << "Failed to invoke test process " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  if ( testFile.empty() || baselineFile.empty() ) {
    std::cout << "testFile: " << testFile << std::endl;
    std::cout << "baselineFile: " << baselineFile << std::endl;
    return EXIT_SUCCESS;
  }

  // Read in comparison images and check that they match
  typedef float                              PixelType;
  const unsigned int                         Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ImageReaderType;

  ImageReaderType::Pointer testReader = ImageReaderType::New();
  testReader->SetFileName( testFile.c_str() );

  ImageReaderType::Pointer baselineReader = ImageReaderType::New();
  baselineReader->SetFileName( baselineFile.c_str() );

  typedef itk::Testing::ComparisonImageFilter< ImageType, ImageType > ComparisonFilterType;
  ComparisonFilterType::Pointer comparisonFilter = ComparisonFilterType::New();
  comparisonFilter->SetTestInput( testReader->GetOutput() );
  comparisonFilter->SetValidInput( baselineReader->GetOutput() );
  try {
    comparisonFilter->UpdateLargestPossibleRegion();

    ComparisonFilterType::AccumulateType totalImageDifference =
      comparisonFilter->GetTotalDifference();
    std::cout << "Total image difference: " << totalImageDifference << std::endl;

    if ( totalImageDifference > itk::NumericTraits< ComparisonFilterType::AccumulateType >::ZeroValue() ) {
      std::cerr << "Image mismatch between test image and baseline\n";
      return EXIT_FAILURE;
    }
  } catch ( itk::ExceptionObject & except ) {
    std::cerr << except << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
