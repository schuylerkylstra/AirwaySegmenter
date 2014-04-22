project(AirwaySegmenter)
cmake_minimum_required(VERSION 2.8)

#Find packages
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

find_package(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

option( SUPPORT_UCHAR_PIXEL "Build with support for input images with unsigned char pixels." OFF )
option( SUPPORT_CHAR_PIXEL "Build with support for input images with char pixels." OFF )
option( SUPPORT_USHORT_PIXEL "Build with support for input images with unsigned short pixels." ON )
option( SUPPORT_SHORT_PIXEL "Build with support for input images with short pixels." ON )
option( SUPPORT_UINT_PIXEL "Build with support for input images with unsigned int pixels." OFF )
option( SUPPORT_INT_PIXEL "Build with support for input images with int pixels." OFF )
option( SUPPORT_ULONG_PIXEL "Build with support for input images with unsigned long pixels." OFF )
option( SUPPORT_LONG_PIXEL "Build with support for input images with long pixels." OFF )
option( SUPPORT_FLOAT_PIXEL "Build with support for input images with float pixels." OFF )
option( SUPPORT_DOUBLE_PIXEL "Build with support for input images with double pixels." OFF )

IF (ITK_VERSION_MAJOR VERSION_LESS "4")
  FATAL_ERROR("ITK version must be at least 4.0 ")
ENDIF()

#Compatibility between slicer and the executable
SET(SLICER_ITK_DIFFERENT true CACHE BOOL "Check if compiling with an ITK build different from the one used by Slicer")
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/CMake/AirwaySegmenterConfig.h.in
               ${CMAKE_BINARY_DIR}/AirwaySegmenterConfig.h)

include_directories(Code)

#Finally the project
SEMMacroBuildCLI(
  NAME ${PROJECT_NAME}
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRAIRIES}
  EXECUTABLE_ONLY
  )
  
#Testing
OPTION(BUILD_TESTING "Build tests." ON)

IF(BUILD_TESTING)
  enable_testing()
  add_subdirectory(Testing)
ENDIF()
