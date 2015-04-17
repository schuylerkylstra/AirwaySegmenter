cmake_minimum_required(VERSION 2.8.9)
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH OFF)
endif()

project(AirwaySegmenter)

# Find packages
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)

# Needed because ITK will clobber the VTK_LIBRARIES and VTK_INCLUDE_DIRS
vtk_module_config(ADDITIONAL_VTK vtkIOXML)
include(${VTK_USE_FILE})

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

#Compatibility between slicer and the executable
SET(SLICER_ITK_DIFFERENT true CACHE BOOL "Check if compiling with an ITK build different from the one used by Slicer")
CONFIGURE_FILE(${AirwaySegmenter_SOURCE_DIR}/CMake/AirwaySegmenterConfig.h.in
               ${AirwaySegmenter_BINARY_DIR}/AirwaySegmenterConfig.h)

include_directories(
  ${AirwaySegmenter_SOURCE_DIR}
  ${AirwaySegmenter_SOURCE_DIR}/Code
  ${ADDITIONAL_VTK_INCLUDE_DIRS}
)

set( AirwaySegmenterCoreLib_SRCS
  Code/AirwaySegmenterCore.h
  Code/AirwaySegmenterCore.cxx
)

add_library( AirwaySegmenterCore ${AirwaySegmenterCoreLib_SRCS} )
target_link_libraries(AirwaySegmenterCore 
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ADDITIONAL_VTK_LIBRARIES}
)
set(CMAKE_SHARED_LINKER_FLAGS -DDLL_EXPORT)

#Finally the projects
SEMMacroBuildCLI(
  NAME ${PROJECT_NAME}
  TARGET_LIBRARIES AirwaySegmenterCore
                   ${ITK_LIBRARIES}
                   ${VTK_LIBRARIES}
  EXECUTABLE_ONLY
)

SEMMacroBuildCLI(
  NAME AirwayGeometryCreator
  TARGET_LIBRARIES AirwaySegmenterCore
                   ${ITK_LIBRARIES}
                   ${ItkVtkGlue_LIBRARIES}
                   ${VTK_LIBRARIES}
  EXECUTABLE_ONLY
)

SEMMacroBuildCLI(
  NAME DrawLines
  TARGET_LIBRARIES AirwaySegmenterCore
                   ${ITK_LIBRARIES}
  EXECUTABLE_ONLY
)

SEMMacroBuildCLI(
  NAME AirwaySegmenterArgs
  TARGET_LIBRARIES AirwaySegmenterCore
                   ${ITK_LIBRARIES}
                   ${ItkVtkGlue_LIBRARIES}
                   ${VTK_LIBRARIES}
  EXECUTABLE_ONLY
)
  
#Testing
OPTION(BUILD_TESTING "Build tests." ON)

IF(BUILD_TESTING)
  enable_testing()
  add_subdirectory(Testing)
ENDIF()
