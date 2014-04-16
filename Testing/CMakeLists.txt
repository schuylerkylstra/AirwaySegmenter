find_path( TEST_DATA_DIR 2005.nrrd )
set( TEST_BASELINES_DIR "${TEST_DATA_DIR}/Baselines" )

set( TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/Temporary )
file( MAKE_DIRECTORY ${TEST_OUTPUT_DIR} )

add_executable( AirwaySegmenterTestDriver AirwaySegmenterTestDriver.cxx )
target_link_libraries( AirwaySegmenterTestDriver ${ITK_LIBRARIES} )

get_property( AIRWAY_SEGMENTER_LOCATION
  TARGET AirwaySegmenter
  PROPERTY LOCATION
)
string( REPLACE "$(OutDir)" "Release" AIRWAY_SEGMENTER_LOCATION_RELEASE
  ${AIRWAY_SEGMENTER_LOCATION}
)
#string( REPLACE ".exe" "" AIRWAY_SEGMENTER_LOCATION_RELEASE
#  ${AIRWAY_SEGMENTER_LOCATION_RELEASE}
#)
message( ${AIRWAY_SEGMENTER_LOCATION_RELEASE} )

add_test( NAME SegmentationTest1
          CONFIGURATIONS Release
          COMMAND AirwaySegmenterTestDriver
            "${AIRWAY_SEGMENTER_LOCATION_RELEASE}"
            --lowerSeed 5.22461,98.0924,-167.275
            --upperSeed 14.8379,201.318,-88.9896
            --RAIImage
            --RAIImagePath "${TEST_OUTPUT_DIR}/2005_RAIImage.nrrd"
            "${TEST_DATA_DIR}/2005.nrrd"
            "${TEST_OUTPUT_DIR}/2005_Segmented.nrrd"
            "${TEST_OUTPUT_DIR}/2005_Segmented.vtp"
          --compare
            "${TEST_OUTPUT_DIR}/2005_Segmented.nrrd"
            "${TEST_BASELINES_DIR}/2005.nrrd"
)
