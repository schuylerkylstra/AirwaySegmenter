# How do I specify a directory without specifying a file?
find_path( TEST_DATA_DIR 1003.nrrd )
set( TEST_BASELINES_DIR "${TEST_DATA_DIR}/Baseline" )

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

add_test( NAME SegmentationTest_1003
          CONFIGURATIONS Release
          COMMAND AirwaySegmenterTestDriver
            "${AIRWAY_SEGMENTER_LOCATION_RELEASE}"
            --lowerSeed 15.2148,163.774,-50.7766
            --upperSeed 34.5117,267.316,67.2616
            --RAIImage
            --RAIImagePath "${TEST_OUTPUT_DIR}/1003_RAIImage.nrrd"
            "${TEST_DATA_DIR}/Input/1003.nrrd"
            "${TEST_OUTPUT_DIR}/1003_Segmented.nrrd"
            "${TEST_OUTPUT_DIR}/1003_Segmented.vtp"
          --compare
            "${TEST_OUTPUT_DIR}/1003_Segmented.nrrd"
            "${TEST_BASELINES_DIR}/1003.nrrd"
)

# add_test( NAME SegmentationTest1
          # CONFIGURATIONS Release
          # COMMAND AirwaySegmenterTestDriver
            # "${AIRWAY_SEGMENTER_LOCATION_RELEASE}"
            # --lowerSeed 5.22461,98.0924,-167.275
            # --upperSeed 14.8379,201.318,-88.9896
            # --RAIImage
            # --RAIImagePath "${TEST_OUTPUT_DIR}/2005_RAIImage.nrrd"
            # "${TEST_DATA_DIR}/Input/2005.nrrd"
            # "${TEST_OUTPUT_DIR}/2005_Segmented.nrrd"
            # "${TEST_OUTPUT_DIR}/2005_Segmented.vtp"
          # --compare
            # "${TEST_OUTPUT_DIR}/2005_Segmented.nrrd"
            # "${TEST_BASELINES_DIR}/2005.nrrd"
# )
