# FIXME - How do I specify a directory without specifying a file?
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
            --debug
            --debugFolder "${TEST_OUTPUT_DIR}"
            --addAirwayFragments
            --airwayFragmentSeed 12.988,203.453,66.529
            --airwayFragmentSeed 12.988,180.216,-18.903
            --breathingMaskThickness 0
            --RAIImage
            --RAIImagePath "${TEST_OUTPUT_DIR}/1003_RAIImage.nrrd"
            "${TEST_DATA_DIR}/Input/1003.nrrd"
            "${TEST_OUTPUT_DIR}/1003_Segmented.nrrd"
            "${TEST_OUTPUT_DIR}/1003_dummy.vtp"
          --compare
            "${TEST_OUTPUT_DIR}/1003_Segmented.nrrd"
            "${TEST_BASELINES_DIR}/1003.nrrd"
)
