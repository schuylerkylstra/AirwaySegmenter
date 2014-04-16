find_path( TEST_DATA_DIR 2005.nrrd )

set( TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/Temporary )
file( MAKE_DIRECTORY ${TEST_OUTPUT_DIR} )

add_test( NAME SegmentationTest1
          COMMAND AirwaySegmenter
          --lowerSeed 5.22461,98.0924,-167.275
          --upperSeed 14.8379,201.318,-88.9896
          --RAIImage
          --RAIImagePath "${TEST_OUTPUT_DIRECTORY}/2005_RAIImage.nrrd"
          "${TEST_DATA_DIR}/2005.nrrd"
          "${TEST_OUTPUT_DIR}/2005_Segmented.nrrd"
          "${TEST_OUTPUT_DIR}/2005_Segmented.vtp"
)
