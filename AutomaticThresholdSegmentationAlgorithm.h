#ifndef vpaw_Segmentation_AutomaticThresholdSegmentationAlgorithm_h_included
#define vpaw_Segmentation_AutomaticThresholdSegmentationAlgorithm_h_included

#include <vector>

#include "Segmentation.h"
#include "SegmentationAlgorithm.h"

#include "AirwaySegmenter.hxx"

#include <itkChangeInformationImageFilter.h>


namespace vpaw {

class DoubleParameter;
class StringParameter;

class Segmentation;


class AutomaticThresholdSegmentationAlgorithm : public SegmentationAlgorithm {
public:
  /** Strings containing parameter names. */
  const static std::string LOWER_SEED_NAME;
  const static std::string LOWER_SEED_RADIUS_NAME;
  const static std::string UPPER_SEED_NAME;
  const static std::string UPPER_SEED_RADIUS_NAME;
  const static std::string MAXIMUM_AIRWAY_RADIUS_NAME;
  const static std::string ERODE_DISTANCE_NAME;
  const static std::string EROSION_PERCENTAGE_NAME;

  AutomaticThresholdSegmentationAlgorithm();
  virtual ~AutomaticThresholdSegmentationAlgorithm();

  void SetLowerSeedPosition( const double position[3] );
  void SetUpperSeedPosition( const double position[3] );

  template< class TInput >
  bool Segment( const TInput * inputImage, Segmentation * segmentation )
  {
    // Assemble program arguments
    AirwaySegmenter::ProgramArguments args;

    // Have to convert from LPS to RAS coordinate system the automatic
    // segmentation code expects.
    args.lowerSeed.resize( 3, 0.0f );
    args.lowerSeed[0]    = -m_LowerSeedPosition[0];
    args.lowerSeed[1]    = -m_LowerSeedPosition[1];
    args.lowerSeed[2]    =  m_LowerSeedPosition[2];
    args.lowerSeedRadius = 20.0;

    args.upperSeed.resize( 3, 0.0f );
    args.upperSeed[0]    = -m_UpperSeedPosition[0];
    args.upperSeed[1]    = -m_UpperSeedPosition[1];
    args.upperSeed[2]    =  m_UpperSeedPosition[2];
    args.upperSeedRadius = 20.0;

    args.dMaxAirwayRadius = 9.0;
    args.dErodeDistance   = 2.0;
    args.iComponent       = -1;

    args.maxillarySinusesSeeds.clear();
    args.maxillarySinusesSeedsRadius = 5.0;

    args.erosionPercentage       = 20.0;
    args.bRemoveMaxillarySinuses = false;

    args.bNoWarning = false;
    args.bDebug     = false;
    args.sDebugFolder = ".";

    args.bRAIImage = false;

    // Change image units from meters to millimeters
    typedef typename itk::ChangeInformationImageFilter< TInput > ChangeFilterType;
    typename ChangeFilterType::Pointer inputInfoChanger = ChangeFilterType::New();
    inputInfoChanger->ChangeSpacingOn();
    inputInfoChanger->ChangeOriginOn();

    typename TInput::SpacingType newSpacing( inputImage->GetSpacing() );
    typename TInput::PointType newOrigin( inputImage->GetOrigin() );
    for ( unsigned int i = 0; i < TInput::ImageDimension; ++i ) {
      newSpacing[i] *= 1e3;
      newOrigin[i]  *= 1e3;
    }
    inputInfoChanger->SetOutputSpacing( newSpacing );
    inputInfoChanger->SetOutputOrigin( newOrigin );
    inputInfoChanger->SetInput( inputImage );

    typename TInput::Pointer outputImage;
    typename TInput::PixelType airwayThreshold;
    int result = AirwaySegmenter::Execute( args,
                                           inputInfoChanger->GetOutput(),
                                           outputImage,
                                           airwayThreshold );

    typename ChangeFilterType::Pointer outputInfoChanger = ChangeFilterType::New();
    outputInfoChanger->ChangeSpacingOn();
    outputInfoChanger->ChangeOriginOn();

    newSpacing = outputImage->GetSpacing();
    newOrigin = outputImage->GetOrigin();
    for ( unsigned int i = 0; i < TInput::ImageDimension; ++i ) {
      newSpacing[i] /= 1e3;
      newOrigin[i]  /= 1e3;
    }
    outputInfoChanger->SetOutputSpacing( newSpacing );
    outputInfoChanger->SetOutputOrigin( newOrigin );
    outputInfoChanger->SetInput( outputImage );
    try {
      outputInfoChanger->Update();
    } catch ( itk::ExceptionObject & excep ) {
      std::cerr << "Exception when updating the output image information\n";
      std::cerr << excep << std::endl;
    }

    segmentation->SetSegmentationImage( outputInfoChanger->GetOutput() );
    segmentation->SetAirwayThreshold( airwayThreshold );

    return ( result == EXIT_SUCCESS );
  }

protected:
  StringParameter * m_LowerSeedName;

  std::vector< float > m_LowerSeedPosition;

  DoubleParameter * m_LowerSeedRadius;

  StringParameter * m_UpperSeedName;

  std::vector< float > m_UpperSeedPosition;

  DoubleParameter * m_UpperSeedRadius;

  /** Set/get the maximal radius for morphological closing (in mm).
   * Should be set roughly a little larger than the maximal
   * expected radius for the airway. */
  DoubleParameter * m_MaximumAirwayRadius;

  /** Erosion distance from estimate of the outer skin layer
   * (in mm) to prevent leaking of the segmentation out of the nose. */
  DoubleParameter * m_ErodeDistance;

  DoubleParameter * m_ErosionPercentage;

};

} // end namespace vpaw

#endif // vpaw_Segmentation_AutomaticThresholdSegmentationAlgorithm_h_included
