/*=============================================================================
//  --- Airway Segmenter ---+
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Authors: Marc Niethammer, Yi Hong, Johan Andruejol, Cory Quammen
=============================================================================*/
#ifndef AirwaySegmenter_FastMarchIt_hxx_included
#define AirwaySegmenter_FastMarchIt_hxx_included

#include "DebugMacros.h"

namespace AirwaySegmenter {

  /*******************************************************************/
  /* This Fast marching is factorised for preventing the memory to be
   * overwhelmed by unused nodes. Now that this function is for
   * internal use, it assumes all its parameters are given correctly.
   */
  /*******************************************************************/
  template <class T>
  itk::Image<float, 3>::Pointer FastMarchIt( typename itk::Image<T, 3>::Pointer image,
                                             std::string type,
                                             double erodedDistance,
                                             double airwayRadius )
  {
    /* Necessaries typedefs */

    typedef itk::Image< float, 3 >                         FloatImageType;
    typedef itk::Image< T, 3 >                             LabelImageType;
    typedef itk::FastMarchingImageFilter< FloatImageType,
                                          FloatImageType > FastMarchingFilterType;
    typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
    typedef typename FastMarchingFilterType::NodeType      NodeType;
    typedef itk::ImageRegionConstIterator<LabelImageType>  ConstIteratorType;

    /* Instantiations */
    typename FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
    typename NodeContainer::Pointer seeds = NodeContainer::New();

    seeds->Initialize();
    NodeType node; //Nodes are created as stack variables and  initialized with a value and an itk::Index position. NodeType node;
    node.SetValue( 0.0 ); // Seed value is 0 for all of them, because these are all starting nodes
    ConstIteratorType it( image, image->GetLargestPossibleRegion() ); // Loop through the output image and set all voxels to 0 seed voxels

    unsigned int uiNumberOfSeeds = 0;

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {
      if (type.compare("Out") == 0) { // Dilation
        if ( it.Get() > 0 ) {
          node.SetIndex( it.GetIndex() );
          seeds->InsertElement( uiNumberOfSeeds++, node );
        }
      }
      else if (type.compare("In") == 0) { // Erosion
        if ( it.Get() == 0 ) {
          node.SetIndex( it.GetIndex() );
          seeds->InsertElement( uiNumberOfSeeds++, node );
        }
      }
    }

    fastMarching->SetTrialPoints( seeds ); // The set of seed nodes is now passed to the FastMarchingImageFilter with the method SetTrialPoints()

    /*
     * The FastMarchingImageFilter requires the user to specify
     * the size of the image to be produced as output
     * This is done using the SetOutputSize()
     */
    fastMarching->SetInput( NULL );
    fastMarching->SetSpeedConstant( 1.0 );  // To solve a simple Eikonal equation
    fastMarching->SetOutputSize( image->GetBufferedRegion().GetSize() );
    fastMarching->SetOutputRegion( image->GetBufferedRegion() );
    fastMarching->SetOutputSpacing( image->GetSpacing() );
    fastMarching->SetOutputOrigin( image->GetOrigin() );
    fastMarching->SetStoppingValue( airwayRadius + erodedDistance + 1 );

    TRY_UPDATE( fastMarching );

    return fastMarching->GetOutput();
  }

} // end namespace AirwaySegmenter

#endif
