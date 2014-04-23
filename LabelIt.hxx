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
#ifndef AirwaySegmenter_LabelIt_hxx_included
#define AirwaySegmenter_LabelIt_hxx_included

namespace AirwaySegmenter {

  /*******************************************************************/
  /* Look for the labels in the spherical region within radius of the
   * ball center and return the most represented label.
   */
  /*******************************************************************/
  template <class T>
  int LabelIt( typename itk::Image<T, 3>::Pointer image,
               std::vector<float> ballCenter,
               double radius,
               bool printLabels )
  {
    typedef itk::Image< T, 3 >                       TImage;
    typedef typename itk::Image< T, 3 >::SizeType    TSize;
    typedef typename itk::Image< T, 3 >::SpacingType TSpacing;
    typedef typename itk::Image< T, 3 >::PointType   TOrigin;
    typedef typename itk::Image< T, 3 >::IndexType   TIndex;

    TOrigin imageOrigin = image->GetOrigin();
    TSpacing imageSpacing = image->GetSpacing();
    TSize imageSize = image->GetBufferedRegion().GetSize();

    /* Convert to LPS system */
    float x, y, z;

    x = -ballCenter[0];
    y = -ballCenter[1];
    z =  ballCenter[2];

    /* Get the bounding box around the sphere. */
    int region[6];

    region[0] = int( floor( ( x - radius - imageOrigin[0] ) / imageSpacing[0] ) );
    region[1] = int( floor( ( y - radius - imageOrigin[1] ) / imageSpacing[1] ) );
    region[2] = int( floor( ( z - radius - imageOrigin[2] ) / imageSpacing[2] ) );
    region[3] = int(  ceil( ( x + radius - imageOrigin[0] ) / imageSpacing[0] ) );
    region[4] = int(  ceil( ( y + radius - imageOrigin[1] ) / imageSpacing[1] ) );
    region[5] = int(  ceil( ( z + radius - imageOrigin[2] ) / imageSpacing[2] ) );

    region[0] = region[0] > 0 ? region[0] : 0;
    region[1] = region[1] > 0 ? region[1] : 0;
    region[2] = region[2] > 0 ? region[2] : 0;

    region[3] = region[3] < (imageSize[0]-1) ? region[3] : (imageSize[0]-1);
    region[4] = region[4] < (imageSize[1]-1) ? region[4] : (imageSize[1]-1);
    region[5] = region[5] < (imageSize[2]-1) ? region[5] : (imageSize[2]-1);

    if (printLabels) {
      std::cout << "Region: " << region[0] << ", " << region[1] << ", " << region[2] << ", " << region[3] << ", " << region[4] << ", "  << region[5] << std::endl;
    }

    std::map<int, int> labels;
    for ( int iI = region[0]; iI <= region[3]; iI++ ) {
      for( int iJ = region[1]; iJ <= region[4]; iJ++ ) {
        for( int iK = region[2]; iK <= region[5]; iK++ ) {

          /* Get real space position */
          double iX = iI * imageSpacing[0] + imageOrigin[0] - x;
          double iY = iJ * imageSpacing[1] + imageOrigin[1] - y;
          double iZ = iK * imageSpacing[2] + imageOrigin[2] - z;

          /* If within the nose ball region */
          if( iX * iX + iY * iY + iZ * iZ <= radius * radius )
          {
            TIndex pixelIndex;
            pixelIndex[0] = iI;
            pixelIndex[1] = iJ;
            pixelIndex[2] = iK;

            int label = image->GetPixel( pixelIndex );

            // Ignore label 0, it's always the background
            if( label != 0 ) {
              labels[label] += 1;
            }
          }
        }
      }
    }

    int finalLabel = 0;
    int labelCount = 0;

    for  (std::map<int, int>::const_iterator it = labels.begin(); it != labels.end(); ++it) {
      if ( printLabels ) {
        std::cout << "Labels " << it->first << " :   " << it->second << std::endl;
      }

      if (it->second > labelCount) {
        labelCount = it->second;
        finalLabel = it->first;
      }
    }

    return finalLabel;
  }

} // end namespace AirwaySegmenter

#endif
