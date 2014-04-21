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

#ifndef __itkPhysicalSpaceBinaryDilateImageFilter_h
#define __itkPhysicalSpaceBinaryDilateImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class PhysicalSpaceBinaryDilateImageFilter
 * \brief Performs a binary dilation in physical coordinate space rather than
 * index space. This is particularly useful for dilations performed on images
 * with non-uniform spacing in each dimension. */
template< typename TInputImage, typename TOutputImage >
class ITK_EXPORT PhysicalSpaceBinaryDilateImageFilter :
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard typedefs. */
  typedef PhysicalSpaceBinaryDilateImageFilter            Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( PhysicalSpaceBinaryDilateImageFilter, ImageToImageFilter );

  typedef TInputImage                            InputImageType;
  typedef typename TInputImage::PixelType        InputPixelType;
  typedef typename TInputImage::Pointer          InputImagePointer;
  typedef typename TInputImage::SizeType         InputSizeType;
  typedef typename TInputImage::IndexType        InputIndexType;
  typedef typename TInputImage::RegionType       InputRegionType;
  typedef typename TInputImage::SpacingType      InputSpacingType;
  itkStaticConstMacro( InputImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  typedef TOutputImage                       OutputImageType;
  typedef typename TOutputImage::PixelType   OutputPixelType;
  typedef typename TOutputImage::Pointer     OutputImagePointer;
  typedef typename TOutputImage::SizeType    OutputSizeType;
  typedef typename TOutputImage::IndexType   OutputIndexType;
  typedef typename TOutputImage::RegionType  OutputRegionType;
  typedef typename TOutputImage::SpacingType OutputSpacingType;
  itkStaticConstMacro( OutputImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  /** Get/set the dilation distance. */
  itkGetConstMacro( DilationDistance, double );
  itkSetMacro( DilationDistance, double );

protected:
  PhysicalSpaceBinaryDilateImageFilter();
  ~PhysicalSpaceBinaryDilateImageFilter();
  void PrintSelf( std::ostream & os, Indent indent ) const;

  void GenerateData();

private:
  PhysicalSpaceBinaryDilateImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented

  double m_DilationDistance;

}; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPhysicalSpaceBinaryDilateImageFilter.hxx"
#endif

#endif // __itkPhysicalSpaceBinaryDilateImageFilter_h
