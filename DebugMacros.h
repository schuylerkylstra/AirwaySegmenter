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
#ifndef AirwaySegmenter_DebugMacros_h_included
#define AirwaySegmenter_DebugMacros_h_included

#define DEBUG_WRITE_LABEL_IMAGE_ONLY( image )                           \
    if ( args.bDebug ) {                                                \
      typename WriterLabelType::Pointer writer =                        \
        WriterLabelType::New();                                         \
                                                                        \
      writer->SetInput( image );                                        \
      std::string filename = args.sDebugFolder;                         \
      filename += "/"#image".nrrd";                                     \
      writer->SetFileName( filename );                                  \
                                                                        \
      try {                                                             \
        writer->Update();                                               \
      } catch ( itk::ExceptionObject & excep ) {                        \
        std::cerr << "Exception caught "                                \
          << "when writing " << image                                   \
          << ".nrrd!" << std::endl;                                     \
        std::cerr << excep << std::endl;                                \
      }                                                                 \
    }                                                                   \

#define DEBUG_WRITE_IMAGE( filter )                                     \
    if ( args.bDebug ) {                                                \
      typedef itk::ImageFileWriter< FloatImageType > WriterFloatType;   \
      typename WriterFloatType::Pointer writer =                        \
        WriterFloatType::New();                                         \
                                                                        \
      writer->SetInput( filter->GetOutput() );                          \
      std::string filename = args.sDebugFolder;                         \
      filename += "/"#filter".nrrd";                                    \
      writer->SetFileName( filename );                                  \
                                                                        \
      try {                                                             \
        writer->Update();                                               \
      } catch ( itk::ExceptionObject & excep ) {                        \
        std::cerr << "Exception caught "                                \
          << "when writing " << filter                                  \
          << ".nrrd!" << std::endl;                                     \
        std::cerr << excep << std::endl;                                \
      }                                                                 \
    }                                                                   \

#define DEBUG_WRITE_LABEL_IMAGE( filter )                               \
    if ( args.bDebug ) {                                                \
      typename WriterLabelType::Pointer writer =                        \
        WriterLabelType::New();                                         \
                                                                        \
      writer->SetInput( filter->GetOutput() );                          \
      std::string filename = args.sDebugFolder;                         \
      filename += "/"#filter".nrrd";                                    \
      writer->SetFileName( filename );                                  \
                                                                        \
      try {                                                             \
        writer->Update();                                               \
      } catch ( itk::ExceptionObject & excep ) {                        \
        std::cerr << "Exception caught "                                \
          << "when writing " << filter                                  \
          << ".nrrd!" << std::endl;                                     \
        std::cerr << excep << std::endl;                                \
      }                                                                 \
    }                                                                   \

#define TRY_UPDATE( filter )                                            \
    try {                                                               \
      filter->Update();                                                 \
    } catch ( itk::ExceptionObject & exception ) {                      \
      std::cerr << "Exception caught when updating filter " #filter "!" \
        << std::endl;                                                   \
      std::cerr << exception << std::endl;                              \
    }                                                                   \

#endif
