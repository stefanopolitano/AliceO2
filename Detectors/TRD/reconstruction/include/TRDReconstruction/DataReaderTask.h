// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   DataReaderTask.h
/// @author Sean Murray
/// @brief  TRD epn task to read incoming data

#ifndef O2_TRD_DATAREADERTASK
#define O2_TRD_DATAREADERTASK

#include "Framework/Task.h"
#include "Framework/DataProcessorSpec.h"
#include "TRDReconstruction/CruRawReader.h"
#include "TRDReconstruction/CompressedRawReader.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/Digit.h"
//#include "DataFormatsTRD/FlpStats.h"

#include <fstream>

using namespace o2::framework;

namespace o2::trd
{

class DataReaderTask : public Task
{
 public:
  DataReaderTask(bool compresseddata, bool byteswap, bool fixdigitendcorruption, int tracklethcheader, bool verbose, bool headerverbose, bool dataverbose) : mCompressedData(compresseddata), mByteSwap(byteswap), mFixDigitEndCorruption(fixdigitendcorruption), mTrackletHCHeaderState(tracklethcheader), mVerbose(verbose), mHeaderVerbose(headerverbose), mDataVerbose(dataverbose) {}
  ~DataReaderTask() override = default;
  void init(InitContext& ic) final;
  void sendData(ProcessingContext& pc, bool blankframe = false);
  void run(ProcessingContext& pc) final;
  bool isTimeFrameEmpty(ProcessingContext& pc);

 private:
  CruRawReader mReader;                  // this will do the parsing, of raw data passed directly through the flp(no compression)
  CompressedRawReader mCompressedReader; //this will handle the incoming compressed data from the flp
                                         // in both cases we pull the data from the vectors build message and pass on.
                                         // they will internally produce a vector of digits and a vector tracklets and associated indexing.
                                         // TODO templatise this and 2 versions of datareadertask, instantiated with the relevant parser.

  bool mVerbose{false};          // verbos output general debuggign and info output.
  bool mDataVerbose{false};      // verbose output of data unpacking
  bool mHeaderVerbose{false};    // verbose output of headers
  bool mCompressedData{false};   // are we dealing with the compressed data from the flp (send via option)
  bool mByteSwap{true};          // whether we are to byteswap the incoming data, mc is not byteswapped, raw data is (too be changed in cru at some point)
                                 //  o2::header::DataDescription mDataDesc; // Data description of the incoming data
  int mTrackletHCHeaderState{0}; // what to do about tracklethcheader, 0 never there, 2 always there, 1 there iff tracklet data, i.e. only there if next word is *not* endmarker 10001000.

  std::string mDataDesc;
  o2::header::DataDescription mUserDataDescription = o2::header::gDataDescriptionInvalid; // alternative user-provided description to pick
  bool mFixDigitEndCorruption{false};                                                     // fix the parsing of corrupt end of digit data. bounce over it.
};

} // namespace o2::trd

#endif // O2_TRD_DATAREADERTASK
