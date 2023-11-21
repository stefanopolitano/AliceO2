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

/// \file Occupancy.cxx
/// \brief Study to measure the occupancy of the ITS layers
/// \author Stefano Politanò (INFN Torino)

#include "ITSStudies/Occupancy.h"
#include "ITSStudies/ITSStudiesConfigParam.h"

#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "Steer/MCKinematicsReader.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "ITStracking/IOUtils.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DetectorsCommonDataFormats/DetID.h"

#include <numeric>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TNtuple.h>

namespace o2
{
namespace its
{
namespace study
{
using namespace o2::framework;
using namespace o2::globaltracking;

using GTrackID = o2::dataformats::GlobalTrackID;
using PVertex = o2::dataformats::PrimaryVertex;
using V0 = o2::dataformats::V0;
using ITSCluster = o2::BaseCluster<float>;
using mask_t = o2::dataformats::GlobalTrackID::mask_t;
using Track = o2::track::TrackParCov;
using TrackITS = o2::its::TrackITS;
using DCA = o2::dataformats::DCA;
using PID = o2::track::PID;

class Occupancy : public Task
{
 public:
  Occupancy(std::shared_ptr<DataRequest> dr,
            std::shared_ptr<o2::base::GRPGeomRequest> gr,
            bool isMC,
            std::shared_ptr<o2::steer::MCKinematicsReader> kineReader) : mDataRequest{dr}, mGGCCDBRequest(gr), mUseMC(isMC), mKineReader(kineReader){};
  ~Occupancy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDict = d; }

 private:
  // Other functions
  void process(o2::globaltracking::RecoContainer&);

  // Helper functions
  void prepareOutput();
  void updateTimeDependentParams(ProcessingContext& pc);
  void fillIBmap(TH2F *histo, int sta, int chipInMod, float weight);
  void fillOBmap(TH2F *histo, int sta, int chipInMod, float weight, int ssta);
  void getClusters(const gsl::span<const o2::itsmft::CompClusterExt>, float weight);
  void saveHistograms();

  // Running options
  bool mUseMC;

  // Data
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::shared_ptr<DataRequest> mDataRequest;
  std::vector<int> mClusterSizes;
  gsl::span<const int> mInputITSidxs;
  std::vector<o2::MCTrack> mMCTracks;
  const o2::itsmft::TopologyDictionary* mDict = nullptr;

  // Output plots
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  std::vector<TH2F*> mOccupancyHistos{};
  std::string mOutName;
  std::shared_ptr<o2::steer::MCKinematicsReader> mKineReader;
};

void Occupancy::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  prepareOutput();
}

void Occupancy::prepareOutput()
{
    auto& params = o2::its::study::ITSAvgClusSizeParamConfig::Instance();
    mOutName = params.outFileName;
    mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>(mOutName.c_str(), "recreate");

    std::vector<int> nStaves{12, 16, 20, 24, 30, 42, 48};

    for (int layer{0}; layer < 7; layer++) {
        if (layer < 3)
        {
            mOccupancyHistos.push_back(new TH2F(Form("Occupancy chip map L%i", layer), "; Chip ID; Stave ID; # Hits / # PVs", 9, -0.5, 8.5, nStaves[layer], -0.5, nStaves[layer] - 0.5));
        }
        else {
            mOccupancyHistos.push:back(new TH2F(Form("Occupancy chip map L%i", layer), "; Chip ID; Stave ID; #LT Cluster size #GT", 49, -0.5, 48.5, 4 * nStaves[layer], -0.5, 4 * nStaves[layer] - 0.5));
        }
        mOccupancyHistos[layer]->SetDirectory(nullptr);
    }
}


void Occupancy::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  updateTimeDependentParams(pc); // Make sure this is called after recoData.collectData, which may load some conditions
  process(recoData);
}

void Occupancy::getClusters(const gsl::span<const o2::itsmft::CompClusterExt> ITSclus, float weight)
{
  for (unsigned int iClus{0}; iClus < ITSclus.size(); ++iClus) {
    auto& clus = ITSclus[iClus];
    o2::itsmft::ChipMappingITS chipMapping;
    int layer, sta, ssta, mod, chipInMod;
    chipMapping.expandChipInfoHW(chipID, layer, sta, ssta, mod, chipInMod);

    if (layer < 3) {
      fillIBmap(mOccupancyHistos[layer], sta, chipInMod, weight);
    } else {
      fillOBmap(mOccupancyHistos[layer], sta, chipInMod, weight, ssta);
    }
  }
}

void Occupancy::fillIBmap(TH2F *histo, int sta, int chipInMod, float weight)
{
    histo->Fill(chipInMod, sta, weight);
}

void Occupancy::fillOBmap(TH2F *histo,  int sta, int chipInMod, float weight, int ssta)
{
    auto xCoord = chipInMod < 7 ? (mod - 1) * 7 + chipInMod : (mod - 1) * 7 + 14 - chipInMod;
    auto yCoord = 4 * sta + ssta * 2 + 1 * (chipInMod < 7);
    histo->Fill(xCoord, yCoord, weight);
}

void Occupancy::process(o2::globaltracking::RecoContainer& recoData)
{
  auto& params = o2::its::study::ITSAvgClusSizeParamConfig::Instance();
  bool isMCTarget = false;
  PVertex pv;

  auto compClus = recoData.getITSClusters();
  getClusters(compClus, 1.); // Default weight is 1
  saveHistograms();

}

void Occupancy::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this param need to be queried only once
    initOnceDone = true;
    o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  }
}

void Occupancy::saveHistograms()
{
  mDBGOut.reset();
  TFile fout(mOutName.c_str(), "RECREATE");

  for (auto& histo : mOccupancyHistos) {
    histo->Write();
  }

  fout.Close();
}

void Occupancy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  // o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj);
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
    return;
  }
}

DataProcessorSpec getOccupancy(mask_t srcTracksMask, mask_t srcClustersMask, bool useMC, std::shared_ptr<o2::steer::MCKinematicsReader> kineReader)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestTracks(srcTracksMask, useMC);
  dataRequest->requestClusters(srcClustersMask, useMC);
  // dataRequest->requestPrimaryVerterticesTMP(useMC);

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              false,                             // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);
  return DataProcessorSpec{
    "its-study-Occupancy",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<Occupancy>(dataRequest, ggRequest, useMC, kineReader)},
    Options{}};
}
} // namespace study
} // namespace its
} // namespace o2
