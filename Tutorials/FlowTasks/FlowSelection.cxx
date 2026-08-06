// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General
// Public License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// O2 includes
#include "FlowSkimTracks.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/PID/PIDTOFParamService.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FlowSelection {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> passName{"passName", "apass5", "Reconstruction pass for TOF"};

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event level

  Configurable<bool> cfgkIsTriggerTVX{"cfgkIsTriggerTVX", true, "FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
  Configurable<bool> cfgkNoTimeFrameBorder{"cfgkNoTimeFrameBorder", true, "Reject collisions close to Time Frame borders"};
  Configurable<bool> cgfkNoITSROFrameBorder{"kNoITSROFrameBorder", true, "Reject events affected by the ITS ROF border"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 1e4, "Maximum occupancy of selected events"};
  Configurable<bool> cfgTVXinTRD{"cfgTVXinTRD", true, "Use kTVXinTRD (reject TRD triggered events)"};
  Configurable<bool> cfgNoSameBunchPileupCut{"cfgNoSameBunchPileupCut", true, "kNoSameBunchPileupCut"};
  Configurable<bool> cfgIsGoodZvtxFT0vsPV{"cfgIsGoodZvtxFT0vsPV", true, "kIsGoodZvtxFT0vsPV"};
  Configurable<bool> cfgNoCollInTimeRangeStandard{"cfgNoCollInTimeRangeStandard", true, "kNoCollInTimeRangeStandard"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "Selects collisions with at least one ITS-TPC track"};
  Configurable<float> cfgVtxZ{"cfgVtxZ", 10, "vertex cut (cm)"};
  Configurable<float> cfgCentrality{"cfgCentrality", 80, "Centrality [0,80] (%)"};

  // Multiplicity
  Configurable<float> pvUp_p0{"pvUp_p0", 4182.12, "PV +3sigma p0"};
  Configurable<float> pvUp_p1{"pvUp_p1", -152.459, "PV +3sigma p1"};
  Configurable<float> pvUp_p2{"pvUp_p2", 2.41955, "PV +3sigma p2"};
  Configurable<float> pvUp_p3{"pvUp_p3", -0.0199481, "PV +3sigma p3"};
  Configurable<float> pvUp_p4{"pvUp_p4", 6.93894e-05, "PV +3sigma p4"};

  Configurable<float> pvLow_p0{"pvLow_p0", 1837.75, "PV -3sigma p0"};
  Configurable<float> pvLow_p1{"pvLow_p1", -60.852, "PV -3sigma p1"};
  Configurable<float> pvLow_p2{"pvLow_p2", 0.724331, "PV -3sigma p2"};
  Configurable<float> pvLow_p3{"pvLow_p3", -0.00366975, "PV -3sigma p3"};
  Configurable<float> pvLow_p4{"pvLow_p4", 6.47562e-06, "PV -3sigma p4"};

  Configurable<float> glUp_p0{"glUp_p0", 3384.43, "Global +3sigma p0"};
  Configurable<float> glUp_p1{"glUp_p1", -118.377, "Global +3sigma p1"};
  Configurable<float> glUp_p2{"glUp_p2", 1.72823, "Global +3sigma p2"};
  Configurable<float> glUp_p3{"glUp_p3", -0.0127887, "Global +3sigma p3"};
  Configurable<float> glUp_p4{"glUp_p4", 4.03432e-05, "Global +3sigma p4"};

  Configurable<float> glLow_p0{"glLow_p0", 885.976, "Global -3sigma p0"};
  Configurable<float> glLow_p1{"glLow_p1", -26.3397, "Global -3sigma p1"};
  Configurable<float> glLow_p2{"glLow_p2", 0.240114, "Global -3sigma p2"};
  Configurable<float> glLow_p3{"glLow_p3", -0.000496168, "Global -3sigma p3"};
  Configurable<float> glLow_p4{"glLow_p4", -1.82704e-06, "Global -3sigma p4"};

  // Track level
  Configurable<float> cfgEta{"cfgEta", 0.8, "eta cut"};
  Configurable<float> cfgPtmin{"cfgPtmin", 0.2, "minimum pt (GeV/c)"};
  Configurable<float> cfgPtmax{"cfgPtmax", 10.0, "maximum pt (GeV/c)"};
  Configurable<bool> cfgIsGlobalTrackSDD{"cfgIsGlobalTrackSDD", true, "Crossed row, DCA, TPC/ITS cuts"};
  Configurable<float> cfgDCAxy{"cfgDCAxy", 0.2, "max DCAxy (in cm)"};
  Configurable<float> cfgDCAz{"cfgDCAz", 2.0, "max DCAz (in cm)"};

  float evalPoly4(float x,
                  float p0, float p1, float p2,
                  float p3, float p4)
  {
    return p0 + p1 * x + p2 * x * x + p3 * x * x * x + p4 * x * x * x * x;
  }

  // Event selection filter
  /*
    Filter VtxZFilter = nabs(aod::collisions::posZ) < cfgVtxZ;
  Filter OccupancyFilter = aod::collision::trackOccupancyInTimeRange > 0 && aod::collision::trackOccupancyInTimeRange < cfgMaxOccupancy;
  Filter noSameBunchPileupFilter = (!cfgNoSameBunchPileupCut) || (aod::evsel::selection_bit(aod::evsel::kNoSameBunchPileup) == true);
  Filter Sel8Filter = (!cfgkIsTriggerTVX || aod::evsel::selection_bit(aod::evsel::kIsTriggerTVX)) &&
  (!cfgkNoTimeFrameBorder ||
   aod::evsel::selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
  (!cgfkNoITSROFrameBorder ||
   aod::evsel::selection_bit(aod::evsel::kNoITSROFrameBorder));
  */

  // Track selection filters
  /*
  Filter PtFilter = (aod::track::pt > cfgPtmin) && (aod::track::pt < cfgPtmax);
  Filter EtaFilter = nabs(track::eta) < cfgEta;
  Filter GlobalTrackInFilter = (!cfgIsGlobalTrackSDD) || (aod::track::isGlobalTrackSDD == true);
  Filter DCAxyFilter = nabs(aod::track::dcaXY) < cfgDCAxy;
  Filter DCAzFilter  = nabs(aod::track::dcaZ)  < cfgDCAz;
  */

  // This is an example of a convenient declaration of "using"
  using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  // using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisDCAxy{100, -10.0f, +10.0f, "DCA_{XY}"};
    const AxisSpec axisDCAz{100, -10.0f, 10.0f, "DCA_{Z}"};
    // create histograms
    histos.add("hist_dcaXY_before", "dhist_dcaXY_before", kTH1F, {axisDCAxy});
    histos.add("hist_dcaXY_after", "hist_dcaXY_after", kTH1F, {axisDCAxy});
    histos.add("hist_dcaZ_before", "hist_dcaZ_before", kTH1F, {axisDCAz});
    histos.add("hist_dcaZ_after", "hist_dcaZ_after", kTH1F, {axisDCAz});

    const AxisSpec axisCent{100, 0, 100, "Centrality (%)"};
    const AxisSpec axisOcc{200, 0, 2000, "Track occupancy"};
    const AxisSpec axisPV{200, 0, 10000, "N_{ch}(PV)"};
    const AxisSpec axisGlobal{200, 0, 10000, "N_{ch}(Global)"};

    histos.add("hist_cent_before", "hist_cent_before", kTH1F, {axisCent});
    histos.add("hist_cent_after", "hist_cent_after", kTH1F, {axisCent});

    histos.add("hist_occ_before", "hist_occ_before", kTH1F, {axisOcc});
    histos.add("hist_occ_after", "hist_occ_after", kTH1F, {axisOcc});

    histos.add("hist_pv_before", "hist_pv_before", kTH1F, {axisPV});
    histos.add("hist_pv_after", "hist_pv_after", kTH1F, {axisPV});

    histos.add("hist_gl_before", "hist_gl_before", kTH1F, {axisGlobal});
    histos.add("hist_gl_after", "hist_gl_after", kTH1F, {axisGlobal});
  }

  void process(CollisionsFull const& collisions, myCompleteTracks const& tracks)
  {
    // ===============================
    // EVENT LOOP
    // ===============================
    for (auto& collision : collisions) {

      float occupancy = collision.trackOccupancyInTimeRange();
      float centrality = collision.centFT0C();
      float PVNch = collision.multNTracksPV(); // PV contributors
      // Global track multiplicity
      auto tracksThisCollision =
        tracks.sliceBy(perCollision, collision.globalIndex());

      int multNTracksGlobal = 0;

      for (auto& track : tracksThisCollision) {
        if (!track.isGlobalTrack())
          continue;
        multNTracksGlobal++;
      }

      float globalNch = multNTracksGlobal;

      // BEFORE event cuts
      histos.fill(HIST("hist_cent_before"), centrality);
      histos.fill(HIST("hist_occ_before"), occupancy);
      histos.fill(HIST("hist_pv_before"), PVNch);
      histos.fill(HIST("hist_gl_before"), globalNch);

      // PV
      float pvUpper = evalPoly4(centrality, pvUp_p0, pvUp_p1, pvUp_p2, pvUp_p3, pvUp_p4);

      float pvLower = evalPoly4(centrality, pvLow_p0, pvLow_p1, pvLow_p2, pvLow_p3, pvLow_p4);

      bool passPV = (PVNch >= pvLower && PVNch <= pvUpper);

      // Global
      float glUpper = evalPoly4(centrality, glUp_p0, glUp_p1, glUp_p2, glUp_p3, glUp_p4);

      float glLower = evalPoly4(centrality, glLow_p0, glLow_p1, glLow_p2, glLow_p3, glLow_p4);

      bool passGlobal = (globalNch >= glLower && globalNch <= glUpper);

      // Example event selection
      if (centrality > 0 && centrality < cfgCentrality) {
        histos.fill(HIST("hist_cent_after"), centrality);
      }

      if (occupancy > 0 && occupancy < cfgMaxOccupancy) {
        histos.fill(HIST("hist_occ_after"), occupancy);
      }

      if (passPV && passGlobal) {
        histos.fill(HIST("hist_pv_after"), PVNch);
        histos.fill(HIST("hist_gl_after"), globalNch);
      }
    }
    for (auto& track : tracks) {

      float dcaXY = track.dcaXY();
      float dcaZ = track.dcaZ();
      // DCAxy without cut
      histos.fill(HIST("hist_dcaXY_before"), dcaXY);
      // DCAxy with DCAxy cut
      if (std::abs(dcaXY) < cfgDCAxy) {
        histos.fill(HIST("hist_dcaXY_after"), dcaXY);
      }

      // DCAz without cut
      histos.fill(HIST("hist_dcaZ_before"), dcaZ);

      // DCAz with DCAz cut
      if (std::abs(dcaZ) < cfgDCAz) {
        histos.fill(HIST("hist_dcaZ_after"), dcaZ);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FlowSelection>(cfgc, TaskName{"flow-selection"})};
  return workflow;
}
