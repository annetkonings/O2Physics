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

struct FlowSkimTracks {

  // O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  // O2_DEFINE_CONFIGURABLE(cfgMpar, int, 8, "Highest order of pt-pt correlations")
  // O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  // O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  // O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  // O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
  // O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, true, "Use additional track cut on phi")
  // O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  // O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  // O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  //  added
  O2_DEFINE_CONFIGURABLE(cfgDoubleTrackFunction, bool, true, "Include track cut at low pt");
  O2_DEFINE_CONFIGURABLE(cfgTrackCutSize, float, 0.06, "Spread of track cut");
  O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 500, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  // O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  // O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");

  // Tables
  Produces<aod::DrCollisions> outputCollisions;
  Produces<aod::DrTracks> outputTracks;

  // TOF & CCDB services
  // Maybe without o2 in front
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // Service<o2::pid::tof::TOFResponse> tofResponse;
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> passName{"passName", "apass5", "Reconstruction pass for TOF"};

  // Using complete tracks including DCA, PID, etc.
  // Collisions with selection + centrality + multiplicity
  using GeneralCollisions =
    soa::Join<
      aod::Collisions,
      aod::BCs,
      aod::EvSels,
      aod::Timestamps,
      aod::Mults,
      aod::CentFT0Cs,
      aod::MultsGlobal,
      aod::Zdcs>;

  // Tracks
  using UnfilteredTracks =
    soa::Join<
      aod::Tracks,
      aod::TracksExtra,
      aod::TrackSelection,
      aod::TracksDCA,
      aod::pidTPCFullPi,
      aod::pidTOFFullPi,
      aod::pidTPCFullKa,
      aod::pidTOFFullKa,
      aod::pidTPCFullPr,
      aod::pidTOFFullPr>;

  // using UsedTracks = soa::Filtered<UnfilteredTracks>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(o2::framework::InitContext&)
  {
    LOG(debug) << "Initializing FlowSkimTracks task";

    // Connect to CCDB
    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);

    // Initialize TOF PID service
    // LOG(debug) << "Initializing the tofSignal task";
    // tofResponse->initSetup(ccdb, initContext);
  }

  void process(GeneralCollisions::iterator const& collision,
               UnfilteredTracks const& tracks)
  {
    // Write collision output
    outputCollisions(
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      collision.flags(),
      collision.timestamp(),
      collision.numContrib(),
      collision.runNumber(),
      collision.multNTracksGlobal(),
      collision.multNGlobalTracksPV(),
      collision.multNTracksPV(),
      collision.centFT0C(),
      collision.energyCommonZNA(),
      collision.energyCommonZNC());

    // Get index of the collision we just stored
    auto drCollisionIndex = outputCollisions.lastIndex();

    // Slice tracks belonging only to this collision
    auto tracksThisCollision = tracks.sliceBy(perCollision, collision.globalIndex());

    for (const auto& track : tracksThisCollision) {
      outputTracks(
        drCollisionIndex,
        track.trackType(),
        track.px(),
        track.py(),
        track.pz(),
        track.eta(),
        track.phi(),
        track.sign(),
        track.isGlobalTrack(),
        track.dcaXY(),
        track.dcaZ(),
        track.tpcSignal(),
        track.tpcNClsFindable(),
        track.tpcNClsFound(),
        track.tpcNClsCrossedRows(),
        track.tpcChi2NCl(),
        track.tpcFractionSharedCls(),
        track.itsClusterSizes(),
        track.itsNCls(),
        track.itsChi2NCl(),
        track.isGlobalTrackSDD(),
        track.flags(),
        track.tpcNSigmaPi(),
        track.tpcNSigmaKa(),
        track.tpcNSigmaPr(),
        track.tofNSigmaPi(),
        track.tofNSigmaKa(),
        track.tofNSigmaPr());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FlowSkimTracks>(cfgc, TaskName{"flow-skim-tracks"})};
  return workflow;
}
