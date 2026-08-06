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

// O2 includes
#include "SkimTracks.h"

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct DerivedFlowTable {
  // Histogram registry: an object to hold your histograms
  // HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // This marks that this task produces a standard derived table
  Produces<aod::DrCollisions> outputCollisions;
  Produces<aod::DrTracks> outputTracks;

  Produces<aod::DrTracksExtra> outputTracksExtra;

  // Look at primary tracks only
  // Filter trackFilter = nabs(aod::track::dcaXY) < maxDCA && nabs(aod::track::eta) < etaWindow && aod::track::pt > minPt;

  // This is an example of a convenient declaration of "using"
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>; // do not forget this!

  // void init(InitContext const&)
  //{
  //  define axes you want to use
  //}

  void process(aod::Collision const& collision, myFilteredTracks const& tracks)
  {

    outputCollisions(collision.posZ());
    for (const auto& track : tracks) {
      outputTracks(
        outputCollisions.lastIndex(), track.trackType(), track.x(), track.alpha(),
        track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
        track.isWithinBeamPipe(), track.px(), track.py(), track.pz(),
        track.sign(), track.eta(), track.phi());

      outputTracksExtra(
        outputTracks.lastIndex(), track.tpcInnerParam());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedFlowTable>(cfgc, TaskName{"derived-flow-table"})};
  return workflow;
}
