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
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.2, "max DCAxy (in cm)"};
  Configurable<float> minTPCCrossedRows{"minTPCCrossedRows", 70, "min crossed rows in the TPC"};
  Configurable<float> maxPVz{"maxPVz", 100.0f, "max value of PVz (cm)"};

  void init(InitContext const&)
  {
    LOG(info) << "Task initialized"; // <- always printed
    // define axes you want to use
    const AxisSpec axisCounter{1, 0.0f, +1.0f, ""};
    const AxisSpec axisEta{30, -1.5f, +1.5f, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0.0f, 10.0f, "p_{T}"};
    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
  }

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {
    if (fabs(collision.posZ()) > maxPVz) {
      return; // skip outside a certain region
    }
    histos.fill(HIST("eventCounter"), 0.5f);
    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < minTPCCrossedRows)
        continue; // badly tracked
      if (fabs(track.dcaXY()) > maxDCAxy)
        continue; // doesn’t point to primary vertex
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}

/*

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct myExampleTask {

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //LOG(info) << "Task initialized";   // <- always printed

    // define axes you want to use
    const AxisSpec axiseta{100, -2.2f, +2.2f, "#eta"};
    const AxisSpec axisEvents{1, 0.0f, 1.0f, ""};

    // create histograms
    histos.add("dNdeta", "dN_{ch}/d#eta;#eta;dN_{ch}/d#eta", kTH1F, {axiseta});
    histos.add("eventCounter", "Event counter", kTH1F, {axisEvents});
  }

  void process(aod::Collision const&, aod::Tracks const& tracks)
  {
    // Count this event
    histos.fill(HIST("eventCounter"), 0.5);

    for (auto& track : tracks) {

        // ---- compute pt ----
        float pt = std::sqrt(track.px() * track.px() + track.py() * track.py());
        if (pt < 0.15) continue;

        // ---- compute eta ----
        float p = std::sqrt(track.px()*track.px() + track.py()*track.py() + track.pz()*track.pz());
        float eta = 0.5f * std::log((p + track.pz()) / (p - track.pz()));
        if (std::abs(eta) > 2.2) continue;

        // ---- compute TPC crossed rows ----
        //auto crossedRows = track.tpcNClsFindable() - track.tpcNClsFindableMinusCrossedRows();
        //if (crossedRows < 70) continue;

        // ---- fill histogram ----
        histos.fill(HIST("dNdeta"), eta);
    }
  }

};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}




#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myExampleTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.2, "max DCAxy (in cm)"};
  Configurable<float> minTPCCrossedRows{"minTPCCrossedRows", 70, "min crossed rows in the TPC"};
  Configurable<float> maxPVz{"maxPVz", 100.0f, "max value of PVz (cm)"};


  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;
  //This is an example of a convenient declaration of "using"
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(InitContext const&)
  {
    //LOG(info) << "Task initialized";   // <- always printed

    // define axes you want to use
    const AxisSpec axisCounter{1, 0.0f, +1.0f, ""};
    const AxisSpec axisEta{30, -1.5f, +1.5f, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0.0f, 10.0f, "p_{T}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    // create histograms
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt});

    histos.add("ptHistogramPion", "ptHistogramPion", kTH1F, {axisPt});
    histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1F, {axisPt});
    histos.add("ptHistogramProton", "ptHistogramProton", kTH1F, {axisPt});
    histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt});
    histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt});
    histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});

    histos.add("numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{10,-0.5f, 9.5f}});
    histos.add("multiplicityCorrelation", "multiplicityCorrelations", kTH2F, {{100, -0.5f, 99.5f}, {100,-0.5f, 99.5f}});
  }


  void processReco(aod::Collision const& collision, myFilteredTracks const& tracks, aod::McParticles const&)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    for (const auto& track : tracks) {
      if( track.tpcNClsCrossedRows() < 70 ) continue; //badly tracked
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());

      if(track.has_mcParticle()){
        auto mcParticle = track.mcParticle();
        histos.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());
        if(mcParticle.isPhysicalPrimary() && fabs(mcParticle.y())<0.5){ // do this in the context of the track ! (context matters!!!)
          if(abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptHistogramPion"), mcParticle.pt());
          if(abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptHistogramKaon"), mcParticle.pt());
          if(abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptHistogramProton"), mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(myExampleTask, processReco, "process reconstructed information", true);

  void processSim(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels,aod::Collisions>> const& collisions, aod::McParticles const& mcParticles, myFilteredTracks const& tracks)
  {
    histos.fill(HIST("numberOfRecoCollisions"), collisions.size()); // number of times coll was reco-ed

    //Now loop over each time this collision has been reconstructed and aggregate tracks
    std::vector<int> numberOfTracks;
    for (auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      // size of grouped tracks may help in understanding why event was split!
      numberOfTracks.emplace_back(groupedTracks.size());
    }
    if( collisions.size() == 2 ) histos.fill(HIST("multiplicityCorrelation"), numberOfTracks[0], numberOfTracks[1]);

    for (const auto& mcParticle : mcParticles) {
      if(mcParticle.isPhysicalPrimary() && fabs(mcParticle.y())<0.5){ // watch out for context!!!
      if(abs(mcParticle.pdgCode())==211) histos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
      if(abs(mcParticle.pdgCode())==321) histos.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
      if(abs(mcParticle.pdgCode())==2212) histos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(myExampleTask, processSim, "process pure simulation information", true);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}
*/
