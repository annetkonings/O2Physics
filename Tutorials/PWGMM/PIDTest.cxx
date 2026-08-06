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

#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <cmath>

using namespace o2;
using namespace o2::framework;
// using namespace o2::pid;

struct PIDTest {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisTPCnSigmaPi{200, -10.0, 10.0, "n#sigma Pi"};
    const AxisSpec axisTPCnSigmaKa{200, -10.0, 10.0, "n#sigma Ka"};
    const AxisSpec axisTPCnSigmaPr{200, -10.0, 10.0, "n#sigma Pr"};

    const AxisSpec axispTPC{200, 0.0f, +10.0f, "p (GeV/c)"};
    const AxisSpec axispTPCSignal{200, 0.0f, +10.0f, "p/|z| (GeV/c)"};
    const AxisSpec axisTPCSignal{200, 0, 1000, "TPC dE/dx (a.u.)"};

    // create histograms
    histos.add("TPCnSigmaPi", "TPCnSigmaPi", kTH2F, {axispTPC, axisTPCnSigmaPi});
    histos.add("TPCnSigmaKa", "TPCnSigmaKa", kTH2F, {axispTPC, axisTPCnSigmaKa});
    histos.add("TPCnSigmaPr", "TPCnSigmaPr", kTH2F, {axispTPC, axisTPCnSigmaPr});

    histos.add("TPCSignal", "TPC Signal", kTH2F, {axispTPCSignal, axisTPCSignal});

    /*
    const AxisSpec axisTOFnSigmaPi{200, -10.0, 10.0, "n#sigma Pi"};
    const AxisSpec axisTOFnSigmaKa{200, -10.0, 10.0, "n#sigma Ka"};
    const AxisSpec axisTOFnSigmaPr{200, -10.0, 10.0, "n#sigma Pr"};

    const AxisSpec axispTOF{200, 0.0f, +10.0f, "p (GeV/c)"};
    const AxisSpec axispTOFbeta{200, 0.0f, +5.0f, "p (GeV/c)"};

    // create histograms
    histos.add("TOFnSigmaPi", "TOFnSigmaPi", kTH2F, {axispTOF, axisTOFnSigmaPi});
    histos.add("TOFnSigmaKa", "TOFnSigmaKa", kTH2F, {axispTOF, axisTOFnSigmaKa});
    histos.add("TOFnSigmaPr", "TOFnSigmaPr", kTH2F, {axispTOF, axisTOFnSigmaPr});


    histos.add("TOFbeta", "TOF #beta", kTH2F, {axispTOFbeta, {100, 0.2f, 1.1f, "#beta"}});
    */
    LOG(info) << "Histos created"; // <- always printed
  }

  void process(soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::TrackSelection> const& tracks)
  // void process(soa::Join<aod::TracksIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr> const& tracks)
  // void process(soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta> const& tracks)
  {
    for (const auto& track : tracks) {
      // if (track.tpcNClsFound() < 70)
      //   continue;

      // if (std::abs(track.eta()) > 0.8)
      //   continue;

      // if (!track.isGlobalTrack())
      //   continue;

      float tpcNSigmaPi = o2::aod::pidutils::tpcNSigma<2>(track);
      float tpcNSigmaKa = o2::aod::pidutils::tpcNSigma<3>(track);
      float tpcNSigmaPr = o2::aod::pidutils::tpcNSigma<4>(track);

      // if (std::abs(tpcNSigmaPi) < 3.0)
      histos.fill(HIST("TPCnSigmaPi"), track.p(), tpcNSigmaPi);
      // if (std::abs(tpcNSigmaKa) < 3.0)
      histos.fill(HIST("TPCnSigmaKa"), track.p(), tpcNSigmaKa);
      // if (std::abs(tpcNSigmaPr) < 3.0)
      histos.fill(HIST("TPCnSigmaPr"), track.p(), tpcNSigmaPr);

      histos.fill(HIST("TPCSignal"), track.p() / std::abs(track.sign()), track.tpcSignal());

      /*
      if (!track.hasTOF()) {
        continue;
      }

      float tofNSigmaPi = o2::aod::pidutils::tofNSigma<2>(track);
      float tofNSigmaKa = o2::aod::pidutils::tofNSigma<3>(track);
      float tofNSigmaPr = o2::aod::pidutils::tofNSigma<4>(track);

      if (std::abs(tofNSigmaPi) < 3.0)
        histos.fill(HIST("TOFnSigmaPi"), track.p(), tofNSigmaPi);
      if (std::abs(tofNSigmaKa) < 3.0)
        histos.fill(HIST("TOFnSigmaKa"), track.p(), tofNSigmaKa);
      if (std::abs(tofNSigmaPr) < 3.0)
        histos.fill(HIST("TOFnSigmaPr"), track.p(), tofNSigmaPr);
      histos.fill(HIST("TOFbeta"), track.p(), track.beta());
      */
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PIDTest>(cfgc)};
}
