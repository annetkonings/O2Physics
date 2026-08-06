#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PIDTest {

  using PiTracks = soa::Join<aod::Tracks, aod::pidTOFPi, aod::pidTPCPi>;

  void init(InitContext const&)
  {
    // LOG(info) << "Task initialized";   // <- always printed

    // define axes you want to use
    const AxisSpec axisTOFnSigmaPi{100, 0.0f, +4.0f, "n#sigma Pi"};
    const AxisSpec axisTOFnSigmaKa{100, 0.0f, +4.0f, "n#sigma Ka"};
    const AxisSpec axisTOFnSigmaPr{100, 0.0f, +4.0f, "n#sigma Pr"};

    // create histograms
    histos.add("TOFnSigmaPi", "TOFnSigmaPi", kTH1F, {axisTOFnSigmaPi});
    histos.add("TOFnSigmaKa", "TOFnSigmaKa", kTH1F, {axisTOFnSigmaKa});
    histos.add("TOFnSigmaPr", "TOFnSigmaPr", kTH1F, {axisTOFnSigmaPr});
  }

  void processPi(PiTracks const& track)
  {
    float tofNSigmaPi = track.tofNSigmaPi();
    float tpcNSigmaPi = track.tpcNSigmaPi();

    float tofNSigmaKa = track.tofNSigmaKa();
    float tpcNSigmaKa = track.tpcNSigmaKa();

    float tofNSigmaPr = track.tofNSigmaPr();
    float tpcNSigmaPr = track.tpcNSigmaPr();

    // Combined nSigma calculation
    float combNSigmaPi = std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi);

    // Individual cuts
    bool passTOFPi = std::abs(tofNSigmaPi) < 3.0;
    bool passTPCPi = std::abs(tpcNSigmaPi) < 3.0;
    bool passAllPi = passTOFPi && passTPCPi;

    bool passTOFKa = std::abs(TOFnSigmaKa) < 3.0;
    bool passTPCKa = std::abs(tpcNSigmaKa) < 3.0;
    bool passAllKa = passTOFKa && passTPCKa;

    bool passTOFPr = std::abs(TOFnSigmaPr) < 3.0;
    bool passTPCPr = std::abs(tpcNSigmaPr) < 3.0;
    bool passAllPr = passTOFPr && passTPCPr;

    if (passTOFPi) {
      histos.fill(HIST("TOFnSigmaPi"), tofNSigmaPi);
    }
    if (passTOFKa) {
      histos.fill(HIST("TOFnSigmaKa"), TOFnSigmaKa);
    }
    if (passTOFPr) {
      histos.fill(HIST("TOFnSigmaPr"), TOFnSigmaPr);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PIDTest>(cfgc)};
}
