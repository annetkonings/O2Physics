#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct FlowSpectatorPlane {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  int collisions = 0;

  void init(InitContext const&)
  {
    LOG(info) << "Task initialized";
    const AxisSpec axisCounter{1, 0.0f, 1.0f, "Counts"};
    histos.add("hCollisions", "Number of collisions", kTH1F, {axisCounter});
  }

  void process(aod::Collision const&)
  {
    histos.fill(HIST("hCollisions"), 0.5f);
    collisions++;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSpectatorPlane>(cfgc)};
}
