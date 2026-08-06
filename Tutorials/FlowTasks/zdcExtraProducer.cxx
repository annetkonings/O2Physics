#include "ZDCExtraTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using namespace o2::aod::evsel;

struct ZDCExtraProducer {
  Produces<aod::ZDCExtra> zdcExtra;

  using ColEvSels = soa::Join<aod::Collisions,
                              aod::EvSels,
                              o2::aod::CentFT0Cs>;
  using BCsRun3 = soa::Join<aod::BCs,
                            aod::Timestamps,
                            aod::BcSels,
                            aod::Run3MatchedToBCSparse>;

  void process(ColEvSels::iterator const& collision, BCsRun3 const& /*bcs*/)
  {
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (!foundBC.has_zdc()) {
      return;
    }
    auto zdc = foundBC.zdc();

    float zna[4] = {0.f};
    float znc[4] = {0.f};
    for (int i = 0; i < 4; ++i) {
      zna[i] = zdc.energySectorZNA()[i];
      znc[i] = zdc.energySectorZNC()[i];
    }

    zdcExtra(zna[0], zna[1], zna[2], zna[3],
             znc[0], znc[1], znc[2], znc[3]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCExtraProducer>(cfgc)};
}
