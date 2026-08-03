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

/// \file   flowSP.cxx
/// \author Noor Koster
/// \since  01/12/2024
/// \brief  task to evaluate flow with respect to spectator plane.

#include "PWGCF/DataModel/SPTableZDC.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowSPstatPID {

  // event selection configurable group
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cEvtUseRCTFlagChecker, bool, true, "Evt sel: use RCT flag checker");
    O2_DEFINE_CONFIGURABLE(cEvtRCTFlagCheckerLabel, std::string, "CBT_hadronPID", "Evt sel: RCT flag checker label (CBT, CBT_hadronPID)"); // all Labels can be found in Common/CCDB/RCTSelectionFlags.h
    O2_DEFINE_CONFIGURABLE(cEvtRCTFlagCheckerZDCCheck, bool, true, "Evt sel: RCT flag checker ZDC check");
    O2_DEFINE_CONFIGURABLE(cEvtRCTFlagCheckerLimitAcceptAsBad, bool, true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad");
    O2_DEFINE_CONFIGURABLE(cEvSelsUseAdditionalEventCut, bool, true, "Bool to enable Additional Event Cut");
    O2_DEFINE_CONFIGURABLE(cEvSelsMaxOccupancy, int, 10000, "Maximum occupancy of selected events");
    O2_DEFINE_CONFIGURABLE(cEvSelsMinOccupancy, int, 0, "Minimum occupancy of selected events");
    O2_DEFINE_CONFIGURABLE(cEvSelsNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
    O2_DEFINE_CONFIGURABLE(cEvSelsIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
    O2_DEFINE_CONFIGURABLE(cEvSelsNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
    O2_DEFINE_CONFIGURABLE(cEvSelsNoCollInTimeRangeNarrow, bool, false, "kNoCollInTimeRangeNarrow");
    O2_DEFINE_CONFIGURABLE(cEvSelsDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
    O2_DEFINE_CONFIGURABLE(cEvSelsIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
    O2_DEFINE_CONFIGURABLE(cEvSelsIsGoodITSLayersAll, bool, true, "Cut time intervals with dead ITS staves");
    O2_DEFINE_CONFIGURABLE(cEvSelsIsGoodITSLayer0123, bool, false, "Cut time intervals with dead ITS staves");

  // QA Plots
    O2_DEFINE_CONFIGURABLE(cFillEventQA, bool, false, "Fill histograms for event QA");
    O2_DEFINE_CONFIGURABLE(cFillTrackQA, bool, false, "Fill histograms for track QA");
    O2_DEFINE_CONFIGURABLE(cFillPIDQA, bool, true, "Fill histograms for PID QA");
    O2_DEFINE_CONFIGURABLE(cFillEventPlaneQA, bool, false, "Fill histograms for Event Plane QA");
    O2_DEFINE_CONFIGURABLE(cFillQABefore, bool, false, "Fill QA histograms before cuts, only for processData");
    O2_DEFINE_CONFIGURABLE(cFillMeanPT, bool, false, "Fill histograms for mean PX/PT");
    O2_DEFINE_CONFIGURABLE(cFillMeanPTextra, bool, false, "Fill histograms for mean PX/PT extra");
    O2_DEFINE_CONFIGURABLE(cUseCentAveragePt, bool, false, "Use <pt> in 1% centrality intervals and not cent average");
    O2_DEFINE_CONFIGURABLE(cFillWithMCParticle, bool, false, "Fill histograms with MCParticle instead of Track");
  // Flags to make and fill histograms
    O2_DEFINE_CONFIGURABLE(cFillGeneralV1Histos, bool, true, "Fill histograms for vn analysis");
    O2_DEFINE_CONFIGURABLE(cFillMixedHarmonics, bool, true, "Flag to make and fill histos for mixed harmonics");
    O2_DEFINE_CONFIGURABLE(cFillEventPlane, bool, false, "Flag to make and fill histos with Event Plane");
    O2_DEFINE_CONFIGURABLE(cFillXandYterms, bool, false, "Flag to make and fill histos for with separate x and y terms for SPM");
    O2_DEFINE_CONFIGURABLE(cFillChargeDependence, bool, true, "Flag to make and fill histos for charge dependent flow");
    O2_DEFINE_CONFIGURABLE(cFillChargeDependenceQA, bool, true, "Flag to make and fill QA histos for charge dependent flow");
    O2_DEFINE_CONFIGURABLE(cFillPID, bool, true, "Flag to make and fill histos for PID flow");
  // Centrality Estimators -> standard is FT0C
    O2_DEFINE_CONFIGURABLE(cCentFT0Cvariant1, bool, false, "Set centrality estimator to CentFT0Cvariant1");
    O2_DEFINE_CONFIGURABLE(cCentFT0M, bool, false, "Set centrality estimator to CentFT0M");
    O2_DEFINE_CONFIGURABLE(cCentFV0A, bool, false, "Set centrality estimator to CentFV0A");
    O2_DEFINE_CONFIGURABLE(cCentNGlobal, bool, false, "Set centrality estimator to CentNGlobal");
  // Standard selections
    O2_DEFINE_CONFIGURABLE(cTrackSelsDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
    O2_DEFINE_CONFIGURABLE(cTrackSelsDCAz, float, 0.2, "Cut on DCA in the longitudinal direction (cm)");
    O2_DEFINE_CONFIGURABLE(cTrackSelsNcls, float, 70, "Cut on number of TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cTrackSelsFshcls, float, 0.4, "Cut on fraction of shared TPC clusters found");
    O2_DEFINE_CONFIGURABLE(cTrackSelsPtmin, float, 0.2, "minimum pt (GeV/c)");
    O2_DEFINE_CONFIGURABLE(cTrackSelsPtmax, float, 10, "maximum pt (GeV/c)");
    O2_DEFINE_CONFIGURABLE(cTrackSelsEta, float, 0.8, "eta cut");
    O2_DEFINE_CONFIGURABLE(cIsMCReco, bool, false, "Is MC Reco");
    O2_DEFINE_CONFIGURABLE(cEvSelsVtxZ, float, 10, "vertex cut (cm)");
    O2_DEFINE_CONFIGURABLE(cMagField, float, 99999, "Configurable magnetic field;default CCDB will be queried");
    O2_DEFINE_CONFIGURABLE(cCentMin, float, 0, "Minimum cenrality for selected events");
    O2_DEFINE_CONFIGURABLE(cCentMax, float, 90, "Maximum cenrality for selected events");
    O2_DEFINE_CONFIGURABLE(cFilterLeptons, bool, true, "Filter out leptons from MCGenerated by requiring |pdgCode| > 100");
  // NUA and NUE weights
    O2_DEFINE_CONFIGURABLE(cFillWeights, bool, true, "Fill NUA weights");
    O2_DEFINE_CONFIGURABLE(cFillWeightsPOS, bool, true, "Fill NUA weights only for positive charges");
    O2_DEFINE_CONFIGURABLE(cFillWeightsNEG, bool, true, "Fill NUA weights only for negative charges");
    O2_DEFINE_CONFIGURABLE(cUseNUA1D, bool, true, "Use 1D NUA weights (only phi)");
    O2_DEFINE_CONFIGURABLE(cUseNUA2D, bool, false, "Use 2D NUA weights (phi and eta)");
    O2_DEFINE_CONFIGURABLE(cUseNUE2D, bool, false, "Use 2D NUE weights");
    O2_DEFINE_CONFIGURABLE(cUseNUE3D, bool, true, "Use 3D NUE weights (pt, eta, centrality)");
    O2_DEFINE_CONFIGURABLE(cUseNUE2Deta, bool, false, "Use 2D NUE weights TRUE: (pt and eta) FALSE: (pt and centrality)");
  // Additional track Selections
    O2_DEFINE_CONFIGURABLE(cTrackSelsUseAdditionalTrackCut, bool, false, "Bool to enable Additional Track Cut");
    O2_DEFINE_CONFIGURABLE(cTrackSelsDoDCApt, bool, true, "Apply Pt dependent DCAz cut");
    O2_DEFINE_CONFIGURABLE(cTrackSelsDCApt1, float, 0.1, "DcaZ < const + (a * b) / pt^1.1 -> this sets a");
    O2_DEFINE_CONFIGURABLE(cTrackSelsDCApt2, float, 0.035, "DcaZ < const + (a * b) / pt^1.1 -> this sets b");
    O2_DEFINE_CONFIGURABLE(cTrackSelsDCAptConsMin, float, 0.1, "DcaZ < const + (a * b) / pt^1.1 -> this sets const");
    O2_DEFINE_CONFIGURABLE(cTrackSelsPIDNsigma, float, 2.0, "nSigma cut for PID");
    O2_DEFINE_CONFIGURABLE(cTrackSelDoTrackQAvsCent, bool, true, "Do track selection QA plots as function of centrality");
  // harmonics for v coefficients
    O2_DEFINE_CONFIGURABLE(cHarm, int, 1, "Flow harmonic n for ux and uy: (Cos(n*phi), Sin(n*phi))");
    O2_DEFINE_CONFIGURABLE(cHarmMixed, int, 2, "Flow harmonic n for ux and uy in mixed harmonics (MH): (Cos(n*phi), Sin(n*phi))");
  // settings for CCDB data
    O2_DEFINE_CONFIGURABLE(cCCDBdir_QQ, std::string, "Users/c/ckoster/ZDC/LHC23_PbPb_pass5/meanQQ/Default", "ccdb dir for average QQ values in 1% centrality bins");
    O2_DEFINE_CONFIGURABLE(cCCDBdir_SP, std::string, "", "ccdb dir for average event plane resolution in 1% centrality bins");
    O2_DEFINE_CONFIGURABLE(cCCDB_NUA, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/Default", "ccdb dir for NUA corrections");
    O2_DEFINE_CONFIGURABLE(cCCDB_NUE, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/NUE/Default", "ccdb dir for NUE corrections (pt)");
    O2_DEFINE_CONFIGURABLE(cCCDB_NUE2D, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/NUE/2D", "ccdb dir for NUE 2D corrections (pt, eta)");
    O2_DEFINE_CONFIGURABLE(cCCDB_NUE3D, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/NUE3D/Default", "ccdb dir for NUE 3D corrections (pt, eta, centrality)");
    O2_DEFINE_CONFIGURABLE(cCCDBdir_centrality, std::string, "", "ccdb dir for Centrality corrections");
    O2_DEFINE_CONFIGURABLE(cCCDBdir_meanPt, std::string, "", "ccdb dir for Mean Pt corrections");

    //Choose either standard PID cut or Bayesian PID cut
    O2_DEFINE_CONFIGURABLE(cUseBayesianPID, bool, false, "Use Bayesian PID instead of standard PID");
    O2_DEFINE_CONFIGURABLE(cUseStatisticalPID, bool, true, "Use Statistical PID instead of standard PID");
    

  // Confogirable axis
  // ConfigurableAxis axisCentrality{"axisCentrality", {20, 0, 100}, "Centrality bins for vn "};
  // ConfigurableAxis axisMomentum{"axisMomentum", {20, 0, 10}, "Momentum bins for vn"};
  //ConfigurableAxis axisNch = {"axisNch", {400, 0, 4000}, "Global N_{ch}"};
  ConfigurableAxis axisMultpv = {"axisMultpv", {400, 0, 4000}, "N_{ch} (PV)"};
  // ConfigurableAxis axisEtaVn{"axisEtaVn", {8, -0.8, 0.8}, "Eta bins for vn"};

  // Added myself
  //ConfigurableAxis axisCent100 = {"axisCent100", {100, 0, 100}, "Centrality (%)"};
  ConfigurableAxis axisOccupancy = {"axisOccupancy", {400, 0, 10000}, "Occupancy"};
  //ConfigurableAxis axisVzcfg = {"axisVzcfg", {40, -15, 15}, "v_{z}"};

  // Configurables containing vector
    O2_DEFINE_CONFIGURABLE(cUsePredeFinedSigma, bool, true, "Use one of the pre-defines settings for the multiplicity vs. centrality plots");
    O2_DEFINE_CONFIGURABLE(cUsePredeFinedSigmaYear, int, 2023, "Predifine what year you want to use (2023/2024)");
    O2_DEFINE_CONFIGURABLE(cUsePredeFinedSigmaNsigma, int, 2, "Sigma used for cuts (1,2,3)");
    Configurable<std::vector<double>> cEvSelsMultPv{"cEvSelsMultPv", std::vector<double>{2223.49, -75.1444, 0.963572, -0.00570399, 1.34877e-05, 3790.99, -137.064, 2.13044, -0.017122, 5.82834e-05}, "Multiplicity cuts (PV) first 5 parameters cutLOW last 5 cutHIGH (Default is +-2sigma pass5) "};
    Configurable<std::vector<double>> cEvSelsMult{"cEvSelsMult", std::vector<double>{1301.56, -41.4615, 0.478224, -0.00239449, 4.46966e-06, 2967.6, -102.927, 1.47488, -0.0106534, 3.28622e-05}, "Multiplicity cuts (Global) first 5 parameters cutLOW last 5 cutHIGH (Default is +-2sigma pass5) "};
    Configurable<std::vector<double>> cPtBinning{"cPtBinning", std::vector<double>{0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pT binning for vn"};

  } cfg;

  RCTFlagsChecker rctChecker;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfg.cEvSelsVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfg.cTrackSelsEta && aod::track::pt > cfg.cTrackSelsPtmin&& aod::track::pt < cfg.cTrackSelsPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true) || cfg.cIsMCReco) && nabs(aod::track::dcaXY) < cfg.cTrackSelsDCAxy&& nabs(aod::track::dcaZ) < cfg.cTrackSelsDCAz;
  //Filter trackFilterMC = nabs(aod::mcparticle::eta) < cfg.cTrackSelsEta && aod::mcparticle::pt > cfg.cTrackSelsPtmin&& aod::mcparticle::pt < cfg.cTrackSelsPtmax;
  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNGlobals>;
  using UnfilteredTracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using UnfilteredTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;

  using UsedTracks = soa::Filtered<UnfilteredTracks>;
  using UsedTracksPID = soa::Filtered<UnfilteredTracksPID>;
  using ZDCCollisions = soa::Filtered<soa::Join<GeneralCollisions, aod::SPTableZDC>>;


  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // struct to hold the correction histos/
  struct Config {
    std::vector<TH1D*> mEfficiency = {};
    std::vector<TH2D*> mEfficiency2D = {};
    std::vector<TH3D*> mEfficiency3D = {};
    std::vector<GFWWeights*> mAcceptance = {};
    std::vector<TH3D*> mAcceptance2D = {};
    bool correctionsLoaded = false;
    int lastRunNumber = 0;

    TProfile* hcorrQQ = nullptr;
    TProfile* hcorrQQx = nullptr;
    TProfile* hcorrQQy = nullptr;
    TProfile* hEvPlaneRes = nullptr;
    TH1D* hCentrality = nullptr;
    TProfile2D* hMeanPt = nullptr;

    bool clQQ = false;
    bool clEvPlaneRes = false;
    bool clCentrality = false;
    bool clMeanPt = false;

  } conf;

  struct SPMvars {
    std::vector<std::map<int, float>> wacc = {{{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}, {{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}, {{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}}; // int for part species, float for weight vector for kIncl, kPos, kNeg
    std::vector<std::map<int, float>> weff = {{{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}, {{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}, {{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}}}; // int for part species, float for weight vector for kIncl, kPos, kNeg
    double centWeight = 1.0;
    double meanPtWeight = 1.0;
    double ux = 0;
    double uy = 0;
    double uxMH = 0;
    double uyMH = 0;
    double qxA = 0;
    double qyA = 0;
    double qxC = 0;
    double qyC = 0;
    double corrQQx = 1;
    double corrQQy = 1;
    double corrQQ = 1;
    double vnA = 0;
    double vnC = 0;
    double vnFull = 0;
    float centrality = 0;
    float vtxz = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    int charge = 0;
    double relPt = 1.;
    double psiA = 0;
    double psiC = 0;
    double psiFull = 0;
    double trackPxA = 0;
    double trackPxC = 0;
    double meanPxA = 0;
    double meanPxC = 0;
  } spm;

  struct PtMaps {
    std::unique_ptr<TProfile> meanPTMap = std::make_unique<TProfile>("meanPTMap", "meanPTMap", 8, -0.8, 0.8);
    std::unique_ptr<TProfile> meanPTMapPos = std::make_unique<TProfile>("meanPTMapPos", "meanPTMapPos", 8, -0.8, 0.8);
    std::unique_ptr<TProfile> meanPTMapNeg = std::make_unique<TProfile>("meanPTMapNeg", "meanPTMapNeg", 8, -0.8, 0.8);

    std::unique_ptr<TProfile> relPxA = std::make_unique<TProfile>("relPxA", "relPxA", 8, -0.8, 0.8);
    std::unique_ptr<TProfile> relPxC = std::make_unique<TProfile>("relPxC", "relPxC", 8, -0.8, 0.8);

    std::unique_ptr<TProfile> relPxANeg = std::make_unique<TProfile>("relPxANeg", "relPxANeg", 8, -0.8, 0.8);
    std::unique_ptr<TProfile> relPxAPos = std::make_unique<TProfile>("relPxAPos", "relPxAPos", 8, -0.8, 0.8);

    std::unique_ptr<TProfile> relPxCNeg = std::make_unique<TProfile>("relPxCNeg", "relPxCNeg", 8, -0.8, 0.8);
    std::unique_ptr<TProfile> relPxCPos = std::make_unique<TProfile>("relPxCPos", "relPxCPos", 8, -0.8, 0.8);
  } ptmaps;

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  OutputObj<GFWWeights> fWeightsPOS{GFWWeights("weights_positive")};
  OutputObj<GFWWeights> fWeightsNEG{GFWWeights("weights_negative")};

  HistogramRegistry registry{"registry"};
  HistogramRegistry histos{"QAhistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Event selection cuts
  std::unique_ptr<TF1> fPhiCutLow = nullptr;
  std::unique_ptr<TF1> fPhiCutHigh = nullptr;
  std::unique_ptr<TF1> fMultPVCutLow = nullptr;
  std::unique_ptr<TF1> fMultPVCutHigh = nullptr;
  std::unique_ptr<TF1> fMultCutLow = nullptr;
  std::unique_ptr<TF1> fMultCutHigh = nullptr;
  std::unique_ptr<TF1> fMultMultPVCut = nullptr;

  enum SelectionCriteria {
    evSel_FilteredEvent,
    evSel_sel8,
    evSel_RCTFlagsZDC,
    evSel_occupancy,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kNoCollInTimeRangeNarrow,
    evSel_kIsVertexITSTPC,
    evSel_kIsGoodITSLayersAll,
    evSel_kIsGoodITSLayer0123,
    evSel_MultCuts,
    evSel_isSelectedZDC,
    evSel_CentCuts,
    nEventSelections
  };

  enum TrackSelections {
    trackSel_ZeroCharge,
    trackSel_Eta,
    trackSel_Pt,
    trackSel_DCAxy,
    trackSel_DCAz,
    trackSel_GlobalTracks,
    trackSel_NCls,
    trackSel_FshCls,
    trackSel_TPCBoundary,
    trackSel_ParticleWeights,
    nTrackSelections
  };

  enum ChargeType {
    kInclusive,
    kPositive,
    kNegative,
    nChargeTypes
  };

  enum FillType {
    kBefore,
    kAfter,
    nFillTypes
  };

  enum ModeType {
    kGen,
    kReco,
    nModeTypes
  };

  enum ParticleType {
    kUnidentified,
    kPions,
    kKaons,
    kProtons,
    nParticleTypes
  };

  enum PIDType {
    kTPC,
    kTOF,
    nPIDTypes
  };

  static constexpr std::string_view Charge[] = {"incl/", "pos/", "neg/"};
  //static constexpr std::string_view Species[] = {"", "pion/", "kaon/", "proton/"};
  static constexpr std::string_view PIDObservable[] = {"TPC/", "TOF/"};
  static constexpr std::string_view Time[] = {"before/", "after/"};

  void init(InitContext const&)
  {

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec axisDCAz = {100, -.5, .5, "DCA_{z} (cm)"};
    AxisSpec axisDCAxy = {100, -.5, .5, "DCA_{xy} (cm)"};
    AxisSpec axisPhiMod = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec axisPhi = {60, 0, constants::math::TwoPI, "#varphi"};
    AxisSpec axisEta = {64, -1.6, 1.6, "#eta"};
    AxisSpec axisEtaVn = {8, -.8, .8, "#eta"};
    AxisSpec axisCentrality = {20, 0, 100, "Centrality (%)"};
    AxisSpec axisVx = {40, -0.01, 0.01, "v_{x}"};
    AxisSpec axisVy = {40, -0.01, 0.01, "v_{y}"};
    AxisSpec axisVz = {40, -10, 10, "v_{z}"};
    AxisSpec axisCent = {90, 0, 90, "Centrality(%)"};
    AxisSpec axisPhiPlane = {100, -constants::math::PI, constants::math::PI, "#Psi"};
    AxisSpec axisT0c = {70, 0, 100000, "N_{ch} (T0C)"};
    AxisSpec axisT0a = {70, 0, 200000, "N_{ch} (T0A)"};
    AxisSpec axisV0a = {70, 0, 200000, "N_{ch} (V0A)"};
    AxisSpec axisShCl = {40, 0, 1, "Fraction shared cl. TPC"};
    AxisSpec axisCl = {80, 0, 160, "Number of cl. TPC"};
    AxisSpec axisNsigma = {100, -10, 10, "Nsigma for TPC and TOF"};
    AxisSpec axisdEdx = {300, 0, 300, "dEdx for PID"};
    AxisSpec axisBeta = {150, 0, 1.5, "Beta for PID"};
    AxisSpec axisCharge = {3, 0, 3, "Charge: 0 = inclusive, 1 = positive, 2 = negative"};
    AxisSpec axisPx = {100, -0.05, 0.05, "p_{x} (GeV/c)"};
    AxisSpec axisNch = {400, 0, 4000, "Global N_{ch}"};
    AxisSpec axisMultpv = {400, 0, 4000, "N_{ch} (PV)"};
    AxisSpec axisOccupancy = {400, 0, 10000, "Occupancy"};
    AxisSpec axisPt = {cfg.cPtBinning, "#it{p}_{T} GeV/#it{c}"};
    int ptbins = cfg.cPtBinning->size() - 1;

    rctChecker.init(cfg.cEvtRCTFlagCheckerLabel, cfg.cEvtRCTFlagCheckerZDCCheck, cfg.cEvtRCTFlagCheckerLimitAcceptAsBad);

    histos.add("hCentrality", "Centrality; Centrality (%); ", {HistType::kTH1D, {axisCent}});
    histos.add("hCentralityCW", "Centrality; Centrality (%); Weighted counts", {HistType::kTH1D, {axisCent}});

    histos.add("hEventCount", "Number of Event; Cut; #Events Passed Cut", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_FilteredEvent + 1, "Filtered events");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_RCTFlagsZDC + 1, "RCT Flags ZDC");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeNarrow + 1, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_CentCuts + 1, "Cenrality range");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kkIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayer0123 + 1, "kkIsGoodITSLayer0123");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_MultCuts + 1, "Multiplicity Cuts Pilup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_isSelectedZDC + 1, "isSelected");

    histos.add("hTrackCount", "Number of Tracks; Cut; #Tracks Passed Cut", {HistType::kTH1D, {{nTrackSelections, 0, nTrackSelections}}});
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_Eta + 1, "Eta");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_Pt + 1, "Pt");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_DCAxy + 1, "DCAxy");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_DCAz + 1, "DCAz");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_GlobalTracks + 1, "GlobalTracks");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_NCls + 1, "nClusters TPC");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_FshCls + 1, "Frac. sh. Cls TPC");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_TPCBoundary + 1, "TPC Boundary");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_ZeroCharge + 1, "Only charged");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_ParticleWeights + 1, "Apply weights");

    if (cfg.cFillWeights) {
      registry.add<TH3>("weights2D/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
      registry.add<TH3>("weights2D/hPhi_Eta_vz_positive", "", kTH3D, {axisPhi, axisEta, axisVz});
      registry.add<TH3>("weights2D/hPhi_Eta_vz_negative", "", kTH3D, {axisPhi, axisEta, axisVz});

      // define output objects
      fWeights->setPtBins(ptbins, &cfg.cPtBinning.value[0]);
      fWeights->init(true, false);

      fWeightsPOS->setPtBins(ptbins, &cfg.cPtBinning.value[0]);
      fWeightsPOS->init(true, false);

      fWeightsNEG->setPtBins(ptbins, &cfg.cPtBinning.value[0]);
      fWeightsNEG->init(true, false);
    }

    if (cfg.cFillEventQA) {
      histos.add("QA/after/hCentFT0C", " ; Cent FT0C (%); ", {HistType::kTH1D, {axisCent}});

      if (cfg.cFillQABefore) {
        histos.addClone("QA/after/", "QA/before/");
      }
    }

    if (doprocessStatisticalPID) {

      if (cfg.cFillTrackQA) {
        histos.add("incl/QA/after/pt_phi", "", {HistType::kTH2D, {axisPt, axisPhiMod}});
        histos.add<TH3>("incl/QA/after/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH3>("incl/QA/after/hPhi_Eta_vz_corrected", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH2>("incl/QA/after/hDCAxy_pt", "", kTH2D, {axisPt, axisDCAxy});
        histos.add<TH2>("incl/QA/after/hDCAz_pt", "", kTH2D, {axisPt, axisDCAz});
        histos.add("incl/QA/after/hSharedClusters_pt", "", {HistType::kTH2D, {axisPt, axisShCl}});
        histos.add("incl/QA/after/hCrossedRows_pt", "", {HistType::kTH2D, {axisPt, axisCl}});
        histos.add("incl/QA/after/hCrossedRows_vs_SharedClusters", "", {HistType::kTH2D, {axisCl, axisShCl}});
        histos.add("incl/QA/after/hMeanPtEta", "", {HistType::kTProfile2D, {axisEta, axisCent}});

        if (cfg.cTrackSelDoTrackQAvsCent) {
          histos.add<TH3>("incl/QA/after/hPt_Eta", "", kTH3D, {axisPt, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPt_Eta_uncorrected", "", kTH3D, {axisPt, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPhi_Eta", "", kTH3D, {axisPhi, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPhi_Eta_uncorrected", "", kTH3D, {axisPhi, axisEta, axisCent});
        } else {
          histos.add<TH3>("incl/QA/after/hPhi_Eta_Pt", "", kTH3D, {axisPhi, axisEta, axisPt});
          histos.add<TH3>("incl/QA/after/hPhi_Eta_Pt_corrected", "", kTH3D, {axisPhi, axisEta, axisPt});
        }

        // Charge dependence done later

        histos.add("incl/QA/after/hdEdxTPC_pt", "", {HistType::kTH2D, {axisPt, axisdEdx}});
        histos.add("incl/QA/after/hBetaTOF_pt", "", {HistType::kTH2D, {axisPt, axisBeta}});
        histos.add("incl/QA/after/hNsigmaTPC_pt", "", {HistType::kTH2D, {axisPt, axisNsigma}});
        histos.add("incl/QA/after/hNsigmaTOF_pt", "", {HistType::kTH2D, {axisPt, axisNsigma}});


        if (cfg.cFillQABefore) {
          histos.addClone("incl/QA/after/", "incl/QA/before/");
        }

        //Charge dependence done later


      if (cfg.cFillEventQA) {
        histos.add("QA/hCentFull", " ; Centrality (%); ", {HistType::kTH1D, {axisCent}});
      }

        if (cfg.cFillGeneralV1Histos) {
          registry.add<TProfile3D>("incl/TPC/vnCodd", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnAodd", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnC", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnA", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnCSetPlane", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnASetPlane", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnOdd", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnEven", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          //Added
          registry.add<TProfile3D>("incl/TPC/vnAvnC", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnA_nw", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnC_nw", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
        }
        
        if (cfg.cFillEventPlane) {
          registry.add<TProfile3D>("incl/TPC/vnA_EP", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnC_EP", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnFull_EP", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
        }
        if (cfg.cFillXandYterms) {
          registry.add<TProfile3D>("incl/TPC/vnAx", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnAy", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnCx", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/vnCy", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
        }
        if (cfg.cFillMixedHarmonics) {
          registry.add<TProfile3D>("incl/TPC/MH/vnAxCxUx_MH", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/MH/vnAyCyUx_MH", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/MH/vnAxCyUy_MH", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
          registry.add<TProfile3D>("incl/TPC/MH/vnAyCxUy_MH", "", kTProfile3D, {axisPt, axisEtaVn, axisCentrality});
        }
        registry.addClone("incl/TPC/", "incl/TOF/");
      }

      if (cfg.cFillChargeDependence || cfg.cFillChargeDependenceQA) {
        LOGF(info, "Cloning charge dependence histograms");
        registry.addClone("incl/", "pos/");
        registry.addClone("incl/", "neg/");

        histos.addClone("incl/", "pos/");
        histos.addClone("incl/", "neg/");
      }

      // Statistical PID histograms

    histos.add<THnSparse>("StatPID/hTPCYieldNeg","TPC distribution;pT;eta;centrality;TPC",kTHnSparseF,{axisPt,axisEta,axisCent,axisdEdx});
    histos.add<THnSparse>("StatPID/hTPCSPOddNeg","SP weighted TPC;pT;eta;centrality;TPC",kTHnSparseF,{axisPt,axisEta,axisCent,axisdEdx});
    histos.add<THnSparse>("StatPID/hTPCYieldPos","TPC distribution;pT;eta;centrality;TPC",kTHnSparseF,{axisPt,axisEta,axisCent,axisdEdx});
    histos.add<THnSparse>("StatPID/hTPCSPOddPos","SP weighted TPC;pT;eta;centrality;TPC",kTHnSparseF,{axisPt,axisEta,axisCent,axisdEdx});

    // Basic QA for statistical PID
    histos.add("StatPID/hEventCount","Event counter",{HistType::kTH1D, {axisEvent}});
    histos.add("StatPID/hTrackCount","Track counter",{HistType::kTH1D, {axisTrack}});
    }

    if (cfg.cEvSelsUseAdditionalEventCut) {
      fMultPVCutLow = std::make_unique<TF1>("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultPVCutHigh = std::make_unique<TF1>("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultCutLow = std::make_unique<TF1>("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultCutHigh = std::make_unique<TF1>("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);

      std::vector<double> paramsMultPVCut;
      std::vector<double> paramsMultCut;
      int y2023 = 2023;
      int y2024 = 2024;
      std::array<int, 3> nSigma = {1, 2, 3};

      if (cfg.cUsePredeFinedSigma) {
        if (cfg.cUsePredeFinedSigmaYear == y2023) {
          if (cfg.cUsePredeFinedSigmaNsigma == nSigma[0]) {
            paramsMultPVCut = {2615.47, -90.5747, 1.25125, -0.00847075, 2.41183e-05, 3399.72, -121.652, 1.84077, -0.0142886, 4.71449e-05};
            paramsMultCut = {1716.84, -56.5663, 0.715202, -0.00426007, 1.05075e-05, 2550.82, -87.4873, 1.22205, -0.00852644, 2.54248e-05};
          } else if (cfg.cUsePredeFinedSigmaNsigma == nSigma[1]) {
            paramsMultPVCut = {2223.49, -75.1444, 0.963572, -0.00570399, 1.34877e-05, 3790.99, -137.064, 2.13044, -0.017122, 5.82834e-05};
            paramsMultCut = {1301.56, -41.4615, 0.478224, -0.00239449, 4.46966e-06, 2967.6, -102.927, 1.47488, -0.0106534, 3.28622e-05};
          } else if (cfg.cUsePredeFinedSigmaNsigma == nSigma[2]) {
            paramsMultPVCut = {1837.75, -60.852, 0.724331, -0.00366975, 6.47562e-06, 4182.12, -152.459, 2.41955, -0.0199481, 6.93894e-05};
            paramsMultCut = {885.976, -26.3397, 0.240114, -0.000496168, -1.82704e-06, 3384.43, -118.377, 1.72823, -0.0127887, 4.03432e-05};
          } else {
            LOGF(fatal, "nSigma can only be 1-3 please reset the variable or give the parameters manually and set cfg.cUsePredeFinedSigma to FALSE");
          }
        } else if (cfg.cUsePredeFinedSigmaYear == y2024) {
          if (cfg.cUsePredeFinedSigmaNsigma == nSigma[0]) {
            paramsMultPVCut = {2726.93, -100.128, 1.45046, -0.0099354, 2.71182e-05, 3404.72, -126.569, 1.92500, -0.0142653, 4.31645e-05};
            paramsMultCut = {1858.77, -66.6070, 0.929146, -0.00606961, 1.57639e-05, 2672.43, -96.7708, 1.39109, -0.00942498, 2.54268e-05};
          } else if (cfg.cUsePredeFinedSigmaNsigma == nSigma[1]) {
            paramsMultPVCut = {2390.04, -87.3154, 1.23176, -0.00806869, 2.06624e-05, 3744.26, -139.927, 2.16863, -0.0165329, 5.17269e-05};
            paramsMultCut = {1451.23, -51.4314, 0.694609, -0.00433959, 1.06698e-05, 3080.42, -112.071, 1.63166, -0.0112533, 3.10348e-05};
          } else if (cfg.cUsePredeFinedSigmaNsigma == nSigma[2]) {
            paramsMultPVCut = {2053.64, -74.5950, 1.01563, -0.00621473, 1.41276e-05, 4083.79, -153.304, 2.41333, -0.0188198, 6.03974e-05};
            paramsMultCut = {1042.50, -35.9374, 0.440681, -0.00222218, 3.20643e-06, 3488.53, -127.396, 1.87339, -0.0131007, 3.67434e-05};
          } else {
            LOGF(fatal, "cUsePredeFinedSigmaNsigma can only be 1-3 please reset the variable or give the parameters manually and set cUsePredeFinedSigma to FALSE");
          }
        } else {
          LOGF(fatal, "cUsePredeFinedSigmaYear can only be 2023/2024 please reset the variable or give the parameters manually and set cUsePredeFinedSigma to FALSE");
        }
      } else {
        paramsMultPVCut = cfg.cEvSelsMultPv;
        paramsMultCut = cfg.cEvSelsMult;
      }

      // number of parameters required in cfg.cEvSelsMultPv and cfg.cEvSelsMult.  (5 Low + 5 High)
      uint64_t nParams = 10;

      if (paramsMultPVCut.size() < nParams) {
        LOGF(fatal, "cfg.cEvSelsMultPv not set properly.. size = %d (should be 10) --> Check your config files!", paramsMultPVCut.size());
      } else if (paramsMultCut.size() < nParams) {
        LOGF(fatal, "cfg.cEvSelsMult not set properly.. size = %d (should be 10) --> Check your config files!", paramsMultCut.size());
      } else {
        fMultPVCutLow->SetParameters(paramsMultPVCut[0], paramsMultPVCut[1], paramsMultPVCut[2], paramsMultPVCut[3], paramsMultPVCut[4]);
        fMultPVCutHigh->SetParameters(paramsMultPVCut[5], paramsMultPVCut[6], paramsMultPVCut[7], paramsMultPVCut[8], paramsMultPVCut[9]);
        fMultCutLow->SetParameters(paramsMultCut[0], paramsMultCut[1], paramsMultCut[2], paramsMultCut[3], paramsMultCut[4]);
        fMultCutHigh->SetParameters(paramsMultCut[5], paramsMultCut[6], paramsMultCut[7], paramsMultCut[8], paramsMultCut[9]);
      }
    }

    if (cfg.cTrackSelsUseAdditionalTrackCut) {
      fPhiCutLow = std::make_unique<TF1>("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = std::make_unique<TF1>("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
  } // end of init

  float getNUA2D(TH3D* hNUA, float eta, float phi, float vtxz)
  {
    int xind = hNUA->GetXaxis()->FindBin(phi);
    int etaind = hNUA->GetYaxis()->FindBin(eta);
    int vzind = hNUA->GetZaxis()->FindBin(vtxz);
    float weight = hNUA->GetBinContent(xind, etaind, vzind);
    if (weight != 0)
      return 1. / weight;
    return 1;
  }

  template <typename TrackObject>
  PIDType getTrackPIDsignal(TrackObject track)
  {
    float TOFsignal = -1.f;
    float TPCsignal = -1.f;

    if track.hasTOF() {
      TOFsignal = std::abs(track.tofSignal());
    }
    if track.hasTPC() {
      TPCsignal = std::abs(track.tpcSignal());
    }
    
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  std::pair<float, uint16_t> getCrossingAngleCCDB(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    auto grpo = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object for Crossing Angle not found for timestamp %llu", timestamp);
      return {0, 0};
    }
    float crossingAngle = grpo->getCrossingAngle();
    uint16_t crossingAngleTime = grpo->getCrossingAngleTime();
    return {crossingAngle, crossingAngleTime};
  }

  // From Generic Framework
  void loadCorrections(uint64_t timestamp)
  {
    // corrections saved on CCDB as TList {incl, pos, neg} of GFWWeights (acc) TH1D (eff) objects!
    if (conf.correctionsLoaded)
      return;

    int nWeights = 3;

    if (cfg.cUseNUA1D) {
      if (cfg.cCCDB_NUA.value.empty() == false) {
        TList* listCorrections = ccdb->getForTimeStamp<TList>(cfg.cCCDB_NUA, timestamp);
        conf.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights")));
        conf.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_positive")));
        conf.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_negative")));
        int sizeAcc = conf.mAcceptance.size();
        if (sizeAcc < nWeights)
          LOGF(fatal, "Could not load acceptance weights from %s", cfg.cCCDB_NUA.value.c_str());
        else
          LOGF(info, "Loaded acceptance weights from %s", cfg.cCCDB_NUA.value.c_str());
      } else {
        LOGF(info, "cfg.cCCDB_NUA empty! No corrections loaded");
      }
    } else if (cfg.cUseNUA2D) {
      if (cfg.cCCDB_NUA.value.empty() == false) {
        TH3D* hNUA2D = ccdb->getForTimeStamp<TH3D>(cfg.cCCDB_NUA, timestamp);
        if (!hNUA2D) {
          LOGF(fatal, "Could not load acceptance weights from %s", cfg.cCCDB_NUA.value.c_str());
        } else {
          LOGF(info, "Loaded acceptance weights from %s", cfg.cCCDB_NUA.value.c_str());
          conf.mAcceptance2D.push_back(hNUA2D);
        }
      } else {
        LOGF(info, "cfg.cCCDB_NUA empty! No corrections loaded");
      }
    }
    // Get Efficiency correction
    if (cfg.cCCDB_NUE.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfg.cCCDB_NUE, timestamp);
      conf.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency")));
      conf.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_pos")));
      conf.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_neg")));
      int sizeEff = conf.mEfficiency.size();
      if (sizeEff < nWeights)
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfg.cCCDB_NUE.value.c_str());
      else
        LOGF(info, "Loaded efficiency histogram from %s", cfg.cCCDB_NUE.value.c_str());
    } else {
      LOGF(info, "cfg.cCCDB_NUE empty! No corrections loaded");
    }
    // Get Efficiency correction
    if (cfg.cCCDB_NUE2D.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfg.cCCDB_NUE2D, timestamp);
      conf.mEfficiency2D.push_back(reinterpret_cast<TH2D*>(listCorrections->FindObject("Efficiency")));
      conf.mEfficiency2D.push_back(reinterpret_cast<TH2D*>(listCorrections->FindObject("Efficiency_pos")));
      conf.mEfficiency2D.push_back(reinterpret_cast<TH2D*>(listCorrections->FindObject("Efficiency_neg")));
      int sizeEff = conf.mEfficiency2D.size();
      if (sizeEff < nWeights)
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfg.cCCDB_NUE.value.c_str());
      else
        LOGF(info, "Loaded efficiency histogram from %s", cfg.cCCDB_NUE.value.c_str());
    } else {
      LOGF(info, "cfg.cCCDB_NUE2 empty! No corrections loaded");
    }

    if (cfg.cCCDB_NUE3D.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfg.cCCDB_NUE3D, timestamp);
      conf.mEfficiency3D.push_back(reinterpret_cast<TH3D*>(listCorrections->FindObject("Efficiency")));
      conf.mEfficiency3D.push_back(reinterpret_cast<TH3D*>(listCorrections->FindObject("Efficiency_pos")));
      conf.mEfficiency3D.push_back(reinterpret_cast<TH3D*>(listCorrections->FindObject("Efficiency_neg")));
      int sizeEff = conf.mEfficiency3D.size();
      if (sizeEff < nWeights)
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfg.cCCDB_NUE.value.c_str());
      else
        LOGF(info, "Loaded efficiency histogram from %s", cfg.cCCDB_NUE.value.c_str());
    } else {
      LOGF(info, "cfg.cCCDB_NUE2 empty! No corrections loaded");
    }
    conf.correctionsLoaded = true;
  }


  //Same for every particle species here
  bool setCurrentTrackWeights(const float& phi, const float& eta, const float& pt, const float& vtxz, const float& centrality)
  {
    float eff = 1.;
    if (cfg.cUseNUE3D) {
        int binx =
          conf.mEfficiency3D->GetXaxis()->FindBin(pt);
        int biny =
          conf.mEfficiency3D->GetYaxis()->FindBin(eta);
        int binz =
          conf.mEfficiency3D->GetZaxis()->FindBin(centrality);
        eff =
          conf.mEfficiency3D
          ->GetBinContent(binx,biny,binz);
    }
    if (eff <= 0)
        return false;
    spm.weff = 1./eff;
    if (cfg.cUseNUA1D) {
        spm.wacc =
          conf.mAcceptance->getNUA(
              phi,
              eta,
              vtxz);
    } else {
        spm.wacc = 1.;
    }
    return true;
  }

template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk)
  {
    if (!collision.sel8())
      return 0;
    histos.fill(HIST("hEventCount"), evSel_sel8);

    if (cfg.cEvtUseRCTFlagChecker && !rctChecker(collision))
      return 0;
    histos.fill(HIST("hEventCount"), evSel_RCTFlagsZDC);

    // Occupancy
    if (cfg.cEvSelsDoOccupancySel) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy > cfg.cEvSelsMaxOccupancy || occupancy < cfg.cEvSelsMinOccupancy) {
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_occupancy);
    }

    if (cfg.cEvSelsNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kNoSameBunchPileup);
    }
    if (cfg.cEvSelsIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsGoodZvtxFT0vsPV);
    }
    if (cfg.cEvSelsNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeStandard);
    }
    if (cfg.cEvSelsNoCollInTimeRangeNarrow) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        // Rejection of the collisions which have other events nearby
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeNarrow);
    }
    if (cfg.cEvSelsIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsVertexITSTPC);
    }

    if (cfg.cEvSelsIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        // New event selection bits to cut time intervals with dead ITS staves
        // https://indico.cern.ch/event/1493023/ (09-01-2025)
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }
    if (cfg.cEvSelsIsGoodITSLayer0123) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsGoodITSLayer0123);
    }

    if (cfg.cEvSelsUseAdditionalEventCut) {
      float vtxz = -999;
      if (collision.numContrib() > 1) {
        vtxz = collision.posZ();
        float zRes = std::sqrt(collision.covZZ());
        float minzRes = 0.25;
        int maxNumContrib = 20;
        if (zRes > minzRes && collision.numContrib() < maxNumContrib)
          vtxz = -999;
      }

      auto multNTracksPV = collision.multNTracksPV();

      if (vtxz > cfg.cEvSelsVtxZ || vtxz < -cfg.cEvSelsVtxZ)
        return 0;
      if (multNTracksPV < fMultPVCutLow->Eval(collision.centFT0C()))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(collision.centFT0C()))
        return 0;
      if (multTrk < fMultCutLow->Eval(collision.centFT0C()))
        return 0;
      if (multTrk > fMultCutHigh->Eval(collision.centFT0C()))
        return 0;

      histos.fill(HIST("hEventCount"), evSel_MultCuts);
    }

    return 1;
  }

  template <typename TrackObject>
  bool trackSelected(TrackObject track, const int& field)
  {
    if (std::fabs(track.eta()) > cfg.cTrackSelsEta)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_Eta);

    if (track.pt() < cfg.cTrackSelsPtmin || track.pt() > cfg.cTrackSelsPtmax)
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_Pt);

    // Edited myself: fabs
    if (std::fabs(track.dcaXY()) > cfg.cTrackSelsDCAxy)
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_DCAxy);

    if (std::fabs(track.dcaZ()) > cfg.cTrackSelsDCAz)
      return false;

    if (cfg.cTrackSelsDoDCApt && std::fabs(track.dcaZ()) > (cfg.cTrackSelsDCAptConsMin + (cfg.cTrackSelsDCApt1 * cfg.cTrackSelsDCApt2) / (std::pow(track.pt(), 1.1))))
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_DCAz);

    if (track.tpcNClsCrossedRows() < cfg.cTrackSelsNcls)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_NCls);

    if (track.tpcFractionSharedCls() > cfg.cTrackSelsFshcls)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_FshCls);

    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
    if (cfg.cFillTrackQA && cfg.cFillQABefore)
      histos.fill(HIST("incl/QA/before/pt_phi"), track.pt(), phimodn);

    if (cfg.cTrackSelsUseAdditionalTrackCut) {
      if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
        return false; // reject track
    }
    if (cfg.cFillTrackQA)
      histos.fill(HIST("incl/QA/after/pt_phi"), track.pt(), phimodn);
    histos.fill(HIST("hTrackCount"), trackSel_TPCBoundary);
    return true;
  }

  template <FillType ft, typename CollisionObject, typename TracksObject>
  inline void fillEventQA(CollisionObject collision, TracksObject tracks)
  {
    if (!cfg.cFillEventQA)
      return;

    static constexpr std::string_view Time[] = {"before", "after"};

    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hCentFT0C"), collision.centFT0C(), spm.centWeight);

    return;
  }

  template <ChargeType ct, PIDType pit, typename TrackObject>
  inline void fillHistograms(TrackObject track)
  {
    double weight = spm.wacc[ct] * spm.weff[ct] * spm.centWeight;
    float scale = 1.0;
    float minusQ = -1.0;
    if (track.eta() < 0)
      scale = -1.0;

    const double invSqrtQQ = 1.0 / std::sqrt(std::fabs(spm.corrQQ));
    const double invSqrtQQx = 1.0 / std::sqrt(std::fabs(spm.corrQQx));
    const double invSqrtQQy = 1.0 / std::sqrt(std::fabs(spm.corrQQy));
    const double uqA = spm.uy * spm.qyA + spm.ux * spm.qxA;
    const double uqC = spm.uy * spm.qyC + spm.ux * spm.qxC;
    const double invMeanPtQQ = 1.0 / spm.meanPtWeight;

    if (cfg.cFillGeneralV1Histos) {
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnAodd"), track.pt(), track.eta(), spm.centrality, scale * (uqA)*invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnCodd"), track.pt(), track.eta(), spm.centrality, scale * (uqC)*invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnOdd"), track.pt(), track.eta(), spm.centrality, scale * 0.5 * ((uqA) - (uqC)) * invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnEven"), track.pt(), track.eta(), spm.centrality, 0.5 * ((uqA) + (uqC)) * invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnA"), track.pt(), track.eta(), spm.centrality, (uqA)*invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnC"), track.pt(), track.eta(), spm.centrality, (uqC)*invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnCSetPlane"), track.pt(), track.eta(), spm.centrality, (spm.uy + spm.ux) * invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnASetPlane"), track.pt(), track.eta(), spm.centrality, (minusQ * spm.ux - spm.uy) * invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnAvnC"), track.pt(), track.eta(), spm.centrality, (uqA)*invSqrtQQ*(uqC)*invSqrtQQ, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnA_nw"), track.pt(), track.eta(), spm.centrality, (uqA)*invSqrtQQ);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnC_nw"), track.pt(), track.eta(), spm.centrality, (uqC)*invSqrtQQ);
    }

    if (cfg.cFillMixedHarmonics) {
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("MH/vnAxCxUx_MH"), track.pt(), track.eta(), spm.centrality, (spm.uxMH * spm.qxA * spm.qxC) / spm.corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("MH/vnAyCyUx_MH"), track.pt(), track.eta(), spm.centrality, (spm.uxMH * spm.qyA * spm.qyC) / spm.corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("MH/vnAxCyUy_MH"), track.pt(), track.eta(), spm.centrality, (spm.uyMH * spm.qxA * spm.qyC) / spm.corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("MH/vnAyCxUy_MH"), track.pt(), track.eta(), spm.centrality, (spm.uyMH * spm.qyA * spm.qxC) / spm.corrQQy, weight);
    }

    if (cfg.cFillXandYterms) {
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnAx"), track.pt(), track.eta(), spm.centrality, (spm.ux * spm.qxA) * invSqrtQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnAy"), track.pt(), track.eta(), spm.centrality, (spm.uy * spm.qyA) * invSqrtQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnCx"), track.pt(), track.eta(), spm.centrality, (spm.ux * spm.qxC) * invSqrtQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnCy"), track.pt(), track.eta(), spm.centrality, (spm.uy * spm.qyC) * invSqrtQQy, weight);
    }

    if (cfg.cFillEventPlane) { // only fill for inclusive!
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnA_EP"), track.pt(), track.eta(), spm.centrality, spm.vnA, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnC_EP"), track.pt(), track.eta(), spm.centrality, spm.vnC, weight);
      registry.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("vnFull_EP"), track.pt(), track.eta(), spm.centrality, spm.vnFull, weight);
    }

    
  }

  template <FillType ft, ChargeType ct, PIDType pit, typename TrackObject>
  inline void fillTrackQA(TrackObject track)
  {
    if (!cfg.cFillTrackQA)
      return;

    double weight = spm.wacc[ct] * spm.weff[ct] * spm.centWeight;

    static constexpr std::string_view Time[] = {"before/", "after/"};
    // NOTE: species[kUnidentified] = "" (when nocfg.cTrackSelDo) {
    if (cfg.cTrackSelDoTrackQAvsCent) {
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPt_Eta"), track.pt(), track.eta(), spm.centrality, weight);
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPt_Eta_uncorrected"), track.pt(), track.eta(), spm.centrality);
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta"), track.phi(), track.eta(), spm.centrality, weight);
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_uncorrected"), track.phi(), track.eta(), spm.centrality);
    } else {
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_Pt"), track.phi(), track.eta(), track.pt());
      histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_Pt_corrected"), track.phi(), track.eta(), track.pt(), weight);
    }

    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz"), track.phi(), track.eta(), spm.vz);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz_corrected"), track.phi(), track.eta(), spm.vz, spm.wacc[ct]);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAxy_pt"), track.pt(), track.dcaXY(), weight);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAz_pt"), track.pt(), track.dcaZ(), weight);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hSharedClusters_pt"), track.pt(), track.tpcFractionSharedCls(), weight);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_pt"), track.pt(), track.tpcNClsCrossedRows(), weight);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_vs_SharedClusters"), track.tpcNClsCrossedRows(), track.tpcFractionSharedCls(), weight);
    histos.fill(HIST(Charge[ct]) + HIST(PIDObservable[pit]) + HIST("QA/") + HIST(Time[ft]) + HIST("hMeanPtEta"), track.eta(), spm.centrality, track.pt(), weight);
  }

  template <FillType ft, ChargeType ct, PIDType pit, typename TrackObject>
  inline void fillPIDQA(TrackObject track)
  {
    if (!cfg.cFillTrackQA)
      return;

    if constexpr (framework::has_type_v<aod::pidtof::TOFNSigmaPi, typename TrackObject::all_columns>) {

      // Always fill these global QA plots for any processed track
      histos.fill(HIST(Charge[ct]) + HIST("QA/") + HIST(Time[ft]) + HIST("hdEdxTPC_pt"), track.pt(), track.tpcSignal());
      histos.fill(HIST(Charge[ct]) + HIST("QA/") + HIST(Time[ft]) + HIST("hBetaTOF_pt"), track.pt(), track.beta());
    }
  }

  template <FillType ft, PIDType pit, typename TrackObject>
  void fillAllQA(TrackObject track)
  {
    fillTrackQA<ft, kInclusive, pit>(track);
    fillPIDQA<ft, kInclusive, pit>(track); // <-- Added pit here

    if (cfg.cFillChargeDependenceQA) {
      switch (spm.charge) {
        case kPositive: {
          fillTrackQA<ft, kPositive, pit>(track);
          fillPIDQA<ft, kPositive, pit>(track); // <-- Added pit here
          break;
        }
        case kNegative: {
          fillTrackQA<ft, kNegative, pit>(track);
          fillPIDQA<ft, kNegative, pit>(track); // <-- Added pit here
          break;
        }
      }
    }
  }


  void processStatisticalPID(ZDCCollisions::iterator const& collision, aod::BCsWithTimestamps const&, UsedTracksPID const& tracks)
  {

    histos.fill(HIST("hEventCount"), evSel_FilteredEvent);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int standardMagField = 99999;
    auto field = (cfg.cMagField == standardMagField) ? getMagneticField(bc.timestamp()) : cfg.cMagField;

    if (bc.runNumber() != conf.lastRunNumber) {
      conf.correctionsLoaded = false;
      conf.clCentrality = false;
      conf.lastRunNumber = bc.runNumber();
      conf.mAcceptance.clear();
      LOGF(info, "Size of mAcceptance: %i (should be 0)", (int)conf.mAcceptance.size());
    }

    if (cfg.cFillQABefore)
      fillEventQA<kBefore>(collision, tracks);

    loadCorrections(bc.timestamp());

    spm.centrality = collision.centFT0C();

    if (cfg.cCentFT0Cvariant1)
      spm.centrality = collision.centFT0CVariant1();
    if (cfg.cCentFT0M)
      spm.centrality = collision.centFT0M();
    if (cfg.cCentFV0A)
      spm.centrality = collision.centFV0A();
    if (cfg.cCentNGlobal)
      spm.centrality = collision.centNGlobal();

    if (!eventSelected(collision, tracks.size()))
      return;

    if (!collision.isSelected()) // selected by ZDCQVectors task (checks signal in ZDC) --> only possible in data not MC
      return;
    histos.fill(HIST("hEventCount"), evSel_isSelectedZDC);

    // Always fill centrality histogram after event selections!
    histos.fill(HIST("hCentrality"), spm.centrality);

    spm.qxA = collision.qxA();
    spm.qyA = collision.qyA();
    spm.qxC = collision.qxC();
    spm.qyC = collision.qyC();

    spm.vz = collision.posZ();

    spm.psiA = 1.0 * std::atan2(spm.qyA, spm.qxA);
    spm.psiC = 1.0 * std::atan2(spm.qyC, spm.qxC);

    // https://twiki.cern.ch/twiki/pub/ALICE/DirectedFlowAnalysisNote/vn_ZDC_ALICE_INT_NOTE_version02.pdf
    spm.psiFull = 1.0 * std::atan2(spm.qyA + spm.qyC, spm.qxA + spm.qxC);

    if (spm.centrality > cfg.cCentMax || spm.centrality < cfg.cCentMin)
      return;

    // Load correlations and SP resolution needed for Scalar Product and event plane methods.
    // Only load once!
    // If not loaded set to 1

    if (cfg.cCCDBdir_QQ.value.empty() == false) {
      if (!conf.clQQ) {
        TList* hcorrList = ccdb->getForTimeStamp<TList>(cfg.cCCDBdir_QQ.value, bc.timestamp());
        conf.hcorrQQ = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCXY"));
        conf.hcorrQQx = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCX"));
        conf.hcorrQQy = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCY"));
        conf.clQQ = true;
      }
      spm.corrQQ = conf.hcorrQQ->GetBinContent(conf.hcorrQQ->FindBin(spm.centrality));
      spm.corrQQx = conf.hcorrQQx->GetBinContent(conf.hcorrQQx->FindBin(spm.centrality));
      spm.corrQQy = conf.hcorrQQy->GetBinContent(conf.hcorrQQy->FindBin(spm.centrality));
    }

    double evPlaneRes = 1.;
    if (cfg.cCCDBdir_SP.value.empty() == false) {
      if (!conf.clEvPlaneRes) {
        conf.hEvPlaneRes = ccdb->getForTimeStamp<TProfile>(cfg.cCCDBdir_SP.value, bc.timestamp());
        conf.clEvPlaneRes = true;
      }
      evPlaneRes = conf.hEvPlaneRes->GetBinContent(conf.hEvPlaneRes->FindBin(spm.centrality));
      if (evPlaneRes < 0)
        LOGF(fatal, "<Cos(PsiA-PsiC)> > 0 for centrality %.2f! Cannot determine resolution.. Change centrality ranges!!!", spm.centrality);
      evPlaneRes = std::sqrt(evPlaneRes);
    }

    spm.centWeight = 1.;
    if (cfg.cCCDBdir_centrality.value.empty() == false) {
      if (!conf.clCentrality) {
        conf.hCentrality = ccdb->getForTimeStamp<TH1D>(cfg.cCCDBdir_centrality.value, bc.timestamp());
        conf.clCentrality = true;
      }
      double centW = conf.hCentrality->GetBinContent(conf.hCentrality->FindBin(spm.centrality));

      if (centW <= 0) {
        spm.centWeight = 0;
        LOGF(fatal, "Centrality weight cannot be negative .. setting to 0. for (%.2f)", spm.centrality);
      }
    }

    fillEventQA<kAfter>(collision, tracks);

    for (const auto& track : tracks) {

      //histos.fill(HIST("hPIDcounts"), trackPID, track.pt());

      if (track.sign() == 0)
        continue;

      histos.fill(HIST("hTrackCount"), trackSel_ZeroCharge);

      spm.charge = ((track.sign() > 0)) ? kPositive : kNegative;

      if (cfg.cFillQABefore) {
        fillTrackQA<kBefore, kInclusive, kNSigmaTPC>(track);

        if (track.hasTOF()) {
            fillPIDQA<kBefore, kInclusive, kNSigmaTOF>(track);
        }

        if (track.hasTPC()) {
            fillPIDQA<kBefore, kInclusive, kTPC>(track);
        }
      }

      if (!trackSelected(track, field))
        continue;



      // constrain angle to 0 -> [0,0+2pi]
      auto phi = RecoDecay::constrainAngle(track.phi(), 0);
      
      // efficiency + acceptance
      if(!setCurrentTrackWeights(phi,eta,pt,spm.vz,spm.centrality))
        continue;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      spm.ux = std::cos(cfg.cHarm * phi);
      spm.uy = std::sin(cfg.cHarm * phi);

      spm.uxMH = std::cos(cfg.cHarmMixed * phi);
      spm.uyMH = std::sin(cfg.cHarmMixed * phi);
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      spm.vnA = std::cos(cfg.cHarm * (phi - spm.psiA)) / evPlaneRes;
      spm.vnC = std::cos(cfg.cHarm * (phi - spm.psiC)) / evPlaneRes;
      spm.vnFull = std::cos(cfg.cHarm * (phi - spm.psiFull)) / evPlaneRes;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      

      float scale = 1.0;
      float minusQ = -1.0;
      if (track.eta() < 0)
        scale = -1.0;

      const double invSqrtQQ = 1.0 / std::sqrt(std::fabs(spm.corrQQ));
      const double invSqrtQQx = 1.0 / std::sqrt(std::fabs(spm.corrQQx));
      const double invSqrtQQy = 1.0 / std::sqrt(std::fabs(spm.corrQQy));
      const double uqA = spm.uy * spm.qyA + spm.ux * spm.qxA;
      const double uqC = spm.uy * spm.qyC + spm.ux * spm.qxC;
      float vnOdd = scale*0.5*(uqA-uqC)*invSqrtQQ;
      double weight = spm.wacc * spm.weff * spm.centWeight;

      
      // PID variable, start with TPC only
      if(!track.hasTPC())
          continue;

      float TPCsignal = track.tpcSignal();

      // histogram coordinates

      std::vector<double> values = {pt,eta,spm.centrality,TPCsignal};

      // Yield and flow distribution 
      if (spm.charge == kPositive) {
        histos.fill(HIST("StatPID/hTPCYieldPos"),values,weight);
        histos.fill(HIST("StatPID/hTPCSPOddPos"),values,weight*vnOdd);
      } else if (spm.charge == kNegative) {
        histos.fill(HIST("StatPID/hTPCYieldNeg"),values,weight);
        histos.fill(HIST("StatPID/hTPCSPOddNeg"),values,weight*vnOdd);
      }    

    } // end of track loop
  } // end of processStatisticalPID

  PROCESS_SWITCH(FlowSPstatPID, processStatisticalPID, "Process statistical PID determination", false);

}

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSPstatPID>(cfgc)};
}
