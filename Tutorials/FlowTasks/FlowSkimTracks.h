#ifndef O2_ANALYSIS_FLOWSKIMTRACKS_H
#define O2_ANALYSIS_FLOWSKIMTRACKS_H

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <vector>

namespace o2::aod
{

// General
// Event / Collision Table, BCs

namespace drcollision
{
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint32_t);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);

DECLARE_SOA_COLUMN(MultNTracksGlobal, multNTracksGlobal, int);
DECLARE_SOA_COLUMN(MultNGlobalTracksPV, multNGlobalTracksPV, int);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);

// Centrality (FT0C)
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);

// ZDC spectator neutrons
DECLARE_SOA_COLUMN(DyEnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(DyEnergyCommonZNC, energyCommonZNC, float);
} // namespace drcollision

DECLARE_SOA_TABLE(DrCollisions, "AOD", "DRCOLLISION",
                  o2::soa::Index<>,
                  drcollision::PosX,
                  drcollision::PosY,
                  drcollision::PosZ,
                  drcollision::Flags,
                  drcollision::Timestamp,
                  drcollision::NumContrib,
                  drcollision::RunNumber,
                  drcollision::MultNTracksGlobal,
                  drcollision::MultNGlobalTracksPV,
                  drcollision::MultNTracksPV,
                  drcollision::CentFT0C,
                  drcollision::DyEnergyCommonZNA,
                  drcollision::DyEnergyCommonZNC);

using DrCollision = DrCollisions::iterator;

// Track Table

namespace drtrack
{

DECLARE_SOA_INDEX_COLUMN(DrCollision, drCollision);

// basic kinematics
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, short);

// track selection flags
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);

// DCA
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);

// TPC quality
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint16_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint16_t);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint16_t);
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);

// ITS quality
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(IsGlobalTrackSDD, isGlobalTrackSDD, uint8_t);

// Flags
DECLARE_SOA_COLUMN(Flags, flags, uint32_t);

// PID
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);

DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);

} // namespace drtrack

DECLARE_SOA_TABLE(DrTracks, "AOD", "DRTRACK",
                  o2::soa::Index<>,
                  drtrack::DrCollisionId,
                  drtrack::TrackType,
                  drtrack::Px,
                  drtrack::Py,
                  drtrack::Pz,
                  drtrack::Eta,
                  drtrack::Phi,
                  drtrack::Sign,
                  drtrack::IsGlobalTrack,
                  drtrack::DcaXY,
                  drtrack::DcaZ,
                  drtrack::TPCSignal,
                  drtrack::TPCNClsFindable,
                  drtrack::TPCNClsFound,
                  drtrack::TPCNClsCrossedRows,
                  drtrack::TPCChi2NCl,
                  drtrack::TPCFractionSharedCls,
                  drtrack::ITSClusterSizes,
                  drtrack::ITSNCls,
                  drtrack::ITSChi2NCl,
                  drtrack::IsGlobalTrackSDD,
                  drtrack::Flags,
                  drtrack::TPCNSigmaPi,
                  drtrack::TPCNSigmaKa,
                  drtrack::TPCNSigmaPr,
                  drtrack::TOFNSigmaPi,
                  drtrack::TOFNSigmaKa,
                  drtrack::TOFNSigmaPr);

using DrTrack = DrTracks::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_FLOWSKIMTRACKS_H
