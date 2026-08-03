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

#ifndef O2_ANALYSIS_SKIMTRACKS_H
#define O2_ANALYSIS_SKIMTRACKS_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(DrCollisions, "AOD", "DRCOLLISION", o2::soa::Index<>,
                  o2::aod::collision::PosX, o2::aod::collision::PosY, o2::aod::collision::PosY);
using DrCollision = DrCollisions::iterator;

namespace exampleTrackSpace
{
DECLARE_SOA_INDEX_COLUMN(DrCollision, drCollision);
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(Snp, snp, float);
DECLARE_SOA_COLUMN(Tgl, tgl, float);
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);
DECLARE_SOA_COLUMN(IsWithinBeamPipe, isWithinBeamPipe, bool);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Sign, sign, short);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace exampleTrackSpace

DECLARE_SOA_TABLE(DrTracks, "AOD", "DRTRACK", o2::soa::Index<>, exampleTrackSpace::DrCollisionId,
                  exampleTrackSpace::TrackType, exampleTrackSpace::X, exampleTrackSpace::Alpha, exampleTrackSpace::Y,
                  exampleTrackSpace::Z, exampleTrackSpace::Snp, exampleTrackSpace::Tgl, exampleTrackSpace::Signed1Pt,
                  exampleTrackSpace::IsWithinBeamPipe, exampleTrackSpace::Px, exampleTrackSpace::Py, exampleTrackSpace::Pz,
                  exampleTrackSpace::Sign, exampleTrackSpace::Eta, exampleTrackSpace::Phi);
using DrTrack = DrTracks::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_SKIMTRACKS_H
