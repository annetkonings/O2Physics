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
#ifndef O2_ANALYSIS_ZDCEXTRATABLE_H
#define O2_ANALYSIS_ZDCEXTRATABLE_H

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <cmath>
#include <vector>

namespace o2::aod
{

namespace zdcextra
{
DECLARE_SOA_COLUMN(EnergyZNA1, energyZNA1, float);
DECLARE_SOA_COLUMN(EnergyZNA2, energyZNA2, float);
DECLARE_SOA_COLUMN(EnergyZNA3, energyZNA3, float);
DECLARE_SOA_COLUMN(EnergyZNA4, energyZNA4, float);

DECLARE_SOA_COLUMN(EnergyZNC1, energyZNC1, float);
DECLARE_SOA_COLUMN(EnergyZNC2, energyZNC2, float);
DECLARE_SOA_COLUMN(EnergyZNC3, energyZNC3, float);
DECLARE_SOA_COLUMN(EnergyZNC4, energyZNC4, float);
} // namespace zdcextra

DECLARE_SOA_TABLE(
  ZDCExtra,
  "AOD",
  "ZDCEXTRA",
  zdcextra::EnergyZNA1,
  zdcextra::EnergyZNA2,
  zdcextra::EnergyZNA3,
  zdcextra::EnergyZNA4,
  zdcextra::EnergyZNC1,
  zdcextra::EnergyZNC2,
  zdcextra::EnergyZNC3,
  zdcextra::EnergyZNC4);

} // namespace o2::aod
#endif // O2_ANALYSIS_ZDCEXTRATABLE_H
