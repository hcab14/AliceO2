// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SAMPAProcessing.cxx
/// \brief Implementation of the SAMPA response
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/CDBInterface.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "FairLogger.h"

using namespace o2::tpc;

SAMPAProcessing::SAMPAProcessing() : mRandomNoiseRing()
{
  updateParameters();
}

SAMPAProcessing::~SAMPAProcessing() = default;

void SAMPAProcessing::updateParameters()
{
  mGasParam = &(ParameterGas::Instance());
  mDetParam = &(ParameterDetector::Instance());
  mEleParam = &(ParameterElectronics::Instance());
  auto& cdb = CDBInterface::instance();
  mPedestalMap = &(cdb.getPedestals());
  mNoiseMap = &(cdb.getNoise());
}

void SAMPAProcessing::getShapedSignal(float ADCsignal, float driftTime, std::vector<float>& signalArray) const
{
  const float timeBinTime = getTimeBinTime(driftTime);
  const float offset = driftTime - timeBinTime;
  Vc::float_v binVector = Vc::float_v::IndexesFromZero() * mEleParam->ZbinWidth;
  const Vc::float_v binInc = Vc::float_v::Size* mEleParam->ZbinWidth;
  Vc::float_v signal;
  float *p = signalArray.data();
  for (int bin = 0; bin < mEleParam->NShapedPoints; bin += Vc::float_v::Size) {
    const Vc::float_v time = timeBinTime + binVector;
    signal = getGamma4(time, Vc::float_v(timeBinTime + offset), Vc::float_v(ADCsignal));
    signal.store(p);
    p += Vc::float_v::Size;
    Vc::prefetchForOneRead(p);
    binVector += binInc;
  }
}
