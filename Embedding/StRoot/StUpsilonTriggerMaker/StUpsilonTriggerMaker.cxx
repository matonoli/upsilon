//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 26 June 2008
//
// StUpsilonTriggerMaker simulates the level-2 trigger offline.
//

// STAR includes
#include "StEventTypes.h"
#include "StDaqLib/EMC/StEmcDecoder.h"
#include "StEmcRawMaker/defines.h"

// Local includes
#include "Quarkonium.hh"
#include "L2Result.h"
#include "StUpsilonTriggerMaker.h"

unsigned short StUpsilonTriggerMaker::ReturnADCvalueFromDAQID(unsigned short DaqID)
{
  if (DaqID<4800){
    return mBemcData[DaqID];
  }
  else {
    return 5000;
  }
}

void StUpsilonTriggerMaker::Clear(Option_t* option)
{
  memset(mTriggerData, 0, sizeof(TrgDataType));
  memset(mBemcData, 0, sizeof(mBemcData));

  StMaker::Clear(option);
}

int StUpsilonTriggerMaker::Init()
{
  mTriggerData = new TrgDataType;
  mEmcDecoder = new StEmcDecoder;

  return StMaker::Init();
}

int StUpsilonTriggerMaker::InitRun(int runNumber)
{
  cout<<"StUpsilonTriggerMaker Init Run is Called !!!"<<endl;
  mEmcDecoder->SetDateTime(GetDate(), GetTime());

  // L2 parameters
  static int userInt[] = { 12, 75, 0, 0, 3 };
  static float userFloat[] = { 4.0, 2.5, 6.0, 15.0, 0.5 };

  // Initialize L2 algorithm
  ups.init(runNumber, userInt, userFloat);

  return StMaker::InitRun(runNumber);
}

int StUpsilonTriggerMaker::Make()
{
  // Get StEvent
  mEvent = (StEvent*)GetDataSet("StEvent");
  if (!mEvent) {
    LOG_WARN << "No StEvent" << endm;
    return kStWarn;
  }

  // Get EMC collection
  StEmcCollection* emc = mEvent->emcCollection();
  if (!emc) {
    LOG_WARN << "No EMC collection" << endm;
    return kStWarn;
  }

  // Get BTOW
  StEmcDetector* det = emc->detector(kBarrelEmcTowerId);
  if (!det) {
    LOG_WARN << "No StEmcDetector(kBarrelEmcTowerId)" << endm;
    return kStWarn;
  }

  // Fill BEMC tower ADC's
  for (size_t m = 1; m <= det->numberOfModules(); ++m) {
    StSPtrVecEmcRawHit& hits = det->module(m)->hits();
    for (size_t i = 0; i < hits.size(); ++i) {
      StEmcRawHit* hit = hits[i];
      int daqId;
      mEmcDecoder->GetDaqIdFromTowerId(hit->softId(BTOW), daqId);
      mBemcData[daqId] = hit->adc();
    }
  }

  // Get L2 decision
  mIsUpsilonTrigger = ups.run(mTriggerData, mBemcData);
  cout<<"Candidates are: "<<endl;
for(int i = 0; i<ups.candidates.size();++i){
  cout<<"candidate "<<i<<" first = "<<ups.candidates[i].first<<" second = "<<ups.candidates[i].second<<endl;}

  return kStOk;
}

int StUpsilonTriggerMaker::Finish()
{
  return kStOk;
}

int StUpsilonTriggerMaker::FinishRun(int runNumber)
{
  return StMaker::FinishRun(runNumber);
}

bool StUpsilonTriggerMaker::isUpsilonTrigger() const
{
  return mIsUpsilonTrigger;
}
