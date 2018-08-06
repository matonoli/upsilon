// -*- mode: C++ -*-
//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 26 June 2008
//
// StUpsilonTriggerMaker simulates the level-2 trigger offline.
//

#ifndef ST_L2_SIMULATOR_MAKER_H
#define ST_L2_SIMULATOR_MAKER_H

// STAR forward declarations
class StEvent;
struct TrgDataType2005;
class StEmcDecoder;

// STAR includes
#include "StMaker.h"

// Local includes
#include "trgStructures.h"

class StUpsilonTriggerMaker : public StMaker {
public:
  StUpsilonTriggerMaker(const char* name = "StUpsilonTriggerMaker") : StMaker(name) {}
  ~StUpsilonTriggerMaker() {}

  void Clear(Option_t* option = "");
  int  Init();
  int  InitRun(int runNumber);
  int  Make();
  int  Finish();
  int  FinishRun(int runNumber);
  bool isUpsilonTrigger() const;
  unsigned short ReturnADCvalueFromDAQID(unsigned short DaqID);

private:
  StEvent* mEvent;
  TrgDataType* mTriggerData;
  unsigned short mBemcData[4800];
  StEmcDecoder* mEmcDecoder;
  bool mIsUpsilonTrigger;

  ClassDef(StUpsilonTriggerMaker, 1);
};

#endif
