//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 5 June 2008
//

#ifndef ST_UPSILON_EVENT_H
#define ST_UPSILON_EVENT_H

// ROOT forward declarations
class TClonesArray;

// STAR forward declarations
class StMuEvent;

// C++ STL includes
#include <set>
using std::set;

// ROOT
#include "TObject.h"
#include "TDatime.h"

class StUpsilonEvent : public TObject {
public:
  StUpsilonEvent();
  ~StUpsilonEvent();

  void Clear(Option_t* opton = "");
  void setEventHeader(StMuEvent* event);

  int runNumber;
  int eventNumber;
  int refMult;
  int numberOfTriggerCandidates;
  TDatime time;
  float magField;
  set<int> triggerIds;
  bool isTrigger(int id) const;
  unsigned int L2Result[64];
  float vx;
  float vy;
  float vz;
  int bbcTimeDiff;
  TClonesArray* candidates;

private:
  StUpsilonEvent(const StUpsilonEvent& event);
  void operator=(const StUpsilonEvent& event);

  ClassDef(StUpsilonEvent,2);
};

#endif
