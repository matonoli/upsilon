//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 5 June 2008
//

// ROOT includes
#include "TClonesArray.h"

// STAR includes
#include "StMuDSTMaker/COMMON/StMuTypes.hh"

// Local includes
#include "StUpsilonCandidate.h"
#include "StUpsilonEvent.h"

ClassImp(StUpsilonEvent);

StUpsilonEvent::StUpsilonEvent()
{
  candidates = new TClonesArray("StUpsilonCandidate");
}

StUpsilonEvent::StUpsilonEvent(const StUpsilonEvent& event)
{
  assert(0);
  runNumber = event.runNumber;
  eventNumber = event.eventNumber;
  refMult = event.refMult;
  time = event.time;
  vx = event.vx;
  vy = event.vy;
  vz = event.vz;
  magField = event.magField;
  triggerIds = event.triggerIds;
  memcpy(L2Result, event.L2Result, sizeof(L2Result));
  bbcTimeDiff = event.bbcTimeDiff;
  candidates = new TClonesArray("StUpsilonCandidate");
  for (int i = 0; i < event.candidates->GetEntriesFast(); ++i) {
    new ((*candidates)[candidates->GetEntriesFast()]) StUpsilonCandidate(*(StUpsilonCandidate*)event.candidates->At(i));
  }
}

StUpsilonEvent::~StUpsilonEvent()
{
  Clear();
  delete candidates;
}

void StUpsilonEvent::Clear(Option_t* option)
{
  candidates->Clear(option);
}

void StUpsilonEvent::setEventHeader(StMuEvent* event)
{
  runNumber = event->runNumber();
  eventNumber = event->eventNumber();
  refMult = event->refMult();
  time.Set(event->eventInfo().time());
  StThreeVectorF v = event->primaryVertexPosition();
  vx = v.x();
  vy = v.y();
  vz = v.z();
  magField = event->magneticField();
  vector<unsigned int> triggerIds = event->triggerIdCollection().nominal().triggerIds();
  copy(triggerIds.begin(), triggerIds.end(), inserter(this->triggerIds, this->triggerIds.begin()));
  memcpy(L2Result, event->L2Result().GetArray(), sizeof(L2Result));
  bbcTimeDiff = event->bbcTriggerDetector().onlineTimeDifference();
}

bool StUpsilonEvent::isTrigger(int id) const
{
  return triggerIds.find(id) != triggerIds.end();
}
