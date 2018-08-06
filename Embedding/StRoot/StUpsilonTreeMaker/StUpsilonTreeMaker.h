// -*- C++ -*-

//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 4 October 2006
//
// Revised 3 June 2008
//

#ifndef StUpsilonTreeMaker_h
#define StUpsilonTreeMaker_h

// ROOT
class TFile;
class TTree;
class TH1;

// STAR
class StMuEvent;
class StMuTrack;
class StMuEmcHit;
class StEmcDecoder;
class StEmcGeom;
class StBemcTables;
class StEmcPosition;

// Local forward declarations
class StUpsilonEvent;

// STAR includes
#include "StMaker.h"

class StUpsilonTreeMaker : public StMaker {
public:
  StUpsilonTreeMaker(const char* name = "StUpsilonTreeMaker") : StMaker(name) {}
  ~StUpsilonTreeMaker() {}

  void Clear(Option_t* option = "");
  int Init();
  int InitRun(int runNumber);
  int Make();
  int Finish();
  int makeBemcStatusTable(const char* filename);

private:
  // bool accept(StMuEvent* event) const;
  //bool accept(StMuTrack* track) const;
  StMuEmcHit* getEmcHit(int id) const;
  StMuEmcHit* addEmcHit(int id);
  //int makeBemcStatusTable(const char* filename);
  //void getTowerNeighbors(int id, vector<int>& neighbors);

  TFile* mFile;
  TTree* mTree;
  StUpsilonEvent* mEvent;
  StEmcDecoder* mEmcDecoder;
  StEmcGeom* mEmcGeom;
  StBemcTables* mBemcTables;
  StEmcPosition* mEmcPosition;
  TH1* hNumberOfCandidates;
  TClonesArray* mEmcHits;

  ClassDef(StUpsilonTreeMaker, 1);
};

#endif
