// -*- mode: C+++ -*- Put Emacs editor in C++ mode
//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 24 June 2008
//
// StUpsilonEmbedMaker is used to run over the Upsilon embedding and
// calculate acceptance, L0 and L2 trigger efficiencies.
//

#ifndef ST_UPSILON_EMBED_MAKER_H
#define ST_UPSILON_EMBED_MAKER_H

// ROOT forward declarations
class TNtuple;
class TH1;
class TH2;
class TList;

// STAR forward declarations
class StMcEvent;
class StMcTrack;
class StEvent;
class StTrack;
class StEmcDecoder;
class StBemcTables;
class StUpsilonTriggerMaker;
class StGlobalTrack;

// STAR includes
#include "StMaker.h"
#include "StThreeVectorF.hh"

class StUpsilonEmbedMaker : public StMaker {
public:
  StUpsilonEmbedMaker(const char* name = "StUpsilonEmbedMaker") : StMaker(name) {}
  ~StUpsilonEmbedMaker() {}

  int    Init();
  int    Make();
  TList* GetTreeList() { return mTreeList; }

  void SetIndex(char* var){mIndex = var;}
  void SetRunNumber(char* var){mRunNumber = var;}
  int SetHighTowerVar(StMcTrack* mcTrack, bool isele);
  //StTrack* findPartner(StMcTrack* mcTrack);
  int getTPADC(int id);
  int InitRun(int runNumber);

private:
  StTrack* findPartner(StMcTrack* mcTrack);
  void getClusterEtaPhiE(int& id, pair<float,float>& etaphi, float& Ecluster, float& E2, float& E3, int& shape);
  void TrackEtaPhiId(const StTrack* mcTrack,pair<float,float>& etaphi, pair<float,float>& etaphiSmd, int& softid);
  void getZPosPhiPos(const StTrack* mcTrack,float* d);
  int FindShape(int& id1, int& id2, int& id3);
  void getTowerNeighbors(int id, vector<int>& neighbors);
  int getHighestNeighbor(int myId);

  Int_t getgRefMult(StEvent *event, StThreeVectorF vtxPos);
  Double_t getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos);

  ofstream ofs;

  StEmcDecoder* mEmcDecoder;
  StBemcTables* mBemcTables;

  
  float bemcEnergy[4801];
  int bemcADC[4801];
  //float bemczPos[4801];
  //float bemcphiPos[4801];
  TNtuple* mTuple;
  StMcEvent* mMcEvent;
  StEvent* mEvent;
  TList* mTreeList;
  char* mIndex;
  char* mRunNumber;
  StUpsilonTriggerMaker* mUpsTrig;
  short unsigned int eleADC;
  short unsigned int posADC;
  double eleE;
  double posE;
  double eleEsum;
  double posEsum;



  ClassDef(StUpsilonEmbedMaker,1);
};

#endif
