// -*- mode: C++ -*-
//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 3 June 2008
//

#ifndef ST_UPSILON_CANDIDATE_H
#define ST_UPSILON_CANDIDATE_H

// ROOT
#include "TObject.h"

class StUpsilonCandidate : public TObject {
public:
  StUpsilonCandidate() {}
  StUpsilonCandidate(const StUpsilonCandidate& ups);
  ~StUpsilonCandidate() {}

  int charge;
  float p;
  float pt;
  float pz;
  float eta;
  float phi;
  float y;
  float m;
  float cosTheta;

  int charge1;
  float p1;
  float pt1;
  float pz1;
  float eta1;
  float phi1;
  float nHitsFit1;
  float dEdx1;
  float Etrig1;
  float Ecluster1;
  float nSigmaElectron1;
  float nSigmaPion1;
  float dcaGlobal1;

  int charge2;
  float p2;
  float pt2;
  float pz2;
  float eta2;
  float phi2;
  float nHitsFit2;
  float dEdx2;
  float Etrig2;
  float Ecluster2;
  float nSigmaElectron2;
  float nSigmaPion2;
  float dcaGlobal2;

private:
  ClassDef(StUpsilonCandidate,1);
};

#endif
