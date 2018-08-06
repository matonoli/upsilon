#ifndef STMTDJPSIUTIL_HH
#define STMTDJPSIUTIL_HH

/***************************************************************************
 *
 * $Id: StMtdJpsiUtil.h,v 1.42 2018/01/18 16:20:32 marr Exp $ 
 * StMtdJpsiUtil
 * Author: Rongrong Ma
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1;
class TGraphErrors;

class StEvent;
class StTrack;
class StGlobalTrack;
class StMtdHit;
class StVertex;
class StMtdPidTraits;

class StMuDst;
class StMuTrack;
class StMuMtdHit;
class StMtdGeometry;

class StMtdTrigger;

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoMtdHit;
class StRefMultCorr;

#include "StPhysicalHelixD.hh"
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "mtdMaps.h"

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<TLorentzVector> LorentzVec;
typedef vector<Double_t> DoubleVec;
typedef vector<StThreeVectorD > PointVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<StThreeVectorD, allocator<StThreeVectorD>> PointVec;
#endif

enum { kNtrig = 4 };  // di-muon -- single-muon -- e-mu -- mb -- vpdNoVtx
enum { kNTrgSetup = 4 };
static const char *trigName[kNtrig] = {"di_mu","si_mu","e_mu","mb"};


class StMtdJpsiUtil : public TNamed {
 public:
  StMtdJpsiUtil(const Char_t *name  = "MtdJpsiUtil",
		const Char_t *title = "Utility class for j/psi analysis using MTD");
  StMtdJpsiUtil(const StMtdJpsiUtil &mtdUtil);
  ~StMtdJpsiUtil();

  void     Init();

  void     setMCflag(Int_t flag = 0);
  Int_t    getMCflag() { return mMCflag; }

  //cuts
  Bool_t   passRun(const int runId);

  Bool_t   passTrigger(StEvent *event);
  Bool_t   passTrigger(StMuDst *event);
  Bool_t   passTrigger(StPicoDst *event);

  Bool_t   passEvent(StEvent *event); // vertex cuts
  Bool_t   passEvent(StMuDst *event); // vertex cuts
  Bool_t   passEvent(StPicoEvent *event);

  Bool_t   passTrack(StMuTrack *t) const ;
  Bool_t   passTrack(StTrack *t, StVertex *vtx) const ;
  Bool_t   passTrack(StPicoTrack *t, StThreeVectorF pVtx, float B) const;

  // trigger selection
  void     setYear(const Int_t year);
  void     setTrigger(Bool_t a, Bool_t b, Bool_t c, Bool_t d);
  void     setTriggerIDs(const IntVec id[kNtrig]); // di-muon -- single-muon -- e-mu
  Int_t    getTriggerType(StEvent *event);
  Int_t    getTriggerType(StMuDst *event);
  Int_t    getTriggerType(StPicoDst *event);
  Int_t    getTrgSetup(const int runId);
  Bool_t   isTriggerOn(Int_t i)   { return mTrigger[i]; }
  Bool_t   isDimuonTrigger(StMuDst *event);
  Bool_t   isSimuonTrigger(StMuDst *event);
  Bool_t   isEmuonTrigger(StMuDst *event);
  Bool_t   isVpdNoVtxTrigger(StMuDst *event);
  Bool_t   isDimuonTrigger(StPicoDst *pico);
  Bool_t   isSimuonTrigger(StPicoDst *pico);
  Bool_t   isEmuonTrigger(StPicoDst *pico);
  Bool_t   isVpdNoVtxTrigger(StPicoDst *pico);
  void     printTriggerIDs();

  // centrality 
  void     setRefMultCorr(StRefMultCorr *relMultCorr)           { mRefMultCorr = relMultCorr;     }
  void     setVzCorrMethod(const Int_t m)                       { mVzCorrMethod = m;              }
  void     setVzCorrFileName(const char* name)                  { mVzCorrFileName = name;         }
  Double_t getRefMultCorr(Int_t year, Int_t runId, Double_t gRefMult, Double_t z_val, Double_t zdc_val);
  Double_t getWeight(Int_t year, Double_t zdc_val, Double_t gRefMultCorr);
  Double_t getCentralityBin16(Int_t year, Double_t zdc_val, Double_t gRefMultCorr);
  
  // vertex cuts
  void     setRequireVtxRanking(const Bool_t r = kFALSE);
  void     setRequireMtdHitForPrimVtx(const Bool_t r = kFALSE);
  void     setUsePrimVtxClosestToVpd(const Bool_t r = kFALSE);
  void     setMaxVtxZ(const Double_t max);
  void     setMaxVtxR(const Double_t max);
  void     setMaxDiffz(const Double_t max);
  Int_t    selectPrimaryVertex(StEvent *event);
  Int_t    selectPrimaryVertex(StMuDst *event);
  Bool_t   isGoodVertex(const Double_t vpd_vz, const Double_t tpc_vz, const Double_t tpc_vr);
  Double_t getMaxVtxZ()                { return mMaxVtxZ;                 }
  Double_t getMaxVtxR()                { return mMaxVtxR;                 }
  Double_t getMaxDiffz()               { return mMaxDiffz;                }
  Bool_t   isRequireMtdHitForPrimVtx() { return mRequireMtdHitForPrimVtx; }
  Bool_t   isUsePrimVtxClosestToVpd()  { return mUsePrimVtxClosestToVpd;  }
  Bool_t   isRequireVtxRanking()       { return mRequireVtxRanking;       }
  
  // track cuts
  void     setTrackType(Int_t type = 0);
  void     setTrackPtLimits(const Double_t min, const Double_t max);
  void     setTrackPhiLimits(const Double_t min, const Double_t max);
  void     setTrackEtaLimits(const Double_t min, const Double_t max);
  void     setMinNHitsFit(const Int_t min);
  void     setMinNHitsDedx(const Int_t min);
  void     setMinFitHitsFaction(const Double_t min);
  void     setMaxDca(const Double_t max);

  Int_t    getTrackType()         { return mTrackType;         }
  Double_t getMinTrkPt()          { return mMinTrkPt;          }
  Double_t getMaxTrkPt()          { return mMaxTrkPt;          }
  Double_t getMinTrkPhi()         { return mMinTrkPhi;         }
  Double_t getMaxTrkPhi()         { return mMaxTrkPhi;         }
  Double_t getMinTrkEta()         { return mMinTrkEta;         }
  Double_t getMaxTrkEta()         { return mMaxTrkEta;         }
  Int_t    getMinNHitsFit()       { return mMinNHitsFit;       }
  Int_t    getMinNHitsDedx()      { return mMinNHitsDedx;      }
  Double_t getMinFitHitsFaction() { return mMinFitHitsFaction; }
  Double_t getMaxDca()            { return mMaxDca;            }

  double   getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos) const;

  // MTD hit cuts
  void     setTrigTimeCut(const Bool_t cut);
  Bool_t   isTrigTimeCutOn()      { return mTrigTimeCut;        }
  Bool_t   isMtdHitInTrigWin(StMuMtdHit *hit, const Double_t trigger_time) const;
  Bool_t   isMtdHitInTrigWin(StMtdHit *hit, const Double_t trigger_time)   const;
  Bool_t   isMtdHitInTrigWin(Int_t backleg, const Int_t module, const Double_t leading_time, const Double_t trigger_time) const;
  Int_t    getMtdPidTraitsIndex(const StPicoMtdHit *hit, const StPicoDst *pico);
  Int_t    getMtdHitIndex(const StPicoTrack *track, const StPicoDst *pico);
  StMtdPidTraits *getMtdPidTraits(StTrack *track);

  Bool_t isInSameTrigUnit(StMuMtdHit *hit1, StMuMtdHit *hit2) const;
  Bool_t isInSameTrigUnit(StMtdHit *hit1, StMtdHit *hit2) const;
  Bool_t isInSameTrigUnit(StPicoMtdHit *hit1, StPicoMtdHit *hit2) const;
  Bool_t isInSameTrigUnit(Int_t backleg1, Int_t module1, Int_t backleg2, Int_t module2) const;
  Int_t  getTrigUnit(Int_t backleg, Int_t module) const;


  // track pair cuts
  Bool_t   checkTrkPairDca(const Double_t dr, const Double_t dz) const;
  void     setTrkPairDca(const Double_t r, const Double_t z);
  Double_t getMaxTrkPairDcaDr()   { return mMaxTrkPairDcaDr; }
  Double_t getMaxTrkPairDcaDz()   { return mMaxTrkPairDcaDz; }  

  // Tof pid
  double   m2(StTrack *track);
  double   m2(const StMuTrack *track);
  double   m2(StPicoTrack *track, const StPicoDst *picoDst);


  // Muon pid
  Bool_t   isMuonCandidate(const StTrack *track);
  Bool_t   isMuonCandidate(const StMuTrack *track, const StMuDst *muDst, StMtdTrigger *trig = 0x0);
  Bool_t   isMuonCandidate(const StPicoTrack *track, const StPicoDst *picoDst, StMtdTrigger *trig = 0x0);
  Bool_t   isMuonCandidate(const Double_t pt, const Double_t nSigmaPi, const Int_t tofindex, const Double_t dz, const Double_t dy, const Double_t dtof, const Bool_t isTrig);
  Bool_t   checkNSigmaPi(const Double_t nSigmaPi) const;
  Bool_t   checkTrkHitDeltaZ(const Double_t dz, const Double_t pt) const;
  Bool_t   checkTrkHitDeltaY(const Double_t dy, const Double_t pt) const;
  Bool_t   checkTrkHitDeltaTof(const Double_t dtof, const Double_t pt) const;
  void     setNsigmaPiCut(const Double_t min, const Double_t max);
  void     setMtdHitTrigger(const Bool_t require);
  void     setBTofMatch(const Bool_t match);

  void     setApplyDzCut(const bool cut)          { mApplyDzCut = cut;        }
  void     setApplyDyCut(const bool cut)          { mApplyDyCut = cut;        }
  void     setApplyDtofCut(const bool cut)        { mApplyDtofCut = cut;      }
  void     setApplyPtDepCutDz(const Bool_t cut)   { mPtDepCutDz = cut;        }
  void     setApplyPtDepCutDy(const Bool_t cut)   { mPtDepCutDy = cut;        }
  void     setApplyPtDepCutDtof(const Bool_t cut) { mPtDepCutDtof = cut;      }
  void     setSigmaDz(double s1, double s2)       { mSigmaDz1 = s1; mSigmaDz2 = s2; }
  void     setSigmaDy(double s1, double s2)       { mSigmaDy1 = s1; mSigmaDy2 = s2; }
  void     setSigmaDtof(double s)                 { mSigmaDtof = s; }
  void     setMuonDeltaZ(const Double_t min, const Double_t max);
  void     setMuonDeltaY(const Double_t min, const Double_t max);
  void     setMuonDeltaTof(const Double_t min, const Double_t max);
  void     setMaxMuonDca(const double max);
  void     setMuonPt(const Double_t min, const Double_t max);
  double   getPtDepDzCut(double pt);
  double   getPtDepDyCut(double pt);
  Double_t getMuonDeltaZ()        { return mMaxMuonDeltaZ;     }
  Double_t getMuonDeltaY()        { return mMaxMuonDeltaY;     }
  Double_t getMaxMuonDeltaTof()   { return mMaxMuonDeltaTof;   }
  Double_t getMinMuonDeltaTof()   { return mMinMuonDeltaTof;   }
  Double_t getMinNsigmaPi()       { return mMinNsigmaPi;       }
  Double_t getMaxNsigmaPi()       { return mMaxNsigmaPi;       }
  Bool_t   getBTofMatch()         { return mBTofMatch;         }
  Double_t getMaxMuonDca()        { return mMaxMuonDca;        }


  // StPicoTrack utility
  StThreeVectorF gMom(const StPicoTrack *track, const StThreeVectorF pVtx, const float B) const;
  double dca(StPicoTrack *track, StThreeVectorF pVtx) const;
  StPhysicalHelixD helix(const StPicoTrack *track, const float B);

  // multiplicity
  Int_t gRefMult(StEvent *event, StThreeVectorF vtxPos);
  Int_t getTofMult(StMuDst *event);
  Int_t getTofMult(StEvent *event);

 
  // utility functions
  void     addCutToHisto(TH1 *h, const Int_t bin, const char *label, const Float_t value = -999) const;

  Double_t getMtdHitGlobalZ(StMuMtdHit *hit) const;
  Double_t getMtdHitGlobalZ(StMtdHit *hit) const;
  Double_t getMtdHitGlobalZ(StPicoMtdHit *hit) const;
  Double_t getMtdHitGlobalZ(Double_t leadingWestTime, Double_t leadingEastTime, Int_t module) const;

  Double_t getMtdHitLocalZ(StMuMtdHit *hit) const;
  Double_t getMtdHitLocalZ(StMtdHit *hit) const;
  Double_t getMtdHitLocalZ(StPicoMtdHit *hit) const;
  Double_t getMtdHitLocalZ(Double_t leadingWestTime, Double_t leadingEastTime) const;
  double   getMtdHitLocalY(StMuMtdHit *hit) const ;
  double   getMtdHitLocalY(StPicoMtdHit *hit) const ;

  Double_t rotatePhi(Double_t phi) const;
  Bool_t   isValidPosition(const StThreeVectorD pos) const;

  void     printConfig();

 private:
  Int_t            mYear;

  Int_t            mMCflag;                                    // 0 - real data; 1 - MC reco; 2 - MC truth; 3 - embedding

  // event cuts
  Bool_t           mTrigger[kNtrig];                           // 0 - off; 1 - on
  IntVec           mTriggerIDs[kNtrig];                        // Valid trigger id collection that will be tested
  IntVec           mBadRunIDs;

  // centrality
  StRefMultCorr    *mRefMultCorr;                              //
  Int_t            mVzCorrMethod;                              // 0 - Polynomial fit; 1 - interpolation using TGraph::Eval()
  TString          mVzCorrFileName; 
  TGraphErrors     *mgVzCorr;
  TGraphErrors     *mgVzCorr1;

  // vertex cuts
  Bool_t           mRequireVtxRanking;                         //
  Bool_t           mRequireMtdHitForPrimVtx;                   // 
  Bool_t           mUsePrimVtxClosestToVpd;                    //
  Double_t         mMaxVtxZ;                                   // Maximum vertex z
  Double_t         mMaxVtxR;                                   // Maximum vertex r
  Double_t         mMaxDiffz;                                  // Maximum TPC-VPD

  // track cuts
  Bool_t           mTrackType;                                 // 0 - primary tracks; 1 - global tracks
  Double_t         mMinTrkPt;                                  // Minimum track pt
  Double_t         mMaxTrkPt;                                  // Maximum track pt
  Double_t         mMinTrkPhi;                                 // Minimum track phi
  Double_t         mMaxTrkPhi;                                 // Maximum track phi
  Double_t         mMinTrkEta;                                 // Minimum track eta
  Double_t         mMaxTrkEta;                                 // Maximum track eta
  Int_t            mMinNHitsFit;                               // Minimum number of hits used for track fit
  Int_t            mMinNHitsDedx;                              // Minimum number of hits used for de/dx
  Double_t         mMinFitHitsFaction;                         // Minimum fraction of NHitsFit/NHitsPoss
  Double_t         mMaxDca;                                    // Maximum track dc

  // MTD hit cuts
  Bool_t           mTrigTimeCut;                               // flag to apply trigger time window cut
  Double_t         mTrigWinCut_low[kNBackleg][kNModule];
  Double_t         mTrigWinCut_high[kNBackleg][kNModule];

  // muon PID
  Double_t         mMinNsigmaPi;                               // Minimum nsigma for pion assumption
  Double_t         mMaxNsigmaPi;                               // Maximum nsigma for pion assumption
  Bool_t           mBTofMatch;                                 // Flag to require TOF matching
  Bool_t           mMtdHitTrigger;                             // Use hits that fire trigger

  Bool_t           mApplyDzCut; 
  Bool_t           mApplyDyCut;
  Bool_t           mApplyDtofCut;
  Bool_t           mPtDepCutDz;                               // Apply pT-dependent cuts
  Bool_t           mPtDepCutDy;                               // Apply pT-dependent cuts
  Bool_t           mPtDepCutDtof;                             // Apply pT-dependent cuts
  TF1              *fResDzVsPt;
  TF1              *fResDyVsPt;
  TF1              *fResDtofVsPt;
  Double_t         mSigmaDz1;
  Double_t         mSigmaDz2;
  Double_t         mSigmaDy1;
  Double_t         mSigmaDy2;
  Double_t         mSigmaDtof;

  Double_t         mMinMuonDeltaZ;                             // dz cut on matched track-hit pair
  Double_t         mMaxMuonDeltaZ;
  Double_t         mMinMuonDeltaY;
  Double_t         mMaxMuonDeltaY;
  Double_t         mMinMuonDeltaTof;
  Double_t         mMaxMuonDeltaTof;
  Double_t         mMaxMuonDca;
  Double_t         mMinMuonPt;
  Double_t         mMaxMuonPt;

  // track pair cuts
  Double_t         mMaxTrkPairDcaDr;                           // dR of the track pair dca w.r.t. beam line
  Double_t         mMaxTrkPairDcaDz;                           // az of the track pair dca w.r.t. vpd vz


  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $Id: built " __DATE__ " " __TIME__ ; return cvs;
  }
  
  ClassDef(StMtdJpsiUtil, 1)
};

inline void StMtdJpsiUtil::setMCflag(Int_t flag)                      { mMCflag = flag;               }

// event cuts
inline void StMtdJpsiUtil::setYear(const Int_t year)                  { mYear = year;                 }
inline void StMtdJpsiUtil::setTrigger(Bool_t a, Bool_t b, Bool_t c, Bool_t d)   
{
  mTrigger[0] = a; mTrigger[1] = b; mTrigger[2] = c; mTrigger[3] = d;
}

// vertex cuts
inline void StMtdJpsiUtil::setRequireVtxRanking(const Bool_t r)       { mRequireVtxRanking = r;       } 
inline void StMtdJpsiUtil::setRequireMtdHitForPrimVtx(const Bool_t r) { mRequireMtdHitForPrimVtx = r; }
inline void StMtdJpsiUtil::setUsePrimVtxClosestToVpd(const Bool_t r)  { mUsePrimVtxClosestToVpd  = r; }
inline void StMtdJpsiUtil::setMaxVtxZ(const Double_t max)             { mMaxVtxZ = max;               }
inline void StMtdJpsiUtil::setMaxVtxR(const Double_t max)             { mMaxVtxR = max;               }
inline void StMtdJpsiUtil::setMaxDiffz(const Double_t max)            { mMaxDiffz = max;              }

// track cuts
inline void StMtdJpsiUtil::setTrackType(Int_t type)                   { mTrackType = type;            }
inline void StMtdJpsiUtil::setTrackPtLimits(const Double_t min, const Double_t max){
  mMinTrkPt  = min; mMaxTrkPt  = max;
}
inline void StMtdJpsiUtil::setTrackPhiLimits(const Double_t min, const Double_t max){
  mMinTrkPhi = min; mMaxTrkPhi = max;
}
inline void StMtdJpsiUtil::setTrackEtaLimits(const Double_t min, const Double_t max){
  mMinTrkEta = min; mMaxTrkEta = max;
}
inline void StMtdJpsiUtil::setMinNHitsFit(const Int_t min)            { mMinNHitsFit = min;           }
inline void StMtdJpsiUtil::setMinNHitsDedx(const Int_t min)           { mMinNHitsDedx = min;          }
inline void StMtdJpsiUtil::setMinFitHitsFaction(const Double_t min)   { mMinFitHitsFaction = min;     }
inline void StMtdJpsiUtil::setMaxDca(const Double_t max) { mMaxDca = max; }
inline void StMtdJpsiUtil::setNsigmaPiCut(const Double_t min, const Double_t max){
  mMinNsigmaPi = min; mMaxNsigmaPi = max;
}
inline void StMtdJpsiUtil::setBTofMatch(const Bool_t match)           { mBTofMatch = match;           }
inline void StMtdJpsiUtil::setMtdHitTrigger(const Bool_t require)     { mMtdHitTrigger = require;     }

// MTD hit cuts
inline void StMtdJpsiUtil::setTrigTimeCut(const Bool_t cut)           { mTrigTimeCut = cut;           }

// track pair cuts
inline void StMtdJpsiUtil::setMuonDeltaZ(const Double_t min, const Double_t max){
  mMinMuonDeltaZ = min; mMaxMuonDeltaZ = max; }

inline void StMtdJpsiUtil::setMuonDeltaY(const Double_t min, const Double_t max){
  mMinMuonDeltaY = min; mMaxMuonDeltaY = max; }

inline void StMtdJpsiUtil::setMuonDeltaTof(const Double_t min, const Double_t max){
  mMinMuonDeltaTof = min; mMaxMuonDeltaTof = max; }

inline void StMtdJpsiUtil::setMuonPt(const Double_t min, const Double_t max){
  mMinMuonPt = min; mMaxMuonPt = max; }


inline void StMtdJpsiUtil::setTrkPairDca(const Double_t r, const Double_t z) {
  mMaxTrkPairDcaDr = r;
  mMaxTrkPairDcaDz = z;
}

inline void StMtdJpsiUtil::setMaxMuonDca(const double max) { mMaxMuonDca = max; }

#endif
