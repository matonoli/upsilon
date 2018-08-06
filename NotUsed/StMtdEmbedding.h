#ifndef STMTDEMBEDDING_HH
#define STMTDEMBEDDING_HH

/***************************************************************************
 *
 * $Id: StMtdEmbedding.h,v 1.42 2018/01/22 20:19:56 marr Exp $ 
 * StMtdEmbedding
 * Author: Rongrong Ma
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TRandom3;

class StMcEvent;
class StMcTrack;
class StMtdHit;
class StPrimaryVertex;

class StPicoDst;
class StPicoDstMaker;
class StPicoTrack;

#include "StMaker.h"
#include "StMtdJpsi/StMtdJpsiUtil.h"
#include "StAssociationMaker/StAssociationMaker.h"


class StMtdEmbedding : public StMaker {
 public:
  StMtdEmbedding(const Char_t *name = "MtdEmbedding");
  ~StMtdEmbedding();

  Int_t    Init();
  Int_t    InitRun(const Int_t runNumber);
  Int_t    Make();
  Int_t    Finish();

  void     setUseVpdVtx(const bool use)               { mUseVpdVtx = use;             }
  void     setEmbedParticle(const char *name)         { mEmbedParticle = name;        }
  void     setPrintConfig(const Bool_t print = kTRUE) { mPrintConfig = print;         }
  void     setMtdUtil(StMtdJpsiUtil*util)             { mMtdUtil = util;              }
  void     setTpcSmearing(double smear)               { mTpcTrkSmear = smear;         }
  void     setInputFileName(const char *name)         { mInFileName = name;           }
  void     setRunEmbedHadron(const bool run)          { mRunEmbedHadron = run;        }
  void     setRunProjection(const bool run)           { mRunProjection = run;         }
  void     setJpsiM(double min, double max)           { mMinJpsiM=min; mMaxJpsiM=max; }
  void     setSaveTree(bool save)                     { mSaveTree = save;             }
  void     setOutTreeName(const char *name)           { mOutTreeFileName = name;      }
  void     setRunSystematics(const bool run)          { mRunSystematics = run;        }
  void     setRunRealDataQa(const bool run)           { mRunRealDataQa = run;         }
  

 protected:
  void     printConfig();
  Int_t    processStEvent();
  Int_t    processPicoDst();
  Int_t    runEmbedHadron();
  void     runProjection();
  double   getSmearedTrkPt(double pt);

  void     bookHistos();
  void     fillInvMass(TLorentzVector muon1, TLorentzVector muon2, THnSparse *hn[2], double mc_pt);
  Int_t    findMatchedRcTrack(StMcTrack *track);
  Bool_t   isMcMuon(StMcTrack *track);
  Bool_t   isRcMuon(StTrack *track);
  bool     isTpcResp(double pt, double eta, double phi);
  bool     isMtdResp(double pt, int bl, int mod);
  bool     isMtdTrig(double pt);
  bool     passPidCut(double pt);
  Bool_t   isMcMtdHit(StMtdHit *hit);
  Int_t    getTrackNdedx(StTrack *track);

  struct StMtdEmbedData
  {
    int centrality;
    int nTracks;
    int    geantId[1000];
    double mcpt[1000];
    double mcphi[1000];
    double mceta[1000];
    double rcpt[1000];
    double rcphi[1000];
    double rceta[1000];
    double nSigmaPi[1000];
    double dz[1000];
    double dy[1000];
    double dtof[1000];
  };

 private:
  TString          mEmbedParticle;                             // embed particle name
  StEvent          *mStEvent;
  StMcEvent        *mMcEvent;
  StMuDst          *mMuDst;                                    // Pointer to MuDst event
  StPicoDst        *mPicoDst;
  Int_t            mRunYear;                                   // Run year
  Int_t            mRunId;                                     // Run number
  Int_t            mEventId;
  Int_t            mTrgSetup;                                  //
  double           mVtxWeight;
  StPrimaryVertex  *mPriVtx;                          
  StThreeVectorF   mVerPosition; 
  Int_t            mTriggerType;                               // Trigger type: 0 di-muon -- 1 single-muon -- 2 e-mu
  Bool_t           mPrintConfig;                               // Flag to print out task configuration
  Int_t            mTrackType;
  Int_t            mCentrality;
  Int_t            mCentType;                                  // 0: 0-20%; 1: 20-40%; 2: 40-60%
  Int_t            mgRefMult;
  Double_t         mgRefMultCorr;
  Int_t            mTofMult;                                   //
  Int_t            mTofMthTrk;
  Double_t         mZdcRate;
  Double_t         mBbcRate;
  Bool_t           mUseVpdVtx;

  StAssociationMaker *mAssoMaker;                              //!
  rcTrackMapType   *mRcTrackMap;                               //!
  mcTrackMapType   *mMcTrackMap;                               //!
  map<Int_t, Int_t> mRcIndices;
  map<Int_t, Int_t> mMcIndices;

  StMtdJpsiUtil    *mMtdUtil;
  StMtdGeometry    *mMtdGeom;                                  //!

  double           mTpcTrkSmear;
  TString          mInFileName;			
  TH1F             *mhVtxWeight;
  TH2F             *mhTpcCorr;
  TF1              *mhMtdRespEmb[30][5];
  TF1              *mhMtdRespCosmic[30][5];
  TF1              *mhMtdTrigEff;
  TF1              *mhMtdTrigElecEff;
  TF1              *mhMtdDtofEff;
  TH1F             *mhInputJpsiPt[4];

  TRandom3         *mRandom;
  double           mMinJpsiM;
  double           mMaxJpsiM;

  Bool_t           mRunEmbedHadron;
  Bool_t           mRunProjection;
  Bool_t           mRunSystematics;
  Bool_t           mRunRealDataQa;

  Bool_t           mSaveTree;
  TString          mOutTreeFileName;
  TTree            *mOutTree;
  TFile            *mOutTreeFile;
  StMtdEmbedData   mEmbedData;

  // List of histograms
  TH1F             *mhAnalysisCuts;
  TH1F             *mhEventStat;
  TH1F             *mhDataVtxZ[kNtrig];
  TH2F             *mhMcVtxZVsDataVtxZ[kNtrig];
  TH2F             *mhSetupVsMcVtxZ[kNtrig];
  TH1F             *mhgRefMult[kNtrig];  
  TH1F             *mhgRefMultCorr[kNtrig];
  TH1F             *mhCentrality[kNtrig];
  TH2F             *mhTofMthTrksVsZdcRate[kNtrig];
  TH2F             *mhTofMthTrksVzBbcRate[kNtrig];
  TH2F             *mhTofMultVsgRefMult[4];
  TH2F             *mhTofMultVsgRefMultCorr[4];
  TH2F             *mhNgTrkVsCent[4];
  TH2F             *mhZdcRateVsCent[4];
  THnSparse        *mhZdcRateVsTrigSetup[kNtrig];

  // hadron embed
  TH2F             *mhHadronDecay;
  TH3F             *mhHadronMcPt;
  TH3F             *mhHadronMcPtTpc;
  TH3F             *mhHadronMcPtMtd;
  TH3F             *mhHadronMcPtMuon;
  TH3F             *mhHadronRcPtTpc;
  TH3F             *mhHadronRcPtMtd;
  TH3F             *mhHadronRcPtMuon;
  THnSparse        *mhHadronMatch;
  TH2F             *mhNsigmaPiVsPt;

  // projection QA
  TH2F             *mhProjTrack;
  TH2F             *mhMatchTrack;
  THnSparse        *mhMtdMatchEff;

  // Embedding QA
  TH1F             *mhNEmbedJpsi[kNtrig];
  THnSparse        *mhgTrkPtRes[kNtrig];
  THnSparse        *mhpTrkPtRes[kNtrig];
  THnSparse        *mhMcTrkInfo[kNtrig];
  THnSparse        *mhMcTrkInfoTpc[kNtrig];
  THnSparse        *mhMcTrkInfoMtd[kNtrig];
  THnSparse        *mhMcTrkInfoFinal[kNtrig];
  THnSparse        *mhMcTrkPtEff[kNtrig][6];
  THnSparse        *mhRcTrkNsigmaPi[kNtrig];
  THnSparse        *mhQaNHitsPoss[kNtrig];
  TH3F             *mhTrkPtBlModMtd[kNtrig];
  TH3F             *mhTrkPtBlModTrig[kNtrig];
  TH2F             *mhRealHitMap[kNtrig];
  TH2F             *mhMcHitMap[kNtrig];
  TH2F             *mhMcHitMapWithEff[kNtrig];

  // Compare with real data
  TH2F             *mhTrkDca[kNtrig][4];
  TH2F             *mhTrkNHitsFit[kNtrig][4];
  TH2F             *mhTrkNHitsDedx[kNtrig][4];
  TH2F             *mhTrkNHitsFrac[kNtrig][4];
  TH3F             *mhTrkNHitsPoss[kNtrig][4];
  THnSparse        *mhTrkEtaPhi[kNtrig][4];
  TH2F             *mhTrkNSigmaPi[kNtrig][4];
  TH2F             *mhTrkDyVsPt[kNtrig][4];
  TH2F             *mhTrkDzVsPt[kNtrig][4];
  TH2F             *mhTrkDtofVsPt[kNtrig][4];
  TH2F             *mhTrkDtofVsMod[kNtrig][4];

  // Run QA on real data
  THnSparse        *mhQaTrkEtaPhi[kNtrig][3];
  TH2F             *mhQaTrkNHitsFit[kNtrig][3][2];
  TH2F             *mhQaTrkNHitsDedx[kNtrig][3][2];
  TH2F             *mhQaTrkDca[kNtrig][3][2];
  TH2F             *mhQaTrkNSigmaPi[kNtrig][3][2];
  THnSparse        *mhQaTrkEtaPhiWithCuts[kNtrig][3];

  // Jpsi Efficiencies
  THnSparse        *mhJpsiInfo[kNtrig][8][2];
  THnSparse        *mhJpsiMatch[kNtrig];

  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $Id: built " __DATE__ " " __TIME__ ; return cvs;
  }
  
  ClassDef(StMtdEmbedding, 1)
};

#endif
