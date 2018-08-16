#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TBranch;
class TClonesArray;
class StPicoPrescales;
class StEmcGeom;
class StEmcPosition;

const int maxTracks = 100;


class StMyAnalysisMaker : public StMaker {
public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    bool SelectEvent(StPicoEvent* eve);
    bool SelectTrack(StPicoTrack* t);
    bool SelectTrigger(StPicoEvent* eve);
    Float_t SelectPair(StPicoTrack* t1, StPicoTrack* t2);
    void DoRefMultCorr(StPicoEvent* eve);
    bool FillEMCEFF(StPicoTrack* t, char* signs);
    bool GetBEMCdist(StPicoTrack* t, float* d);
    bool FillTree();
    
    void    DeclareHistograms();
    void    WriteHistograms();
    
private:
    
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mEvent;
    StEmcPosition  *mEmcPos;
    StEmcGeom      *geomBEMC;
    TString    mOutName;

    StPicoPrescales *mPrescales;

    // GLOBAL VARIABLES
    int cent9;
    int cent16;
    float wMcorr;
    float grefMultC;
    int evFlag;

    // TREE VARIABLES
    TTree* upsTree;
    Int_t tEventId;
    Int_t tEventRunId;
    Float_t tEventTpcVx;
    Float_t tEventTpcVy;
    Float_t tEventTpcVz;
    Float_t tEventVpdVz;
    Int_t tEventCent9;
    Int_t tEventCent16;
    Int_t tEventRefMult;
    Int_t tEventRefMultC;
    Int_t tNElectrons;
    Float_t tElePt[maxTracks];
    Float_t tEleP[maxTracks];
    Float_t tElePx[maxTracks];
    Float_t tElePy[maxTracks];
    Float_t tElePz[maxTracks];
    Float_t tEleEta[maxTracks];
    Int_t tEleCharge[maxTracks];
    Int_t tEleNHits[maxTracks];
    Float_t tEleNRat[maxTracks];
    Int_t tEleNDedx[maxTracks];
    Float_t tEleDca[maxTracks];
    Float_t tEleDedx[maxTracks];
    Float_t tEleNSigE[maxTracks];
    Float_t tEleE[maxTracks];
    Float_t tEleEcl[maxTracks];
    Float_t tEleZDist[maxTracks];
    Float_t tElePhiDist[maxTracks];
    Float_t tEleR[maxTracks];
    Int_t tEleTrig[maxTracks];

    //Event Quality Histograms
    TH1F* hEventSelection;
    TH1F* hEventCuts;
    TH1F* hEventTrigger;

    TH2F* hEventVzvsVzvpdO;
    TH1F* hEventdVzO;
    TH1F* hEventVrO;
    TH1F* hEventnEmcPidsO;
    
    TH1F* hEventrunId;
    TH1F* hEventeventId;
    TH1F* hEventfillId;
    TH1F* hEventrefMult;
    TH2F* hEventVzvsVzvpd;
    TH1F* hEventdVz;
    TH1F* hEventVr;
    TH2F* hEventnTrigTowers;
    TH2F* hEventnBtowEmc;
    TH2F* hEventnEmcTracks;
    TH2F* hEventnBtowTracks;
    TH1F* hEventAllTracks;
    TH1F* hEventTriggerTracks;
    TH2F* hEventNPrimaries;
    TH2F* hEventNPrimariesCent;
    TH3F* hEventVzvsNPrim;
    TH3F* hEventVzvsNPrimHard;
    TH1F* hEta;

    TH2F* hCorrgrefMult;
    TH2F* hCorrgrefMultWT;
    TH2F* hCorrCentrality;
    TH2F* hCorrCentrality2;
    TH1F* hCorrWeight;

    TH1F* hTrigFlag;
    TH2F* hTrigAdcId;

    TH2F* hBtowAdcId;

    TH2F* hEmcAdcId;
    TH1F* hEmcE;
    TH1F* hEmcE0vE;
    TH2F* hEmczDistphiDist;
    TH1F* hEmcTrackIndex;

    TH1F* hTracksSelection;
    TH2F* hTrackEoverpvsp;
    TH2F* hTrackdEdxvsp;
    TH2F* hTracknSigmaElectronvsp;
    TH2F* hTracknSigmaPionvsp;
    TH3F* hTrackEtaPhiPtG;
    TH3F* hTrackEtaPhiPtP;
    TH1F* hTrackDca;
    TH1F* hTracknHitsRatio;
    TH2F* hTrackAdc0vsTower;
    TH2F* hTrackzDistphiDist;
    TH1F* hTrackE0vE;
    TH1F* hTracknHitsFit;
    TH2F* hTrackDistvE0E;
    TH2F* hTrackEtaPhiD;
    TH1F* hTrackR;
    TH2F* hTrackRvDist;
    TH2F* hTracketaDistzDist;
    TH2F* hTrackphiDistphiDist;

    TH2F* hTrackphiDist2;

    TH2F* hTrackSigmavsEoverp;
    TH2F* hTrackSigmavsE0vE;
    TH2F* hTrackpMomgMom;

    TH2F* hElectronPhivdEdx;
    TH2F* hElectronTowervp;
    TH1F* hElectronPt;
    TH1F* hElectronP;
    TH2F* hElectronzDistphiDist;
    TH1F* hElectronDca;
    TH2F* hElectronDistvE0E;
    TH2F* hElectronEtaPhiD;
    TH1F* hElectronR;
    TH2F* hElectronRvDist;
    TH2F* hElectronetaDistzDist;
    TH2F* hElectronphiDistphiDist;

    TH1F* hPairCuts;
    TH3F* hMotherPtEtaPhi;
    TH2F* hMotherYvEta;
    TH3F* hUpsPtYPhi;

    //Kuba
    TH1F* hFillTree;
    TH1F* hFillTreeElectrons;
    TH1F* hFillTreeNElectrons;
    TH1F* hPidTraitsIndex;
    TH1F* hPidTraitsIndexTree;
    TH2F* hElectronTrigAdcId;
    TH2F* hTrigEtaPhi;
    TH2F* hElectronTrigEtaPhi;
    TH3F* hTrackEtaPhiPtPrimOnly;
    TH2F* hTrigElectronEvsP;
    TH2F* hTrigElectronE0vsP;
    TH2F* hTrigElectronGlEvsP;
    TH2F* hTrigElectronGlE0vsP;
    //;

    TH3F* hIMpp;
    TH3F* hIMmm;
    TH3F* hIMpm;
    TH1F* hIMun;
    TH1F* hIMli;
    TH1F* hIMmix;
    TH2F* hIMmixC;

    //FOR NSIGMAEFF
    TH3F* hElectronnSigmavsP;
    TH3F* hElectrondEdxvsP;
    TH2F* hPionnSigmaEvsP;
    TH2F* hKaonnSigmaEvsP;
    TH2F* hProtonnSigmaEvsP;

    //FOR TRIGEFF
    TH1F* hElectronPtMB;
    TH1F* hElectronPMB;
    TH1F* hElectronPtMBps;
    TH1F* hElectronPMBps;
    TH1F* hElectronPtHT;
    TH1F* hElectronPHT;
    TH1F* hPS9;
    TH1F* hPS10;

    //FOR EMCEFF
    TH3F* hEMCEFF0matchUS;
    TH3F* hEMCEFF0epUS;
    TH3F* hEMCEFF0peUS;
    TH3F* hEMCEFF0zdistUS;
    TH3F* hEMCEFF0phidistUS;

    TH1F* hEMCEFFcUS;
    TH3F* hEMCEFFmatchUS;
    TH3F* hEMCEFFepUS;
    TH3F* hEMCEFFpeUS;
    TH3F* hEMCEFFzdistUS;
    TH3F* hEMCEFFphidistUS;

    TH1F* hEMCEFFcLS;
    TH3F* hEMCEFFmatchLS;
    TH3F* hEMCEFFepLS;
    TH3F* hEMCEFFpeLS;
    TH3F* hEMCEFFzdistLS;
    TH3F* hEMCEFFphidistLS;

    TH2F* hMotherDPhi;
    TH2F* hMotherDPhivM;

    TH3F* hPureDca;

    
    ClassDef(StMyAnalysisMaker, 1)
};

#endif
