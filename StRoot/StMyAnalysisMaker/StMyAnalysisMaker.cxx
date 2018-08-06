// uncomment to select picos version
//#define VERS_P15_LOWMID
//#define VERS_P16_LOWMID
//#define VERS_P15_HIGH
#define VERS_P17

#include "StMyAnalysisMaker.h"
#ifndef VERS_P17
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoEmcTrigger.h"
#include "StRoot/StPicoDstMaker/StPicoBTOWHit.h"
#endif
#ifdef VERS_P17
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#endif

#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TParameter.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StPicoPrescales/StPicoPrescales.h"

#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"

#include <iostream>

//#define MIXEVENTS
//#define NSIGMAEFF
//#define TRIGEFF
//#define WRITETREE
//#define EMCEFF
#define FULLTREE

#ifdef MIXEVENTS
    TTree* mixTree;
    int mtreeEntries = 800;
    int mtreeCounter;
    Float_t mtreeVz;
    Int_t mtreeCentrality;
    TClonesArray mtreeTrigP("TLorentzVector",mtreeEntries);
    TClonesArray mtreeAllP("TLorentzVector",mtreeEntries);

    #include <vector>
    int mixSaved = 400;
    int mixCounter = 0;
    std::vector<TLorentzVector*> allElectrons;
    std::vector<TLorentzVector*> trigElectrons;
    std::vector<Float_t> allElectronsVz;
    std::vector<Float_t> trigElectronsVz;
    std::vector<Int_t> allElectronsC;
    std::vector<Int_t> trigElectronsC;
    #ifdef __MAKECINT__
    #pragma link C++ class vector<TLorentzVector>+;
    #endif
#endif

#ifdef WRITETREE
    TTree* upsMass;
    int utreeEntries = 200;     // max entries per one chain
    int utreeCounter;
    Int_t utreeEventID;
    Int_t utreeCentrality;
    TClonesArray utreeUpsMass("TParameter<float>",utreeEntries);
    TClonesArray utreeUpsSigns("TParameter<int>",utreeEntries);
#endif

const float elmass = 0.000510998910;
StThreeVectorF primVpos;

ClassImp(StMyAnalysisMaker)

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mOutName = outName;
}

//-----------------------------------------------------------------------------
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Init() {

    DeclareHistograms();

    #ifdef FULLTREE
    upsTree = new TTree("upsTree","upsilon event tree");
    upsTree->Branch("tEventId",&tEventId,"tEventId/I");
    upsTree->Branch("tEventRunId",&tEventRunId,"tEventRunId/I");
    upsTree->Branch("tEventTpcVx",&tEventTpcVx,"tEventTpcVx/F");
    upsTree->Branch("tEventTpcVy",&tEventTpcVy,"tEventTpcVy/F");
    upsTree->Branch("tEventTpcVz",&tEventTpcVz,"tEventTpcVz/F");
    upsTree->Branch("tEventVpdVz",&tEventVpdVz,"tEventVpdVz/F");
    upsTree->Branch("tEventCent9",&tEventCent9,"tEventCent9/I");
    upsTree->Branch("tEventCent16",&tEventCent16,"tEventCent16/I");
    upsTree->Branch("tEventRefMult",&tEventRefMult,"tEventRefMult/I");
    upsTree->Branch("tEventRefMultC",&tEventRefMultC,"tEventRefMultC/I");
    upsTree->Branch("tNElectrons",&tNElectrons,"tNElectrons/I");
    upsTree->Branch("tElePt",tElePt,"tElePt[tNElectrons]/F");
    upsTree->Branch("tEleP",tEleP,"tEleP[tNElectrons]/F");
    upsTree->Branch("tElePx",tElePx,"tElePx[tNElectrons]/F");
    upsTree->Branch("tElePy",tElePy,"tElePy[tNElectrons]/F");
    upsTree->Branch("tElePz",tElePz,"tElePz[tNElectrons]/F");
    upsTree->Branch("tEleEta",tEleEta,"tEleEta[tNElectrons]/F");
    upsTree->Branch("tEleCharge",tEleCharge,"tEleCharge[tNElectrons]/I");
    upsTree->Branch("tEleNHits",tEleNHits,"tEleNHits[tNElectrons]/I");
    upsTree->Branch("tEleNRat",tEleNRat,"tEleNRat[tNElectrons]/F");
    upsTree->Branch("tEleNDedx",tEleNDedx,"tEleNDedx[tNElectrons]/I");
    upsTree->Branch("tEleDca",tEleDca,"tEleDca[tNElectrons]/F");
    upsTree->Branch("tEleDedx",tEleDedx,"tEleDedx[tNElectrons]/F");
    upsTree->Branch("tEleNSigE",tEleNSigE,"tEleNSigE[tNElectrons]/F");
    upsTree->Branch("tEleE",tEleE,"tEleE[tNElectrons]/F");
    upsTree->Branch("tEleEcl",tEleEcl,"tEleEcl[tNElectrons]/F");
    upsTree->Branch("tEleZDist",tEleZDist,"tEleZDist[tNElectrons]/F");
    upsTree->Branch("tElePhiDist",tElePhiDist,"tElePhiDist[tNElectrons]/F");
    upsTree->Branch("tEleR",tEleR,"tEleR[tNElectrons]/F");
    upsTree->Branch("tEleTrig",tEleTrig,"tEleTrig[tNElectrons]/I");
    #endif

    #ifdef TRIGEFF
    string psdir("StRoot/run14AuAu200GeVPrescales");
    mPrescales = new StPicoPrescales(psdir);
    hPS9  = new TH1F("hPS9","",mPrescales->numberOfRuns(),0,mPrescales->numberOfRuns());
    hPS10 = new TH1F("hPS10","",mPrescales->numberOfRuns(),0,mPrescales->numberOfRuns());
    mPrescales->fillPrescalesHist(hPS9,9);
    mPrescales->fillPrescalesHist(hPS10,10);
    #endif

    #ifdef MIXEVENTS
    mixTree = new TTree("mixTree","");
    mixTree->Branch("Vz",&mtreeVz,"Vz/f");
    mixTree->Branch("Centrality",&mtreeCentrality,"Centrality/i");
    mixTree->Branch("TrigP",&mtreeTrigP);
    mixTree->Branch("AllP",&mtreeAllP);
    #endif

    #ifdef WRITETREE
    upsMass = new TTree("upsMass","");
    upsMass->Branch("EventID",&utreeEventID,"EventID/i");
    upsMass->Branch("Centrality",&utreeCentrality,"Centrality/i");
    upsMass->Branch("UpsMass",&utreeUpsMass);
    upsMass->Branch("UpsSigns",&utreeUpsSigns);
    #endif

    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Finish() {
    
    if(mOutName!="") {
        TString histFileName = "hist_";
        histFileName += mOutName.Data();
        TFile *fout = new TFile(histFileName,"RECREATE");
        fout->cd();
        WriteHistograms();
        fout->Close();
    }

    #ifdef MIXEVENTS    
    if(mOutName!="") {
        TString mixeFileName = "mixe_";
        mixeFileName += mOutName.Data();
        TFile *mixeFile = new TFile(mixeFileName,"RECREATE");
        mixeFile->cd();
        mixTree->Write();
        mixeFile->Close();  }
    #endif

    #ifdef WRITETREE
    if(mOutName!="") {
        TString massFileName = "upsMass_";
        massFileName += mOutName.Data();
        TFile *massFile = new TFile(massFileName,"RECREATE");
        massFile->cd();
        upsMass->Write();
        massFile->Close();  }
    #endif

    #ifdef FULLTREE
    if(mOutName!="") {
        TString treeFileName = "upsTree_";
        treeFileName += mOutName.Data();
        TFile *treeFile = new TFile(treeFileName,"RECREATE");
        treeFile->cd();
        upsTree->Write();
        treeFile->Close();  }
    #endif

    #ifdef TRIGEFF
    delete mPrescales;
    #endif

    return kStOK;
}
//-----------------------------------------------------------------------------

void StMyAnalysisMaker::DeclareHistograms() {


    hEventSelection         = new TH1F("hEventSelection","",12,-0.5,11.5);
    hEventCuts              = new TH1F("hEventCuts","",8,-0.5,7.5);
    hEventTrigger           = new TH1F("hEventTrigger","",10,-0.5,9.5);

    hEventVzvsVzvpdO        = new TH2F("hEventVzvsVzvpdO","hEventVzvsVzvpdO",100,-100,100,100,-100,100);
    hEventdVzO              = new TH1F("hEventdVzO","hEventdVzO",120,-6,6);
    hEventVrO               = new TH1F("hEventVrO","hEventVrO",80,-4,4);
    hEventnEmcPidsO         = new TH1F("hEventnEmcPidsO","hEventnEmcPidsO",1000,-1,999);
    
    hEventrunId             = new TH1F("hEventrunId","hEventrunId",10000,15070000,15170000);
    hEventeventId           = new TH1F("hEventeventId","hEventeventId",100,0,10000000);
    hEventfillId            = new TH1F("hEventfillId","hEventfillId",100,0,100000); 
    hEventrefMult           = new TH1F("hEventrefMult","hEventrefMult",800,0,800);
    hEventVzvsVzvpd         = new TH2F("hEventVzvsVzvpd","hEventVzvsVzvpd",100,-100,100,100,-100,100);
    hEventdVz               = new TH1F("hEventdVz","hEventdVz",120,-6,6);
    hEventVr                = new TH1F("hEventVr","hEventVr",80,-4,4);
    hEventnTrigTowers       = new TH2F("hEventnTrigTowers","hEventnTrigTowers",20,0,20,20,0,20);
    hEventnBtowEmc          = new TH2F("hEventnBtowEmc","hEventnBtowEmc",200,0,1000,200,0,1000);
    hEventnEmcTracks        = new TH2F("hEventnEmcTracks","hEventnEmcTracks",200,0,1000,200,0,1000);
    hEventnBtowTracks       = new TH2F("hEventnBtowTracks","hEventnBtowTracks",200,0,1000,200,0,1000);
    hEventAllTracks         = new TH1F("hEventAllTracks","hEventAllTracks",20,0,20);
    hEventTriggerTracks     = new TH1F("hEventTriggerTracks","",10,0,10);

    hEventNPrimaries        = new TH2F("hEventNPrimaries","",2000,0,2000,2000,0,2000);
    hEventNPrimariesCent    = new TH2F("hEventNPrimariesCent","",10,0,10,2000,0,2000);

    hEventVzvsNPrim         = new TH3F("hEventVzvsNPrim","",400,-50,50,2000,0,2000,10,0,10);
    hEventVzvsNPrimHard     = new TH3F("hEventVzvsNPrimHard","",400,-50,50,2000,0,2000,10,0,10);
    hEta                    = new TH1F("hEta","",300,-3,3);

    hCorrgrefMult           = new TH2F("hCorrgrefMult","", 400,0,800,400,0,800);
    hCorrgrefMultWT         = new TH2F("hCorrgrefMultWT","", 400,0,800,400,0,800);
    hCorrCentrality         = new TH2F("hCorrCentrality","",9,-0.5,8.5,9,-0.5,8.5);
    hCorrCentrality2        = new TH2F("hCorrCentrality2","",9,-0.5,8.5,9,-0.5,8.5);
    hCorrWeight				= new TH1F("hCorrWeight","",100,0.9,1.1);

    hTrigFlag               = new TH1F("hTrigFlag","hTrigFlag",20,0,20);
    hTrigAdcId              = new TH2F("hTrigAdcId","hTrigAdcId",100,0,100,5000,0,5000);

    hBtowAdcId              = new TH2F("hBtowAdcId","hBtowAdcId",200,0,1000,5000,0,5000);

    hEmcAdcId               = new TH2F("hEmcAdcId","hEmcAdcId",200,0,1000,5000,0,5000);
    hEmcE                   = new TH1F("hEmcE","hEmcE",150,0,15);
    hEmcE0vE                = new TH1F("hEmcE0vE","hEmcE0vE",150,0,1.5);
    hEmczDistphiDist        = new TH2F("hEmczDistphiDist","hEmczDistphiDist",200,-10,10,200,-0.1,0.1);
    hEmcTrackIndex          = new TH1F("hEmcTrackIndex","hEmcTrackIndex",5000,-10,4990);

    hTracksSelection        = new TH1F("hTracksSelection","hTracksSelection",21,-0.5,20.5);
    hTrackEoverpvsp         = new TH2F("hTrackEoverpvsp","hTrackEoverpvsp",450,0,15,200,0,4);
    hTrackdEdxvsp           = new TH2F("hTrackdEdxvsp","hTrackdEdxvsp",450,0,15,250,0,10);
    hTracknSigmaElectronvsp = new TH2F("hTracknSigmaElectronvsp","hTracknSigmaElectronvsp",450,0,15,100,-10,10);
    hTracknSigmaPionvsp     = new TH2F("hTracknSigmaPionvsp","hTracknSigmaPionvsp",450,0,15,100,-10,10);
    hTrackEtaPhiPtG         = new TH3F("hTrackEtaPhiPtG","hTrackEtaPhiPtG",450,0,15,160,-1.3,1.3,100,-3.2,3.2);
    hTrackEtaPhiPtP         = new TH3F("hTrackEtaPhiPtP","hTrackEtaPhiPtP",450,0,15,160,-1.3,1.3,100,-3.2,3.2);
    hTrackDca               = new TH1F("hTrackDca","hTrackDca",100,0,20);
    hTracknHitsRatio        = new TH1F("hTracknHitsRatio","hTracknHitsRatio",100,0,1);
    hTrackAdc0vsTower       = new TH2F("hTrackAdc0vsTower","hTrackAdc0vsTower",200,0,1000,5000,0,5000);
    hTrackzDistphiDist      = new TH2F("hTrackzDistphiDist","hTrackzDistphiDist",200,-10,10,200,-0.1,0.1);
    hTrackE0vE              = new TH1F("hTrackE0vE","hTrackE0vE",150,0,1.5);
    hTracknHitsFit          = new TH1F("hTracknHitsFit","hTracknHitsFit",50,0,50);
    hTrackDistvE0E          = new TH2F("hTrackDistvE0E","",200,0,20,150,0,1.5);
    hTrackEtaPhiD           = new TH2F("hTrackEtaPhiD","",400,-0.1,0.1,400,-0.1,0.1);
    hTrackR                 = new TH1F("hTrackR","",400,0,0.2);
    hTrackRvDist            = new TH2F("hTrackRvDist","",100,0,0.1,100,0,10.0);
    hTracketaDistzDist      = new TH2F("hTracketaDistzDist","",400,-0.1,0.1,400,-10.0,10.0);
    hTrackphiDistphiDist    = new TH2F("hTrackphiDistphiDist","",400,-0.1,0.1,400,-0.1,0.1);

    hTrackphiDist2          = new TH2F("hTrackphiDist2","hTrackphiDist2",200,-0.1,0.1,200,-0.1,0.1);

    hTrackSigmavsEoverp     = new TH2F("hTrackSigmavsEoverp","",100,-10,10,200,0,4);
    hTrackSigmavsE0vE       = new TH2F("hTrackSigmavsE0vE","",100,-10,10,150,0,1.5);
    hTrackpMomgMom          = new TH2F("hTrackpMomgMom","",450,0,15,450,0,15);

    hElectronPhivdEdx		= new TH2F("hElectronPhivdEdx","",600,-3.2,3.2,600,1,7);
    hElectronTowervp		= new TH2F("hElectronTowervp","",150,0,15,5000,0,5000);
    hElectronPt             = new TH1F("hElectronPt","",150,0,15);
    hElectronP              = new TH1F("hElectronP","",150,0,15);
    hElectronzDistphiDist   = new TH2F("hElectronzDistphiDist","",200,-10,10,200,-0.1,0.1);
    hElectronDca            = new TH1F("hElectronDca","hElectronDca",100,0,10);
    hElectronDistvE0E       = new TH2F("hElectronDistvE0E","",200,0,20,150,0,1.5);

    hElectronEtaPhiD        = new TH2F("hElectronEtaPhiD","",200,-0.1,0.1,200,-0.1,0.1);
    hElectronR              = new TH1F("hElectronR","",200,0,0.2);
    hElectronRvDist         = new TH2F("hElectronRvDist","",100,0,0.1,100,0,10.0);
    hElectronetaDistzDist   = new TH2F("hElectronetaDistzDist","",400,-0.1,0.1,400,-10.0,10.0);
    hElectronphiDistphiDist = new TH2F("hElectronphiDistphiDist","",400,-0.1,0.1,400,-0.1,0.1);

    hPairCuts               = new TH1F("hPairCuts","",12,-0.5,11.5);
    hMotherPtEtaPhi         = new TH3F("hMotherPtEtaPhi","",450,0,15,160,-1.3,1.3,100,-3.2,3.2);
    hMotherYvEta            = new TH2F("hMotherYvEta","",160,-1.3,1.3,160,-1.3,1.3);
    hUpsPtYPhi              = new TH3F("hUpsPtYPhi","",200,0,20,160,-1.3,1.3,100,-3.2,3.2);

    hIMpp           		= new TH3F("hIMpp","",200,0,20,10,1,21,11,-1.5,9.5);
    hIMpp->Sumw2();
    hIMmm       		    = new TH3F("hIMmm","",200,0,20,10,1,21,11,-1.5,9.5);
    hIMmm->Sumw2();
    hIMpm    		        = new TH3F("hIMpm","",200,0,20,10,1,21,11,-1.5,9.5);
    hIMpm->Sumw2();
    hIMun   		        = new TH1F("hIMun","",200,0,20);
    hIMun->Sumw2();
    hIMli    		        = new TH1F("hIMli","",200,0,20);
    hIMli->Sumw2();
    hIMmix                  = new TH1F("hIMmix","",200,0,20);
    hIMmix->Sumw2();
    hIMmixC                 = new TH2F("hIMmixC","",200,0,20,10,0,10);

    #ifdef NSIGMAEFF
    //Int_t kmaxT = 5000;
    /*mTree = new TTree("mTree","Electrons");
    TnSigmaE = new TClonesArray("TParameter<float>",kmaxT);
    mTree->Branch("TnSigmaE",&TnSigmaE,64000,99);
    TPt = new TClonesArray("TParameter<float>",kmaxT);
    mTree->Branch("TPt",&TPt,64000,99);*/
    hElectronnSigmavsP      = new TH3F("hElectronnSigmavsP","",150,0,15,300,-10,5,10,0,10);
    hElectrondEdxvsP        = new TH3F("hElectrondEdxvsP","",150,0,15,200,1,5,10,0,10);
    hPionnSigmaEvsP         = new TH2F("hPionnSigmaEvsP","",150,0,15,300,-10,5);
    hKaonnSigmaEvsP         = new TH2F("hKaonnSigmaEvsP","",150,0,15,300,-10,5);
    hProtonnSigmaEvsP         = new TH2F("hProtonnSigmaEvsP","",150,0,15,300,-10,5);
    #endif

    #ifdef TRIGEFF
    hElectronPtMB           = new TH1F("hElectronPtMB","",150,0,15);
    hElectronPMB            = new TH1F("hElectronPMB","",150,0,15);
    hElectronPtMBps         = new TH1F("hElectronPtMBps","",150,0,15);
    hElectronPMBps          = new TH1F("hElectronPMBps","",150,0,15);
    hElectronPtHT           = new TH1F("hElectronPtHT","",150,0,15);
    hElectronPHT            = new TH1F("hElectronPHT","",150,0,15);
    #endif

    #ifdef EMCEFF
    hEMCEFF0matchUS         = new TH3F("hEMCEFF0matchUS","",150,0,15,3,0,3,10,0,10);
    hEMCEFF0epUS            = new TH3F("hEMCEFF0epUS","",150,0,15,200,0,4,10,0,10);
    hEMCEFF0peUS            = new TH3F("hEMCEFF0peUS","",150,0,15,200,0,4,10,0,10);
    hEMCEFF0zdistUS         = new TH3F("hEMCEFF0zdistUS","",150,0,15,200,-10,10,10,0,10);
    hEMCEFF0phidistUS       = new TH3F("hEMCEFF0phidistUS","",150,0,15,200,-0.1,0.1,10,0,10);
    hElectronnSigmavsP      = new TH3F("hElectronnSigmavsP","",150,0,15,300,-10,5,10,0,10);

    hEMCEFFcUS              = new TH1F("hEMCEFFcUS","",15,0,15);
    hEMCEFFmatchUS          = new TH3F("hEMCEFFmatchUS","",150,0,15,3,0,3,10,0,10);
    hEMCEFFepUS             = new TH3F("hEMCEFFepUS","",150,0,15,200,0,4,10,0,10);
    hEMCEFFpeUS             = new TH3F("hEMCEFFpeUS","",150,0,15,200,0,4,10,0,10);
    hEMCEFFzdistUS          = new TH3F("hEMCEFFzdistUS","",150,0,15,200,-10,10,10,0,10);
    hEMCEFFphidistUS        = new TH3F("hEMCEFFphidistUS","",150,0,15,200,-0.1,0.1,10,0,10);

    hEMCEFFcLS              = new TH1F("hEMCEFFcLS","",15,0,15);
    hEMCEFFmatchLS          = new TH3F("hEMCEFFmatchLS","",150,0,15,3,0,3,10,0,10);
    hEMCEFFepLS             = new TH3F("hEMCEFFepLS","",150,0,15,200,0,4,10,0,10);
    hEMCEFFpeLS             = new TH3F("hEMCEFFpeLS","",150,0,15,200,0,4,10,0,10);
    hEMCEFFzdistLS          = new TH3F("hEMCEFFzdistLS","",150,0,15,200,-10,10,10,0,10);
    hEMCEFFphidistLS        = new TH3F("hEMCEFFphidistLS","",150,0,15,200,-0.1,0.1,10,0,10);

    hMotherDPhi             = new TH2F("hMotherDPhi","",150,0,15,100,0,3.2);
    hMotherDPhivM           = new TH2F("hMotherDPhivM","",150,0,15,100,0,3.2);
    hPureDca                = new TH3F("hPureDca","",150,0,15,100,0,10,10,0,10);
    #endif
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
    
    // WRITE HISTOGRAMS    
    hEventSelection->Write();
    hEventCuts->Write();
    hEventTrigger->Write();

    hEventVzvsVzvpdO->Write();
    hEventdVzO->Write();
    hEventVrO->Write();
    hEventnEmcPidsO->Write();
    
    hEventrunId->Write();
    hEventeventId->Write();
    hEventfillId->Write();
    hEventrefMult->Write();
    hEventVzvsVzvpd->Write();
    hEventdVz->Write();
    hEventVr->Write();
    hEventnTrigTowers->Write();
    hEventnBtowEmc->Write();
    hEventnEmcTracks->Write();
    hEventnBtowTracks->Write();
    hEventAllTracks->Write();
    hEventTriggerTracks->Write();
    hEventNPrimaries->Write();
    hEventNPrimariesCent->Write();
    hEventVzvsNPrim->Write();
    hEventVzvsNPrimHard->Write();
    hEta->Write();

    hCorrgrefMult->Write();
    hCorrgrefMultWT->Write();
    hCorrCentrality->Write();
    hCorrCentrality2->Write();
    hCorrWeight->Write();

    hTrigFlag->Write();
    hTrigAdcId->Write();

    hBtowAdcId->Write();

    hEmcAdcId->Write();
    hEmcE->Write();
    hEmcE0vE->Write();
    hEmczDistphiDist->Write();
    hEmcTrackIndex->Write();

    hTracksSelection->Write();
    hTrackEoverpvsp->Write();
    hTrackdEdxvsp->Write();
    hTracknSigmaElectronvsp->Write();
    hTracknSigmaPionvsp->Write();
    hTrackEtaPhiPtG->Write();
    hTrackEtaPhiPtP->Write();
    hTrackDca->Write();
    hTracknHitsRatio->Write();
    hTrackAdc0vsTower->Write();
    hTrackzDistphiDist->Write();
    hTrackE0vE->Write();
    hTracknHitsFit->Write();
    hTrackDistvE0E->Write();
    hTracketaDistzDist->Write();
    hTrackphiDistphiDist->Write();

    hTrackEtaPhiD->Write();
    hTrackR->Write();
    hTrackRvDist->Write();

    hTrackphiDist2->Write();

    hTrackSigmavsEoverp->Write();
    hTrackSigmavsE0vE->Write();
    hTrackpMomgMom->Write();

    hElectronPhivdEdx->Write();
    hElectronTowervp->Write();
    hElectronPt->Write();
    hElectronP->Write();
    hElectronzDistphiDist->Write();
    hElectronDca->Write();
    hElectronDistvE0E->Write();

    hElectronEtaPhiD->Write();
    hElectronR->Write();
    hElectronRvDist->Write();
    hElectronetaDistzDist->Write();
    hElectronphiDistphiDist->Write();

    hPairCuts->Write();
    hMotherPtEtaPhi->Write();
    hMotherYvEta->Write();
    hUpsPtYPhi->Write();

    hIMpp->Write();
    hIMmm->Write();
    hIMpm->Write();
    hIMun->Write();
    hIMli->Write();
    hIMmix->Write();
    hIMmixC->Write();

    #ifdef NSIGMAEFF
    //mTree->Write();
    hElectronnSigmavsP->Write();
    hElectrondEdxvsP->Write();
    hPionnSigmaEvsP->Write();
    hKaonnSigmaEvsP->Write();
    hProtonnSigmaEvsP->Write();
    #endif

    #ifdef TRIGEFF
    hElectronPtMB->Write();
    hElectronPMB->Write();
    hElectronPtMBps->Write();
    hElectronPMBps->Write();
    hElectronPtHT->Write();
    hElectronPHT->Write();
    hPS9->Write();
    hPS10->Write();
    #endif

    #ifdef EMCEFF
    hEMCEFF0matchUS->Write();
    hEMCEFF0epUS->Write();
    hEMCEFF0peUS->Write();
    hEMCEFF0zdistUS->Write();
    hEMCEFF0phidistUS->Write();
    hElectronnSigmavsP->Write();

    hEMCEFFcUS->Write();
    hEMCEFFmatchUS->Write();
    hEMCEFFepUS->Write();
    hEMCEFFpeUS->Write();
    hEMCEFFzdistUS->Write();
    hEMCEFFphidistUS->Write();

    hEMCEFFcLS->Write();
    hEMCEFFmatchLS->Write();
    hEMCEFFepLS->Write();
    hEMCEFFpeLS->Write();
    hEMCEFFzdistLS->Write();
    hEMCEFFphidistLS->Write();
    hMotherDPhi->Write();
    hMotherDPhivM->Write();

    hPureDca->Write();
    #endif

    
}

//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::SelectEvent(StPicoEvent* eve){		// some cuts are already in picos

	hEventCuts->Fill(0);

    //remove this
    //if (eve->runId() != 15094014) return false;
    //if (eve->runId() != 15138012) return false;
    //if (eve->runId() != 15167014) return false;
    //

    if (! SelectTrigger(eve)) return false;
    hEventCuts->Fill(1);

    //if ( fabs( eve->vzVpd() ) > 30) return false;      // NEVER USE THIS
    hEventCuts->Fill(2);

    if ( fabs( eve->primaryVertex().z() ) > 30) return false;
    hEventCuts->Fill(3);

    if ( fabs( eve->primaryVertex().z() - eve->vzVpd() ) > 4) return false;
    hEventCuts->Fill(4);

    return true;
}

//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::SelectTrack(StPicoTrack* t){ 	// some cuts already in picos

	// TPC cuts, dca cut is in Make()
	if (t->nHitsFit() < 20) return false;
    hTracksSelection->Fill(3);

    if ((float)t->nHitsFit()/t->nHitsMax() < 0.52) return false;
    hTracksSelection->Fill(4);

    if (t->nHitsDedx() < 10) return false;
    hTracksSelection->Fill(5);

    if (t->pMom().mag() == 0) return false;	// ->isPrimary()
    hTracksSelection->Fill(6);

    #if !defined(NSIGMAEFF) && !defined(EMCEFF)
    if (t->nSigmaElectron() < -1.5) return false;
    hTracksSelection->Fill(7);
        
    if (t->nSigmaElectron() > 3) return false;
    hTracksSelection->Fill(8);
    #endif

    #ifdef EMCEFF
    if (t->nSigmaElectron() < -0.75) return false;
    hTracksSelection->Fill(7);
        
    if (t->nSigmaElectron() > 3.0) return false;
    hTracksSelection->Fill(8);
    #endif

    if (fabs(t->pMom().pseudoRapidity() ) > 1) return false;
    hTracksSelection->Fill(9);

    #ifdef NSIGMAEFF
    if (fabs(t->nSigmaPion()) < 0.005) hPionnSigmaEvsP->Fill(t->pMom().mag(),t->nSigmaElectron());
    if (fabs(t->nSigmaKaon()) < 0.005) hKaonnSigmaEvsP->Fill(t->pMom().mag(),t->nSigmaElectron());
    if (fabs(t->nSigmaProton()) < 0.005) hProtonnSigmaEvsP->Fill(t->pMom().mag(),t->nSigmaElectron());
    #endif

    #ifndef VERS_P17
    Short_t index = t->emcPidTraitsIndex();
    StPicoEmcPidTraits* emctraits = mPicoDst->emcPidTraits(index);
    #endif
    #ifdef VERS_P17
    Short_t index = t->bemcPidTraitsIndex();
    StPicoBEmcPidTraits* emctraits = mPicoDst->bemcPidTraits(index);
    #endif
    
    #ifndef EMCEFF
    // EMC cuts
    //if ( fabs( emctraits->zDist() ) > 5) return false;        //using R now
    hTracksSelection->Fill(10);

    //if ( fabs( emctraits->phiDist() ) > 0.05) return false;
    hTracksSelection->Fill(11);

    if ( (emctraits->e()/t->pMom().mag() < 0.3) || (emctraits->e()/t->pMom().mag() > 1.8)) return false;
    hTracksSelection->Fill(12);

    if (emctraits->e() < 0.1) return false;
    hTracksSelection->Fill(13);

    //if ( (emctraits->e0()/emctraits->e()) < 0.5) return false;
    hTracksSelection->Fill(18);
    #endif
    
    // momentum cut
    #if !defined(NSIGMAEFF) && !defined(TRIGEFF) && !defined(EMCEFF)
    if (t->pMom().mag() < 3.0) return false;
    hTracksSelection->Fill(14);
    #endif

    #if defined(NSIGMAEFF) || defined(EMCEFF)
    if (t->pMom().mag() < 2.0) return false;
    hTracksSelection->Fill(14);
    #endif

    #ifdef TRIGEFF
    if (t->pMom().perp() < 2.0) return false;
    hTracksSelection->Fill(14);
    #endif

    return true;
}

//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::SelectTrigger(StPicoEvent* eve){

	hEventTrigger->Fill(0);

    #ifndef TRIGEFF
	// P15ic lowmid version
	#ifdef VERS_P15_LOWMID
		int tw = eve->triggerWord();
    	if (tw>>19 & 0x1 || tw>>20 & 0x1) hEventTrigger->Fill(1);
    	if (tw>>21 & 0x1 || tw>>22 & 0x1) hEventTrigger->Fill(2);
    	if (tw>>23 & 0x1 || tw>>24 & 0x1) hEventTrigger->Fill(3);
    	if ( !(tw>>19 & 0x1 || tw>>20 & 0x1) &&  
    		 !(tw>>21 & 0x1 || tw>>22 & 0x1) && 
    		 !(tw>>23 & 0x1 || tw>>24 & 0x1) ) hEventTrigger->Fill(4);
    	if (tw>>1 & 0x1) hEventTrigger->Fill(5);
    	if (tw>>18 & 0x1) hEventTrigger->Fill(6);

    	if ( (tw>>21 & 0x1) || (tw>>22 & 0x1) ) return true;
    	return false;
	#endif

    // P15ic highlumi version
	#ifdef VERS_P15_HIGH
		int tw = eve->triggerWord();
    	if (tw>>0 & 0x1 || tw>>1 & 0x1) hEventTrigger->Fill(1);
    	if (tw>>2 & 0x1 || tw>>3 & 0x1) hEventTrigger->Fill(2);
    	if (tw>>4 & 0x1 || tw>>5 & 0x1) hEventTrigger->Fill(3);
    	if ( !(tw>>0 & 0x1 || tw>>1 & 0x1) &&  
    		 !(tw>>2 & 0x1 || tw>>3 & 0x1) && 
    		 !(tw>>4 & 0x1 || tw>>5 & 0x1) ) hEventTrigger->Fill(4);
    	if (tw>>6 & 0x1) hEventTrigger->Fill(5);
    	if (tw>>7 & 0x1) hEventTrigger->Fill(6);

    	if (tw>>2 & 0x1) hEventTrigger->Fill(7);
    	if (tw>>3 & 0x1) hEventTrigger->Fill(8);		// 450212 should be in P15ic_high

    	if ( (tw>>2 & 0x1) || (tw>>3 & 0x1) ) return true;
    	return false;
	#endif

    // P16id lowmid version
	#ifdef VERS_P16_LOWMID
		if ((eve->isTrigger(450201) || eve->isTrigger(450211) )) hEventTrigger->Fill(1);
    	if ((eve->isTrigger(450202) || eve->isTrigger(450212) )) hEventTrigger->Fill(2);
    	if ((eve->isTrigger(450203) || eve->isTrigger(450213) )) hEventTrigger->Fill(3);
    	if ( !(eve->isTrigger(450201) || eve->isTrigger(450202)) && 
    		 !(eve->isTrigger(450203) || eve->isTrigger(450211)) &&
    		 !(eve->isTrigger(450212) || eve->isTrigger(450213)) ) hEventTrigger->Fill(4);
    	if (eve->isTrigger(450060)) hEventTrigger->Fill(5);
    	if (eve->isTrigger(450103)) hEventTrigger->Fill(6);

    	if ((eve->isTrigger(450202) || eve->isTrigger(450212) )) return true;
    	return false;
	#endif

    #ifdef VERS_P17 
        if ((eve->isTrigger(450201) || eve->isTrigger(450211) )) hEventTrigger->Fill(1);
        if ((eve->isTrigger(450202) || eve->isTrigger(450212) )) hEventTrigger->Fill(2);
        if ((eve->isTrigger(450203) || eve->isTrigger(450213) )) hEventTrigger->Fill(3);
        if ( !(eve->isTrigger(450201) || eve->isTrigger(450202)) && 
             !(eve->isTrigger(450203) || eve->isTrigger(450211)) &&
             !(eve->isTrigger(450212) || eve->isTrigger(450213)) ) hEventTrigger->Fill(4);
        if (eve->isTrigger(450060)) hEventTrigger->Fill(5);
        if (eve->isTrigger(450103)) hEventTrigger->Fill(6);

        if ((eve->isTrigger(450202) || eve->isTrigger(450212) )) return true;
        return false;
    #endif
    #endif

    #ifdef TRIGEFF // look at VPDMB-30 instead
    #ifdef VERS_P15_LOWMID
        int tw = eve->triggerWord();
        if (tw>>9 & 0x1 || tw>>10 & 0x1) hEventTrigger->Fill(1);
        if (tw>>9 & 0x1) hEventTrigger->Fill(2);
        if (tw>>10 & 0x1) hEventTrigger->Fill(3);
        if ((tw>>21 & 0x1) || (tw>>22 & 0x1)) hEventTrigger->Fill(4);

        // set event flag
        evFlag = 0;
        if ((tw>>21 & 0x1) || (tw>>22 & 0x1)) evFlag = 1;
        if (tw>>9 & 0x1) evFlag = 2;
        if (tw>>10 & 0x1) evFlag = 3;

        if ( (tw>>9 & 0x1) || (tw>>10 & 0x1) ||
            (tw>>21 & 0x1) || (tw>>22 & 0x1) ) return true;
        return false;
    #endif

    #ifdef VERS_P16_LOWMID
        if ((eve->isTrigger(450010) || eve->isTrigger(450020) )) hEventTrigger->Fill(1);
        if ((eve->isTrigger(450202) || eve->isTrigger(450212) )) hEventTrigger->Fill(2);

        if ((eve->isTrigger(450010) || eve->isTrigger(450020) )) return true;
        return false;
    #endif    
    #endif
}

//-----------------------------------------------------------------------------

Float_t StMyAnalysisMaker::SelectPair(StPicoTrack* t1, StPicoTrack* t2){

	StThreeVectorF p1 = t1->pMom();         // THIS VERSION USES PMOM, NOT GMOM
    StThreeVectorF p2 = t2->pMom();
    
    TLorentzVector* tl1 = new TLorentzVector(p1.x(),p1.y(),p1.z(),sqrt(p1.mag()*p1.mag() + elmass*elmass));
    TLorentzVector* tl2 = new TLorentzVector(p2.x(),p2.y(),p2.z(),sqrt(p2.mag()*p2.mag() + elmass*elmass));
    
    TLorentzVector* q = new TLorentzVector(*tl1+*tl2);

    hMotherPtEtaPhi->Fill(q->Perp(),q->PseudoRapidity(),q->Phi());
    hMotherYvEta->Fill(q->Rapidity(),q->PseudoRapidity());

    // DAUGHTER PAIR CUTS
    hPairCuts->Fill(0);

	if (q->M() < 0.2) return -99;    
    hPairCuts->Fill(1);

    // likely candidates
    if (q->M() > 9 && q->M() < 11 && t1->charge()*t2->charge() < 0) {
		hUpsPtYPhi->Fill(q->Perp(),q->Rapidity(),q->Phi());
		hPairCuts->Fill(2); }

    if (q->Perp() > 10) return -99;
    hPairCuts->Fill(3);

    if (fabs(q->Rapidity()) > 1) return -99;
    hPairCuts->Fill(4);

    if (t1->pMom().mag() < 4.5 && t2->pMom().mag() < 4.5) return -99; // 1 e needs to have 4.5
    hPairCuts->Fill(5);

    return q->M();

}

//-----------------------------------------------------------------------------

void StMyAnalysisMaker::DoRefMultCorr(StPicoEvent* eve) {

	// P15ic version
    #if defined(VERS_P15_LOWMID) || defined(VERS_P15_HIGH)
    	StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr();   // for VPDMB5 this needs weighting
    	grefmultCorrUtil->init(15075008);
    	grefmultCorrUtil->initEvent(eve->grefMult(), eve->primaryVertex().z(), eve->ZDCx());
    #endif
    
    // P16id version
    #ifdef VERS_P16_LOWMID
    	StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_VpdMB30() ;
		grefmultCorrUtil->init(15075008);
		//grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);			//THESE ONLY FOR VPDMB5
		//grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
    	grefmultCorrUtil->initEvent(eve->grefMult(), eve->primaryVertex().z(), eve->ZDCx());
    #endif

    // P17 version -- needs to be defined! using P16v for now
    #ifdef VERS_P17
        StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_VpdMB30() ;
        grefmultCorrUtil->init(15075008);
        //grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);           //THESE ONLY FOR VPDMB5
        //grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
        grefmultCorrUtil->initEvent(eve->grefMult(), eve->primaryVertex().z(), eve->ZDCx());
    #endif

	//     * 5% increment centrality bins (16 bins)	
    cent16 = grefmultCorrUtil->getCentralityBin16();

    //     * 5% increment in 0-10%, and 10% increment in 10-80% (9 bins)      	
    cent9 = grefmultCorrUtil->getCentralityBin9();

    //     weight for histograms in peripheral        	
    wMcorr = grefmultCorrUtil->getWeight();

    //     corrected grefMult          		
    grefMultC = grefmultCorrUtil->getRefMultCorr();  

}

//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::FillEMCEFF(StPicoTrack* t, char* signs) {
#ifdef EMCEFF
    if (signs != "LS" && signs != "US") {
        cout << "Argument signs of FillEMCEFF must be either LS or US" << endl;
        return false;
    }

    StThreeVectorF p1 = t->pMom();
    Short_t in1 = t->emcPidTraitsIndex();
    Int_t match1 = 1;
    if (in1 < 0) match1 = 0;        
    StPicoEmcPidTraits* emc1 = mPicoDst->emcPidTraits(in1);
    if (! emc1) match1 = 0;
    
    if (signs == "US") {
        hEMCEFF0matchUS->Fill(p1.mag(),match1,cent9);
        if(p1.mag() != 0) hEMCEFF0epUS->Fill(p1.mag(),emc1->e()/p1.mag(),cent9);
        if(emc1->e() != 0) hEMCEFF0peUS->Fill(p1.mag(),p1.mag()/emc1->e(),cent9);
        hEMCEFF0zdistUS->Fill(p1.mag(),emc1->zDist(),cent9);
        hEMCEFF0phidistUS->Fill(p1.mag(),emc1->phiDist(),cent9);
        hElectronnSigmavsP->Fill(p1.mag(),t->nSigmaElectron(),cent9);

        Float_t dca = (t->dcaGeometry().helix().origin() - primVpos ).mag();
        hPureDca->Fill(p1.mag(),dca,cent9);
    }

    if (signs == "US") {
    hEMCEFFcUS->Fill(0);
    hEMCEFFmatchUS->Fill(p1.mag(),match1,cent9);
    if (match1 != 1) return false;
    hEMCEFFcUS->Fill(1);

    if(p1.mag() != 0) hEMCEFFepUS->Fill(p1.mag(),emc1->e()/p1.mag(),cent9);
    if(emc1->e() != 0) hEMCEFFpeUS->Fill(p1.mag(),p1.mag()/emc1->e(),cent9);
    if (p1.mag() == 0 || emc1->e()/p1.mag() < 0.3 || emc1->e()/p1.mag() > 1.8) return false;
    hEMCEFFcUS->Fill(2);

    hEMCEFFzdistUS->Fill(p1.mag(),emc1->zDist(),cent9);
    if (emc1->zDist() < -5 || emc1->zDist() > 5) return false;;
    hEMCEFFcUS->Fill(3);

    hEMCEFFphidistUS->Fill(p1.mag(),emc1->phiDist(),cent9);
    if (emc1->phiDist() < -0.05 || emc1->phiDist() > 0.05) return false;
    hEMCEFFcUS->Fill(4);

    return true; }

    if (signs == "LS") {
    hEMCEFFcLS->Fill(0);
    hEMCEFFmatchLS->Fill(p1.mag(),match1,cent9);
    if (match1 != 1) return false;
    hEMCEFFcLS->Fill(1);

    if(p1.mag() != 0) hEMCEFFepLS->Fill(p1.mag(),emc1->e()/p1.mag(),cent9);
    if(emc1->e() != 0) hEMCEFFpeLS->Fill(p1.mag(),p1.mag()/emc1->e(),cent9);
    if (p1.mag() == 0 || emc1->e()/p1.mag() < 0.3 || emc1->e()/p1.mag() > 1.8) return false;
    hEMCEFFcLS->Fill(2);

    hEMCEFFzdistLS->Fill(p1.mag(),emc1->zDist(),cent9);
    if (emc1->zDist() < -5 || emc1->zDist() > 5) return false;;
    hEMCEFFcLS->Fill(3);

    hEMCEFFphidistLS->Fill(p1.mag(),emc1->phiDist(),cent9);
    if (emc1->phiDist() < -0.05 || emc1->phiDist() > 0.05) return false;
    hEMCEFFcLS->Fill(4);

    return true; }
#endif
    return false;
}

//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::GetBEMCdist(StPicoTrack* t, float* d) {

    //THIS CALCULATES DISTANCE OF TRACK PROJECTION FROM CLUSTER CENTER OF MASS
    StThreeVectorD* positionBEMC = new StThreeVectorD;
    StThreeVectorD* momentumBEMC = new StThreeVectorD;
    float etaCluster = 0, phiCluster = 0;
    float etaTrack = 0, phiTrack = 0;

    #ifndef VERS_P17
    Short_t index = t->emcPidTraitsIndex();
    if (index < 0) return false;        
    StPicoEmcPidTraits* emc = mPicoDst->emcPidTraits(index);          //this accesses the cluster 
    #endif

    #ifdef VERS_P17
    Short_t index = t->bemcPidTraitsIndex();
    if (index < 0) return false;        
    StPicoBEmcPidTraits* emc = mPicoDst->bemcPidTraits(index);          //this accesses the cluster 
    #endif
    
    if (! emc) return false;

    int id1 = emc->btowId();
    if (id1 < 0) return false;      //somehow this can be -1???
    // do cluster by yourself:
    std::vector<int> neighbors;
    for (int deta = -1; deta <= 1; ++deta) {
        for (int dphi = -1; dphi <= 1; ++dphi) {
            if (deta==0 && dphi==0) continue;
            int idn = mEmcPos->getNextTowerId(id1, deta, dphi);
            if (idn) neighbors.push_back(idn);
        }
    }
    int id2 = 0;//emc->btowId2();
    int id3 = 0;//emc->btowId3(); 
    float e2 = 0;
    float e3 = 0;

    int nBtow = mPicoDst->numberOfBTOWHits();
    for (int i = 0; i < neighbors.size(); ++i)
    {
        float Ecur = 0;
        for (int j = 0; j < nBtow; ++j)
        {
            StPicoBTOWHit* bhit = mPicoDst->btowHit(j);
            if (bhit->id() != neighbors[i]) continue;
            Ecur = bhit->energy();
            break;
        }
        if (Ecur > e2) {        // form from 2 highest, not closest
            e3 = e2;
            id3 = id2;
            e2 = Ecur;
            id2 = neighbors[i]; }
        else if (Ecur > e3) {
            e3 = Ecur;
            id3 = neighbors[i]; }
    }
    float eta1 = 0, eta2 = 0, eta3 = 0;
    float phi1 = 0, phi2 = 0, phi3 = 0;
    float clx = 0, cly = 0;

    geomBEMC->getEta(id1,eta1);
    geomBEMC->getPhi(id1,phi1);
    if (id2>0) {
        geomBEMC->getEta(id2,eta2);
        geomBEMC->getPhi(id2,phi2); }
    if (id3>0) {
        geomBEMC->getEta(id3,eta3);
        geomBEMC->getPhi(id3,phi3); }

    float ecluster = emc->e1() + e2 + e3;
    if (ecluster==0) return false;
    etaCluster = (eta1*emc->e1()+eta2*e2+eta3*e3)/ecluster;
    clx = (cos(phi1)*emc->e1()+cos(phi2)*e2+cos(phi3)*e3)/ecluster;
    cly = (sin(phi1)*emc->e1()+sin(phi2)*e2+sin(phi3)*e3)/ecluster;
    phiCluster = atan2(cly,clx);

    StPhysicalHelixD th = t->helix();
    bool okBEMC = mEmcPos->projTrack(positionBEMC,momentumBEMC,&th,mEvent->bField(),geomBEMC->Radius());
    if (!okBEMC) return false;
    etaTrack = positionBEMC->pseudoRapidity();
    phiTrack = positionBEMC->phi();

    float dphi = phiTrack - phiCluster;
    if(dphi>=TMath::Pi()) dphi=dphi-TMath::TwoPi();
    if(dphi<-TMath::Pi()) dphi=dphi+TMath::TwoPi();

    d[0] = etaTrack - etaCluster;
    d[1] = dphi;

    /*cout << "= ids are " << id1 << " " << id2 << " " << id3 << endl;
    cout << "= etas are " << eta1 << " " << eta2 << " " << eta3 << endl;
    cout << "= calc: " << eta1 << "*" << emc->e1() << " + " << eta2 << "*" << emc->e2() << " + " << eta3 << "*" << emc->e3() << " / " << emc->e() << endl;
    cout << "= NEW calc: " << eta1 << "*" << emc->e1() << " + " << eta2 << "*" << e2 << " + " << eta3 << "*" << e3 << " / " << ecluster << endl;
    cout << "= center eta is " << etaCluster << " projected eta is " << etaTrack << endl;
    cout << "= ev is " << mEvent->eventId() << endl;
    cout << "= #clusters is " << mPicoDst->numberOfEmcPidTraits() << endl;
    cout << "=1= track id is " << t->id() << " emctraits id is " << emc->bemcId() << endl;  
    cout << "=1= eta is " << d[0] << " phi is " << d[1] << endl;*/

    //cout << "distt " << d[0] << " and " << d[1] << endl;
    /*delete geomBEMC;
    delete mEmcPos;
    delete positionBEMC;
    delete momentumBEMC;*/
    return true;
}

//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::FillTree() {

    // select event
    if (! SelectTrigger(mEvent)) return false;
    if ( fabs( mEvent->primaryVertex().z() ) > 100) return false;
    if ( fabs( mEvent->primaryVertex().z() - mEvent->vzVpd() ) > 4) return false;

    #ifndef VERS_P17
    int nEmcPids = mPicoDst->numberOfEmcPidTraits();
    #endif
    #ifdef VERS_P17
    int nEmcPids = mPicoDst->numberOfBEmcPidTraits();
    #endif

    if ( nEmcPids < 1) return false;

    // find triggered tower
    vector<int> TrigTowers;
    vector<int> TrigTowersAdcOnly;
    int nTrig = mPicoDst->numberOfEmcTriggers();
    int nBtow = mPicoDst->numberOfBTOWHits();
    for (int iTow = 0; iTow < nBtow; iTow++)
    {
        StPicoBTOWHit* bhit = mPicoDst->btowHit(iTow);
        if (bhit->adc()>>4 < 19) continue;
        TrigTowersAdcOnly.push_back(iTow);
        for (int iTrig = 0; iTrig < nTrig; iTrig++)
        {
            StPicoEmcTrigger* trig = mPicoDst->emcTrigger(iTrig);
            if (bhit->id() == trig->id())   TrigTowers.push_back(iTow);
        }
    }
    if (TrigTowers.size() < 1 && TrigTowersAdcOnly.size() < 1) return false;

    //find electrons    
    int nTracks = mPicoDst->numberOfTracks();
    int nElectrons=0;
    for (int iTrk = 0; iTrk < nTracks; iTrk++)
    {
        StPicoTrack* t = mPicoDst->track(iTrk);
        if (! t) continue;
        #ifndef VERS_P17
        Short_t index = t->emcPidTraitsIndex();
        if (index < 0) return false;        
        StPicoEmcPidTraits* emctraits = mPicoDst->emcPidTraits(index);          //this accesses the cluster 
        #endif
        #ifdef VERS_P17
        Short_t index = t->bemcPidTraitsIndex();
        if (index < 0) return false;        
        StPicoBEmcPidTraits* emctraits = mPicoDst->bemcPidTraits(index);          //this accesses the cluster 
        #endif
        if (! emctraits) continue;
        
        float etaphi[2];
        etaphi[0]=999, etaphi[1]=999;
        bool getBEMC = GetBEMCdist(t, etaphi);
        Float_t dca = t->helix().distance(primVpos);

        if (t->nHitsFit() < 10)             continue;
        if ((float)t->nHitsFit()/t->nHitsMax() < 0.52) continue;
        if (t->nHitsDedx() < 10)            continue;
        if (t->pMom().mag() == 0)           continue; // ->isPrimary()
        if (t->nSigmaElectron() < -2)       continue;
        if (t->nSigmaElectron() > 3.5)      continue;
        if (fabs(t->pMom().pseudoRapidity() ) > 1.1) continue;
        if (dca > 3)                        continue;
        if (emctraits->e()/t->pMom().mag() < 0.2) continue;
        if (emctraits->e()/t->pMom().mag() > 1.9) continue;
        if (emctraits->e() < 0.1)           continue;
        if (sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1]) > 0.06) continue;
        if (t->pMom().mag() < 2.0)          continue;
    
        //is trigger?
        int isTrigger = 0;
        for (int iTrg = 0; iTrg < TrigTowers.size(); iTrg++)
        {
            StPicoBTOWHit* bhit = mPicoDst->btowHit(TrigTowers[iTrg]);
            if ( fabs(emctraits->e1() - bhit->energy() ) > 0.01 &&
                 fabs(emctraits->e0() - bhit->energy() ) > 0.01 ) continue;
            isTrigger = 1;  }
        for (int iTrg = 0; iTrg < TrigTowersAdcOnly.size(); iTrg++)
        {
            StPicoBTOWHit* bhit = mPicoDst->btowHit(TrigTowersAdcOnly[iTrg]);
            if ( fabs(emctraits->e1() - bhit->energy() ) > 0.01 &&
                 fabs(emctraits->e0() - bhit->energy() ) > 0.01 ) continue;
            isTrigger = isTrigger ? 3 : 2;  }

        tElePt[nElectrons]          = t->pMom().perp();
        tEleP[nElectrons]           = t->pMom().mag();
        tElePx[nElectrons]          = t->pMom().x();
        tElePy[nElectrons]          = t->pMom().y();
        tElePz[nElectrons]          = t->pMom().z();
        tEleEta[nElectrons]         = t->pMom().pseudoRapidity();
        
        tEleCharge[nElectrons]      = t->charge();
        tEleNHits[nElectrons]       = t->nHitsFit();
        tEleNRat[nElectrons]        = (float)t->nHitsFit()/t->nHitsMax();
        tEleNDedx[nElectrons]       = t->nHitsDedx();
        tEleDca[nElectrons]         = dca;
        tEleDedx[nElectrons]        = t->dEdx();
        tEleNSigE[nElectrons]       = t->nSigmaElectron();
        
        tEleE[nElectrons]           = emctraits->e0();
        tEleEcl[nElectrons]         = emctraits->e();
        tEleZDist[nElectrons]       = emctraits->zDist();
        tElePhiDist[nElectrons]     = emctraits->phiDist();
        tEleR[nElectrons]           = sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1]);
        tEleTrig[nElectrons]        = isTrigger;

        nElectrons++;
    }
    if (nElectrons < 2) return false; 
    
    tEventId        = mEvent->eventId();
    tEventRunId     = mEvent->runId();
    tEventTpcVx     = mEvent->primaryVertex().x();
    tEventTpcVy     = mEvent->primaryVertex().y();
    tEventTpcVz     = mEvent->primaryVertex().z();
    tEventVpdVz     = mEvent->vzVpd();
    tEventCent9     = cent9;
    tEventCent16    = cent16;
    tEventRefMult   = mEvent->grefMult();
    tEventRefMultC  = grefMultC;
    tNElectrons     = nElectrons;

    return true;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::Clear(Option_t *opt) {
}


//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Make() {
    
    hEventSelection->Fill(0); // counts how many times make is called
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    
    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    mEvent = mPicoDst->event();
    if (!mEvent){
        LOG_WARN << "No PicoEvent! Skip!" << endm;
        return kStWarn;
    }
    
    // CENTRALITY CORRECTION AND DEFINITION
    DoRefMultCorr(mEvent);

    // PRESELECTION HISTOGRAMS
    hEventVzvsVzvpdO->Fill(mEvent->vzVpd(),mEvent->primaryVertex().z());
    hEventdVzO->Fill(mEvent->primaryVertex().z() - mEvent->vzVpd());    
    hEventVrO->Fill(mEvent->primaryVertex().perp());
    #ifndef VERS_P17 
    hEventnEmcPidsO->Fill(mPicoDst->numberOfEmcPidTraits());
    #endif
    #ifdef VERS_P17
    hEventnEmcPidsO->Fill(mPicoDst->numberOfBEmcPidTraits());
    #endif

    //  
    primVpos = mEvent->primaryVertex();
    geomBEMC = StEmcGeom::getEmcGeom("bemc");
    mEmcPos = new StEmcPosition;

    #ifdef FULLTREE
    if (FillTree()) upsTree->Fill();
    #endif

    // SELECT GOOD EVENT INCLUDING TRIGGER
    if (! SelectEvent(mEvent) ) return kStOk;
    #ifndef VERS_P17
    int nEmcPids = mPicoDst->numberOfEmcPidTraits();
    #endif
    #ifdef VERS_P17
    int nEmcPids = mPicoDst->numberOfBEmcPidTraits();
    #endif
    #ifndef EMCEFF
    if ( nEmcPids < 1) return kStOk;
    #endif
    hEventCuts->Fill(5);
    hEventSelection->Fill(1);

    // POSTSELECTION HISTOGRAMS
    hEventrunId->Fill(mEvent->runId());
    hEventeventId->Fill(mEvent->eventId());
    hEventfillId->Fill(mEvent->fillId());
    hEventrefMult->Fill(mEvent->grefMult());
    hEventVzvsVzvpd->Fill(mEvent->vzVpd(),mEvent->primaryVertex().z());
    hEventdVz->Fill(mEvent->primaryVertex().z() - mEvent->vzVpd());    
    hEventVr->Fill(mEvent->primaryVertex().perp());
    Int_t nPrimaries = 0;
    for (int iTr = 0; iTr < mPicoDst->numberOfTracks(); ++iTr)
    {
        StPicoTrack* t = mPicoDst->track(iTr);
        if (! t) continue;
        if (t->pMom().mag() == 0) continue;
        nPrimaries++;
    }
    hEventNPrimaries->Fill(mEvent->grefMult(),nPrimaries);
    hEventNPrimariesCent->Fill(cent9,nPrimaries);     

    // DEFINE PRESCALE FACTOR
    #ifdef TRIGEFF
    float psF;
    switch (evFlag) {
        case 1: psF = 1;
        break;
        case 2: psF = mPrescales->prescale(mEvent->runId(),9);
        break;
        case 3: psF = mPrescales->prescale(mEvent->runId(),10);
        break;
        default: psF = 1;
        break; }
    #endif


    // CHECK MULT CORRECTION--------------------------------
    hCorrgrefMult->Fill(mEvent->grefMult(),grefMultC);
    hCorrgrefMultWT->Fill(mEvent->grefMult(),grefMultC,wMcorr);
    hCorrWeight->Fill(wMcorr);

    #if defined(VERS_P15_LOWMID)
    if ( mEvent->grefMult()>441 ) hCorrCentrality->Fill(8,cent9);
    if ( mEvent->grefMult()>374 && mEvent->grefMult()<442 ) hCorrCentrality->Fill(7,cent9);
    if ( mEvent->grefMult()>263 && mEvent->grefMult()<375 ) hCorrCentrality->Fill(6,cent9);
    if ( mEvent->grefMult()>179 && mEvent->grefMult()<264 ) hCorrCentrality->Fill(5,cent9);
    if ( mEvent->grefMult()>116 && mEvent->grefMult()<180 ) hCorrCentrality->Fill(4,cent9);
    if ( mEvent->grefMult()>71 && mEvent->grefMult()<117 ) hCorrCentrality->Fill(3,cent9);
    if ( mEvent->grefMult()>40 && mEvent->grefMult()<72 ) hCorrCentrality->Fill(2,cent9);
    if ( mEvent->grefMult()>21 && mEvent->grefMult()<41 ) hCorrCentrality->Fill(1,cent9);
    if ( mEvent->grefMult()>10 && mEvent->grefMult()<22 ) hCorrCentrality->Fill((double)0,cent9);

    if ( grefMultC>441 ) hCorrCentrality2->Fill(8,cent9);
    if ( grefMultC>374 && grefMultC<442 ) hCorrCentrality2->Fill(7,cent9);
    if ( grefMultC>263 && grefMultC<375 ) hCorrCentrality2->Fill(6,cent9);
    if ( grefMultC>179 && grefMultC<264 ) hCorrCentrality2->Fill(5,cent9);
    if ( grefMultC>116 && grefMultC<180 ) hCorrCentrality2->Fill(4,cent9);
    if ( grefMultC>71 && grefMultC<117 ) hCorrCentrality2->Fill(3,cent9);
    if ( grefMultC>40 && grefMultC<72 ) hCorrCentrality2->Fill(2,cent9);
    if ( grefMultC>21 && grefMultC<41 ) hCorrCentrality2->Fill(1,cent9);
    if ( grefMultC>10 && grefMultC<22 ) hCorrCentrality2->Fill((double)0,cent9);
    #endif

    #ifdef VERS_P15_HIGH
    if ( mEvent->grefMult()>467 ) hCorrCentrality->Fill(8,cent9);
    if ( mEvent->grefMult()>392 && mEvent->grefMult()<468 ) hCorrCentrality->Fill(7,cent9);
    if ( mEvent->grefMult()>275 && mEvent->grefMult()<393 ) hCorrCentrality->Fill(6,cent9);
    if ( mEvent->grefMult()>188 && mEvent->grefMult()<276 ) hCorrCentrality->Fill(5,cent9);
    if ( mEvent->grefMult()>123 && mEvent->grefMult()<189 ) hCorrCentrality->Fill(4,cent9);
    if ( mEvent->grefMult()>75 && mEvent->grefMult()<124 ) hCorrCentrality->Fill(3,cent9);
    if ( mEvent->grefMult()>43 && mEvent->grefMult()<76 ) hCorrCentrality->Fill(2,cent9);
    if ( mEvent->grefMult()>22 && mEvent->grefMult()<44 ) hCorrCentrality->Fill(1,cent9);
    if ( mEvent->grefMult()>10 && mEvent->grefMult()<23 ) hCorrCentrality->Fill((double)0,cent9);

    if ( grefMultC>467 ) hCorrCentrality2->Fill(8,cent9);
    if ( grefMultC>392 && grefMultC<468 ) hCorrCentrality2->Fill(7,cent9);
    if ( grefMultC>275 && grefMultC<393 ) hCorrCentrality2->Fill(6,cent9);
    if ( grefMultC>188 && grefMultC<276 ) hCorrCentrality2->Fill(5,cent9);
    if ( grefMultC>123 && grefMultC<189 ) hCorrCentrality2->Fill(4,cent9);
    if ( grefMultC>75 && grefMultC<124 ) hCorrCentrality2->Fill(3,cent9);
    if ( grefMultC>43 && grefMultC<76 ) hCorrCentrality2->Fill(2,cent9);
    if ( grefMultC>22 && grefMultC<44 ) hCorrCentrality2->Fill(1,cent9);
    if ( grefMultC>10 && grefMultC<23 ) hCorrCentrality2->Fill((double)0,cent9);
    #endif

    #ifdef VERS_P16_LOWMID
    if ( mEvent->grefMult()>447 ) hCorrCentrality->Fill(8,cent9);
    if ( mEvent->grefMult()>379 && mEvent->grefMult()<448 ) hCorrCentrality->Fill(7,cent9);
    if ( mEvent->grefMult()>268 && mEvent->grefMult()<380 ) hCorrCentrality->Fill(6,cent9);
    if ( mEvent->grefMult()>183 && mEvent->grefMult()<269 ) hCorrCentrality->Fill(5,cent9);
    if ( mEvent->grefMult()>120 && mEvent->grefMult()<184 ) hCorrCentrality->Fill(4,cent9);
    if ( mEvent->grefMult()>74 && mEvent->grefMult()<121 ) hCorrCentrality->Fill(3,cent9);
    if ( mEvent->grefMult()>42 && mEvent->grefMult()<75 ) hCorrCentrality->Fill(2,cent9);
    if ( mEvent->grefMult()>22 && mEvent->grefMult()<43 ) hCorrCentrality->Fill(1,cent9);
    if ( mEvent->grefMult()>10 && mEvent->grefMult()<23 ) hCorrCentrality->Fill((double)0,cent9);

    if ( grefMultC>447 ) hCorrCentrality2->Fill(8,cent9);
    if ( grefMultC>379 && grefMultC<448 ) hCorrCentrality2->Fill(7,cent9);
    if ( grefMultC>268 && grefMultC<380 ) hCorrCentrality2->Fill(6,cent9);
    if ( grefMultC>183 && grefMultC<269 ) hCorrCentrality2->Fill(5,cent9);
    if ( grefMultC>120 && grefMultC<184 ) hCorrCentrality2->Fill(4,cent9);
    if ( grefMultC>74 && grefMultC<121 ) hCorrCentrality2->Fill(3,cent9);
    if ( grefMultC>42 && grefMultC<75 ) hCorrCentrality2->Fill(2,cent9);
    if ( grefMultC>22 && grefMultC<43 ) hCorrCentrality2->Fill(1,cent9);
    if ( grefMultC>10 && grefMultC<23 ) hCorrCentrality2->Fill((double)0,cent9);
    #endif

    #ifdef WRITETREE
    // WRITE EVENT INTO TO TREE-----------------------------
    utreeUpsMass.Clear();
    utreeUpsSigns.Clear();
    utreeCounter = 0;
    utreeEventID = mEvent->eventId();
    utreeCentrality = cent9;
    #endif

    #ifdef MIXEVENTS
    mtreeTrigP.Clear();
    mtreeAllP.Clear();
    mtreeCounter = 0;
    mtreeVz = mEvent->primaryVertex().z();
    mtreeCentrality = cent9;
    #endif

    // QUANTIFY TRIGGERS------------------------------------
    int nTrig = mPicoDst->numberOfEmcTriggers();
    for (int i = 0; i < nTrig; i++)
    {
        StPicoEmcTrigger* trig = mPicoDst->emcTrigger(i);
        hTrigFlag->Fill(trig->flag());
        hTrigAdcId->Fill(trig->adc(), trig->id());
    }
    //------------------------------------------------------

    // QUANTIFY BTOWERS-------------------------------------
    int nBtow = mPicoDst->numberOfBTOWHits();
    //cout << "-----EVENT ID " << mEvent->eventId() << endl;
    for (int i = 0; i < nBtow; i++)
    {
        StPicoBTOWHit* bhit = mPicoDst->btowHit(i);
        hBtowAdcId->Fill(bhit->adc(), bhit->id());
        //printf("--BTOW ID %i -- adc  -- %i -- energy -- %f \n",bhit->id(), bhit->adc(), bhit->energy());
    }
    //------------------------------------------------------

    // QUANTIFY EMCPIDS--------------------------------------
    for (int i = 0; i < nEmcPids; i++)
    {
        #ifndef VERS_P17
        StPicoEmcPidTraits* emctraits = mPicoDst->emcPidTraits(i);
        #endif
        #ifdef VERS_P17
        StPicoBEmcPidTraits* emctraits = mPicoDst->bemcPidTraits(i);
        #endif
        hEmcAdcId->Fill(emctraits->adc0(),emctraits->btowId());
        hEmcE->Fill(emctraits->e());
        if (emctraits->e() != 0) hEmcE0vE->Fill(emctraits->e0()/emctraits->e());
        hEmczDistphiDist->Fill(emctraits->zDist(),emctraits->phiDist());
        hEmcTrackIndex->Fill(emctraits->trackIndex());  // how many clusters dont have associated tracks ?
    }
    //-------------------------------------------------------


    // TRIG TOWER SELECTION----------------------------------
    vector<int> TrigTowers1;
    for (int i = 0; i < nBtow; i++)
    {
        StPicoBTOWHit* bhit = mPicoDst->btowHit(i);
        if (bhit->adc()>>4 < 19) continue;
        for (int j = 0; j < nTrig; j++)
        {
            StPicoEmcTrigger* trig = mPicoDst->emcTrigger(j);
            if (bhit->id() == trig->id())
            {
                TrigTowers1.push_back(i);
            }
        }
    }
    hEventnTrigTowers->Fill(nTrig,TrigTowers1.size());	//needs change to 2D
    //-------------------------------------------------------

    #ifndef TRIGEFF
    if ( TrigTowers1.size() < 1 )
    {
        return kStOk;
    }
    #endif
    hEventSelection->Fill(2);
    


    // ELECTRONS SELECTION-----------------------------------
    int nTracks = mPicoDst->numberOfTracks();
    StThreeVectorF pvtx = mEvent->primaryVertex();
    float bfield = mEvent->bField();

    vector<int> AllTracks1;
    vector<int> TriggerTracks1;

    hEventnBtowEmc->Fill(nBtow,nEmcPids);
    hEventnEmcTracks->Fill(nEmcPids,nTracks);
    hEventnBtowTracks->Fill(nBtow,nTracks);

    #ifdef NSIGMAEFF
    //int fkTr = 0;
    #endif

    int nPrim = 0;
    int nPrimHard = 0;

    for (int i = 0; i < nTracks; i++)
    {
        hTracksSelection->Fill(0);
        StPicoTrack* t = mPicoDst->track(i);
        if (! t) continue;
        hTracksSelection->Fill(1);

        if (t->pMom().mag() > 0 && fabs(t->pMom().pseudoRapidity())<2.0 && t->nSigmaElectron()>0) nPrim++;
        if (t->pMom().mag() > 3.5 && fabs(t->pMom().pseudoRapidity())<2.0 && t->nSigmaElectron()>0) nPrimHard++;
        if (t->pMom().mag() > 0) hEta->Fill(t->pMom().pseudoRapidity());


        #ifndef EMCEFF
        #ifndef VERS_P17
        Short_t index = t->emcPidTraitsIndex();
        if (index < 0) continue;        
        StPicoEmcPidTraits* emctraits = mPicoDst->emcPidTraits(index);          //this accesses the cluster 
        #endif
        #ifdef VERS_P17
        Short_t index = t->bemcPidTraitsIndex();
        if (index < 0) continue;        
        StPicoBEmcPidTraits* emctraits = mPicoDst->bemcPidTraits(index);          //this accesses the cluster 
        #endif
        if (! emctraits) continue;
        #endif
        hTracksSelection->Fill(2);

        // TRACK HISTOGRAMS
        Float_t dca = (t->dcaGeometry().helix().origin() - primVpos ).mag();
        hTrackdEdxvsp->Fill(t->pMom().mag(),t->dEdx());
        hTracknSigmaElectronvsp->Fill(t->pMom().mag(),t->nSigmaElectron());
        hTracknSigmaPionvsp->Fill(t->pMom().mag(),t->nSigmaPion());
        hTrackEtaPhiPtG->Fill(t->gPt(),t->gMom(pvtx, bfield).pseudoRapidity(),t->gMom(pvtx, bfield).phi());
        hTrackEtaPhiPtP->Fill(t->pMom().perp(),t->pMom().pseudoRapidity(),t->pMom().phi());
        hTrackDca->Fill(dca);
        hTracknHitsRatio->Fill((float)t->nHitsFit()/t->nHitsMax());
        hTracknHitsFit->Fill(t->nHitsFit());

        #ifndef EMCEFF
        if (t->pMom().mag() != 0) hTrackEoverpvsp->Fill(t->pMom().mag(), emctraits->e()/t->pMom().mag());
        hTrackAdc0vsTower->Fill(emctraits->adc0(),emctraits->btowId());
        hTrackzDistphiDist->Fill(emctraits->zDist(),emctraits->phiDist());
        if (emctraits->e() != 0) hTrackE0vE->Fill(emctraits->e0()/emctraits->e());
        hTrackphiDist2->Fill(emctraits->phiDist(),emctraits->phiTowDist());

        if (t->pMom().mag() != 0) hTrackSigmavsEoverp->Fill(t->nSigmaElectron(),emctraits->e()/t->pMom().mag());
        if (emctraits->e() != 0) hTrackSigmavsE0vE->Fill(t->nSigmaElectron(),emctraits->e0()/emctraits->e());
        if (emctraits->e() != 0) hTrackDistvE0E->Fill(emctraits->zDist()*emctraits->zDist()+10000*emctraits->phiDist()*emctraits->phiDist(),emctraits->e0()/emctraits->e());
        float etaphi[2];
        etaphi[0]=999, etaphi[1]=999;

        /*cout << "============" << endl;
        cout << "= ev is " << mEvent->eventId() << endl;
        cout << "= #clusters is " << mPicoDst->numberOfEmcPidTraits() << endl;
        cout << "=1= track id is " << t->id() << " emctraits id is " << emctraits->bemcId() << " towid " << emctraits->btowId() << endl; 
        cout << "= nbtow is " << nBtow << endl;
        for (int i = 0; i < nBtow; i++)
            {
                StPicoBTOWHit* bhit = mPicoDst->btowHit(i);
                //hBtowAdcId->Fill(bhit->adc(), bhit->id());
                printf("--BTOW ID %i -- adc  -- %i -- energy -- %f \n",bhit->id(), bhit->adc(), bhit->energy());
            } */

        //if (t->pMom().mag() < 4) continue;
        bool getBEMC = GetBEMCdist(t, etaphi);
        if (getBEMC) {
            hTrackEtaPhiD->Fill(etaphi[0],etaphi[1]);
            hTrackR->Fill(sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1])); 
            hTrackRvDist->Fill(sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1]),sqrt(emctraits->zDist()*emctraits->zDist()+10000*emctraits->phiDist()*emctraits->phiDist()));   
            hTracketaDistzDist->Fill(etaphi[0],emctraits->zDist());
            hTrackphiDistphiDist->Fill(etaphi[1],emctraits->phiDist());    }
        #endif  
        hTrackpMomgMom->Fill(t->pMom().mag(),t->gPtot());

        //cout << "=2= eta is " << etaphi[0] << " phi is " << etaphi[1] << endl;
        
        //---
        if (! SelectTrack(t)) continue;

        #ifndef EMCEFF
   		// DCA STUDY
        double primvpos2[] = {primVpos.x(), primVpos.y(), primVpos.z()};
        Float_t dca2 = t->dcaGeometry().thelix().Dca(primvpos2);
        Float_t dca3 = t->helix().distance(primVpos);
        //cout << "moving by " << t->helix().pathLength(primVpos) << endl;
        StPhysicalHelixD trHelix = t->helix();
        trHelix.moveOrigin(trHelix.pathLength(primVpos)); // doesnt do shit
        Float_t dca4 = (trHelix.origin() - primVpos ).mag(); 
        //cout << "dca is " << dca << " vs " << dca2 << " vs " << dca3 << " vs " << dca4 << endl;

        if (dca3 > 1.5) continue;
        #endif
    	hTracksSelection->Fill(15);

    	// ELECTRON HISTOS
    	hElectronPhivdEdx->Fill(t->pMom().phi(),t->dEdx());
        hElectronPt->Fill(t->pMom().perp());
        hElectronP->Fill(t->pMom().mag());
        hElectronzDistphiDist->Fill(emctraits->zDist(),emctraits->phiDist());
        hElectronDca->Fill(dca3);
        if (emctraits->e() != 0) hElectronDistvE0E->Fill(emctraits->zDist()*emctraits->zDist()+10000*emctraits->phiDist()*emctraits->phiDist(),emctraits->e0()/emctraits->e());
        if (getBEMC) {
            hElectronEtaPhiD->Fill(etaphi[0],etaphi[1]);
            hElectronR->Fill(sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1])); 
            hElectronRvDist->Fill(sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1]),sqrt(emctraits->zDist()*emctraits->zDist()+10000*emctraits->phiDist()*emctraits->phiDist()));   
            hElectronetaDistzDist->Fill(etaphi[0],emctraits->zDist());
            hElectronphiDistphiDist->Fill(etaphi[1],emctraits->phiDist());     }


        //cout << "=3= eta is " << etaphi[0] << " phi is " << etaphi[1] << endl;

        if (sqrt(etaphi[0]*etaphi[0]+etaphi[1]*etaphi[1]) > 0.025) continue;

        // SELECT ALL ELECTRONS
        AllTracks1.push_back(i);
        hTracksSelection->Fill(16);

        //#ifdef NSIGMAEFF
        //hElectronnSigmavsPt->Fill(t->pMom().perp(),t->nSigmaElectron());
        //hElectrondEdxvsPt->Fill(t->pMom().perp(),t->dEdx());

        /*TParameter<float> parf;
        TClonesArray &TnSigmaEad = *TnSigmaE;
        parf.SetVal(t->dEdx());
        new (TnSigmaEad[fkTr]) TParameter<float>(parf);

        TClonesArray &TPtad = *TPt;
        parf.SetVal(t->pMom().perp());
        new (TPtad[fkTr]) TParameter<float>(parf);
        // p.mag can be added here
        fkTr++;*/
        //#endif

        #ifdef TRIGEFF
        if (evFlag == 1) {
            hElectronPtHT->Fill(t->pMom().perp());
            hElectronPHT->Fill(t->pMom().mag()); }

        if (evFlag == 2 || evFlag == 3) {
            hElectronPtMB->Fill(t->pMom().perp());
            hElectronPMB->Fill(t->pMom().mag());
            hElectronPtMBps->Fill(t->pMom().perp(),psF);
            hElectronPMBps->Fill(t->pMom().mag(),psF); }
        #endif

        // SELECT TRIGGER ELECTRONS----------
        bool isTrigger1 = true;
        #ifndef EMCEFF
        isTrigger1 = false;
        for (int j = 0; j < TrigTowers1.size(); j++)
        {
            StPicoBTOWHit* bhit = mPicoDst->btowHit(TrigTowers1[j]);
            if ( fabs(emctraits->e1() - bhit->energy() ) > 0.01 &&
            	 fabs(emctraits->e0() - bhit->energy() ) > 0.01 ) continue;
            isTrigger1 = true;
        	hElectronTowervp->Fill(t->pMom().mag(),bhit->id());
        }
        #endif
        if (isTrigger1) 
        {
            TriggerTracks1.push_back(i);
            hTracksSelection->Fill(17);
        }   
    
    }

    hEventAllTracks->Fill(AllTracks1.size());
    hEventTriggerTracks->Fill(TriggerTracks1.size());

    hEventVzvsNPrim->Fill(mEvent->primaryVertex().z(), nPrim, cent9);
    hEventVzvsNPrimHard->Fill(mEvent->primaryVertex().z(), nPrimHard, cent9);
    //-------------------------------------------------------


    if ( AllTracks1.size() < 2 ) return kStOk;
    hEventSelection->Fill(3);

    #ifdef NSIGMAEFF
    for (int i = 0; i < AllTracks1.size()-1; ++i)
    {
        StPicoTrack* t1 = (StPicoTrack*)mPicoDst->track(AllTracks1[i]);
        StThreeVectorF p1 = t1->pMom();
        TLorentzVector* tl1 = new TLorentzVector(p1.x(),p1.y(),p1.z(),sqrt(p1.mag()*p1.mag() + elmass*elmass));
        for (int j = i+1; j < AllTracks1.size(); ++j)
        {
            StPicoTrack* t2 = (StPicoTrack*)mPicoDst->track(AllTracks1[j]);
            if (t1->charge()*t2->charge()>0) continue;
            StThreeVectorF p2 = t2->pMom();
            TLorentzVector* tl2 = new TLorentzVector(p2.x(),p2.y(),p2.z(),sqrt(p2.mag()*p2.mag() + elmass*elmass));

            TLorentzVector* q = new TLorentzVector(*tl1+*tl2);
            //if (q->M() > 0.1) continue;
            if (q->M() > 3.25) continue;
            if (q->M() < 2.9) continue;
            hElectronnSigmavsP->Fill(p1.mag(),t1->nSigmaElectron(),cent9);
            hElectronnSigmavsP->Fill(p2.mag(),t2->nSigmaElectron(),cent9);
            hElectrondEdxvsP->Fill(p1.mag(),t1->dEdx(),cent9);
            hElectrondEdxvsP->Fill(p2.mag(),t2->dEdx(),cent9);
        }
    }
    #endif

    #ifdef EMCEFF
    for (int i = 0; i < AllTracks1.size()-1; ++i)
    {
        StPicoTrack* t1 = (StPicoTrack*)mPicoDst->track(AllTracks1[i]);
        StThreeVectorF p1 = t1->pMom();
        TLorentzVector* tl1 = new TLorentzVector(p1.x(),p1.y(),p1.z(),sqrt(p1.mag()*p1.mag() + elmass*elmass));
        for (int j = i+1; j < AllTracks1.size(); ++j)
        {
            StPicoTrack* t2 = (StPicoTrack*)mPicoDst->track(AllTracks1[j]);
            StThreeVectorF p2 = t2->pMom();
            TLorentzVector* tl2 = new TLorentzVector(p2.x(),p2.y(),p2.z(),sqrt(p2.mag()*p2.mag() + elmass*elmass));

            hMotherDPhi->Fill(0.5*(p1.mag()+p2.mag()),fabs(p1.phi()-p2.phi()));

            TLorentzVector* q = new TLorentzVector(*tl1+*tl2);
            hMotherDPhivM->Fill(q->M(),fabs(p1.phi()-p2.phi()));

            hEMCEFFcUS->Fill(11);
            if (q->M() > 3.25) continue;
            if (q->M() < 2.9) continue;
            hEMCEFFcUS->Fill(12);

            //if (fabs(p1.phi()-p2.phi()) < 0.10) continue;

            if (t1->charge()*t2->charge()<0) {
                FillEMCEFF(t1,"US");
                FillEMCEFF(t2,"US"); }

            if (t1->charge()*t2->charge()>0) {
                FillEMCEFF(t1,"LS");
                FillEMCEFF(t2,"LS"); }
        }
    }
    #endif

    if ( TriggerTracks1.size() < 1 ) return kStOk;
    hEventSelection->Fill(4);

    #ifdef MIXEVENTS
    if (mixCounter < 2*mixSaved)
    {
        if (mixCounter%2==0)
        {
            StPicoTrack* t1 = (StPicoTrack*)mPicoDst->track(TriggerTracks1[0]);
            StThreeVectorF p1 = t1->pMom();
            TLorentzVector* tl1 = new TLorentzVector(p1.x(),p1.y(),p1.z(),sqrt(p1.mag()*p1.mag() + elmass*elmass));

            new (mtreeTrigP[mtreeCounter]) TLorentzVector(*tl1);



            /*trigElectrons.push_back(tl1);
            trigElectronsVz.push_back(mEvent->primaryVertex().z());
            trigElectronsC.push_back(cent9);*/
        }
        if ((mixCounter+1)%2==0)
        {
            StPicoTrack* t1 = (StPicoTrack*)mPicoDst->track(AllTracks1[0]);
            StThreeVectorF p1 = t1->pMom();
            TLorentzVector* tl1 = new TLorentzVector(p1.x(),p1.y(),p1.z(),sqrt(p1.mag()*p1.mag() + elmass*elmass));

            new (mtreeAllP[mtreeCounter]) TLorentzVector(*tl1);

            /*allElectrons.push_back(tl1);
            allElectronsVz.push_back(mEvent->primaryVertex().z());
            allElectronsC.push_back(cent9);*/
        }
        mixCounter++;
        mtreeCounter++;
    }
    #endif

    //-------------------------------------------------------

    vector<int> UsedTracks1;
    for (int i = 0;i < TriggerTracks1.size();i++){            //
        StPicoTrack* t1 = (StPicoTrack*)mPicoDst->track(TriggerTracks1[i]);
        for (int j = 0;j < AllTracks1.size();j++){
            if (TriggerTracks1[i] == AllTracks1[j])
                continue;

            bool isUsed = false;
            for (int k = 0; k < UsedTracks1.size(); k++)
            {
                if (UsedTracks1[k] == AllTracks1[j])
                {
                    isUsed = true;
                }
            }
            if (isUsed) continue;

            StPicoTrack* t2 = (StPicoTrack*)mPicoDst->track(AllTracks1[j]);
            
            float mass = SelectPair(t1,t2);

            if (mass < 0) continue;

            #ifdef WRITETREE
            TParameter<float> paramf;
            TParameter<int> parami;
            #endif
            
            if ((t1->charge() > 0) && (t2->charge() > 0)) {
                hIMpp->Fill(mass,t2->pMom().perp(),cent9);
            }
            if ((t1->charge() < 0) && (t2->charge() < 0)) {
            	hIMmm->Fill(mass,t2->pMom().perp(),cent9);
            }
            if (t1->charge()*t2->charge() < 0){
               	hIMpm->Fill(mass,t2->pMom().perp(),cent9);
               	hIMun->Fill(mass);

                #ifdef WRITETREE 
                paramf.SetVal(mass);
                parami.SetVal(0);
                new (utreeUpsMass[utreeCounter]) TParameter<float> (paramf);
                new (utreeUpsSigns[utreeCounter]) TParameter<int> (parami);
                utreeCounter++;
                #endif
            }
            if (t1->charge()*t2->charge() > 0){
               	hIMli->Fill(mass);

                #ifdef WRITETREE 
                paramf.SetVal(mass);
                parami.SetVal(1);
                new (utreeUpsMass[utreeCounter]) TParameter<float> (paramf);
                new (utreeUpsSigns[utreeCounter]) TParameter<int> (parami);
                utreeCounter++;
                #endif
            }
        }
        UsedTracks1.push_back(TriggerTracks1[i]);
    }

    #ifdef MIXEVENTS
    if (mtreeCounter > 0) mixTree->Fill();
    #endif

    #ifdef WRITETREE
    if (utreeCounter > 0) upsMass->Fill();
    #endif
  	
  	//-------------------------------------------------------
    return kStOK;
}

