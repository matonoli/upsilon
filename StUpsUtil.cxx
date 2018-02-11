#include "headers.h"
#include "StMtdTrigger.h"
#include "StMtdJpsiUtil.h"

ClassImp(StMtdJpsiUtil)

//_____________________________________________________________________________
StMtdJpsiUtil::StMtdJpsiUtil(const Char_t *name, const Char_t *title) : TNamed(name,title)
{
  // default constructor

  mMCflag                  = 0;

  // event cuts
  mTrigger[0]              = kTRUE;
  mTrigger[1]              = kFALSE;
  mTrigger[2]              = kFALSE;
  mTrigger[3]              = kFALSE;
  mYear                    = 2014;
  
  for(Int_t i=0; i<kNtrig; i++)
    mTriggerIDs[i].clear();
  mBadRunIDs.clear();

  // centrality
  mRefMultCorr             = NULL;
  mVzCorrMethod            = 0;
  mVzCorrFileName          = "";
  mgVzCorr                 = NULL;
  mgVzCorr1                = NULL;

  // vertex cuts
  mRequireVtxRanking       = kFALSE;
  mRequireMtdHitForPrimVtx = kFALSE;
  mUsePrimVtxClosestToVpd  = kTRUE;
  mMaxVtxZ                 = 100;
  mMaxVtxR                 = 2;
  mMaxDiffz                = 3;

  // track cuts
  mTrackType               = 0;
  mMinTrkPt                = 1;
  mMaxTrkPt                = 1e4;
  mMinTrkPhi               = 0;
  mMaxTrkPhi               = 2*pi;
  mMinTrkEta               = -0.8;
  mMaxTrkEta               = 0.8;
  mMinNHitsFit             = 15;
  mMinNHitsDedx            = 10;
  mMinFitHitsFaction       = 0;
  mMaxDca                  = 1;

  // MTD hit cuts
  mTrigTimeCut             = kFALSE;
  memset(mTrigWinCut_low,  -1, sizeof(mTrigWinCut_low));
  memset(mTrigWinCut_high, -1, sizeof(mTrigWinCut_high));
  
  // muon PID
  mMinNsigmaPi             = -1;
  mMaxNsigmaPi             = 3;
  mBTofMatch               = kFALSE;
  mMtdHitTrigger           = kFALSE;

  mApplyDzCut              = kTRUE;
  mApplyDyCut              = kTRUE;
  mApplyDtofCut            = kTRUE;
  mPtDepCutDz              = kTRUE;
  mPtDepCutDy              = kTRUE;
  mPtDepCutDtof            = kFALSE;
  fResDzVsPt               = NULL;
  fResDyVsPt               = NULL;
  fResDtofVsPt             = NULL;
  mSigmaDz1                = 2.0;
  mSigmaDz2                = 2.5;
  mSigmaDy1                = 2.0;
  mSigmaDy2                = 2.5;
  mSigmaDtof               = 1.0;
  mMinMuonDeltaZ           = -20.;
  mMaxMuonDeltaZ           = 20.;
  mMinMuonDeltaY           = -20.;
  mMaxMuonDeltaY           = 20.;
  mMinMuonDeltaTof         = -1.e4;
  mMaxMuonDeltaTof         = 1.0;
  mMinMuonPt               = 1.;
  mMaxMuonPt               = 1e4;

  // global track pair cuts
  mMaxTrkPairDcaDr         = 5.;
  mMaxTrkPairDcaDz         = 4.;
}

//_____________________________________________________________________________
StMtdJpsiUtil::StMtdJpsiUtil(const StMtdJpsiUtil &mtdUtil) : TNamed(mtdUtil.fName,mtdUtil.fTitle)
{
  mYear   = mtdUtil.mYear;
  mMCflag = mtdUtil.mMCflag;
  for(int i=0; i<kNtrig; i++)
    {
      mTrigger[i] = mtdUtil.mTrigger[i];
      mTriggerIDs[i].clear();
    }
  mBadRunIDs.clear();

  // vertex cuts
  mRequireVtxRanking       = mtdUtil.mRequireVtxRanking;
  mRequireMtdHitForPrimVtx = mtdUtil.mRequireMtdHitForPrimVtx;
  mUsePrimVtxClosestToVpd  = mtdUtil.mUsePrimVtxClosestToVpd;
  mMaxVtxZ                 = mtdUtil.mMaxVtxZ;
  mMaxVtxR                 = mtdUtil.mMaxVtxR;
  mMaxDiffz                = mtdUtil.mMaxDiffz;

  // track cuts
  mTrackType               = mtdUtil.mTrackType;
  mMinTrkPt                = mtdUtil.mMinTrkPt;
  mMaxTrkPt                = mtdUtil.mMaxTrkPt;
  mMinTrkPhi               = mtdUtil.mMinTrkPhi;
  mMaxTrkPhi               = mtdUtil.mMaxTrkPhi;
  mMinTrkEta               = mtdUtil.mMinTrkEta;
  mMaxTrkEta               = mtdUtil.mMaxTrkEta;
  mMinNHitsFit             = mtdUtil.mMinNHitsFit;
  mMinNHitsDedx            = mtdUtil.mMinNHitsDedx;
  mMinFitHitsFaction       = mtdUtil.mMinFitHitsFaction;
  mMaxDca                  = mtdUtil.mMaxDca;

  // MTD hit cuts
  mTrigTimeCut             = mtdUtil.mTrigTimeCut;
  memset(mTrigWinCut_low,  -1, sizeof(mTrigWinCut_low));
  memset(mTrigWinCut_high, -1, sizeof(mTrigWinCut_high));
  
  // muon PID
  mMinNsigmaPi             = mtdUtil.mMinNsigmaPi;
  mMaxNsigmaPi             = mtdUtil.mMaxNsigmaPi;
  mBTofMatch               = mtdUtil.mBTofMatch;
  mMtdHitTrigger           = mtdUtil.mMtdHitTrigger;

  mApplyDzCut              = mtdUtil.mApplyDzCut;
  mApplyDyCut              = mtdUtil.mApplyDyCut;
  mApplyDtofCut            = mtdUtil.mApplyDtofCut;
  mPtDepCutDz              = mtdUtil.mPtDepCutDz;
  mPtDepCutDy              = mtdUtil.mPtDepCutDy;
  mPtDepCutDtof            = mtdUtil.mPtDepCutDtof;
  fResDzVsPt               = NULL;
  fResDyVsPt               = NULL;
  fResDtofVsPt             = NULL;
  mSigmaDz1                = mtdUtil.mSigmaDz1;
  mSigmaDz2                = mtdUtil.mSigmaDz2;
  mSigmaDy1                = mtdUtil.mSigmaDy1;
  mSigmaDy2                = mtdUtil.mSigmaDy2;
  mSigmaDtof               = mtdUtil.mSigmaDtof;
  mMinMuonDeltaZ           = mtdUtil.mMinMuonDeltaZ;
  mMaxMuonDeltaZ           = mtdUtil.mMaxMuonDeltaZ;
  mMinMuonDeltaY           = mtdUtil.mMinMuonDeltaY;
  mMaxMuonDeltaY           = mtdUtil.mMaxMuonDeltaY;
  mMinMuonDeltaTof         = mtdUtil.mMinMuonDeltaTof;
  mMaxMuonDeltaTof         = mtdUtil.mMaxMuonDeltaTof;
  mMinMuonPt               = mtdUtil.mMinMuonPt;
  mMaxMuonPt               = mtdUtil.mMaxMuonPt;

  // global track pair cuts
  mMaxTrkPairDcaDr         = mtdUtil.mMaxTrkPairDcaDr;
  mMaxTrkPairDcaDz         = mtdUtil.mMaxTrkPairDcaDz;
}
 
//_____________________________________________________________________________
StMtdJpsiUtil::~StMtdJpsiUtil()
{
  // default destructor
  if(fResDzVsPt) delete fResDzVsPt;
  if(fResDyVsPt) delete fResDyVsPt;
  if(fResDtofVsPt) delete fResDtofVsPt;
}

//_____________________________________________________________________________
void StMtdJpsiUtil::Init()
{
  printf("==================================\n");
  LOG_INFO << "Initialize triggers for run " << mYear << endm;

  // trigger id
  if(mYear == 2013)
    {
      // pp 500 GeV
      const Int_t ntrig = 2; 
      Int_t di_muon[ntrig] = {430103, 430113};
      Int_t sg_muon[ntrig] = {430101, 430111};
      Int_t el_muon[ntrig] = {430102, 430112};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[2].push_back(430122);
      mTriggerIDs[3].push_back(430001);
      mTriggerIDs[3].push_back(430011);
      mTriggerIDs[3].push_back(430021);
      mTriggerIDs[3].push_back(430031);
      // mTriggerIDs[3].push_back(430005);
      // mTriggerIDs[3].push_back(430015);
    }
  else if(mYear == 2014)
    {
      // AuAu 200 GeV
      const Int_t ntrig = 5; 
      Int_t di_muon[ntrig] = {450601, 450611, 450621, 450631, 450641};
      Int_t sg_muon[ntrig] = {450600, 450610, 450620, 450630, 450640};
      Int_t el_muon[ntrig] = {450602, 450612, 450622, 450632, 450642};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      //mTriggerIDs[0].push_back(450604);
      //mTriggerIDs[0].push_back(450605);
      //mTriggerIDs[0].push_back(450606);
      //mTriggerIDs[0].push_back(450635); // HLT
      //mTriggerIDs[0].push_back(450645); // HLT
      mTriggerIDs[3].push_back(450013); //VPD-ZDC-novtx-mon
      mTriggerIDs[3].push_back(450023); //VPD-ZDC-novtx-mon
      // mTriggerIDs[3].push_back(450050); //VPDMB-5-p-nobsmd-hlt
      // mTriggerIDs[3].push_back(450060); //VPDMB-5-p-nobsmd-hlt
      // mTriggerIDs[3].push_back(450010); //VPDMB-30
      // mTriggerIDs[3].push_back(450020); //VPDMB-30

      // bad run list
      const int nRuns = 38;
      const int runIDs[nRuns] = {15078103, 15078104, 15078107, 15078108, 15079059, 
				 15079061, 15084002, 15084022, 15084052, 15088003, 
				 15090006, 15097032, 15097034, 15102021, 15104018, 
				 15104039, 15104059, 15106008, 15106009, 15106010, 
				 15106011, 15107077, 15110032, 15110038, 15119021, 
				 15132026, 15142054, 15151035, 15151036, 15151037, 
				 15151038, 15151039, 15151040, 15151041, 15151042, 
				 15151043, 15162019, 15166023};

      for(Int_t i=0; i<nRuns; i++)
	{
	  mBadRunIDs.push_back(runIDs[i]);
	}
    }
  else if(mYear == 2015)
    {
      // pp 200 GeV
      const Int_t ntrig = 5; 
      Int_t di_muon[ntrig] = {470602, 480602, 490602, 500602, 510602};
      Int_t sg_muon[ntrig] = {470600, 480600, 490600, 500600, 510600};
      Int_t el_muon[ntrig] = {470601, 480601, 490601, 500601, 510601};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      // bad run list
      const int nRuns = 194;
      const int runIDs[nRuns] = {16047004, 16047005, 16047008, 16047104, 16052036,
				 16052037, 16052038, 16054059, 16054060, 16054061,
				 16054062, 16054063, 16054064, 16054069, 16054070,
				 16054072, 16054073, 16054074, 16054075, 16054077,
				 16054078, 16054079, 16054080, 16054082, 16054086,
				 16054087, 16055002, 16055003, 16055004, 16055005,
				 16055007, 16055010, 16055011, 16055012, 16055013,
				 16055018, 16055019, 16055021, 16055022, 16055024,
				 16055025, 16055124, 16055127, 16058072, 16060018,
				 16060025, 16060036, 16060053, 16060054, 16060055,
				 16060056, 16060057, 16060058, 16060059, 16060060,
				 16060061, 16060062, 16060063, 16060064, 16060065,
				 16062008, 16062009, 16062011, 16062014, 16062048,
				 16063096, 16063097, 16063099, 16065011, 16065059,
				 16066028, 16069045, 16069050, 16069060, 16071043,
				 16071044, 16072047, 16073004, 16073007, 16073015,
				 16082014, 16083044, 16084017, 16088013, 16089002,
				 16091057, 16091058, 16091059, 16091061, 16092040,
				 16092044, 16094018, 16095031, 16095034, 16098007,
				 16100023, 16100024, 16100025, 16101018, 16104002,
				 16104046, 16105044, 16106001, 16108026, 16108031,
				 16108032, 16109010, 16109011, 16109012, 16114032,
				 16115056,
				 16124017, 16124018, 16124019, 16124033, 16124035,
				 16124037, 16125001, 16125003, 16125014, 16125015,
				 16125016, 16125024, 16125035, 16125038, 16125039,
				 16127048, 16127049, 16128005, 16128006, 16128014,
				 16128056, 16130012, 16130015, 16130016, 16130032,
				 16132021, 16132022, 16132046, 16133085, 16134042,
				 16135047, 16138013, 16140015, 16140016, 16141036,
				 16146002, 16149001, 16149002, 16149003, 16149004,
				 16149005, 16149008, 16149009, 16149010, 16149011,
				 16149013, 16149014, 16150001, 16150042, 16152032,
				 16152034, 16152035, 16154009, 16154010, 16154011,
				 16154021, 16155017, 16155039, 16156010, 16156028,
				 16156030, 16156031, 16156032, 16156033, 16156034, 
				 16157034, 16157047, 16157071, 16158021, 16158032,
				 16158039, 16158042, 16158043, 16158044, 16158045,
				 16159009, 16159019, 16161034};

      for(Int_t i=0; i<nRuns; i++)
	{
	  mBadRunIDs.push_back(runIDs[i]);
	}

    }
  else if (mYear == 2016)
    {
      const Int_t ntrig = 3; 
      Int_t di_muon[ntrig] = {520602, 520612, 520622};
      Int_t sg_muon[ntrig] = {520604, 520614, 520624};
      Int_t el_muon[ntrig] = {520606, 520616, 520626};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[0].push_back(46);
      mTriggerIDs[0].push_back(47);
      mTriggerIDs[0].push_back(520803);

      mTriggerIDs[3].push_back(520005);  //VPD-ZDC-novtx
      mTriggerIDs[3].push_back(520015);  //VPD-ZDC-novtx
      mTriggerIDs[3].push_back(570005);  //VPD-ZDC-novtx
    }
  printTriggerIDs();

  // centrality
  if(mVzCorrMethod==1)
    {
      // use linear extraplation implemented in TGraph::Eval()
      // for Vz correction
      if(mVzCorrFileName.Length()==0)
	{
	  LOG_ERROR << "[e] No input file for VzCorr graph! Please double check! " <<endm;
	}
      else
	{
	  TFile *fVzIn = TFile::Open(mVzCorrFileName.Data(),"read");
	  if(fVzIn->IsOpen())
	    {
	      LOG_INFO << "Open file " << mVzCorrFileName.Data() << endm;
	      mgVzCorr  = (TGraphErrors*)fVzIn->Get("ggRefVsVz_VpdMB_noVtx");
	      mgVzCorr1 = (TGraphErrors*)fVzIn->Get("ggRefVsVz_VpdMB_noVtx_3");
	      fVzIn->Close();
	    }
	  else
	    {
	      LOG_ERROR << "[e] Cannot open input file for VzCorr graph! Please double check! " <<endm;
	    }
	}
    }  


  // dy, dz, dtof resolution vs pt
  fResDtofVsPt = new TF1("fResDtofVsPt","[0]+[1]*exp([2]/x)");
  fResDtofVsPt->SetParameters(0.0817528, 0.0169419, 4.34897);
  fResDzVsPt = new TF1("fResDzVsPt","[0]+[1]*exp([2]/x)");
  fResDyVsPt = new TF1("fResDyVsPt","[0]+[1]*exp([2]/x)");
  fResDzVsPt->SetParameters(-21.04, 21.09, 0.693);
  fResDyVsPt->SetParameters(-12.61, 13.43, 0.889);
  if (mYear == 2013)
   {
     fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.444217);
     fResDyVsPt->SetParameters(-17.6867, 18.4528, 0.637142);
   }
  else if(mYear==2014)
    {
      fResDzVsPt->SetParameters(-21.04, 21.09, 0.693);
      fResDyVsPt->SetParameters(-12.61, 13.43, 0.889);
    }
  else if(mYear==2015)
    {
      fResDzVsPt->SetParameters(-32.6793, 32.6034, 0.44421);
      fResDyVsPt->SetParameters(-17.6867, 18.4528, 0.637142);
    }
  else if(mYear==2016)
    {
      fResDzVsPt->SetParameters(-21.04, 21.09, 0.693);
      fResDyVsPt->SetParameters(-12.61, 13.43, 0.889);
    }
}

//_____________________________________________________________________________
void StMtdJpsiUtil::setTriggerIDs(const IntVec id[kNtrig])
{
  for(Int_t i=0; i<kNtrig; i++)
    mTriggerIDs[i] = id[i];

  printTriggerIDs();
}

//_____________________________________________________________________________
void StMtdJpsiUtil::printTriggerIDs()
{
  for(Int_t j=0; j<kNtrig; j++)
    {
      cout << trigName[j] << ": ";
      for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	{
	  cout << mTriggerIDs[j][i] << " ";
	}
      cout << endl;
    }
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passRun(const int runId)
{ 
  for(UInt_t i=0; i<mBadRunIDs.size(); i++)
    {
      if(runId==mBadRunIDs[i])
	{
	  return kFALSE;
	}
    }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrigger(StEvent *event)
{ 
  if(getTriggerType(event)==-1) return kFALSE;
  else return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrigger(StMuDst *event)
{ 
  if(getTriggerType(event)==-1) return kFALSE;
  else return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrigger(StPicoDst *pico)
{ 
  if(getTriggerType(pico)==-1) return kFALSE;
  else return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passEvent(StEvent *event)
{
  // Cut on vertex z
  Int_t vtx_index = selectPrimaryVertex(event);
  if(vtx_index<0) return kFALSE;
  StPrimaryVertex* priVertex = event->primaryVertex(vtx_index);
  if(!priVertex) return kFALSE;
  StThreeVectorF verPos = priVertex->position();

  StBTofCollection * tofCollection = event->btofCollection(); 
  if(!tofCollection) return kFALSE;
  StBTofHeader *tofHeader = tofCollection->tofHeader();
  if(!tofHeader)     return kFALSE;
  if(mRequireVtxRanking && priVertex->ranking()<0) return kFALSE; 
  if(! (TMath::Abs(verPos.x())>0 || TMath::Abs(verPos.y())>0 || TMath::Abs(verPos.z())>0) ) return kFALSE;
  double vr = TMath::Sqrt(verPos.x()*verPos.x() + verPos.y()*verPos.y());
  if(!isGoodVertex(tofHeader->vpdVz(), verPos.z(), vr)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passEvent(StMuDst *event)
{
  // Cut on vertex z
  Int_t vtx_index = selectPrimaryVertex(event);
  if(vtx_index<0) return kFALSE;
  StMuPrimaryVertex* priVertex = event->primaryVertex(vtx_index);
  if(!priVertex) return kFALSE;
  StThreeVectorF verPos = priVertex->position();

  StBTofHeader *tofHeader = event->btofHeader();
  if(!tofHeader) return kFALSE;
  if(mRequireVtxRanking && priVertex->ranking()<0) return kFALSE; 
  if(! (TMath::Abs(verPos.x())>0 || TMath::Abs(verPos.y())>0 || TMath::Abs(verPos.z())>0) ) return kFALSE;
  double vr = TMath::Sqrt(verPos.x()*verPos.x() + verPos.y()*verPos.y());
  if(!isGoodVertex(tofHeader->vpdVz(), verPos.z(), vr)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passEvent(StPicoEvent *event)
{
  // Cut on vertex z
  StThreeVectorF verPos = event->primaryVertex();
  if(mRequireVtxRanking && event->ranking()<0) return kFALSE; 
  if(! (TMath::Abs(verPos.x())>0 || TMath::Abs(verPos.y())>0 || TMath::Abs(verPos.z())>0) ) return kFALSE;
  double vr = TMath::Sqrt(verPos.x()*verPos.x() + verPos.y()*verPos.y());
  if(!isGoodVertex(event->vzVpd(), verPos.z(), vr)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getRefMultCorr(Int_t year, Int_t runId, Double_t gRefMult, Double_t z_val, Double_t zdc_val)
{
  if(!mRefMultCorr)
    {
      LOG_ERROR << "No RefMult class avaiable" << endm;
      return -1;
    }

  Double_t gRefMult_corr = -1;
  mRefMultCorr->init(runId);
  mRefMultCorr->initEvent(gRefMult,z_val,zdc_val*1e3);
  gRefMult_corr = mRefMultCorr->getRefMultCorr();
 
  if(year==2014 && zdc_val > 55)
    {
      // http://www.star.bnl.gov/protected/heavy/marr/paper/Run14_AuAu200_Jpsi/SupportMaterial/readMiniTree_interpolation.C
      const Double_t par0 =   558.41   ;
      const Double_t par1 =   0.342761   ;
      const Double_t par2 =   -0.0138915   ;
      const Double_t par3 =   -1.19606e-05   ;
      const Double_t par4 =   3.2903e-06   ;
      const Double_t par5 =   -7.61524e-10   ;
      const Double_t par6 =   -1.70019e-10   ;    
      const Double_t par7 =   0.0;
    
      Double_t z = z_val;
      const Double_t  gRefMult_ref = par0; // Reference mean gRefMult at z=0
      const Double_t  gRefMult_z = par0 + par1*z + par2*z*z + par3*z*z*z + par4*z*z*z*z + par5*z*z*z*z*z +  par6*z*z*z*z*z*z;
      Double_t  Hovno = 0;
      if(mVzCorrMethod==0) Hovno = (gRefMult_ref+par7)/gRefMult_z; // Correction factor for gRefMult, takes into account z_vertex dependence
      else if(mVzCorrMethod==1) 
	{
	  if(mgVzCorr) Hovno = mgVzCorr->Eval(0)/mgVzCorr->Eval(z)*mgVzCorr1->Eval(0)/mgVzCorr1->Eval(z);
	  else LOG_ERROR << "[e] Vz correction is not available!" << endm;
	}
      else
	{
	  LOG_WARN << "[w] Unknown scheme for Vz correction. Please choose 0 or 1!" << endm;
	}

    
      Double_t zdc = zdc_val;
      const Int_t nVzBins =1;
      const Double_t par00[nVzBins] ={184.786  }  ;
      const Double_t par11[nVzBins] ={-0.395819 }  ;
      const Double_t parBaseLine=par00[0]+par11[0]*0.;//for Vz [-1,1],Zdc=0kHz
    
      Double_t tmp=0;
      for(int i=0;i<nVzBins;i++) {
        tmp=(par00[i]+par11[i]*zdc);//zdc unit should be kHz,since when i did correction, the int is based on kHz
      }
      Double_t zdcFactor= parBaseLine/tmp;
    
      gRefMult_corr  = (gRefMult+gRandom->Rndm())*Hovno*zdcFactor;
    }

  return gRefMult_corr ;
}

//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getWeight(Int_t year, Double_t zdc_val, Double_t gRefMultCorr)
{
  double weight = 1;
  if(!mRefMultCorr)
    {
      LOG_ERROR << "No RefMult class avaiable" << endm;
      return weight;
    }

  weight = mRefMultCorr->getWeight();
  if(year==2014 && zdc_val > 55)
    {
      Double_t par0 =  1.56227 ;
      Double_t par1 =  -22.4469 ;
      Double_t par2 =  0.694748 ;
      Double_t par3 =  7.00885 ;
      Double_t par4 =  -0.00457379 ;
      Double_t par6 =  579.825 ;
      Double_t par7 =  9.77725e-06  ;
      Double_t A    =  0.00000e+00 ;
    
      Double_t mVz = 100.0; // this has to be modified...
    
      if(gRefMultCorr != -(par3/par2)) // avoid denominator = 0
	{
	  weight = par0 + par1/(par2*gRefMultCorr + par3) + par4*(par2*gRefMultCorr + par3) + par6/pow(par2*gRefMultCorr+par3 ,2) + par7*pow(par2*gRefMultCorr+par3 ,2);
	  weight = weight + (weight-1.0)*(A*mVz*mVz); // z-dependent weight correction
	}
    
      if(gRefMultCorr>=350.) weight = 1.;
    }
  return weight;
}

//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getCentralityBin16(Int_t year, Double_t zdc_val, Double_t gRefMultCorr)
{
  double cent = -1;
  if(!mRefMultCorr)
    {
      LOG_ERROR << "No RefMult class avaiable" << endm;
      return cent;
    }
  cent = mRefMultCorr->getCentralityBin16();


  if(year==2014 && zdc_val > 55)
    {
      cent = -1;
      double boundaries[17] = {10, 15, 22, 31, 43, 57, 75, 97, 123, 153, 188, 229, 275, 329, 392, 467,1000};
      for(int i=0; i<16; i++)
	{
	  if(gRefMultCorr>boundaries[i] && gRefMultCorr<=boundaries[i+1])
	    {
	      cent = i;
	      break;
	    }
	}
    }
  return cent;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMuonCandidate(const StTrack *track)
{
  if(!track) return kFALSE;
  double nSigmaPi = -999.;
  int    tofIndex = -1;
  double dz = -999., dy = -999, dtof = -999.;
  bool isMtdTrig = kFALSE;
 
  // pt
  double pt = track->geometry()->momentum().perp();

  // nSigmaPi cut
  StTpcDedxPidAlgorithm pidAlgorithm;
  const StParticleDefinition *pd = track->pidTraits(pidAlgorithm);
  if(pd && pidAlgorithm.traits())
    {
      static StPionPlus* Pion = StPionPlus::instance();
      nSigmaPi = pidAlgorithm.numberOfSigma(Pion);
    }

  // tof match
  tofIndex = track->isBToFMatched() ? 0 : -1;

  // dz cut
  StMtdPidTraits* mtdpid = 0;
  const StSPtrVecTrackPidTraits& traits = track->pidTraits();
  for(unsigned int it=0; it<traits.size(); it++)
    {
      if (traits[it]->detector() == kMtdId)
        {
          mtdpid = dynamic_cast<StMtdPidTraits*>(traits[it]);
          break;
        }
    }
  if(mtdpid) 
    {
      dz = mtdpid->deltaZ();
      dy = mtdpid->deltaY();
      dtof = mtdpid->timeOfFlight() - mtdpid->expTimeOfFlight();
    }

  return isMuonCandidate(pt, nSigmaPi, tofIndex, dz, dy, dtof, isMtdTrig);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMuonCandidate(const StMuTrack *track, const StMuDst *muDst, StMtdTrigger *trig)
{
  if(!track) return kFALSE;
  double nSigmaPi = -999.;
  int    tofIndex = -1;
  double dz = -999., dy = -999, dtof = -999.;
  bool isMtdTrig = kFALSE;

  // pt
  double pt = track->pt();

  // nSigmaPi cut
  nSigmaPi = track->nSigmaPion();

  // tof match
  tofIndex = track->index2BTofHit();

  // dz cut
  int iMtd = track->index2MtdHit();
  if(iMtd>=0 && iMtd<(Int_t)muDst->numberOfMTDHit())
    {
      const StMuMtdPidTraits mtdPid = track->mtdPidTraits();
      dy     = mtdPid.deltaY();
      dz     = mtdPid.deltaZ();
      dtof   = mtdPid.timeOfFlight() - mtdPid.expTimeOfFlight();
      StMuMtdHit *hit = muDst->mtdHit(iMtd);
      isMtdTrig = trig->isMtdHitFiredTrigger(hit);
    }

  return isMuonCandidate(pt, nSigmaPi, tofIndex, dz, dy, dtof, isMtdTrig);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMuonCandidate(const StPicoTrack *track, const StPicoDst *picoDst, StMtdTrigger *trig)
{
  if(!track) return kFALSE;
  double nSigmaPi = -999.;
  int    tofIndex = -1;
  double dz = -999., dy = -999, dtof = -999.;
  bool isMtdTrig = false;

  // pt
  StThreeVectorF mom;
  if(mTrackType==0) mom = track->pMom();
  if(mTrackType==1) 
    {
      float bField = picoDst->event()->bField();
      StThreeVectorF verPos = picoDst->event()->primaryVertex();
      mom = gMom(track, verPos, bField);
    }
  double pt = mom.perp();

  // nSigmaPi cut
  nSigmaPi = track->nSigmaPion();

  // tof match
  tofIndex = track->bTofPidTraitsIndex();

  // dz cut
  int iMtd = track->mtdPidTraitsIndex();
  if(iMtd>=0)
    {
      StPicoMtdPidTraits *mtdPid = picoDst->mtdPidTraits(iMtd);
      dy = mtdPid->deltaY();
      dz = mtdPid->deltaZ();
      dtof = mtdPid->deltaTimeOfFlight();

      int hitIndex = getMtdHitIndex(track, picoDst);
      StPicoMtdHit *hit = picoDst->mtdHit(hitIndex);
      isMtdTrig = hit->triggerFlag()>0 ? true : false;
    }
	  
  return isMuonCandidate(pt, nSigmaPi, tofIndex, dz, dy, dtof, isMtdTrig);
}


//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMuonCandidate(const Double_t pt, const Double_t nSigmaPi, const Int_t tofindex, const Double_t dz, const Double_t dy, const Double_t dtof, const Bool_t isTrig)
{
  if(pt < mMinMuonPt   || pt > mMaxMuonPt)               return kFALSE;
  if(!checkTrkHitDeltaZ(dz,pt))                          return kFALSE;
  if(!checkTrkHitDeltaY(dy,pt))                          return kFALSE;
  if(!checkTrkHitDeltaTof(dtof,pt))                      return kFALSE;
  if(!checkNSigmaPi(nSigmaPi))                           return kFALSE;
  if(mBTofMatch && tofindex<0)                           return kFALSE;
  if(mMtdHitTrigger && !isTrig)                          return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTriggerType(StMuDst *event)
{
  if(mMCflag>0) return 0;

  Int_t trigger = -1;

  for(Int_t j=0; j<kNtrig; j++)
    {
      if(!mTrigger[j]) continue;
      for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	{
	  if(event->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[j][i]))
	    {
	      trigger = j;
	      break;
	    }
	}
    }
  if(mYear==2014)
    {
      StMuMtdHeader *muMtdHeader = event->mtdHeader();
      if(muMtdHeader && muMtdHeader->shouldHaveRejectEvent()==1 && trigger==0)
      	trigger = -1;
    }
  return trigger;
}


//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTriggerType(StEvent *event)
{
  if(mMCflag>0) return 0;

  Int_t trigger = -1;

  for(Int_t j=0; j<kNtrig; j++)
    {
      if(!mTrigger[j]) continue;
      for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	{
	  if(event->triggerIdCollection()->nominal()->isTrigger(mTriggerIDs[j][i]))
	    {
	      trigger = j;
	      break;
	    }
	}
    }
  return trigger;

}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTriggerType(StPicoDst *pico)
{
  if(mMCflag>0) return 0;

  Int_t trigger = -1;

#if (PICOVERSION == 2016) || (PICOVERSION == 201402)
  for(Int_t j=0; j<kNtrig; j++)
    {
      if(!mTrigger[j]) continue;
      for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	{
	  if(pico->event()->isTrigger(mTriggerIDs[j][i]))
	    {
	      trigger = j;
	      break;
	    }
	}
    }
#else
  for(Int_t j=0; j<kNtrig; j++)
    {
      if(!mTrigger[j]) continue;

      if(j==0 && pico->event()->isDiMuon()) trigger = j;
      else if(j==1 && pico->event()->isSingleMuon()) trigger = j;
      else if(j==2 && pico->event()->isEMuon()) trigger = j;
      else if(j==3)
	{
	  int triggerWord = pico->event()->triggerWord();
	  // select VPD-ZDC-novtx-mon trigger
	  if( (triggerWord & (1<<11)) ||  (triggerWord & (1<<12)) )
	    trigger = j;
	}
      else trigger = -1;
    }
#endif



#if PICOVERSION == 2014
  if(mYear==2014)
    {
      StPicoMtdTrigger *mtdTrigger = (StPicoMtdTrigger*)pico->mtdTrigger(0);
      if(mtdTrigger && mtdTrigger->shouldHaveRejectEvent()==1 && trigger==0)
	trigger = -1;
    }
#endif

  return trigger;
}


//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTrgSetup(const int runId)
{
  if(mYear==2014)
    {
      for(int i=0; i<kNTrgSetup_2014; i++)
	{
	  for(int j=0; j<kRuns_TrgSetup_2014[i]; j++)
	    {
	      if(i==0 && runId == kRuns_AuAu_200_production_2014[j])      return 0;
	      if(i==1 && runId == kRuns_AuAu_200_production_low_2014[j])  return 1;
	      if(i==2 && runId == kRuns_AuAu_200_production_mid_2014[j])  return 2;
	      if(i==3 && runId == kRuns_AuAu_200_production_high_2014[j]) return 3;
	    }
	}
    }
  return 0;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isDimuonTrigger(StMuDst *event)
{
  for(UInt_t i=0; i<mTriggerIDs[0].size(); i++)
    {
      if(event->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[0][i]))
	{
	  return kTRUE;
	}
    }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isSimuonTrigger(StMuDst *event)
{
  for(UInt_t i=0; i<mTriggerIDs[1].size(); i++)
    {
      if(event->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[1][i]))
	{
	  return kTRUE;
	}
    }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isEmuonTrigger(StMuDst *event)
{
  for(UInt_t i=0; i<mTriggerIDs[2].size(); i++)
    {
      if(event->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[2][i]))
	{
	  return kTRUE;
	}
    }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isVpdNoVtxTrigger(StMuDst *event)
{
  for(UInt_t i=0; i<mTriggerIDs[3].size(); i++)
    {
      if(event->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[3][i]))
	{
	  return kTRUE;
	}
    }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isDimuonTrigger(StPicoDst *pico)
{
#if (PICOVERSION == 2016) || (PICOVERSION == 201402)
  for(UInt_t i=0; i<mTriggerIDs[0].size(); i++)
    {
      if(pico->event()->isTrigger(mTriggerIDs[0][i]))
	{
	  return kTRUE;
	}
    }
#else
  if(pico->event()->isDiMuon()) return kTRUE;
#endif

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isSimuonTrigger(StPicoDst *pico)
{
#if (PICOVERSION == 2016) || (PICOVERSION == 201402)
  for(UInt_t i=0; i<mTriggerIDs[1].size(); i++)
    {
      if(pico->event()->isTrigger(mTriggerIDs[1][i]))
	{
	  return kTRUE;
	}
    }
#else
  if(pico->event()->isSingleMuon()) return kTRUE;
#endif

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isEmuonTrigger(StPicoDst *pico)
{
#if (PICOVERSION == 2016) || (PICOVERSION == 201402)
  for(UInt_t i=0; i<mTriggerIDs[2].size(); i++)
    {
      if(pico->event()->isTrigger(mTriggerIDs[2][i]))
	{
	  return kTRUE;
	}
    }
#else
  if(pico->event()->isEMuon()) return kTRUE;
#endif

  return kFALSE;
}


//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isVpdNoVtxTrigger(StPicoDst *pico)
{
#if (PICOVERSION == 2016) || (PICOVERSION == 201402)
  for(UInt_t i=0; i<mTriggerIDs[3].size(); i++)
    {
      if(pico->event()->isTrigger(mTriggerIDs[3][i]))
	{
	  return kTRUE;
	}
    }
#else
  int triggerWord = pico->event()->triggerWord();
  if(mYear==2014)
    {
      if( (triggerWord & (1<<11)) ||  (triggerWord & (1<<12)) )
	return kTRUE;
    }
#endif

  return kFALSE;
}



//_____________________________________________________________________________
Int_t StMtdJpsiUtil::selectPrimaryVertex(StMuDst *event)
{
  Int_t nPrimVtx = event->numberOfPrimaryVertices();

  Int_t vtx_index = 0;

  if(mUsePrimVtxClosestToVpd)
    {
      vtx_index = 0;

      // Get VPD vz
      Double_t mVpdVz = 0;
      StBTofHeader *tofHeader = event->btofHeader();
      if(!tofHeader) return vtx_index;
      mVpdVz = tofHeader->vpdVz();

      //////////////////////////////////////
      // select the right vertex using VPD
      for(int i=0;i<nPrimVtx;i++)
	{
	  StMuPrimaryVertex *vtx = event->primaryVertex(i);
	  if(!vtx) continue;
	  double dz = vtx->position().z() - mVpdVz;
	  if(abs(mVpdVz)<200. && fabs(dz)<3.) 
	    {
	      vtx_index = i;
	      break;
	    }
	}
      return vtx_index;
    }

  if(mRequireMtdHitForPrimVtx)
    {
      // Find the primary vertex matched with MTD hits
      vtx_index = -1;
      for(Int_t i=0; i<nPrimVtx; i++)
	{
	  event->setVertexIndex(i);
	  Int_t nPrimTrk = event->numberOfPrimaryTracks();
	  Int_t nMtdTrk = 0;
	  for(Int_t j=0; j<nPrimTrk; j++)
	    {
	      StMuTrack* pTrack = event->primaryTracks(j);
	      if(!pTrack || !passTrack(pTrack)) continue;
	      Int_t iMtd = pTrack->index2MtdHit();
	      if(iMtd<0) continue;
	      nMtdTrk ++;
	    }

	  if(nMtdTrk>=2)
	    {
	      vtx_index = i;
	      break;
	    }
	}
      return vtx_index;
    }
  return vtx_index;
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::selectPrimaryVertex(StEvent *event)
{
  Int_t nPrimVtx = event->numberOfPrimaryVertices();

  Int_t vtx_index = 0;

  if(mUsePrimVtxClosestToVpd)
    {
      vtx_index = -1;

      // Get VPD vz
      Double_t mVpdVz = 0.;
      StBTofCollection * tofCollection = event->btofCollection(); 
      if(!tofCollection) return vtx_index;
      StBTofHeader *tofHeader = tofCollection->tofHeader();
      if(!tofHeader)     return vtx_index;
      mVpdVz  = tofHeader->vpdVz();

      // find the primary vertex closest to VPD vz
      for(Int_t i=0; i<nPrimVtx; i++)
	{
	  StPrimaryVertex *vertex = event->primaryVertex(i);
	  Double_t dz = vertex->position().z() - mVpdVz;
	  if(abs(mVpdVz)<200. && fabs(dz)<3.)
	    {
	      vtx_index = i;
	      break;
	    }
	}
      return vtx_index;
    }

  if(mRequireMtdHitForPrimVtx)
    {
      // Find the primary vertex matched with MTD hits
      vtx_index = -1;
      for(Int_t i=0; i<nPrimVtx; i++)
	{
	  StSPtrVecPrimaryTrack &  tracks = event->primaryVertex(i)->daughters();
	  Int_t nPrimTrk = tracks.size();
	  Int_t nMtdTrk = 0;
	  for(Int_t j=0; j<nPrimTrk; j++)
	    {
	      StTrack* pTrack = tracks[j];
	      if(!pTrack || !passTrack(pTrack,event->primaryVertex(i))) continue;
	      StSPtrVecTrackPidTraits& traits = pTrack->pidTraits();
	      StMtdPidTraits* mtdpid = NULL;
	      for(UInt_t it=0; it<traits.size(); it++)
		{
		  if (traits[it]->detector() == kMtdId)
		    {
		      mtdpid = dynamic_cast<StMtdPidTraits*>(traits[it]);
		      break;
		    }
		}
	      if(mtdpid)
		{
		  StMtdHit* hit = mtdpid->mtdHit();
		  if(!hit) continue;
		  nMtdTrk ++;
		}
	    }
	  if(nMtdTrk>=2)
	    {
	      vtx_index = i;
	      break;
	    }
	}
      return vtx_index;
    }
  return vtx_index;
}


//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isGoodVertex(const Double_t vpd_vz, const Double_t tpc_vz, const Double_t tpc_vr)
{
  if(mMaxVtxZ<1e4 && TMath::Abs(tpc_vz)>=mMaxVtxZ)          return kFALSE;
  if(mMaxDiffz<1e4 && TMath::Abs(tpc_vz-vpd_vz)>=mMaxDiffz) return kFALSE;
  if(mMaxVtxR<1e4 && tpc_vr>=mMaxVtxR)                      return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
double StMtdJpsiUtil::getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos) const
{
  if(!globalTrack) return 999;
  StDcaGeometry* trDcaGeom = globalTrack->dcaGeometry();
  if(!trDcaGeom) return 999;
  StPhysicalHelixD dcahh = trDcaGeom->helix();
  double pathlength = dcahh.pathLength(vtxPos, false);
  return (dcahh.at(pathlength)-vtxPos).mag();
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrack(StTrack *track, StVertex *vtx) const 
{
  if(!track) return kFALSE;
  StThreeVectorF mom = track->geometry()->momentum();
  Float_t pt = mom.perp();
  Float_t eta = mom.pseudoRapidity();
  Float_t phi = mom.phi();
  if(phi<0) phi += 2*pi;
  Int_t nHitsFit = track->fitTraits().numberOfFitPoints(kTpcId);
  Int_t nHitsPoss = track->numberOfPossiblePoints(kTpcId);

  if(pt < mMinTrkPt   || pt > mMaxTrkPt)       return kFALSE;
  if(eta < mMinTrkEta || eta > mMaxTrkEta)     return kFALSE;
  if(phi < mMinTrkPhi || phi > mMaxTrkPhi)     return kFALSE;
  if(nHitsFit<mMinNHitsFit)                    return kFALSE;
  if(1.*nHitsFit/nHitsPoss<mMinFitHitsFaction) return kFALSE;

  if(mMCflag==0 || mMCflag==3)
    {
      StTpcDedxPidAlgorithm pidAlgorithm;
      const StParticleDefinition *pd = track->pidTraits(pidAlgorithm);
      if(!pd || !pidAlgorithm.traits()) return kFALSE;
      if(pidAlgorithm.traits()->numberOfPoints()<mMinNHitsDedx) return kFALSE;

      StGlobalTrack *globalTrack = 0x0;
      if(mTrackType==0) globalTrack = dynamic_cast<StGlobalTrack*>(track->node()->track(global));
      else              globalTrack = dynamic_cast<StGlobalTrack*>(track);
      double dca = getDca(globalTrack,vtx->position());
      if(mMaxDca<1e4 && dca>mMaxDca) return kFALSE;
    }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrack(StMuTrack *track) const 
{
   StThreeVectorF mom = track->momentum();
   Float_t pt = mom.perp();
   Float_t eta = mom.pseudoRapidity();
   Float_t phi = mom.phi();
   if(phi<0) phi += 2*pi;
   Int_t nHitsFit  = track->nHitsFit(kTpcId);
 
   if(pt < mMinTrkPt   || pt > mMaxTrkPt)             return kFALSE;
   if(eta < mMinTrkEta || eta > mMaxTrkEta)           return kFALSE;
   if(phi < mMinTrkPhi || phi > mMaxTrkPhi)           return kFALSE;
   if(nHitsFit<mMinNHitsFit)                          return kFALSE;
   if(mMCflag==0 || mMCflag==3)
     {
       if(track->nHitsDedx()<mMinNHitsDedx)                   return kFALSE;
       if(mMaxDca<1e4 && track->dcaGlobal().mag()>mMaxDca)  return kFALSE;
     }
   if(nHitsFit/(1.0*track->nHitsPoss(kTpcId))<mMinFitHitsFaction) return kFALSE;
   return kTRUE;
 }


//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::passTrack(StPicoTrack *track, StThreeVectorF pVtx, float B) const 
{
  StThreeVectorF mom;
  if(mTrackType==0) mom = track->pMom();
  if(mTrackType==1) mom = gMom(track, pVtx, B);

  Float_t pt = mom.perp();
  Float_t eta = mom.pseudoRapidity();
  Float_t phi = mom.phi();
  if(phi<0) phi += 2*pi;
  Int_t nHitsFit  = track->nHitsFit();

  //cout << pt << "  " << eta << "  " << phi << "  " << nHitsFit << "  " << track->nHitsDedx() << "  " << track->dca().x() << "  " << track->dca().y() << "  " << track->dca().z() << "  " << track->dca().magnitude()<< "  " << track->dca().mag() << endl;
 
   if(pt < mMinTrkPt   || pt > mMaxTrkPt)             return kFALSE;
   if(eta < mMinTrkEta || eta > mMaxTrkEta)           return kFALSE;
   if(phi < mMinTrkPhi || phi > mMaxTrkPhi)           return kFALSE;
   if(nHitsFit<mMinNHitsFit)                          return kFALSE;
   if(mMCflag==0 || mMCflag==3)
     {
       if(track->nHitsDedx()<mMinNHitsDedx)           return kFALSE;
       if(mMaxDca<1e4 && dca(track,pVtx)>mMaxDca)  return kFALSE;
     }
   if(nHitsFit/(1.0*track->nHitsMax())<mMinFitHitsFaction) return kFALSE;
   return kTRUE;
 }

//_____________________________________________________________________________
StThreeVectorF StMtdJpsiUtil::gMom(const StPicoTrack *track, const StThreeVectorF pVtx, const float B) const
{
#if ( (PICOVERSION == 2013) || (PICOVERSION == 2015) )
  return track->gMom();
#elif ( (PICOVERSION == 2014) || (PICOVERSION == 201402) )
  return track->gMom(pVtx, B);
#elif PICOVERSION == 2016
  return track->gMom();
#endif
}

//_____________________________________________________________________________
double StMtdJpsiUtil::dca(StPicoTrack *track, StThreeVectorF pVtx) const
{
#if ( (PICOVERSION == 2013) || (PICOVERSION == 2015) )
  return track->dca();
#elif (PICOVERSION == 2014) || (PICOVERSION == 201402)
  StPhysicalHelixD dcahh = track->helix();
  return dcahh.distance(pVtx);
#elif PICOVERSION == 2016
  return (track->dca()-pVtx).mag();
#endif
}

//_____________________________________________________________________________
StPhysicalHelixD StMtdJpsiUtil::helix(const StPicoTrack *track, const float B)
{
#if ( (PICOVERSION == 2013) || (PICOVERSION == 2015) )
  return StPhysicalHelixD(track->gMom(), track->origin(), B * kilogauss, track->charge());
#elif (PICOVERSION == 2014) || (PICOVERSION == 201402)
  return track->helix();
#elif PICOVERSION == 2016
  return track->helix(B);
#endif
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::gRefMult(StEvent *event, StThreeVectorF vtxPos)
{
  int grefmult = 0;
  if(!event) return grefmult;

  StSPtrVecTrackNode &trackNodes = event->trackNodes();
  int nNodes = trackNodes.size();
  for(int i=0; i<nNodes; i++)
    {
      StGlobalTrack *pTrack = dynamic_cast<StGlobalTrack*>(trackNodes[i]->track(global));
      if(!pTrack) continue;
      //if(pTrack->idTruth()<200) continue;

      StThreeVectorF momentum = pTrack->geometry()->momentum();
      if (fabs(momentum.pseudoRapidity()) <  0.5 && 
	  fabs(getDca(pTrack,vtxPos)) < 3 &&
	  pTrack->fitTraits().numberOfFitPoints(kTpcId) >= 10)
	grefmult++;  
    }
  return grefmult;

}


//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTofMult(StMuDst *event)
{
  int nTofMatch = 0;
  int nNodes = event->numberOfPrimaryTracks();
  for(int i=0; i<nNodes; i++)
    {
      StMuTrack *pTrack = event->primaryTracks(i);
      if(pTrack->nHitsFit(kTpcId)<10) continue;
      if(pTrack->dcaGlobal().mag()>3) continue;
      if(pTrack->tofHit()) nTofMatch++;
    }
  return nTofMatch;
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTofMult(StEvent *event)
{
  int nTofMatch = 0;
  StPrimaryVertex *vertex = event->primaryVertex();
  int nTracks = vertex->numberOfDaughters();
  for(int i=0; i<nTracks; i++)
    {
      StTrack *pTrack = vertex->daughter(i);
      if(pTrack->type()!=primary) continue;
      if(pTrack->fitTraits().numberOfFitPoints(kTpcId)<10) continue;
      StGlobalTrack *globalTrack = dynamic_cast<StGlobalTrack*>(pTrack->node()->track(global));
      StDcaGeometry* trDcaGeom = globalTrack->dcaGeometry();
      if(!trDcaGeom) continue;
      StPhysicalHelixD dcahh = trDcaGeom->helix();
      double dca = dcahh.distance(vertex->position(),kFALSE);
      if(dca>3) continue;
      if(pTrack->isBToFMatched()) nTofMatch++;
    }
  return nTofMatch;
}


//_____________________________________________________________________________
double StMtdJpsiUtil::getPtDepDzCut(double pt)
{ 
  return fResDzVsPt->Eval(pt);
}

//_____________________________________________________________________________
double StMtdJpsiUtil::getPtDepDyCut(double pt)
{ 
  return fResDyVsPt->Eval(pt);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::checkTrkHitDeltaZ(const Double_t dz, const Double_t pt) const
{
  if(!mApplyDzCut) return true;

  double min = -999, max = -999;
  if(!mPtDepCutDz)
    {
      min = mMinMuonDeltaZ;
      max = mMaxMuonDeltaZ;
    }
  else
    {
      double cut = fResDzVsPt->Eval(pt);
      if(pt<3) max = mSigmaDz1 * cut;
      else     max = mSigmaDz2 * cut;

      min = -1 * max;
    }
  return (dz>=min && dz<=max);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::checkTrkHitDeltaY(const Double_t dy, const Double_t pt) const
{
  if(!mApplyDyCut) return true;

  double min = -999, max = -999;
  if(!mPtDepCutDy)
    {
      min = mMinMuonDeltaY;
      max = mMaxMuonDeltaY;
    }
  else
    {
      double cut = fResDyVsPt->Eval(pt);
      if(pt<3)  max = mSigmaDy1 * cut;
      else      max = mSigmaDy2 * cut;

      min = -1 * max;
    }
  return (dy>=min && dy<=max);
}


//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::checkTrkHitDeltaTof(const Double_t dtof, const Double_t pt) const
{
  if(!mApplyDtofCut) return true;

  double min = -999, max = -999;
  if(!mPtDepCutDtof)
    {
      min = mMinMuonDeltaTof;
      max = mMaxMuonDeltaTof;
    }
  else
    {
      max = mSigmaDtof * fResDtofVsPt->Eval(pt);
    }

  return (dtof<=max && dtof>=min);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::checkNSigmaPi(const Double_t nSigmaPi) const
{
  return (nSigmaPi>=mMinNsigmaPi && nSigmaPi<=mMaxNsigmaPi);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMtdHitInTrigWin(StMtdHit *hit, const Double_t trigger_time) const
{
  return isMtdHitInTrigWin( hit->backleg(), hit->module(), hit->leadingEdgeTime().first, trigger_time);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMtdHitInTrigWin(StMuMtdHit *hit, const Double_t trigger_time) const
{
  return isMtdHitInTrigWin( hit->backleg(), hit->module(), hit->leadingEdgeTime().first, trigger_time);
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isMtdHitInTrigWin(Int_t backleg, const Int_t module, const Double_t leading_time, const Double_t trigger_time) const
{
  if(!mTrigTimeCut) return kTRUE;
  if(mMCflag>0) return kTRUE;

  Double_t tDiff = leading_time - trigger_time;
  while(tDiff<0) tDiff += 51200;
  if(tDiff > mTrigWinCut_low[backleg-1][module-1] && tDiff < mTrigWinCut_high[backleg-1][module-1])
    return kTRUE;
  else
    return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isInSameTrigUnit(StMuMtdHit *hit1, StMuMtdHit *hit2) const
{
  return isInSameTrigUnit(hit1->backleg(), hit1->module(), hit2->backleg(), hit2->module());
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isInSameTrigUnit(StMtdHit *hit1, StMtdHit *hit2) const
{
  return isInSameTrigUnit(hit1->backleg(), hit1->module(), hit2->backleg(), hit2->module());
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isInSameTrigUnit(StPicoMtdHit *hit1, StPicoMtdHit *hit2) const
{
  return isInSameTrigUnit(hit1->backleg(), hit1->module(), hit2->backleg(), hit2->module());
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isInSameTrigUnit(Int_t backleg1, Int_t module1, Int_t backleg2, Int_t module2) const
{
  return (kTrgID[backleg1-1][module1-1] == kTrgID[backleg2-1][module2-1]);
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getTrigUnit(Int_t backleg, Int_t module) const
{
  return kTrgID[backleg-1][module-1];
}


//_____________________________________________________________________________
double StMtdJpsiUtil::m2(StPicoTrack *track, const StPicoDst *picoDst)
{
  double m2 = -1;
  double p = 0;
  if(mTrackType==0) p = track->pMom().mag();
  else              
    {
      StThreeVectorF mom = gMom(track, picoDst->event()->primaryVertex(), picoDst->event()->bField());
      p = mom.mag();
    }
  int index = track->bTofPidTraitsIndex();
  if(index>=0)
    {
      StPicoBTofPidTraits* tofPid = picoDst->btofPidTraits(index);
      double beta = tofPid->btofBeta();
      if(beta!=0) m2 = TMath::Power(p,2) * (1/TMath::Power(beta,2)-1);
    }
  return m2;
}

//_____________________________________________________________________________
double StMtdJpsiUtil::m2(StTrack *track)
{
  double m2 = -1;
  double p = track->geometry()->momentum().mag();

  StBTofPidTraits *tofPid = 0;
  StSPtrVecTrackPidTraits &traits = track->pidTraits();
  for(UInt_t itrait=0; itrait<traits.size(); itrait++)
    {
      if(traits[itrait]->detector() == kTofId)
	{
	  tofPid = dynamic_cast<StBTofPidTraits*>(traits[itrait]);
	  break;
	}
    }

  if(tofPid)
    {
      double beta = tofPid->beta();
      if(beta!=0) m2 = TMath::Power(p,2) * (1/TMath::Power(beta,2)-1);
    }
  return m2;
}

//_____________________________________________________________________________
StMtdPidTraits *StMtdJpsiUtil::getMtdPidTraits(StTrack *track)
{
  StSPtrVecTrackPidTraits &traits = track->pidTraits();
  StMtdPidTraits *mtdPid = 0;
  for(UInt_t itrait=0; itrait<traits.size(); itrait++)
    {
      if(traits[itrait]->detector() == kMtdId)
	{
	  mtdPid = dynamic_cast<StMtdPidTraits*>(traits[itrait]);
	  break;
	}
    }
  return mtdPid;
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getMtdPidTraitsIndex(const StPicoMtdHit *hit, const StPicoDst *pico)
{
  Int_t index = -1;

  Int_t nPidTraits = pico->numberOfMtdPidTraits();
  for(Int_t i=0; i<nPidTraits; i++)
    {
      StPicoMtdPidTraits *mtdPid = pico->mtdPidTraits(i);
      if(mtdPid->backleg()==hit->backleg() &&
	 mtdPid->module()==hit->module() &&
	 mtdPid->cell()==hit->cell())
	{
	  index = i;
	  break;
	}
    }
  return index;
}

//_____________________________________________________________________________
Int_t StMtdJpsiUtil::getMtdHitIndex(const StPicoTrack *track, const StPicoDst *pico)
{
  Int_t index = -1;
  if(track->mtdPidTraitsIndex()>=0)
    {
      StPicoMtdPidTraits *mtdPid = pico->mtdPidTraits(track->mtdPidTraitsIndex());
      Int_t nMtdHits = pico->numberOfMtdHits();
      for(Int_t i=0; i<nMtdHits; i++)
	{
	  StPicoMtdHit *hit = pico->mtdHit(i);
	  if(!hit) continue;
	  if(mtdPid->backleg()==hit->backleg() &&
	     mtdPid->module()==hit->module() &&
	     mtdPid->cell()==hit->cell())
	    {
	      index = i;
	      break;
	    }
	}
    }
  return index;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::checkTrkPairDca(const Double_t dr, const Double_t dz) const
{
  if(TMath::Abs(dz) < mMaxTrkPairDcaDz && TMath::Abs(dr) < mMaxTrkPairDcaDr) return kTRUE;
  else return kFALSE;
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitGlobalZ(StMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitGlobalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second, hit->module());
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitGlobalZ(StMuMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitGlobalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second, hit->module());
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitGlobalZ(StPicoMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitGlobalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second, hit->module());
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitGlobalZ(Double_t leadingWestTime, Double_t leadingEastTime, Int_t module) const
{
  Double_t z = (module-3)*stripLength - (leadingWestTime-leadingEastTime)/2./gMtdCellDriftV*1e3;
  return z;
}

//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitLocalZ(StMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitLocalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second);
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitLocalZ(StMuMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitLocalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second);
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitLocalZ(StPicoMtdHit *hit) const
{
  Int_t backleg = hit->backleg();
  if(backleg<1 || backleg>30)
    {
      LOG_WARN << "Wrong backleg id: " << backleg << endm;
      return -999;
    }
  return getMtdHitLocalZ(hit->leadingEdgeTime().first, hit->leadingEdgeTime().second);
}


//_____________________________________________________________________________
Double_t StMtdJpsiUtil::getMtdHitLocalZ(Double_t leadingWestTime, Double_t leadingEastTime) const
{
  return (leadingEastTime-leadingWestTime)/2./gMtdCellDriftV*1e3;
}

//_____________________________________________________________________________
double StMtdJpsiUtil::getMtdHitLocalY(StMuMtdHit *hit) const
{
  int backleg = hit->backleg();
  int module  = hit->module();
  int cell    = hit->cell();

  double cell_width = gMtdCellWidth+gMtdCellGap;
  double y_center = (module<4? 1 : -1) * (cell-gMtdNCells/2+0.5) * cell_width;

  if(backleg==8)  y_center -= 3 * cell_width;
  if(backleg==24) y_center += 2 * cell_width;

  return y_center;
}

//_____________________________________________________________________________
double StMtdJpsiUtil::getMtdHitLocalY(StPicoMtdHit *hit) const
{
  int backleg = hit->backleg();
  int module  = hit->module();
  int cell    = hit->cell();

  double cell_width = gMtdCellWidth+gMtdCellGap;
  double y_center = (module<4? 1 : -1) * (cell-gMtdNCells/2+0.5) * cell_width;

  if(backleg==8)  y_center -= 3 * cell_width;
  if(backleg==24) y_center += 2 * cell_width;

  return y_center;
}

//_____________________________________________________________________________
Double_t StMtdJpsiUtil::rotatePhi(Double_t phi) const
{
  Double_t outPhi = phi;
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
}

//_____________________________________________________________________________
Bool_t StMtdJpsiUtil::isValidPosition(const StThreeVectorD pos) const
{
  if(TMath::IsNaN(pos.x())||TMath::IsNaN(pos.y())||TMath::IsNaN(pos.z())) return kFALSE;
  if(TMath::Abs(pos.perp())>450.||TMath::Abs(pos.z())>300.) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
void StMtdJpsiUtil::addCutToHisto(TH1 *h, const Int_t bin, const char *label, const Float_t value) const
{
  if(!h) return;
  h->GetXaxis()->SetBinLabel(bin,label);
  if(value!=-999)
    h->SetBinContent(bin,value);
}

//_____________________________________________________________________________
void StMtdJpsiUtil::printConfig()
{
  const char *decision[2] = {"no","yes"};
  const char *track[2] = {"primary","global"};
  const char *data[2] = {"real data","simulation"};
  const char *trigName[4] = {"di-muon","single-muon","e-mu","mb"};
  printf("\n=== Configuration for StMtdJpsiUtil ===\n");
  printf("Run on %s\n",data[mMCflag]);
  for(Int_t i=0; i<4; i++)
    {
      if(mTrigger[i]) printf("Enable %s trigger\n",trigName[i]);
    }
  if(mVzCorrMethod==0) printf("Use standard Vz correction for high luminosity\n");
  if(mVzCorrMethod==1) printf("Use intrapolation Vz correction for high luminosity\n");
  printf("Require vertex rannking: %s\n",decision[mRequireVtxRanking]);
  printf("Use the primary vtx closest to VPD: %s\n",decision[mUsePrimVtxClosestToVpd]);
  printf("Use the primary vtx with MTD hits: %s\n",decision[mRequireMtdHitForPrimVtx]);
  printf("Use %s tracks\n",track[mTrackType]);
  printf("Maximum vertex r: %1.1f cm\n",mMaxVtxR);
  printf("Maximum vertex z: %1.1f cm\n",mMaxVtxZ);
  printf("|TPC-VPD| < %1.1f cm\n",mMaxDiffz);
  printf("Track pt  range: [%1.2f, %1.2f]\n",mMinTrkPt,mMaxTrkPt);
  printf("Track phi range: [%1.2f, %1.2f]\n",mMinTrkPhi,mMaxTrkPhi);
  printf("Track eta range: [%1.2f, %1.2f]\n",mMinTrkEta,mMaxTrkEta);
  printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
  printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
  printf("Minimum fraction of fit hits: %4.2f\n",mMinFitHitsFaction);
  printf("Maximum dca: %1.1f cm\n",mMaxDca);
  printf("MTD hit trig window cut: %s\n",decision[mTrigTimeCut]);
  printf("Muon PID cuts:\n");
  printf("    %1.2f < NsigmaPi < %1.2f\n",mMinNsigmaPi,mMaxNsigmaPi);
  printf("    TOF match: %s\n",decision[mBTofMatch]);
  printf("    MTD hit trigger: %s\n",decision[mMtdHitTrigger]);
  printf("    %1.2f < pt < %1.0f GeV/c\n",mMinMuonPt,mMaxMuonPt);
  if(mApplyDzCut)   
    {
      if(!mPtDepCutDz) 
	printf("    %1.0f < dz < %1.0f cm\n",mMinMuonDeltaZ,mMaxMuonDeltaZ);
      else 
	printf("    (%1.2fsigma, %1.2fsigma)\n",mSigmaDz1,mSigmaDz2);
    }
  else printf("    No dz cut\n");

  if(mApplyDyCut)   
    {
      if(!mPtDepCutDy)
	printf("    %1.0f < dy < %1.0f cm\n",mMinMuonDeltaY,mMaxMuonDeltaY);
      else
	printf("    (%1.2fsigma, %1.2fsigma)\n",mSigmaDy1,mSigmaDy2);
    }
  else  printf("    No dy cut\n");

  if(mApplyDtofCut)   
    {
      if(!mPtDepCutDtof)
	printf("    %1.1f < dtof < %2.2f ns\n",mMinMuonDeltaTof,mMaxMuonDeltaTof);
      else
	printf("    %1.1fsigma\n",mSigmaDtof);
    }
  else
    printf("    No dtof cut\n");


  if(mTrackType==1)
    {
      printf("Track pair DCA: |dca_{r}-beam_{r}|<%1.1f cm, |dca_{z}-vpd_{z}|<%1.1f cm\n",mMaxTrkPairDcaDr,mMaxTrkPairDcaDz);
    }
  printf("=======================================\n\n");
}









