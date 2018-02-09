#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <iterator>

#include "THnSparse.h"
#include "TRandom3.h"
#include "TStreamerInfo.h"
#include "TH3F.h"

#include "SystemOfUnits.h"   // has "tesla" in it

#include "StEventTypes.h"
#include "Stypes.h"
#include "headers.h"
#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcTrack.hh"
#include "StMcEvent/StMcMtdHitCollection.hh"
#include "StMcEvent/StMcMtdHit.hh"
#include "StEventUtilities/StuRefMult.hh"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "StAssociationMaker/StMcParameterDB.h"
#include "StMtdEmbedding.h"

ClassImp(StMtdEmbedding)

//_____________________________________________________________________________
StMtdEmbedding::StMtdEmbedding(const Char_t *name) : StMaker(name)
{
  // default constructor

  mEmbedParticle           = "";
  mStEvent                 = NULL;
  mMcEvent                 = NULL;
  mMuDst                   = NULL;
  mPicoDst                 = NULL;
  mRunYear                 = -1;
  mRunId                   = -1;
  mEventId                 = -1;
  mTrgSetup                = -1;
  mVtxWeight               = 1;
  mPriVtx                  = NULL;
  mTriggerType             = -1;
  mTrackType               = -1;
  mCentrality              = -1;
  mCentType                = -1;
  mgRefMult                = -1;
  mTofMult                 = -1;
  mTofMthTrk               = -1;
  mZdcRate                 = 0;
  mBbcRate                 = 0;
  mUseVpdVtx               = false;
  mVerPosition.set(999,999,999);

  mAssoMaker               = NULL;
  mRcTrackMap              = NULL;
  mMcTrackMap              = NULL;
  mRcIndices.clear();
  mMcIndices.clear();

  mPrintConfig             = kTRUE;
  mMtdUtil                 = NULL;
  mMtdGeom                 = NULL;

  mTpcTrkSmear             = 0;
  mInFileName              = "";
  mhTpcCorr                = 0x0;
  mhMtdTrigEff             = 0x0;
  mhMtdTrigElecEff         = 0x0;
  mhMtdDtofEff             = 0x0;
  memset(&mhMtdRespEmb,    0, sizeof(mhMtdRespCosmic));
  memset(&mhMtdRespCosmic, 0, sizeof(mhMtdRespCosmic));
  memset(&mhInputJpsiPt,   0, sizeof(mhInputJpsiPt));

  mRandom                  = new TRandom3();
  mRandom->SetSeed(0);
  mMinJpsiM                = 3.0;
  mMaxJpsiM                = 3.2;

  mRunEmbedHadron          = kFALSE;
  mRunProjection           = kFALSE;
  mRunSystematics          = kFALSE;
  mRunRealDataQa           = kFALSE;


  mSaveTree                = kFALSE;
  mOutTreeFileName         = "";
  mOutTree                 = NULL;
  mOutTreeFile             = NULL;
  memset(&mEmbedData, 0, sizeof(mEmbedData));
}
 
//_____________________________________________________________________________
StMtdEmbedding::~StMtdEmbedding()
{
  // default destructor
  if(mRandom)    delete mRandom;
  if(mOutTreeFile) delete mOutTreeFile;
  if(mMtdGeom)    delete mMtdGeom;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Init()
{
  if(!mMtdUtil) 
    {
      mMtdUtil = new StMtdJpsiUtil();
      mMtdUtil->printConfig();
    }

  if(mMtdUtil->getTrackType()==0) mTrackType = primary;
  else                            mTrackType = global;

  printf("\n\n+++++++++++++++++++++++++\n");
  printf("+++++ Get efficiency +++++\n");
  if(mInFileName.Length()!=0)
    {
      TFile *finput = TFile::Open(mInFileName.Data(),"read");
      if(finput->IsOpen())
	{
	  mhVtxWeight = (TH1F*)finput->Get("VtxWeight");
	  if(!mhVtxWeight)
	    printf("[e] Can't get vertex weights.\n");


	  mhTpcCorr = (TH2F*)finput->Get("hTpcCorr");
	  if(!mhTpcCorr)
	    printf("[e] Can't get TPC response efficiency.\n");


	  for(int i=0; i<30; i++)
	    {
	      for(int j=0; j<5; j++)
		{
		  mhMtdRespEmb[i][j] = (TF1*)finput->Get(Form("Embed_MtdRespEff_BL%d_Mod%d",i+1,j+1));
		  if(!mhMtdRespEmb[i][j])
		    printf("[e] Can't get MTD response efficiency for embedding (%d,%d).\n",i+1,j+1);

		  mhMtdRespCosmic[i][j] =  (TF1*)finput->Get(Form("Cosmic_MtdRespEff_BL%d_Mod%d",i+1,j+1));
		  if(!mhMtdRespCosmic[i][j])
		    printf("[e] Can't get MTD response efficiency for cosmic (%d,%d).\n",i+1,j+1);
		}
	    }

	  mhMtdTrigEff =  (TF1*)finput->Get("MtdTrigEff_FitFunc");
	  if(!mhMtdTrigEff)
	    printf("[e] Can't get muon efficiency!\n");

	  mhMtdTrigElecEff = (TF1*)finput->Get("TrigElecEff_FitFunc");
	  if(!mhMtdTrigElecEff)
	    printf("[e] Can't get trigger electronics efficiency!\n");

	  mhMtdDtofEff = (TF1*)finput->Get("DtofEff0.75_FitFunc");
	  if(!mhMtdDtofEff)
	    printf("[e] Can't get dtof cut efficiency!\n");

	  for(int j=0; j<4; j++)
	    {
	      mhInputJpsiPt[j] =  (TH1F*)finput->Get(Form("hInputJpsiShape_Cent%d",j));
	      if(!mhInputJpsiPt[j])
		printf("[e] Can't get input Jpsi shape (%d).\n",j);
	    }
	}
      else
	{
	  printf("[e] File %s can't be opened!\n",mInFileName.Data());
	}
    }
  else
    {
      printf("[w] No input file is set for efficiencies\n");
    }
  printf("+++++ Done +++++\n");
  printf("+++++++++++++++++++++++++\n\n\n");

  if(mPrintConfig) printConfig();
  bookHistos();

  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::InitRun(const Int_t runNumber)
{
  if(mRunProjection)
    {
      LOG_INFO << "Initializing MTD Geometry:" << endm;
      if(gGeoManager)
	{
	  LOG_INFO << "TGeoManager (" << gGeoManager->GetName() << "," << gGeoManager->GetTitle() << ") exists" << endm;
	}
      else
	{
	  GetDataBase("VmcGeometry");
	  if(!gGeoManager)
	    {
	      int year = GetDateTime().GetYear();
	      TString geomTag = Form("y%da",year);
	      TString ts = Form("$STAR/StarVMC/Geometry/macros/loadStarGeometry.C(\"%s\",1)",geomTag.Data());
	      int ierr=0;
	      gROOT->Macro(ts.Data(),&ierr);
	      assert(!ierr);
	    }
	}
      assert(gGeoManager);

      if(!gMtdGeometry)
	{
	  mMtdGeom = new StMtdGeometry("mtdGeometry","mtdGeometry in StMtdMatchMaker",gGeoManager);
	  mMtdGeom->Init(this);
	  LOG_INFO<<" Created a new mtdGeometry ..."<<endm;
	}
      assert(gMtdGeometry);
      mMtdGeom = gMtdGeometry;  
    }    

  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Finish()
{  
  if(mSaveTree && mOutTreeFile)
    {
      mOutTreeFile->Write();
      mOutTreeFile->Close();
      LOG_INFO << "StMtdEmbedding::Finish() -> write out tree in " << mOutTreeFileName.Data() << endm;
    }
  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Make()
{
  Int_t iret = -1;

  memset(&mEmbedData, 0, sizeof(mEmbedData));
  mRcIndices.clear();
  mMcIndices.clear();

  mStEvent=(StEvent *) GetInputDS("StEvent");
  StPicoDstMaker *picoDstMaker = (StPicoDstMaker*) GetMaker("picoDst");
  StMuDstMaker *muDstMaker = (StMuDstMaker*) GetMaker("MuDst");

  if(mStEvent)
    {
      mMcEvent = (StMcEvent*) GetDataSet("StMcEvent");
      if(!mMcEvent)
	{
	  LOG_WARN << "No McEvent is available " << endm;
	  return kStWarn;
	}

      if(muDstMaker)
	{
	  mMuDst = muDstMaker->muDst();
	}
      
      mAssoMaker = (StAssociationMaker*) GetMaker("StAssociationMaker");
      if (mAssoMaker)
	{
	  mRcTrackMap = mAssoMaker->rcTrackMap();
	  mMcTrackMap = mAssoMaker->mcTrackMap();
	}
      else
	{
	  LOG_WARN << "No association map" << endm;
	  return kStWarn;
	}
      iret = processStEvent();
    }
  else if(picoDstMaker)
    {
      mPicoDst = picoDstMaker->picoDst();
      iret = processPicoDst();
    }

  if(mSaveTree) mOutTree->Fill();

  return iret;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::processStEvent()
{ 
  mhEventStat->Fill(0.5);
  mRunId = mStEvent->runInfo()->runId();
  mEventId = mStEvent->id();
  mRandom->SetSeed(mEventId*100+mRunId);
  mRunYear = mRunId/1000000 - 1 + 2000;
  if(!mMtdUtil->passRun(mRunId)) return kStOK;
  mTrgSetup = mMtdUtil->getTrgSetup(mRunId);
  if(mEmbedParticle.Contains("Jpsi",TString::kIgnoreCase))
    {
      if(!mMtdUtil->passTrigger(mStEvent)) return kStOK;
      mTriggerType = mMtdUtil->getTriggerType(mStEvent);
    }
  else
    {
      mTriggerType = 0;
    }
  mhEventStat->Fill(1.5);
  mhEventStat->Fill(mTriggerType+2.5);

  // multiplicity & centrality
  mTofMthTrk = mMtdUtil->getTofMult(mStEvent);
  mZdcRate = mStEvent->runInfo()->zdcCoincidenceRate() * 1e-3;
  mBbcRate = mStEvent->runInfo()->bbcCoincidenceRate() * 1e-3;

  mCentrality = -1;
  mgRefMult = -1;
  mgRefMultCorr = -1;

  if(!mMtdUtil->passEvent(mStEvent)) return kStOK;
  int rc_vtx_index = mMtdUtil->selectPrimaryVertex(mStEvent);
  mPriVtx = mStEvent->primaryVertex(rc_vtx_index);
  mVerPosition = mPriVtx->position();
  double tpcVz = mVerPosition.z();
  StMuPrimaryVertex* priVertex = NULL;
  if(mMuDst)
    {
      mMuDst->setVertexIndex(rc_vtx_index);
      priVertex = mMuDst->primaryVertex();
    }

  StBTofCollection * tofCollection = mStEvent->btofCollection(); 
  if(!tofCollection) return kStOK;
  //StBTofHeader *tofHeader = tofCollection->tofHeader();
  //double vpdVz = tofHeader->vpdVz();
  mhDataVtxZ[mTriggerType]->Fill(tpcVz);
  mhTofMthTrksVsZdcRate[mTriggerType]->Fill(mZdcRate,mTofMthTrk);
  mhTofMthTrksVzBbcRate[mTriggerType]->Fill(mBbcRate,mTofMthTrk);

  // vertex weight
  mVtxWeight = mhVtxWeight->GetBinContent(mhVtxWeight->FindFixBin(tpcVz));

  StMcVertex *mcVtx = mMcEvent->primaryVertex();
  if(mcVtx)
    {
      mhMcVtxZVsDataVtxZ[mTriggerType]->Fill(tpcVz, mcVtx->position().z());
      mhSetupVsMcVtxZ[mTriggerType]->Fill(mcVtx->position().z(),mTrgSetup);
    }


  mTofMult = 0;
  StSPtrVecBTofHit& tofVec = tofCollection->tofHits();
  int nTofHits = tofVec.size();
  for(int i=0; i<nTofHits; i++)
    {
      StBTofHit *tofHit = (StBTofHit*)tofVec[i];
      if(tofHit->idTruth()==0) mTofMult++;
    }

  // Get centrality bin
  double eventWeight = 1;
  if(mRunYear==2014 || mRunYear ==2016)
    {
      if(mMuDst)
	{
	  mgRefMult = mMuDst->event()->grefmult();
	}
      else
	{
	  mgRefMult = mMtdUtil->gRefMult(mStEvent,mVerPosition);
	}
      mgRefMultCorr = mMtdUtil->getRefMultCorr(mRunYear, mRunId, mgRefMult, tpcVz, mZdcRate);
      mCentrality = mMtdUtil->getCentralityBin16(mRunYear, mZdcRate, mgRefMultCorr);
      eventWeight = mMtdUtil->getWeight(mRunYear, mZdcRate, mgRefMultCorr);
      if(mCentrality>=12 && mCentrality<=15) mCentType = 0;
      if(mCentrality>=8 && mCentrality<=11)  mCentType = 1;
      if(mCentrality>=4 && mCentrality<=7)   mCentType = 2;
      if(mCentrality>=0 && mCentrality<=3)   mCentType = 3;
    }
  else
    {
      mCentType = 0;
    }

  // apply additional event weight for 2014 analysis
  mVtxWeight *= eventWeight;

  StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();
  int nGlobalTracks = 0;
  Int_t nNodes = trackNodes.size();
  for(Int_t i=0; i<nNodes; i++)
    {
      StTrack *pTrack = trackNodes[i]->track(global);
      if(!pTrack) continue;
      if(mRcIndices.count(i)>0) continue;
      nGlobalTracks++;
    }

  if(mTriggerType==0 && mTrgSetup>-1)
    {
      mhTofMultVsgRefMult[mTrgSetup]->Fill(mgRefMult, mTofMult);
      mhTofMultVsgRefMultCorr[mTrgSetup]->Fill(mgRefMultCorr, mTofMult);
      mhNgTrkVsCent[mTrgSetup]->Fill(mCentrality, nGlobalTracks);
      mhZdcRateVsCent[mTrgSetup]->Fill(mCentrality, mZdcRate);
    }


  // reject events with extremely large or small
  // TofMult/gRefMultCorr values
  if(mRunYear==2014)
    {
      if(mTofMult>30*mgRefMultCorr || mTofMult<0.5*mgRefMultCorr) return kStOK;
    }

  //printf("[i] Run = %d, event = %d, gRefMult = %d, vz = %4.2f, ZDC = %4.2f, gRefMultCorr = %4.2f, cent = %d\n",
  //	 mRunId,  mStEvent->id(), mgRefMult, tpcVz, mZdcRate, mgRefMultCorr, mCentrality);
  // cout << mcVtx->position().z() << " =? " << tpcVz << " =? " << priVertex->position().z() << endl;
  // cout << mgRefMult << " =? " << mMuDst->event()->grefmult() << endl;

  mhCentrality[mTriggerType]->Fill(mCentrality);
  mhgRefMult[mTriggerType]->Fill(mgRefMult);
  mhgRefMultCorr[mTriggerType]->Fill(mgRefMultCorr);

  mhEventStat->Fill(2.5+3);
  mhEventStat->Fill(mTriggerType+3.5+3);
  double fill[] = {mZdcRate, mCentrality*1., mTrgSetup*1.};
  mhZdcRateVsTrigSetup[mTriggerType]->Fill(fill);

  //-------------------------------------------------------------
  // build map between mc and reconstructed tracks
  StSPtrVecMcTrack mctracks = mMcEvent->tracks();
  Int_t nMcTracks = mctracks.size();
  LOG_DEBUG << "# of mc tracks: " << nMcTracks << endm;
  for(Int_t i=0; i<nMcTracks; i++) 
    {
      StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
      if(!mctrack) continue;
      Int_t rcIndex = findMatchedRcTrack(mctrack);
      mMcIndices[i] = rcIndex;
      if(rcIndex>=0) 
	{
	  LOG_DEBUG << "mc track " << i << " is matched to rc track " << rcIndex << endm;
	  mRcIndices[rcIndex] = i;
	}
    }

  if(mRunEmbedHadron) 
    {
      int iret = runEmbedHadron();
      if(!mEmbedParticle.Contains("Jpsi",TString::kIgnoreCase))
	return iret;
    }

  if(mRunProjection) runProjection();

  StMtdCollection* mtdCollection = mStEvent->mtdCollection();
  StSPtrVecMtdHit& mtdVec = mtdCollection->mtdHits();  
  Int_t nHits = mtdVec.size();
  for(int i=0; i<nHits; i++)
    {
      StMtdHit *hit = (StMtdHit*)mtdVec[i];
      int backleg = hit->backleg();
      int gchannel = 12 * (hit->module()-1) + hit->cell();
      if(isMcMtdHit(hit))  mhMcHitMap[mTriggerType]->Fill(backleg,gchannel);
      else mhRealHitMap[mTriggerType]->Fill(backleg,gchannel);
    }

  //-------------------------------------------------------------
  // Reset the dy in the PidTraits for MC
  const double cell_width = 4.4; //cm
  for(Int_t i=0; i<nNodes; i++)
    {
      StTrack *pTrack = trackNodes[i]->track(mTrackType);
      if(!pTrack) continue;
      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(pTrack);
      if(!mtdPid) continue;
      StMtdHit* hit = mtdPid->mtdHit();
      if(!isMcMtdHit(hit)) continue;
      int backleg = hit->backleg();
      double dy = mtdPid->deltaY();
      if(backleg==8)  dy -= 3 * cell_width;
      if(backleg==24) dy += 2 * cell_width;
      mtdPid->setDeltaY(dy);
    }
  //-------------------------------------------------------------
  

  // MC truth
  for(Int_t i=0; i<nMcTracks; i++) 
    {
      StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
      if(!mctrack) continue;
      double pt = mctrack->pt();
      double eta = mctrack->pseudoRapidity();
      double phi = mMtdUtil->rotatePhi(mctrack->momentum().phi());
      int rcIndex = mMcIndices[i];
      if(!isMcMuon(mctrack)) continue;

      double charge = mctrack->particleDefinition()->charge();
      double mc_fill[] = {pt, eta, phi, charge, mCentrality*1.};
      double mc_fill_2[] = {pt, eta, mCentrality*1., mZdcRate, mTrgSetup*1.};
      mhMcTrkInfo[mTriggerType]->Fill(mc_fill, mVtxWeight);
      mhMcTrkPtEff[mTriggerType][0]->Fill(mc_fill_2, mVtxWeight);
      if(rcIndex<0) continue;
	    
      StTrack *rcTrack = trackNodes[rcIndex]->track(mTrackType);
      if(!mMtdUtil->passTrack(rcTrack,mPriVtx)) continue;
      StThreeVectorF momentum = rcTrack->geometry()->momentum();

      // get momentum resolution
      double rc_pt_1 = momentum.perp();
      double rc_pt_2 = 0;
      if(mTrackType==primary)
	{
	  double fill_p[] = {(pt-rc_pt_1)/pt, rc_pt_1, pt, mCentrality*1.};
	  mhpTrkPtRes[mTriggerType]->Fill(fill_p);
	  StTrack *rcTrack2 = trackNodes[rcIndex]->track(global);
	  if(rcTrack2) 
	    {
	      rc_pt_2 = rcTrack2->geometry()->momentum().perp();
	      double fill_g[] = {(pt-rc_pt_2)/pt, rc_pt_2, pt, mCentrality*1.};
	      mhgTrkPtRes[mTriggerType]->Fill(fill_g);
	    }
	}
      else
	{
	  double fill_g[] = {(pt-rc_pt_1)/pt, rc_pt_1, pt, mCentrality*1.};
	  mhgTrkPtRes[mTriggerType]->Fill(fill_g);
	  StTrack *rcTrack2 = trackNodes[rcIndex]->track(primary);
	  if(rcTrack2) 
	    {
	      rc_pt_2 = rcTrack2->geometry()->momentum().perp();
	      double fill_p[] = {(pt-rc_pt_2)/pt, rc_pt_2, pt, mCentrality*1.};
	      mhpTrkPtRes[mTriggerType]->Fill(fill_p);
	    }
	}

      double rc_pt = getSmearedTrkPt(momentum.perp());
      if(mRunSystematics) rc_pt = momentum.perp();
      double rc_eta       = momentum.pseudoRapidity();
      double rc_phi       = mMtdUtil->rotatePhi(momentum.phi());
      if(!mRunSystematics && !isTpcResp(rc_pt, rc_eta, rc_phi)) continue;
      //double rc_fill[] = {rc_pt, rc_eta, rc_phi, charge, mCentrality*1.};
      double rc_fill_2[] = {rc_pt, rc_eta, mCentrality*1., mZdcRate, mTrgSetup*1.};
      mhMcTrkInfoTpc[mTriggerType]->Fill(mc_fill, mVtxWeight);
      mhMcTrkPtEff[mTriggerType][1]->Fill(rc_fill_2, mVtxWeight);

      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(rcTrack);
      if(!mtdPid) continue;
      mhMcTrkInfoMtd[mTriggerType]->Fill(mc_fill, mVtxWeight);

      StMtdHit* hit = mtdPid->mtdHit();
      int backleg = hit->backleg();
      int module = hit->module();
      mhTrkPtBlModMtd[mTriggerType]->Fill(rc_pt, backleg, module);
      
      if(!mRunSystematics && !isMtdResp(rc_pt, backleg, module)) continue;
      mhMcTrkPtEff[mTriggerType][2]->Fill(rc_fill_2, mVtxWeight);

      bool mcMtdHit = isMcMtdHit(hit);
      if(!mcMtdHit) mhMcTrkPtEff[mTriggerType][3]->Fill(rc_fill_2, mVtxWeight);

      if(!isRcMuon(rcTrack)) continue;
      if(!mRunSystematics && !passPidCut(rc_pt)) continue;
      mhMcTrkPtEff[mTriggerType][4]->Fill(rc_fill_2, mVtxWeight);

      if(!isMtdTrig(rc_pt)) continue;
      mhMcTrkInfoFinal[mTriggerType]->Fill(mc_fill, mVtxWeight);
      mhMcTrkPtEff[mTriggerType][5]->Fill(rc_fill_2, mVtxWeight);
      mhMcHitMapWithEff[mTriggerType]->Fill(hit->backleg(),12 * (hit->module()-1) + hit->cell());
      mhTrkPtBlModTrig[mTriggerType]->Fill(rc_pt, backleg, module);
    }


  int nEmbedJpsi = 0;
  for(Int_t i=0; i<nMcTracks; i++) 
    {
      StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
      if(!mctrack || !isMcMuon(mctrack)) continue;
      const StThreeVectorF mom = mctrack->momentum();
      TLorentzVector muon1(mom.x(),mom.y(),mom.z(),mctrack->energy());
      Int_t gq = mctrack->particleDefinition()->charge();
      for(Int_t j=i+1; j<nMcTracks; j++) 
	{
	  StMcTrack *mctrack2 = dynamic_cast<StMcTrack *>(mctracks[j]);
	  if(!mctrack2 || !isMcMuon(mctrack2)) continue;
	  if(!mctrack->parent() || !mctrack2->parent() ||
	     mctrack->parent()->key() != mctrack2->parent()->key()) continue;

	  const StThreeVectorF mom2 = mctrack2->momentum();
	  TLorentzVector muon2(mom2.x(),mom2.y(),mom2.z(),mctrack2->energy());
	  TLorentzVector parent = muon1 + muon2;
	  double mc_pt = parent.Pt();
	  Int_t gq2 = mctrack2->particleDefinition()->charge();
	  if(gq*gq2>0) continue;
	  nEmbedJpsi++;
	  fillInvMass(muon1,muon2,mhJpsiInfo[mTriggerType][0],mc_pt);

	  if(mMcIndices[i]<0 || mMcIndices[j]<0) continue;
	  // matched to recontructed TPC tracks
	  StTrack *rcTrack1 = trackNodes[mMcIndices[i]]->track(mTrackType);
	  StTrack *rcTrack2 = trackNodes[mMcIndices[j]]->track(mTrackType);
	  if(!mMtdUtil->passTrack(rcTrack1,mPriVtx) || !mMtdUtil->passTrack(rcTrack2,mPriVtx)) continue;
	  const StThreeVectorF imom = rcTrack1->geometry()->momentum();
	  TLorentzVector muoni;
	  muoni.SetXYZM(imom.x(),imom.y(),imom.z(),muMass);
		      
	  const StThreeVectorF jmom = rcTrack2->geometry()->momentum();
	  TLorentzVector muonj;
	  muonj.SetXYZM(jmom.x(),jmom.y(),jmom.z(),muMass);
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][7],mc_pt);
	  muoni.SetPtEtaPhiM(getSmearedTrkPt(muoni.Pt()), muoni.Eta(), muoni.Phi(), muMass);
	  muonj.SetPtEtaPhiM(getSmearedTrkPt(muonj.Pt()), muonj.Eta(), muonj.Phi(), muMass);

	  double rc_pt1        = muoni.Pt();
	  double rc_eta1       = muoni.Eta();
	  double rc_phi1       = mMtdUtil->rotatePhi(muoni.Phi());
	  double rc_pt2        = muonj.Pt();
	  double rc_eta2       = muonj.Eta();
	  double rc_phi2       = mMtdUtil->rotatePhi(muonj.Phi());
	  if(!isTpcResp(rc_pt1, rc_eta1, rc_phi1) || !isTpcResp(rc_pt2, rc_eta2, rc_phi2)) continue;

	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][1],mc_pt);

	  StMtdPidTraits *mtdPid1 = mMtdUtil->getMtdPidTraits(rcTrack1);
	  StMtdPidTraits *mtdPid2 = mMtdUtil->getMtdPidTraits(rcTrack2);
	  if(!mtdPid1 || !mtdPid2) continue;
	  StMtdHit* hit1 = mtdPid1->mtdHit();
	  StMtdHit* hit2 = mtdPid2->mtdHit();
	  if(!isMtdResp(rc_pt1, hit1->backleg(), hit1->module()) || !isMtdResp(rc_pt2, hit2->backleg(), hit2->module())) continue;
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][2],mc_pt);
	  
	  if(!isMcMtdHit(hit1) || !isMcMtdHit(hit2)) continue;
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][3],mc_pt);

	  if(!isRcMuon(rcTrack1) || !passPidCut(rc_pt1) || !isRcMuon(rcTrack2) || !passPidCut(rc_pt2)) continue;
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][4],mc_pt);

	  if(!isMtdTrig(imom.perp()) || !isMtdTrig(jmom.perp())) continue;
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][5],mc_pt);

	  if(mMtdUtil->isInSameTrigUnit(hit1, hit2)) continue;
	  fillInvMass(muoni,muonj,mhJpsiInfo[mTriggerType][6],mc_pt);

	  // match the mc and reco jpsi for spectrum weighting
	  TLorentzVector jpsi_mc = muon1 + muon2;
	  double jpsi_pt_mc = jpsi_mc.Pt();
	  double jpsi_y_mc  = jpsi_mc.Rapidity();
	  TLorentzVector jpsi_rc = muoni + muonj;
	  double jpsi_pt_rc = jpsi_rc.Pt();
	  double jpsi_m_rc  = jpsi_rc.M();
	  double jpsi_y_rc  = jpsi_rc.Rapidity();
	  double pt1_rc = (muoni.Pt() > muonj.Pt()) ? muoni.Pt() : muonj.Pt();
	  double pt2_rc = (muoni.Pt() < muonj.Pt()) ? muoni.Pt() : muonj.Pt();
	  if(jpsi_pt_mc>0.15 && jpsi_pt_rc>0.15 &&
	     fabs(jpsi_y_mc) <= 0.5 &&fabs(jpsi_y_rc) <= 0.5 &&
	     jpsi_m_rc <= 3.4 && jpsi_m_rc >= 2.8 &&
	     pt1_rc >= 1.5)
	    {
	      double fill_mth[] = {jpsi_pt_mc, jpsi_pt_rc, pt2_rc, mCentrality*1., mZdcRate};
	      mhJpsiMatch[mTriggerType]->Fill(fill_mth, mVtxWeight); 
	    }

	}
    }
  mhNEmbedJpsi[mTriggerType]->Fill(nEmbedJpsi);
 
  // Find reconstructed mc tracks
  IntVec trkId;
  trkId.clear();
  LOG_DEBUG << "Number of nodes: " << nNodes << endm;
  for(Int_t i=0; i<nNodes; i++)
    {
      StTrack *pTrack = trackNodes[i]->track(mTrackType);
      if(!pTrack) continue;
      if(mRcIndices.count(i)==0) continue;

      StGlobalTrack  *globalTrack = 0;
      if(mTrackType==primary) 
	{
	  globalTrack = dynamic_cast<StGlobalTrack*>(pTrack->node()->track(global));
	}
      else 
	{
	  globalTrack = dynamic_cast<StGlobalTrack*>(pTrack);
	}
      StThreeVectorF momentum = pTrack->geometry()->momentum();
      double pt        = momentum.perp();
      double eta       = momentum.pseudoRapidity();
      double phi       = mMtdUtil->rotatePhi(momentum.phi());
      double dca       = mMtdUtil->getDca(globalTrack,mPriVtx->position());
      int    nHitsFit  = pTrack->fitTraits().numberOfFitPoints(kTpcId);
      int    nHitsPoss = pTrack->numberOfPossiblePoints(kTpcId);
      int    nHitsDedx = 0;
      double nSigmaPi  = -999.;
      static StTpcDedxPidAlgorithm pidAlgorithm;
      static StPionPlus* Pion = StPionPlus::instance();
      const StParticleDefinition *pd = pTrack->pidTraits(pidAlgorithm);
      if(pd && pidAlgorithm.traits())
	{
	  nHitsDedx = pidAlgorithm.traits()->numberOfPoints();
	  nSigmaPi  = pidAlgorithm.numberOfSigma(Pion);
	}
      if(!mMtdUtil->passTrack(pTrack,mPriVtx)) continue;

      mhTrkDca[mTriggerType][0]->Fill(pt,dca);
      mhTrkNHitsFit[mTriggerType][0]->Fill(pt,nHitsFit);
      mhTrkNHitsDedx[mTriggerType][0]->Fill(pt,nHitsDedx);
      mhTrkNHitsPoss[mTriggerType][0]->Fill(pt,nHitsPoss,mTrgSetup);
      mhTrkNHitsFrac[mTriggerType][0]->Fill(pt,nHitsFit*1./nHitsPoss);
      double fill[] = {pt, eta, phi};
      mhTrkEtaPhi[mTriggerType][0]->Fill(fill);
      mhTrkNSigmaPi[mTriggerType][0]->Fill(pt,nSigmaPi);

      double fillNsigmaPi[] = {pt, nSigmaPi, eta, mVerPosition.z(), mCentrality*1., mTrgSetup*1.};
      mhRcTrkNsigmaPi[mTriggerType]->Fill(fillNsigmaPi);

      double fillNHitsPoss[] = {pt,nHitsPoss*1.,mVerPosition.z(), mTrgSetup*1.};
      mhQaNHitsPoss[mTriggerType]->Fill(fillNHitsPoss);

      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(pTrack);
      if(!mtdPid) continue;
      StMtdHit* hit = mtdPid->mtdHit();
      Bool_t isMcHit = isMcMtdHit(hit);
      if(isMcHit)
	{
	  double dy = mtdPid->deltaY();
	  double dz = mtdPid->deltaZ();
	  mhTrkDyVsPt[mTriggerType][0]->Fill(pt,dy);
	  mhTrkDzVsPt[mTriggerType][0]->Fill(pt,dz);
	  double tof_mc = mtdPid->timeOfFlight();
	  double tof_exp = mtdPid->expTimeOfFlight();
	  mhTrkDtofVsPt[mTriggerType][0]->Fill(pt,tof_mc-tof_exp);
	  mhTrkDtofVsMod[mTriggerType][0]->Fill((hit->backleg()-1)*5.+hit->module(),tof_mc-tof_exp);
	}
    }

  // Find pion candidates from real data
  int nTracks = mPriVtx->numberOfDaughters();
  for(int i=0; i<nTracks; i++)
    {
      StTrack *pTrack = mPriVtx->daughter(i);
      if(!pTrack) continue;
      if(pTrack->type()!=mTrackType) continue;
      if(pTrack->idTruth()<200) continue;

      StGlobalTrack *globalTrack = 0x0;
      if(mTrackType==primary) globalTrack = dynamic_cast<StGlobalTrack*>(pTrack->node()->track(global));
      else                    globalTrack = dynamic_cast<StGlobalTrack*>(pTrack);

      StThreeVectorF momentum = pTrack->geometry()->momentum();
      double pt        = momentum.perp();
      double eta       = momentum.pseudoRapidity();
      double phi       = mMtdUtil->rotatePhi(momentum.phi());
      double dca       = mMtdUtil->getDca(globalTrack,mPriVtx->position());
      int    nHitsFit  = pTrack->fitTraits().numberOfFitPoints(kTpcId);
      int    nHitsPoss = pTrack->numberOfPossiblePoints(kTpcId);
      int    nHitsDedx = 0;
      double nSigmaPi  = -999.;
      static StTpcDedxPidAlgorithm pidAlgorithm;
      static StPionPlus* Pion = StPionPlus::instance();
      const StParticleDefinition *pd = pTrack->pidTraits(pidAlgorithm);
      if(pd && pidAlgorithm.traits())
	{
	  nHitsDedx = pidAlgorithm.traits()->numberOfPoints();
	  nSigmaPi  = pidAlgorithm.numberOfSigma(Pion);
	}
      double fill[] = {pt, eta, phi};
      
      if(mRunRealDataQa)
	{
	  int trg_index = mTrgSetup - 1;
	  int tpc_index = 0;
	  if(phi<=6.1 && phi>=5.6 && eta <=0.2) tpc_index = 1;
	  if(trg_index>-1)
	    {
	      mhQaTrkEtaPhi[mTriggerType][trg_index]->Fill(fill);
	      mhQaTrkNHitsFit[mTriggerType][trg_index][tpc_index]->Fill(pt, nHitsFit);
	      mhQaTrkNHitsDedx[mTriggerType][trg_index][tpc_index]->Fill(pt, nHitsDedx);
	      mhQaTrkDca[mTriggerType][trg_index][tpc_index]->Fill(pt, dca);
	      mhQaTrkNSigmaPi[mTriggerType][trg_index][tpc_index]->Fill(pt, nSigmaPi);
	    }
	}
      if(!mMtdUtil->passTrack(pTrack,mPriVtx)) continue;

      if(mRunRealDataQa)
	{
	  int trg_index = mTrgSetup - 1;
	  if(trg_index>-1)
	    {
	      mhQaTrkEtaPhiWithCuts[mTriggerType][trg_index]->Fill(fill);
	    }
	}
      
      if(fabs(nSigmaPi)>2) continue;

      mhTrkDca[mTriggerType][1]->Fill(pt,dca);
      mhTrkNHitsFit[mTriggerType][1]->Fill(pt,nHitsFit);
      mhTrkNHitsDedx[mTriggerType][1]->Fill(pt,nHitsDedx);
      mhTrkNHitsPoss[mTriggerType][1]->Fill(pt,nHitsPoss,mTrgSetup);
      mhTrkNHitsFrac[mTriggerType][1]->Fill(pt,nHitsFit*1./nHitsPoss);

      mhTrkEtaPhi[mTriggerType][1]->Fill(fill);
      mhTrkNSigmaPi[mTriggerType][1]->Fill(pt,nSigmaPi);

      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(pTrack);
      if(!mtdPid) continue;
      StMtdHit* hit = mtdPid->mtdHit();
      if(isMcMtdHit(hit)) continue;
      double dy = mtdPid->deltaY();
      double dz = mtdPid->deltaZ();
      double dtof = mtdPid->timeOfFlight()-mtdPid->expTimeOfFlight();
      mhTrkDyVsPt[mTriggerType][1]->Fill(pt,dy);
      mhTrkDzVsPt[mTriggerType][1]->Fill(pt,dz);
      mhTrkDtofVsPt[mTriggerType][1]->Fill(pt,dtof);
      mhTrkDtofVsMod[mTriggerType][1]->Fill((hit->backleg()-1)*5.+hit->module(),dtof);
    }

  if(mRunSystematics)
    {
      //-------------------------------------------------------------
      // move back dy in the PidTraits when running systematic uncertainty study
      for(Int_t i=0; i<nNodes; i++)
	{
	  StTrack *pTrack = trackNodes[i]->track(mTrackType);
	  if(!pTrack) continue;
	  StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(pTrack);
	  if(!mtdPid) continue;
	  StMtdHit* hit = mtdPid->mtdHit();
	  if(!isMcMtdHit(hit)) continue;
	  int backleg = hit->backleg();
	  double dy = mtdPid->deltaY();
	  if(backleg==8)  dy += 3 * cell_width;
	  if(backleg==24) dy -= 2 * cell_width;
	  mtdPid->setDeltaY(dy);
	}
      //-------------------------------------------------------------
    }
  
  return kStOK;

}

//_____________________________________________________________________________
Int_t StMtdEmbedding::processPicoDst()
{
  mhEventStat->Fill(0.5);

  StPicoEvent *picoEvent = mPicoDst->event();
  if(!picoEvent) return kStWarn;
  float bField = picoEvent->bField();
  if(!mMtdUtil->passTrigger(mPicoDst)) return kStOK;
  mTriggerType = mMtdUtil->getTriggerType(mPicoDst);
  mhEventStat->Fill(1.5);
  mhEventStat->Fill(mTriggerType+2.5);

  if(!mMtdUtil->passEvent(picoEvent)) return kStOK;
  mhEventStat->Fill(2.5+3);
  mhEventStat->Fill(mTriggerType+3.5+3);
  StThreeVectorF verPos = picoEvent->primaryVertex();

  // TPC tracks
  IntVec trkId;
  trkId.clear();
  Int_t nNodes = mPicoDst->numberOfTracks();
  for(Int_t i=0; i<nNodes; i++)
    {
      StPicoTrack *pTrack = mPicoDst->track(i);
      if(!pTrack) continue;
      StThreeVectorF momentum;
      if(mTrackType==global) 
	momentum = mMtdUtil->gMom(pTrack,verPos, bField);
      else if(mTrackType==primary)
	{
	  momentum = pTrack->pMom();
	  if(momentum.mag()<=0) continue;
	}
      if(!mMtdUtil->passTrack(pTrack,verPos, bField)) continue;

      double pt        = momentum.perp();
      double eta       = momentum.pseudoRapidity();
      double phi       = mMtdUtil->rotatePhi(momentum.phi());
      int    nHitsFit  = pTrack->nHitsFit();
      int    nHitsPoss = pTrack->nHitsMax();
      int    nHitsDedx = pTrack->nHitsDedx();
      double nSigmaPi  = pTrack->nSigmaPion();
      double m2        = mMtdUtil->m2(pTrack, mPicoDst);
      double dca       = mMtdUtil->dca(pTrack, verPos);

      if(m2>0.016 && m2<0.021 && abs(nSigmaPi)<2)
	{
	  mhTrkDca[mTriggerType][1]->Fill(pt,dca);
	  mhTrkNHitsFit[mTriggerType][1]->Fill(pt,nHitsFit);
	  mhTrkNHitsDedx[mTriggerType][1]->Fill(pt,nHitsDedx);
	  mhTrkNHitsFrac[mTriggerType][1]->Fill(pt,nHitsFit*1./nHitsPoss);
	  double fill[] = {pt, eta, phi};
	  mhTrkEtaPhi[mTriggerType][1]->Fill(fill);
	  mhTrkNSigmaPi[mTriggerType][1]->Fill(pt,nSigmaPi);
	}

      int iMtdPid = pTrack->mtdPidTraitsIndex();
      if(iMtdPid<0) continue;
      StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(iMtdPid);
      Double_t dy = mtdPid->deltaY();
      Double_t dz = mtdPid->deltaZ();

      if(m2>0.016 && m2<0.021 && abs(nSigmaPi)<2)
	{
	  mhTrkDyVsPt[mTriggerType][1]->Fill(pt,dy);
	  mhTrkDzVsPt[mTriggerType][1]->Fill(pt,dz);
	}

      if(!mMtdUtil->isMuonCandidate(pTrack, mPicoDst)) continue;
      trkId.push_back(i);
    }

  // Look at Jpsi muons
  UInt_t nMuon = trkId.size();
  LOG_DEBUG << nMuon << " muon candidates." << endm;
  StPicoTrack *ipTrack = 0x0, *jpTrack = 0x0;
  for(UInt_t i=0; i<nMuon; i++)
    {
      ipTrack = mPicoDst->track(trkId[i]);
      StThreeVectorF mom1;
      if(mTrackType==primary)       mom1 = ipTrack->pMom();
      else if (mTrackType==global)  mom1 = mMtdUtil->gMom(ipTrack,verPos, bField);
      Int_t q1 = ipTrack->charge();
      double pt1 = mom1.perp();
      TLorentzVector muon1;
      muon1.SetXYZM(mom1.x(),mom1.y(),mom1.z(),muMass);
      for(UInt_t j=i+1; j<nMuon; j++)
	{
	  jpTrack = mPicoDst->track(trkId[j]);
	  StThreeVectorF mom2;;
	  if(mTrackType==primary)       mom2 = jpTrack->pMom();
	  else if (mTrackType==global)  mom2 = mMtdUtil->gMom(jpTrack,verPos, bField);
	  Int_t q2 = jpTrack->charge();
	  double pt2 = mom2.perp();
	  TLorentzVector muon2;
	  muon2.SetXYZM(mom2.x(),mom2.y(),mom2.z(),muMass);

	  TLorentzVector jpsi = muon1 + muon2;
	  if(jpsi.M()>mMaxJpsiM || jpsi.M()<mMinJpsiM) continue;

	  int type = 2;
	  if(q1*q2>0) type = 3;
	  mhTrkDca[mTriggerType][type]->Fill(pt1,mMtdUtil->dca(ipTrack, verPos));
	  mhTrkNHitsFit[mTriggerType][type]->Fill(pt1,ipTrack->nHitsFit());
	  mhTrkNHitsDedx[mTriggerType][type]->Fill(pt1,ipTrack->nHitsDedx());
	  mhTrkNHitsFrac[mTriggerType][type]->Fill(pt1,ipTrack->nHitsFit()*1./ipTrack->nHitsMax());
	  mhTrkNHitsPoss[mTriggerType][type]->Fill(pt1,ipTrack->nHitsMax(),mTrgSetup);
	  double fill1[] = {pt1, mom1.pseudoRapidity(), mMtdUtil->rotatePhi(mom1.phi())};
	  mhTrkEtaPhi[mTriggerType][type]->Fill(fill1);
	  mhTrkNSigmaPi[mTriggerType][type]->Fill(pt1,ipTrack->nSigmaPion());
	  mhTrkDyVsPt[mTriggerType][type]->Fill(pt1,mPicoDst->mtdPidTraits(ipTrack->mtdPidTraitsIndex())->deltaY());
	  mhTrkDzVsPt[mTriggerType][type]->Fill(pt1,mPicoDst->mtdPidTraits(ipTrack->mtdPidTraitsIndex())->deltaZ());
	      
	  mhTrkDca[mTriggerType][type]->Fill(pt2,mMtdUtil->dca(jpTrack, verPos));
	  mhTrkNHitsFit[mTriggerType][type]->Fill(pt2,jpTrack->nHitsFit());
	  mhTrkNHitsDedx[mTriggerType][type]->Fill(pt2,jpTrack->nHitsDedx());
	  mhTrkNHitsFrac[mTriggerType][type]->Fill(pt2,jpTrack->nHitsFit()*1./jpTrack->nHitsMax());
	  mhTrkNHitsPoss[mTriggerType][type]->Fill(pt2,jpTrack->nHitsMax(),mTrgSetup);
	  double fill2[] = {pt2, mom2.pseudoRapidity(), mMtdUtil->rotatePhi(mom2.phi())};
	  mhTrkEtaPhi[mTriggerType][type]->Fill(fill2);
	  mhTrkNSigmaPi[mTriggerType][type]->Fill(pt2,jpTrack->nSigmaPion());
	  mhTrkDyVsPt[mTriggerType][type]->Fill(pt2,mPicoDst->mtdPidTraits(jpTrack->mtdPidTraitsIndex())->deltaY());
	  mhTrkDzVsPt[mTriggerType][type]->Fill(pt2,mPicoDst->mtdPidTraits(jpTrack->mtdPidTraitsIndex())->deltaZ());
	}
    }
  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::runEmbedHadron()
{
  StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();

  StSPtrVecMcTrack mctracks = mMcEvent->tracks();
  Int_t nMcTracks = mctracks.size();
  LOG_DEBUG << "# of mc tracks: " << nMcTracks << endm;
  int trackCount = 0;
  for(Int_t i=0; i<nMcTracks; i++) 
    {
      StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
      if(!mctrack) continue;
      int geantId = mctrack->geantId();
      double pt = mctrack->pt();
      double eta = mctrack->pseudoRapidity();
      StMcTrack *parent = mctrack->parent();
      int parent_geant_id = -1;
      if(parent) parent_geant_id = parent->geantId();
      LOG_DEBUG  << "MC track " << i << " with id = " << geantId 
		 << ", parent = " << parent_geant_id << endm;
      bool pass = false;
      if(mEmbedParticle.Contains("Phi",TString::kIgnoreCase))
      	{
      	  if((geantId==11 || geantId==12) && parent_geant_id==10151) pass = true;
      	}
      else if(mEmbedParticle.Contains("Hadron",TString::kIgnoreCase))
	{
	  if( (geantId==8 || geantId==9 || geantId==11 || geantId==12)
	      && parent_geant_id==-1) pass = true;
	}
      else if(mEmbedParticle.Contains("Jpsi",TString::kIgnoreCase))
	{
	  if( (geantId==5 || geantId==6) && parent_geant_id==168) pass = true;
	}
      if(!pass) continue;
      if(fabs(eta)<0.5) mhHadronMcPt->Fill(pt,geantId,mCentrality);
      StMcVertex *stop_vtx  = mctrack->stopVertex();
      double stop_z = -1, stop_r = -1;
      int geantProcess = 0;
      if(stop_vtx) 
	{
	  StThreeVectorF pos = stop_vtx->position();
	  stop_z = pos.z();
	  stop_r = sqrt(pos.x()*pos.x()+pos.y()*pos.y());
	  geantProcess = stop_vtx->geantProcess();
	  bool ismuon = false;
	  int nDaug = stop_vtx->numberOfDaughters();
	  for(int iDaug=0; iDaug<nDaug; iDaug++)
	    {
	      StMcTrack *daug = stop_vtx->daughter(iDaug);
	      int daug_id = daug->geantId();
	      if(daug_id==5 || daug_id==6)
		{
		  ismuon = true;
		  break;
		}
	    }
	  if(!ismuon && geantProcess==5) geantProcess = 40;
	}
      mhHadronDecay->Fill(geantProcess, geantId);
      LOG_DEBUG << "Stop: z = " << stop_z << " r = " << stop_r << endm;

      if(0)
	{
	  if(stop_vtx)
	    {
	      cout << "Geant process: " << stop_vtx->geantProcess() << "  ->  " << geantProcess << endl;
	      cout << "Generator process: " << stop_vtx->generatorProcess() << endl;
	      cout << "Decay daughters: " << endl;
	      int nDaug = stop_vtx->numberOfDaughters();
	      for(int iDaug=0; iDaug<nDaug; iDaug++)
		{
		  StMcTrack *daug = stop_vtx->daughter(iDaug);
		  cout << "  " << iDaug << " is " <<  daug->particleDefinition()->name() << endl;
		}
	      cout << endl;
	    }
	}


      int rcIndex = mMcIndices[i];
      if(rcIndex<0) continue;
      StTrack *rcTrack = trackNodes[rcIndex]->track(mTrackType);
      StThreeVectorF momentum = rcTrack->geometry()->momentum();
      double rc_pt        = momentum.perp();
      double rc_eta       = momentum.pseudoRapidity();
      double rc_phi       = momentum.phi();
      if(!mMtdUtil->passTrack(rcTrack,mPriVtx)) continue;
      if(fabs(eta)<0.5) 
	{
	  mhHadronMcPtTpc->Fill(pt,geantId,mCentrality);
	  mhHadronRcPtTpc->Fill(rc_pt,geantId,mCentrality);
	}

      double nSigmaPi  = -999.;
      static StTpcDedxPidAlgorithm pidAlgorithm;
      static StPionPlus* Pion = StPionPlus::instance();
      const StParticleDefinition *pd = rcTrack->pidTraits(pidAlgorithm);
      if(pd && pidAlgorithm.traits())
	{
	  nSigmaPi  = pidAlgorithm.numberOfSigma(Pion);
	  if(fabs(eta)<0.5)
	    {
	      mhNsigmaPiVsPt->Fill(pt, nSigmaPi);
	    }
	}

      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(rcTrack);
      if(!mtdPid) continue;
      StMtdHit* hit = mtdPid->mtdHit();
      if(!hit) continue;
      if(!isMcMtdHit(hit)) continue;
      if(fabs(eta)<0.5) 
	{
	  mhHadronMcPtMtd->Fill(pt,geantId,mCentrality);
	  mhHadronRcPtMtd->Fill(rc_pt,geantId,mCentrality);
	}
      cout << "backleg = " << hit->backleg() << endl;

      double dy      = mtdPid->deltaY();
      double dz      = mtdPid->deltaZ();
      double tof_mc  = mtdPid->timeOfFlight();
      double tof_exp = mtdPid->expTimeOfFlight();
      double dtof    = tof_mc-tof_exp;
      cout << tof_mc << endl;

      double fill[] = {geantId*1., rc_pt, dz, dy, dtof, stop_r, geantProcess*1.};
      mhHadronMatch->Fill(fill);

      if(mSaveTree)
	{
	  mEmbedData.geantId[trackCount]  = geantId;
	  mEmbedData.mcpt[trackCount]     = pt;
	  mEmbedData.mcphi[trackCount]    = mctrack->momentum().phi();
	  mEmbedData.mceta[trackCount]    = eta;
	  mEmbedData.rcpt[trackCount]     = rc_pt;
	  mEmbedData.rcphi[trackCount]    = rc_phi;
	  mEmbedData.rceta[trackCount]    = rc_eta;
	  mEmbedData.nSigmaPi[trackCount] = nSigmaPi;
	  mEmbedData.dz[trackCount]       = dz;
	  mEmbedData.dy[trackCount]       = dy;
	  mEmbedData.dtof[trackCount]     = dtof;
	  trackCount++;
	}
      if(!isRcMuon(rcTrack)) continue;
      if(fabs(eta)<0.5) 
	{
	  mhHadronMcPtMuon->Fill(pt,geantId,mCentrality);
	  mhHadronRcPtMuon->Fill(rc_pt,geantId,mCentrality);
	}
    }
  if(mSaveTree) 
    {
      mEmbedData.centrality = mCentrality;
      mEmbedData.nTracks = trackCount;
    }
  return kStOK;
}

//_____________________________________________________________________________
void StMtdEmbedding::runProjection()
{
  StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();
  StSPtrVecMcTrack mctracks = mMcEvent->tracks();
  Int_t nMcTracks = mctracks.size();
  IntVec idVec;
  DoubleVec pathVec,  tofVec;
  PointVec crossVec;
  int bl, mod, cell;
  for(Int_t i=0; i<nMcTracks; i++) 
    {
      StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
      if(!mctrack) continue;
      int rcIndex = mMcIndices[i];
      if(rcIndex<0) continue;
      StTrack *rcTrack = trackNodes[rcIndex]->track(mTrackType);
      if(!mMtdUtil->passTrack(rcTrack,mPriVtx)) continue;
      double pt = rcTrack->geometry()->momentum().perp();
      double phi = mMtdUtil->rotatePhi(rcTrack->geometry()->momentum().phi());
      double eta = rcTrack->geometry()->momentum().pseudoRapidity();

      StPhysicalHelixD helix = rcTrack->outerGeometry()->helix();
      mMtdGeom->HelixCrossCellIds(helix,mVerPosition,idVec,pathVec,crossVec,tofVec);
      int nMatch = idVec.size();
      if(nMatch<=0) continue;
      mMtdGeom->DecodeCellId(idVec[0],bl,mod,cell);
      StMtdGeoModule *mMtdGeoModule = mMtdGeom->GetGeoModule(bl,mod);
      if(!mMtdGeoModule) continue;
      double local[3]={0,0,0};
      double global[3]={crossVec[0].x(),crossVec[0].y(),crossVec[0].z()};
      mMtdGeoModule->MasterToLocal(global,local);

      // MTD matching efficiency
      double isMatch = 0;
      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(rcTrack);
      if(mtdPid && isMcMtdHit(mtdPid->mtdHit()) && isMtdResp(pt, mtdPid->mtdHit()->backleg(), mtdPid->mtdHit()->module()))
	{
	  isMatch = 1;
	}
      double fill[] = {pt, bl*1.0, mVerPosition.z(), isMatch, phi, eta};
      mhMtdMatchEff->Fill(fill);

      // MTD response efficiency
      if(nMatch==1)
	{
	  bool isproj = false;
	  double dzsig = fabs(mMtdUtil->getPtDepDzCut(pt));
	  Float_t zcut = 3*dzsig < 15 ? 3*dzsig : 15;
	  Int_t cellshift = 0;
	  if(bl==8 ) cellshift = +3;
	  if(bl==24) cellshift = -2;
	  if(fabs(local[2])<= (stripLength/2.-zcut) && cell>=3+cellshift && cell<=8+cellshift)  isproj = true;
	  if(isproj)
	    {
	      int gMod = (bl-1)*5 + mod;
	      mhProjTrack->Fill(pt,gMod);
	      StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(rcTrack);
	      if(mtdPid && isMcMtdHit(mtdPid->mtdHit()) )
		{
		  mhMatchTrack->Fill(pt,gMod);
		}
	    }
	}
    }
}

//_____________________________________________________________________________
void StMtdEmbedding::fillInvMass(TLorentzVector muon1, TLorentzVector muon2, THnSparse *hn[2], double mc_pt)
{
  TLorentzVector jpsi = muon1 + muon2;
  Double_t pt1 = muon1.Pt();
  Double_t pt2 = muon2.Pt();
  if(pt1<pt2)
    {
      pt1 = pt2;
      pt2 = muon1.Pt();
    }

  if(jpsi.Pt()>0.15)
    {
      Double_t fill[] = {jpsi.M(),jpsi.Pt(),jpsi.Rapidity(),pt1,pt2, 1.*mCentrality, mZdcRate};
      hn[0]->Fill(fill, mVtxWeight);
    }
  // double weight = 1;
  // if(mhInputJpsiPt[mCentType]) 
  //   {
  //     weight = mhInputJpsiPt[mCentType]->GetBinContent(mhInputJpsiPt[mCentType]->FindFixBin(mc_pt));
  //   }
  //hn[1]->Fill(fill,weight);
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::findMatchedRcTrack(StMcTrack *track)
{
  pair<mcTrackMapIter,mcTrackMapIter> mcBounds = mMcTrackMap->equal_range(track);
  Int_t maxCommonHits = 0;
  const StGlobalTrack *rcCandTrack = 0;
  for(mcTrackMapIter mcMapIter = mcBounds.first; mcMapIter != mcBounds.second; mcMapIter ++)
    {
      StTrackPairInfo *pair = mcMapIter->second;
      const StGlobalTrack *rcTrack = pair->partnerTrack();
      Int_t commonHits = pair->commonTpcHits();
      if(commonHits > maxCommonHits)
	{
	  maxCommonHits = commonHits;
	  rcCandTrack = rcTrack;
	}
    }

  Int_t rcIndex = -1;
  if(maxCommonHits>=10)
    {
      StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();
      Int_t nNodes = trackNodes.size();
      for(Int_t i=0; i<nNodes; i++)
	{
	  StTrack *rcTrack = trackNodes[i]->track(mTrackType);
	  if(!rcTrack) continue;
	  if(rcTrack->key()==rcCandTrack->key())
	    {
	      rcIndex = i;
	      break;
	    }
	}
    }
  return rcIndex;
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMcMuon(StMcTrack *track)
{
  Int_t geantId = track->geantId();
  if(geantId==5 || geantId==6) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isRcMuon(StTrack *track)
{
  StMtdPidTraits *mtdPid = mMtdUtil->getMtdPidTraits(track);
  if(!mtdPid) return kFALSE; 
  if(!isMcMtdHit(mtdPid->mtdHit())) return kFALSE;

  return mMtdUtil->isMuonCandidate(track);
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMcMtdHit(StMtdHit *hit)
{
  return (hit && hit->idTruth()>0);
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isTpcResp(double pt, double eta, double phi)
{
  if(!mhTpcCorr) return true;
  if(eta>0.2) return true;
  if(pt<1) pt = 1;
  if(pt>=20) pt = 10;
  int binx = mhTpcCorr->GetXaxis()->FindFixBin(pt);
  int biny = mhTpcCorr->GetYaxis()->FindFixBin(phi);
  if(binx==0 || binx>mhTpcCorr->GetNbinsX()) return true;
  if(biny==0 || biny>mhTpcCorr->GetNbinsY()) return true;
  double eff = mhTpcCorr->GetBinContent(binx,biny);
  double pro = mRandom->Rndm();
  if(pro<eff) return true;
  else        return false;		       
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMtdResp(double pt, int bl, int mod)
{
  if(!mhMtdRespEmb[bl-1][mod-1] || !mhMtdRespCosmic[bl-1][mod-1]) return kTRUE;
  double pt_tmp = pt;
  if(pt>10) pt_tmp = 10;
  double eff = mhMtdRespCosmic[bl-1][mod-1]->Eval(pt_tmp)/mhMtdRespEmb[bl-1][mod-1]->Eval(pt_tmp);
  double pro = mRandom->Rndm();
  if(pro<eff) return kTRUE;
  else        return kFALSE;			       
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMtdTrig(double pt)
{
  if(!mhMtdTrigElecEff || !mhMtdTrigEff) return kTRUE;

  double pt_tmp = pt;
  if(pt>10) pt_tmp = 10;
  double eff = mhMtdTrigElecEff->Eval(pt_tmp) * mhMtdTrigEff->Eval(pt_tmp);
  double pro = mRandom->Rndm();
  if(pro<eff) return kTRUE;
  else        return kFALSE;			       
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::passPidCut(double pt)
{
  if(!mhMtdDtofEff) return kTRUE;

  double pt_tmp = pt;
  if(pt>10) pt_tmp = 10;
  double eff = mhMtdDtofEff->Eval(pt_tmp);
  double pro = mRandom->Rndm();
  if(pro<eff) return kTRUE;
  else        return kFALSE;			       
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::getTrackNdedx(StTrack *track)
{
  int nDedxHits = -1;
  StTpcDedxPidAlgorithm pidAlgorithm;
  const StParticleDefinition *pd = track->pidTraits(pidAlgorithm);
  if(pd && pidAlgorithm.traits()) nDedxHits = pidAlgorithm.traits()->numberOfPoints();
  return nDedxHits;
}

//_____________________________________________________________________________
double StMtdEmbedding::getSmearedTrkPt(double pt)
{
  double smear_pt = pt;
  if(mRunYear==2014)
    {
      smear_pt = pt * mRandom->Gaus(1,sqrt(pt)*mTpcTrkSmear);
    }
  //cout << pt << " -> " << smear_pt << " with " << sqrt(pt)*mTpcTrkSmear*pt << endl;

  return smear_pt;
}

//_____________________________________________________________________________
void StMtdEmbedding::bookHistos()
{
  // event histograms
  const char *tracks[2] = {"global","primary"};
  const char *qaTrackType[4] = {"MCreco","DataPion","DataMuonUL","DataMuonLS"};
  const char *trkEffType[8] = {"MC","Tpc","MtdMth","Fake","MuonPid","MtdTrig","TrigUnit","Embed"};
  const char *weight_name[2] = {"","_w"};
  const char *trgSetupName[3] = {"prod_low","prod_mid","prod_high"};
  const char *trgSetupName2[4] = {"prod","prod_low","prod_mid","prod_high"};
  const char *tpcStatus[2] = {"good","bad"};

  mhAnalysisCuts = new TH1F("hAnalysisCuts","Cuts used for analysis",20,0,20);
  AddHist(mhAnalysisCuts);
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 1,  "|vtx_z|",            mMtdUtil->getMaxVtxZ());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 2,  "trk_pt_min",         mMtdUtil->getMinTrkPt());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 3,  "trk_pt_max",         mMtdUtil->getMaxTrkPt());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 4,  "trk_eta",            mMtdUtil->getMaxTrkEta());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 5,  "MinNHitsFit",        mMtdUtil->getMinNHitsFit());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 6,  "MinNHitsDedx",       mMtdUtil->getMinNHitsDedx());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 7,  "mMaxDca",            mMtdUtil->getMaxDca());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 8,  "mMinNsigmaPi",       mMtdUtil->getMinNsigmaPi());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 9,  "mMaxNsigmaPi",       mMtdUtil->getMaxNsigmaPi());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 10, "mMinFitHitsFaction", mMtdUtil->getMinFitHitsFaction());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 11, "mMaxDiffz",          mMtdUtil->getMaxDiffz());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 12, "mMaxTrkPairDcaDr",   mMtdUtil->getMaxTrkPairDcaDr());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 13, "mMaxTrkPairDcaDz",   mMtdUtil->getMaxTrkPairDcaDz());
  mMtdUtil->addCutToHisto(mhAnalysisCuts, 14, "mMuonDeltaZ",        mMtdUtil->getMuonDeltaZ());

  mhEventStat = new TH1F("hEventStat","Event statistics",10,0,10);
  mhEventStat->GetXaxis()->SetBinLabel(1,"All events");
  mhEventStat->GetXaxis()->SetBinLabel(2,"Good trigger");
  for(Int_t i=0; i<3; i++)
    mhEventStat->GetXaxis()->SetBinLabel(3+i,Form("%s_raw",trigName[i]));
  mhEventStat->GetXaxis()->SetBinLabel(3+3,"Good vertex");
  for(Int_t i=0; i<3; i++)
    mhEventStat->GetXaxis()->SetBinLabel(4+3+i,Form("%s_acc",trigName[i]));
  AddHist(mhEventStat);

  const int nTrkPtBin = 200;
  const double lowTrkPtBin = 0, hiTrkPtBin = 20.;

  const int dimZdc = 3;
  const int nBinsZdc[dimZdc] = {120, 20, 5};
  const double lowBinZdc[dimZdc] = {0, 0, 0};
  const double upBinZdc[dimZdc] = {120, 20, 5};

  const Int_t dimJpsi = 7;
  const Int_t nBinsJpsi[dimJpsi] =     {100,  80,  20, 20, 100, 20, 12};
  const Double_t lowBinJpsi[dimJpsi] = {2,    0,   -1, 0,  0,   0,  0};
  const Double_t upBinJpsi[dimJpsi]  = {4,    20,  1,  10, 10,  20, 120};

  const int dimJpsiMatch = 5;
  const int nBinsJpsiMatch[dimJpsiMatch] =     {200, 80,  100, 20, 12};
  const double lowBinJpsiMatch[dimJpsiMatch] = {0,   0,   0,   0,  0};
  const double upBinJpsiMatch[dimJpsiMatch] =  {20,  20,  10,  20, 120};

  const int dimTrkInfo = 5;
  const int nBinsTrkInfo[dimTrkInfo] = {nTrkPtBin, 40, 360, 3, 16};
  const double lowBinTrkInfo[dimTrkInfo] = {lowTrkPtBin, -2, 0, -1, 0};
  const double upBinTrkInfo[dimTrkInfo] = {hiTrkPtBin, 2, 2*pi, 2, 16};

  const int dimTrkEff = 5;
  const int nBinsTrkEff[dimTrkEff]     = {nTrkPtBin,  40,    16, 20,  4};
  const double lowBinTrkEff[dimTrkEff] = {lowTrkPtBin, -2, 0,  0,   0};
  const double upBinTrkEff[dimTrkEff]  = {hiTrkPtBin, 2,   16, 200, 4};

  const int dimTrkNsigmaPi = 6;
  const int nBinsTrkNsigmaPi[dimTrkNsigmaPi]     = {10, 60,  4,  20,   16, 4};
  const double lowBinTrkNsigmaPi[dimTrkNsigmaPi] = {0,  -3,  -1, -100, 0,  0};
  const double upBinTrkNsigmaPi[dimTrkNsigmaPi]  = {10, 3,   1,  100,  16, 4};

  const Int_t    dimQaNHitsPoss = 4;
  const Int_t    nBinsQaNHitsPoss[dimQaNHitsPoss]  = {10, 50, 40,  5};
  const Double_t lowBinQaNHitsPoss[dimQaNHitsPoss] = {0,  0,  -100, 0};
  const Double_t upBinQaNHitsPoss[dimQaNHitsPoss]  = {10, 50, 100,  5};

  const Int_t dimQaTrk = 3;
  const Int_t nBinsQaTrk[dimQaTrk] = {nTrkPtBin , 20, 360};
  const Double_t lowBinQaTrk[dimQaTrk] = {lowTrkPtBin, -1, 0};
  const Double_t upBinQaTrk[dimQaTrk]  = {hiTrkPtBin, 1, 2*pi};

  const int dimTrkRes = 4;
  const int nBinsTrkRes[dimTrkRes] =     {200,   nTrkPtBin,   nTrkPtBin,   16};
  const double lowBinTrkRes[dimTrkRes] = {-0.95, lowTrkPtBin, lowTrkPtBin, 0};
  const double upBinTrkRes[dimTrkRes] =  {1.05,  hiTrkPtBin,   hiTrkPtBin, 16};

  const Int_t dimHadron = 7;
  const Int_t nBinsHadron[dimHadron] = {20, 200, 400, 200, 150, 450, 41};
  const Double_t lowBinHadron[dimHadron] = {0, 0, -199.5, -99.5, -5, 0, 0};
  const Double_t upBinHadron[dimHadron]  = {20, 20, 200.5, 100.5, 10, 450, 41};

  const int dimMthEff = 6;
  const int nBinsMthEff[dimMthEff]     = {150, 30, 100,    2, 90, 10};
  const double lowBinMthEff[dimMthEff] = {0,   1,   -100,  0, 0, -1.0};
  const double upBinMthEff[dimMthEff]  = {15,  31,  100,   2, 2*TMath::Pi(), 1.0};

  if(mRunEmbedHadron)
    {
      mhHadronDecay = new TH2F("mhHadronDecay","MC particle decay mode;geant process;geant id",41,0,41,20,0,20);
      AddHist(mhHadronDecay); 

      mhHadronMatch = new THnSparseF(Form("mhHadronMatch"),Form("geant id vs p_{T,rec} vs #Deltaz vs #Deltay vs #Deltatof vs decay radius vs geant process;geant id;p_{T,rec} (GeV/c);#Deltaz (cm);#Deltay (cm);#Deltatof (ns);r (cm);process"),dimHadron,nBinsHadron,lowBinHadron,upBinHadron);
      AddHist((TH1*)mhHadronMatch); 

      mhHadronMcPt = new TH3F("mhHadronMcPt","p_{T} of MC tracks (|#eta|<0.5);p_{T} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronMcPt); 

      mhHadronMcPtTpc = new TH3F("mhHadronMcPtTpc","p_{T} of MC tracks recontructed in TPC (|#eta|<0.5);p_{T} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronMcPtTpc); 

      mhHadronMcPtMtd = new TH3F("mhHadronMcPtMtd","p_{T} of MC tracks matched to MTD (|#eta|<0.5);p_{T} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronMcPtMtd); 

      mhHadronMcPtMuon = new TH3F("mhHadronMcPtMuon","p_{T} of MC tracks with muon PID (|#eta|<0.5);p_{T} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronMcPtMuon); 

      mhHadronRcPtTpc = new TH3F("mhHadronRcPtTpc","p_{T} of MC tracks recontructed in TPC (|#eta|<0.5);p_{T,rc} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronRcPtTpc); 

      mhHadronRcPtMtd = new TH3F("mhHadronRcPtMtd","p_{T} of MC tracks matched to MTD (|#eta|<0.5);p_{T,mc} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronRcPtMtd); 

      mhHadronRcPtMuon = new TH3F("mhHadronRcPtMuon","p_{T} of MC tracks with muon PID (|#eta|<0.5);p_{T,rc} (GeV/c);geantId;centrality",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,20,0,20,16,0,16);
      AddHist(mhHadronRcPtMuon); 

      mhNsigmaPiVsPt = new TH2F("mhNsigmaPiVsPt","nSigmaPi vs pt distribution;p_{T} (GeV/c);n#sigma_{#pi}",nTrkPtBin,lowTrkPtBin,hiTrkPtBin,200,-5,5);
      AddHist(mhNsigmaPiVsPt); 
    }

  if(mRunProjection)
    {
      mhProjTrack = new TH2F("mhProjTrack","Backleg vs pt of tracks projected to MTD;p_{T} (GeV/c);BL",100,0,10,150,0.5,150.5);
      AddHist(mhProjTrack); 

      mhMatchTrack = new TH2F("mhMatchTrack","Backleg vs pt of tracks matched to MTD;p_{T} (GeV/c);BL",100,0,10,150,0.5,150.5);
      AddHist(mhMatchTrack); 

      mhMtdMatchEff = new THnSparseF("mhMtdMatchEff","p_{T} vs BL vs v_{z} vs IsMatch vs #varphi vs #eta;p_{T} (GeV/c);BL;v_{z} (cm);IsMatch;#varphi;#eta",dimMthEff,nBinsMthEff,lowBinMthEff,upBinMthEff);
      AddHist((TH1*)mhMtdMatchEff); 
    }

  for(Int_t i=0; i<kNtrig; i++)
    {
      if(!mMtdUtil->isTriggerOn(i)) continue;

      if(mRunRealDataQa)
	{
	  for(int j=0; j<3; j++)
	    {
	      mhQaTrkEtaPhi[i][j] = new THnSparseF(Form("mhQaTrkEtaPhi_%s_%s",trigName[i],trgSetupName[j]),Form("%s: p_{T} vs #eta vs #varphi;p_{T} (GeV/c);#eta;#varphi",trigName[i]),dimQaTrk,nBinsQaTrk,lowBinQaTrk,upBinQaTrk);
	      AddHist((TH1*)mhQaTrkEtaPhi[i][j]); 

	      mhQaTrkEtaPhiWithCuts[i][j] = new THnSparseF(Form("mhQaTrkEtaPhiWithCuts_%s_%s",trigName[i],trgSetupName[j]),Form("%s: p_{T} vs #eta vs #varphi;p_{T} (GeV/c);#eta;#varphi",trigName[i]),dimQaTrk,nBinsQaTrk,lowBinQaTrk,upBinQaTrk);
	      AddHist((TH1*)mhQaTrkEtaPhiWithCuts[i][j]); 

	      for(int k=0; k<2; k++)
		{
		  mhQaTrkNHitsFit[i][j][k] = new TH2F(Form("mhQaTrkNHitsFit_%s_%s_%s",trigName[i],trgSetupName[j],tpcStatus[k]),Form("%s: NHitsFit vs p_{T} of tracks;p_{T} (GeV/c);NHitsFit",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,50,0,50);
		  AddHist(mhQaTrkNHitsFit[i][j][k]); 

		  mhQaTrkNHitsDedx[i][j][k] = new TH2F(Form("mhQaTrkNHitsDedx_%s_%s_%s",trigName[i],trgSetupName[j],tpcStatus[k]),Form("%s: NHitsDedx vs p_{T} of tracks;p_{T} (GeV/c);NHitsDedx",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,50,0,50);
		  AddHist(mhQaTrkNHitsDedx[i][j][k]); 

		  mhQaTrkDca[i][j][k] = new TH2F(Form("mhQaTrkDca_%s_%s_%s",trigName[i],trgSetupName[j],tpcStatus[k]),Form("%s: Dca vs p_{T} of tracks;p_{T} (GeV/c);DCA (cm)",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,30,0,3);
		  AddHist(mhQaTrkDca[i][j][k]); 

		  mhQaTrkNSigmaPi[i][j][k] = new TH2F(Form("mhQaTrkNSigmaPi_%s_%s_%s",trigName[i],trgSetupName[j],tpcStatus[k]),Form("%s: n#sigma_{#pi} vs p_{T} of tracks;p_{T} (GeV/c);n#sigma_{#pi}",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,100,-5,5);
		  AddHist(mhQaTrkNSigmaPi[i][j][k]); 
		}
	    }
	}

      mhZdcRateVsTrigSetup[i] = new THnSparseF(Form("mhZdcRateVsTrigSetup_%s",trigName[i]),Form("%s: ZdcRate vs. centrality vs. trigSetup;ZdcRate (kHz);centrality;trigSetup",trigName[i]),dimZdc,nBinsZdc,lowBinZdc,upBinZdc);

      mhTofMthTrksVsZdcRate[i] = new TH2F(Form("mhTofMthTrksVsZdcRate_%s",trigName[i]),Form("%s: mTofMthTrks vs zdc rate;zdc rate (kHz);mTofMthTrks",trigName[i]),1000,0,1000,500,0,500);

      mhTofMthTrksVzBbcRate[i] = new TH2F(Form("mhTofMthTrksVzBbcRate_%s",trigName[i]),Form("%s: mTofMthTrks vs bbc rate;bbc rate (kHz);mTofMthTrks",trigName[i]),1000,0,10000,500,0,500);

      if(i==0)
	{
	  for(int j=0; j<4; j++)
	    {
	      mhTofMultVsgRefMult[j] = new TH2F(Form("mhTofMultVsgRefMult_%s_%s",trigName[i],trgSetupName2[j]),Form("%s: mTofMult vs gRefMult;gRefMult;mTofMult",trigName[i]),100,0,1000,100,0,5000);
	      
	      mhTofMultVsgRefMultCorr[j] = new TH2F(Form("mhTofMultVsgRefMultCorr_%s_%s",trigName[i],trgSetupName2[j]),Form("%s: mTofMult vs gRefMultCorr;gRefMultCorr;mTofMult",trigName[i]),100,0,1000,100,0,5000);

	      mhNgTrkVsCent[j] = new TH2F(Form("mhNgTrkVsCent_%s_%s",trigName[i],trgSetupName2[j]),Form("%s: nGlobalTrack vs cent;cent;ngTrk",trigName[i]),20,0,20,300,0,30000);

	      mhZdcRateVsCent[j] = new TH2F(Form("mhZdcRateVsCent_%s_%s",trigName[i],trgSetupName2[j]),Form("%s: ZDCrate vs cent;cent;ZDCx (kHz)",trigName[i]),20,0,20,120,0,120);
	    }
	}

      mhgRefMult[i]  = new TH1F(Form("mhgRefMult_%s",trigName[i]),Form("%s: gRefMult distribution;gRefMult",trigName[i]),1000,0,1000); 

      mhgRefMultCorr[i]  = new TH1F(Form("mhgRefMultCorr_%s",trigName[i]),Form("%s: corrected gRefMult distribution;gRefMultCorr",trigName[i]),1000,0,1000);

      mhCentrality[i]  = new TH1F(Form("mhCentrality_%s",trigName[i]),Form("%s: centrality distribution;centrality",trigName[i]),20,0,20);

      mhDataVtxZ[i]  = new TH1F(Form("mhDataVtxZ_%s",trigName[i]),Form("%s: z distribution of data vertex;vz (cm)",trigName[i]),200,-200,200);

      mhMcVtxZVsDataVtxZ[i]  = new TH2F(Form("mhMcVtxZVsDataVtxZ_%s",trigName[i]),Form("%s: z distribution of MC vertex vs. Data;vz_{data} (cm);vz_{mc} (cm)",trigName[i]),200,-200,200,200,-200,200);

      mhSetupVsMcVtxZ[i]  = new TH2F(Form("mhSetupVsMcVtxZ_%s",trigName[i]),Form("%s: z distribution of MC vertex;vz (cm);TrgSetup",trigName[i]),200,-200,200,4,0,4);

      mhNEmbedJpsi[i] = new TH1F(Form("hNEmbedJpsi_%s",trigName[i]),Form("%s: # of embedded j/psi per event;N",trigName[i]),50,0,50);

      mhgTrkPtRes[i] = new THnSparseF(Form("mhgTrkPtRes_%s",trigName[i]),Form("%s: (p_{T,true}-p_{T,reco})/p_{T,true} vs p_{T,reco} vs p_{T,true} vs centrality of global muon tracks;(p_{T,true}-p_{T,reco})/p_{T,true};p_{T,reco} (GeV/c);p_{T,true} (GeV/c);centrality",trigName[i]),dimTrkRes,nBinsTrkRes,lowBinTrkRes,upBinTrkRes);

      mhpTrkPtRes[i] = new THnSparseF(Form("mhpTrkPtRes_%s",trigName[i]),Form("%s: (p_{T,true}-p_{T,reco})/p_{T,true} vs p_{T,reco} vs p_{T,true} vs centrality of primary muon tracks;(p_{T,true}-p_{T,reco})/p_{T,true};p_{T,reco} (GeV/c);p_{T,true} (GeV/c);centrality",trigName[i]),dimTrkRes,nBinsTrkRes,lowBinTrkRes,upBinTrkRes);

      mhMcTrkInfo[i] = new THnSparseF(Form("mhMcTrkInfo_%s",trigName[i]),Form("%s: p_{T} vs #eta vs #varphi vs charge vs centrality of MC tracks;p_{T}^{mc} (GeV/c);#eta_{mc};#varphi_{mc};centrality",trigName[i]),dimTrkInfo,nBinsTrkInfo,lowBinTrkInfo,upBinTrkInfo);

      mhMcTrkInfoTpc[i] = new THnSparseF(Form("mhMcTrkInfoTpc_%s",trigName[i]),Form("%s: p_{T} vs #eta vs #varphi vs charge vs centrality of MC tracks matched to TPC;p_{T}^{mc} (GeV/c);#eta_{mc};#varphi_{mc};charge;centrality",trigName[i]),dimTrkInfo,nBinsTrkInfo,lowBinTrkInfo,upBinTrkInfo);

      mhMcTrkInfoMtd[i] = new THnSparseF(Form("mhMcTrkInfoMtd_%s",trigName[i]),Form("%s: p_{T} vs #eta vs #varphi vs charge vs centrality of MC tracks matched to MTD;p_{T}^{mc} (GeV/c);#eta_{mc};#varphi_{mc};charge;centrality",trigName[i]),dimTrkInfo,nBinsTrkInfo,lowBinTrkInfo,upBinTrkInfo);

      mhMcTrkInfoFinal[i] = new THnSparseF(Form("mhMcTrkInfoFinal_%s",trigName[i]),Form("%s: p_{T} vs #eta vs #varphi vs charge vs centrality of MC tracks after all cuts;p_{T}^{mc} (GeV/c);#eta_{mc};#varphi_{mc};charge;centrality",trigName[i]),dimTrkInfo,nBinsTrkInfo,lowBinTrkInfo,upBinTrkInfo);

      if(!mRunSystematics)
	{
	  AddHist((TH1*)mhZdcRateVsTrigSetup[i]);
	  AddHist(mhTofMthTrksVsZdcRate[i]); 
	  AddHist(mhTofMthTrksVzBbcRate[i]); 
	  if(i==0)
	    {
	      for(int j=0; j<4; j++)
		{
		  AddHist(mhTofMultVsgRefMult[j]);
		  AddHist(mhTofMultVsgRefMultCorr[j]);
		  AddHist(mhNgTrkVsCent[j]);
		  AddHist(mhZdcRateVsCent[j]);
		}
	    }
	  AddHist(mhgRefMult[i]);
	  AddHist(mhgRefMultCorr[i]); 
	  AddHist(mhCentrality[i]);
	  AddHist(mhDataVtxZ[i]);
	  AddHist(mhMcVtxZVsDataVtxZ[i]);
	  AddHist(mhSetupVsMcVtxZ[i]);
	  AddHist(mhNEmbedJpsi[i]); 
	  AddHist((TH1*)mhgTrkPtRes[i]); 
	  AddHist((TH1*)mhpTrkPtRes[i]); 
	  AddHist((TH1*)mhMcTrkInfo[i]); 
	  AddHist((TH1*)mhMcTrkInfoTpc[i]); 
	  AddHist((TH1*)mhMcTrkInfoMtd[i]);
	  AddHist((TH1*)mhMcTrkInfoFinal[i]);
	}

      for(int j=0; j<6; j++)
	{
	  mhMcTrkPtEff[i][j] = new THnSparseF(Form("mhMcTrkPtEff_%s_%s",trkEffType[j],trigName[i]),Form("%s: p_{T} vs centrality vs ZdcRate vs TrgSetup;p_{T}^{mc} (GeV/c);centrality;ZdcRate (kHz);Setup",trigName[i]),dimTrkEff,nBinsTrkEff,lowBinTrkEff,upBinTrkEff);
	  if(!mRunSystematics) AddHist((TH1*)mhMcTrkPtEff[i][j]); 
	  else
	    {
	      if(j==0 || j==1 || j==4) AddHist((TH1*)mhMcTrkPtEff[i][j]); 
	    }
	}

      mhRcTrkNsigmaPi[i] = new THnSparseF(Form("mhRcTrkNsigmaPi_%s",trigName[i]),Form("%s: p_{T} vs n#sigma_{#pi} vs #eta vs vz vs centrality vs TrgSetup;p_{T}^{mc} (GeV/c);n#sigma_{#pi};#eta;vz (cm);centrality;Setup",trigName[i]),dimTrkNsigmaPi,nBinsTrkNsigmaPi,lowBinTrkNsigmaPi,upBinTrkNsigmaPi);

      mhQaNHitsPoss[i] = new THnSparseF(Form("mhQaNHitsPoss_%s",trigName[i]),"Track p_{T} vs. NHitsPoss vs. vz vs. luminosity;p_{T} (GeV/c);NHitsPoss;vz (cm);TrigSetup",dimQaNHitsPoss, nBinsQaNHitsPoss, lowBinQaNHitsPoss, upBinQaNHitsPoss);

      mhTrkPtBlModMtd[i] = new TH3F(Form("mhTrkPtBlModMtd_%s",trigName[i]),Form("%s: reco track p_{T} vs backleg vs module;p_{T} (GeV/c);Backleg;Module",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,30,0.5,30.5,5,0.5,5.5);

      mhTrkPtBlModTrig[i] = new TH3F(Form("mhTrkPtBlModTrig_%s",trigName[i]),Form("%s: reco track p_{T} vs backleg vs module;p_{T} (GeV/c);Backleg;Module",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,30,0.5,30.5,5,0.5,5.5);

      mhRealHitMap[i] = new TH2F(Form("mhRealHitMap_%s",trigName[i]),Form("%s: real hit map;Backleg;Channel",trigName[i]),30,0.5,30.5,60,0,60);

      mhMcHitMap[i] = new TH2F(Form("mhMcHitMap_%s",trigName[i]),Form("%s: MC hit map;Backleg;Channel",trigName[i]),30,0.5,30.5,60,0,60);

      mhMcHitMapWithEff[i] = new TH2F(Form("mhMcHitMapWithEff_%s",trigName[i]),Form("%s: MC hit map with efficiency;Backleg;Channel",trigName[i]),30,0.5,30.5,60,0,60);

      if(!mRunSystematics)
	{
	  AddHist((TH1*)mhRcTrkNsigmaPi[i]); 
	  AddHist((TH1*)mhQaNHitsPoss[i]);
	  AddHist(mhTrkPtBlModMtd[i]); 
	  AddHist(mhTrkPtBlModTrig[i]); 
	  AddHist(mhRealHitMap[i]); 
	  AddHist(mhMcHitMap[i]); 
	  AddHist(mhMcHitMapWithEff[i]); 
	}
 
      for(int j=0; j<4; j++)
	{
	  mhTrkDca[i][j] = new TH2F(Form("hTrkDca_%s_%s",qaTrackType[j],trigName[i]),Form("%s: dca vs p_{T} of tracks;p_{T} (GeV/c);dca (cm)",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,100,0,10);

	  mhTrkNHitsFit[i][j] = new TH2F(Form("hTrkNHitsFit_%s_%s",qaTrackType[j],trigName[i]),Form("%s: NHitsFit vs p_{T} of tracks;p_{T} (GeV/c);NHitsFit",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,50,0,50);

	  mhTrkNHitsDedx[i][j] = new TH2F(Form("hTrkNHitsDedx_%s_%s",qaTrackType[j],trigName[i]),Form("%s: NHitsDedx vs p_{T} of tracks;p_{T} (GeV/c);NHitsDedx",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,50,0,50);

	  mhTrkNHitsFrac[i][j] = new TH2F(Form("hTrkNHitsFrac_%s_%s",qaTrackType[j],trigName[i]),Form("%s: NHitsFrac vs p_{T} of tracks;p_{T} (GeV/c);NHitsFrac",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,100,0,1);

	  mhTrkNHitsPoss[i][j] = new TH3F(Form("hTrkNHitsPoss_%s_%s",qaTrackType[j],trigName[i]),Form("%s: NHitsPoss vs p_{T} of tracks;p_{T} (GeV/c);NHitsPoss;TrigSetup",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,50,0,50,5,0,5);

	  mhTrkEtaPhi[i][j] = new THnSparseF(Form("hTrkEtaPhi_%s_%s",qaTrackType[j],trigName[i]),Form("%s: p_{T} vs #eta vs #varphi;p_{T} (GeV/c);#eta;#varphi",trigName[i]),dimQaTrk,nBinsQaTrk,lowBinQaTrk,upBinQaTrk);

	  mhTrkNSigmaPi[i][j] = new TH2F(Form("hTrkNSigmaPi_%s_%s",qaTrackType[j],trigName[i]),Form("%s: n#sigma_{#pi} vs p_{T} of tracks;p_{T} (GeV/c);n#sigma_{#pi}",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,100,-10,10);

	  mhTrkDyVsPt[i][j] = new TH2F(Form("hTrkDyVsPt_%s_%s",qaTrackType[j],trigName[i]),Form("%s: #Deltay vs p_{T} of tracks;p_{T} (GeV/c);#Deltay (cm)",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,200,-100,100);

	  mhTrkDzVsPt[i][j] = new TH2F(Form("hTrkDzVsPt_%s_%s",qaTrackType[j],trigName[i]),Form("%s: #Deltaz vs p_{T} of tracks;p_{T} (GeV/c);#Deltaz (cm)",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,400,-200,200);

	  mhTrkDtofVsPt[i][j] = new TH2F(Form("hTrkDtofVsPt_%s_%s",qaTrackType[j],trigName[i]),Form("%s: #Deltatof vs p_{T} of tracks;p_{T} (GeV/c);#Deltatof (ns)",trigName[i]),nTrkPtBin,lowTrkPtBin,hiTrkPtBin,400,-5,5);

	  mhTrkDtofVsMod[i][j] = new TH2F(Form("hTrkDtofVsMod_%s_%s",qaTrackType[j],trigName[i]),Form("%s: #Deltatof vs module of tracks;mod;#Deltatof (ns)",trigName[i]),150,0.5,150.5,400,-5,5);

	  if(!mRunSystematics)
	    {
	      AddHist(mhTrkDca[i][j]);
	      AddHist(mhTrkNHitsFit[i][j]);
	      AddHist(mhTrkNHitsDedx[i][j]);
	      AddHist(mhTrkNHitsFrac[i][j]);
	      AddHist(mhTrkNHitsPoss[i][j]);
	      AddHist((TH1*)mhTrkEtaPhi[i][j]); 
	      AddHist(mhTrkNSigmaPi[i][j]);
	      AddHist(mhTrkDyVsPt[i][j]);
	      AddHist(mhTrkDzVsPt[i][j]);
	      AddHist(mhTrkDtofVsPt[i][j]);
	      AddHist(mhTrkDtofVsMod[i][j]);
	    }
	}

      mhJpsiMatch[i] = new THnSparseF(Form("mhJpsiMatch_%s",trigName[i]),Form("%s: p_{T,true} vs p_{T,reco} vs p_{T2,rc} vs centrality vs ZdcRate (Unlike-sign,%s);p_{T}^{mc} (GeV/c);p_{T}^{rec} (GeV/c);p_{T2,rc} (GeV/c);cent;ZdcRate (kHz)",trigName[i],tracks[mTrackType]),dimJpsiMatch,nBinsJpsiMatch,lowBinJpsiMatch,upBinJpsiMatch);
      mhJpsiMatch[i]->Sumw2();
      if(!mRunSystematics)
	{
	  AddHist((TH1*)mhJpsiMatch[i]); 
	}

      for(int j=0; j<8; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      mhJpsiInfo[i][j][k] = new THnSparseF(Form("hJpsiInfo_%s_%s%s",trkEffType[j],trigName[i],weight_name[k]),Form("%s: invariant mass vs p_{T} vs rapidity vs pt1 vs pt2 vs centrality vs ZDCrate (US,%s);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);y;p_{T,1} (GeV/c);p_{T,2} (GeV/c);cent;ZDCrate (kHz)",trigName[i],tracks[mTrackType]),dimJpsi,nBinsJpsi,lowBinJpsi,upBinJpsi);
	      mhJpsiInfo[i][j][k]->Sumw2();
	      if(!mRunSystematics)
		{
		  AddHist((TH1*)mhJpsiInfo[i][j][k]); 
		}
	    }
	}
    }


  if(mSaveTree)
    {
      if(mOutTreeFileName.Length()==0) mOutTreeFileName = "jpsi.embed.tree.root";
      mOutTreeFile = new TFile(mOutTreeFileName.Data(),"recreate");
      LOG_INFO << "Create the output to store the hadron embed tree: " << mOutTreeFileName.Data() << endm;
      
      mOutTree = new TTree("EmbedTree","Embedding tree");
      mOutTree->SetAutoSave(100000); // 100 MB

      mOutTree->Branch("centrality",&mEmbedData.centrality,"centrality/I");
      mOutTree->Branch("nTracks",   &mEmbedData.nTracks,   "nTracks/I");
      mOutTree->Branch("geantId",   &mEmbedData.geantId,   "geantId[nTracks]/I");
      mOutTree->Branch("mcpt",      &mEmbedData.mcpt,      "mcpt[nTracks]/D");
      mOutTree->Branch("mcphi",     &mEmbedData.mcphi,     "mcphi[nTracks]/D");
      mOutTree->Branch("mceta",     &mEmbedData.mceta,     "mceta[nTracks]/D");
      mOutTree->Branch("rcpt",      &mEmbedData.rcpt,      "rcpt[nTracks]/D");
      mOutTree->Branch("rcphi",     &mEmbedData.rcphi,     "rcphi[nTracks]/D");
      mOutTree->Branch("rceta",     &mEmbedData.rceta,     "rceta[nTracks]/D");
      mOutTree->Branch("nSigmaPi",  &mEmbedData.nSigmaPi,  "nSigmaPi[nTracks]/D");
      mOutTree->Branch("dz",        &mEmbedData.dz,        "dz[nTracks]/D");
      mOutTree->Branch("dy",        &mEmbedData.dy,        "dy[nTracks]/D");
      mOutTree->Branch("dtof",      &mEmbedData.dtof,      "dtof[nTracks]/D");
    }
}

//_____________________________________________________________________________
void StMtdEmbedding::printConfig()
{
  const char *decision[2] = {"no","yes"};
  const char *trigName[3] = {"di-muon","single-muon","e-mu"};
  printf("\n=== Configuration for StMtdEmbedding ===\n");
  for(Int_t i=0; i<3; i++)
    {
      if(mMtdUtil->isTriggerOn(i)) printf("Enable %s trigger\n",trigName[i]);
    }
  printf("Use the primary vtx closest to VPD: %s\n",decision[mMtdUtil->isUsePrimVtxClosestToVpd()]);
  printf("Use the primary vtx with MTD hits: %s\n",decision[mMtdUtil->isRequireMtdHitForPrimVtx()]);
  printf("Use VPD vertex: %s\n",decision[mUseVpdVtx]);
  printf("TPC track smear: %4.2f%%\n",mTpcTrkSmear*100);
  printf("=======================================\n\n");
}
