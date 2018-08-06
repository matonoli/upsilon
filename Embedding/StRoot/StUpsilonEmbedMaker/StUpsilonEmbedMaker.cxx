// Modified by Oliver Matonoha
// CTU
// 1 February 2018

// Rosi Reed <rjreed@ucdavis.edu>
// UC Davis
// 1 July 2009
//
// Editted from:
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 24 June 2008
//
// StUpsilonEmbedMaker runs over the Upsilon embedding and
// calculates acceptance, L0 and L2 trigger efficiencies
// for the Upsilon.
//
#include <iostream>
#include <fstream>

#include "StUpsilonTriggerMaker/StUpsilonTriggerMaker.h"
#include "StUpsilonTriggerMaker/Quarkonium.hh"
#include "StDaqLib/EMC/StEmcDecoder.h"

// ROOT includes 
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"

// STAR includes
#include "StEventTypes.h"
#include "StMcEvent/StMcEventTypes.hh"
#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "StEmcTriggerMaker/StEmcTriggerMaker.h"
#include "tables/St_g2t_event_Table.h"
#include "tables/St_g2t_pythia_Table.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StContainers.h"
#include "StEmcUtil/database/StBemcTables.h"


// Local includes

#include "StUpsilonEmbedMaker.h"

ClassImp(StUpsilonEmbedMaker);

float scaleFactor(double Eta, int hitType=0) {
  // This function is used for determining the acceptance
  //     hitType = 0 - towers
  //     hitType = 1 - pre shower
  //     hitType = 2 - shower max Eta
  //     hitType = 3 - shower max Phi
  //function calculates the proper scale factor given a hit type and an
  //eta so that energy can be calculated
    float P0[]={14.69,559.7,0.1185e6,0.1260e6};
    float P1[]={-0.1022,-109.9,-0.3292e5,-0.1395e5};
    float P2[]={0.7484,-97.81,0.3113e5,0.1971e5};
    
    float x=fabs(Eta);
    return P0[hitType]+P1[hitType]*x+P2[hitType]*x*x;
}

int StUpsilonEmbedMaker::SetHighTowerVar(StMcTrack* mcTrack, bool isele){
  //This function finds the soft tower id of the highest towers assocated with
  //a particular mcTrack.  It also sets class variables eleEsum and posEsum, 
  //which are the sum of the energies left in the BEMC by the electron and 
  //positron.

  //Assumes that the mcTracks are the upsilon daughters

  StMcCalorimeterHit* MaxCalHit = 0;
  double maxEnergy = 0;
  double sumE = 0;
  int MaxSoftId = 0;
 
  for (size_t i=0; i<mcTrack->bemcHits().size(); ++i) {
    //loops over all the bemc hits and sums the energy and determines
    //the maximum hit
    StMcCalorimeterHit* calHit = mcTrack->bemcHits()[i];
    int modN = calHit->module();
    int etaN = calHit->eta();
    int subN = calHit->sub();
    float dE = calHit->dE();
    Float_t eta = 0;
    Float_t phi = 0;	
    StEmcGeom *geomBemc = StEmcGeom::getEmcGeom(1);
    geomBemc->getEta(modN,etaN,eta);
    geomBemc->getPhi(modN,subN,phi);
    sumE = sumE + dE*scaleFactor(eta);
    if (dE*scaleFactor(eta) > maxEnergy){
      maxEnergy = dE*scaleFactor(eta);
      MaxCalHit = calHit;
      geomBemc->getId(modN, etaN, subN, MaxSoftId);
    }
  }
 
  cout<<"MaxCalHit = "<<MaxCalHit<<endl;
  //set class variables
  if (isele)
    eleEsum = sumE;
  else
    posEsum = sumE;
  
  return MaxSoftId;
}

int StUpsilonEmbedMaker::Init()
{
  mEmcDecoder = new StEmcDecoder;
  mBemcTables = new StBemcTables;
  // Create tree list
  mTreeList = new TList;

  // Create ntuples
  TString varlist = "eventId:runId:vx:vy:vz:vxRc:vyRc:vzRc:pid:ntracks:ntracksRec:grefmult:upsGeantId:upsM:upsCosTheta:upsRap:upsP:upsPt:upsPz:upsEta:upsPhi:eleP:elePt:elePz:eleEta:elePhi:posP:posPt:posPz:posEta:posPhi:eleEsum:eleAdc:eleEcluster:eleEtower:eleEtower2:eleEtower3:eleShape:eleEtaTower:elePhiTower:eleRcSoftID:eleRcAdc:eleRcEnergy:eleRczD:eleRcphiD:eleMaxNeighborId:eleMaxNeighborAdc:eleMaxNeighborEnergy:posEsum:posAdc:posEcluster:posEtower:posEtower2:posEtower3:posShape:posEtaTower:posPhiTower:posRcSoftID:posRcAdc:posRcEnergy:posRczD:posRcphiD:posMaxNeighborId:posMaxNeighborAdc:posMaxNeighborEnergy:L2cosTheta:L2invM:isL0:isL2:eleReco:eleRcDca:elePRc:elePtRc:elePxRc:elePyRc:elePzRc:eleEtaRc:elePhiRc:eleR:eleRsmd:eleAcceptTrack:eleMcFirstTpcP:eleMcLastTpcP:eleMcFirstTpcX:eleMcFirstTpcY:eleMcFirstTpcZ:eleMcLastTpcX:eleMcLastTpcY:eleMcLastTpcZ:eleRcFirstTpcX:eleRcFirstTpcY:eleRcFirstTpcZ:eleRcLastTpcX:eleRcLastTpcY:eleRcLastTpcZ:eleTpcHits:eleTpcPos:posReco:posRcDca:posPRc:posPtRc:posPxRc:posPyRc:posPzRc:posEtaRc:posPhiRc:posR:posRsmd:posAcceptTrack:posMcFirstTpcP:posMcLastTpcP:posMcFirstTpcX:posMcFirstTpcY:posMcFirstTpcZ:posMcLastTpcX:posMcLastTpcY:posMcLastTpcZ:posRcFirstTpcX:posRcFirstTpcY:posRcFirstTpcZ:posRcLastTpcX:posRcLastTpcY:posRcLastTpcZ:posTpcHits:posTpcPos:upsMRc:upsPRc:upsPtRc:upsPzRc:cosThetaRc";
  mTuple = new TNtuple("ups", "Upsion ntuple", varlist);
  mTreeList->Add(mTuple);

  return StMaker::Init();
}

int StUpsilonEmbedMaker::InitRun(int runNumber)
{
  mEmcDecoder->SetDateTime(GetDate(), GetTime());
  mBemcTables->loadTables(this);
  return StMaker::InitRun(runNumber);
}

int StUpsilonEmbedMaker::Make()
{
  const double emass = 0.511e-3;
  int L0HTthreshold = 18;
  int L2HTthreshold = 80;
  float E1clusterthreshold = 4.5; //gev
  float E2clusterthreshold = 3.0; //gev
  float cosThetathreshold = 0.0;
  float invMassHighthreshold = 25;
  float invMassLowthreshold = 6.5;

  // Get StMcEvent
  mMcEvent = (StMcEvent*)GetDataSet("StMcEvent");
  if (!mMcEvent) {
    LOG_WARN << "No StMcEvent" << endm;
    return kStWarn;
  }
 
  // Print StMcEvent
  LOG_INFO << "eventGeneratorEventLabel = " << mMcEvent->eventGeneratorEventLabel() << endm;
  LOG_INFO << "eventNumber = " << mMcEvent->eventNumber() << endm;
  /*LOG_INFO << "runNumber = " << mMcEvent->runNumber() << endm;
  LOG_INFO << "type = " << mMcEvent->type() << endm;
  LOG_INFO << "zWest = " << mMcEvent->zWest() << endm;
  LOG_INFO << "nWest = " << mMcEvent->nWest() << endm;
  LOG_INFO << "zEast = " << mMcEvent->zEast() << endm;
  LOG_INFO << "nEast = " << mMcEvent->nEast() << endm;
  LOG_INFO << "eventGeneratorFinalStateTracks = " << mMcEvent->eventGeneratorFinalStateTracks() << endm;
  LOG_INFO << "numberOfPrimaryTracks = " << mMcEvent->numberOfPrimaryTracks() << endm;
  LOG_INFO << "subProcessId = " << mMcEvent->subProcessId() << endm;
  LOG_INFO << "impactParameter = " << mMcEvent->impactParameter() << endm;
  LOG_INFO << "phiReactionPlane = " << mMcEvent->phiReactionPlane() << endm;
  LOG_INFO << "triggerTimeOffset = " << mMcEvent->triggerTimeOffset() << endm;
  LOG_INFO << "nBinary = " << mMcEvent->nBinary() << endm;
  LOG_INFO << "nWoundedEast = " << mMcEvent->nWoundedEast() << endm;
  LOG_INFO << "nWoundedWest = " << mMcEvent->nWoundedWest() << endm;
  LOG_INFO << "nJets = " << mMcEvent->nJets() << endm;*/

  // Get Pythia record
  TDataSet* geant = GetDataSet("geant");
  TDataSetIter geantIter(geant);

  St_g2t_pythia* pythia = (St_g2t_pythia*)geantIter("g2t_pythia");

  if (pythia) {
    g2t_pythia_st* g2t_pythia = pythia->GetTable();
    assert(g2t_pythia);
    LOG_INFO << "Pythia subprocess_id = " << g2t_pythia->subprocess_id << endm;
  }
  else {
    LOG_WARN << "No Pythia record" << endm;
  }

  // Get StEvent
  mEvent = (StEvent*)GetDataSet("StEvent");
  if (!mEvent) {
    LOG_WARN << "No StEvent" << endm;
    return kStWarn;
  }

  Int_t mgrefmult = getgRefMult(mEvent,mMcEvent->primaryVertex()->position());
  cout << "HEYYY grefmult is " << mgrefmult << endl;

  // Get BBC trigger maker
  // StBbcTriggerMaker* bbcTrig = (StBbcTriggerMaker*)GetMakerInheritsFrom("StBbcTriggerMaker");
  // assert(bbcTrig);

  // LOG_INFO << "BBC = " << bbcTrig->isTrigger() << endm;

  // Get EMC trigger maker
  StEmcTriggerMaker* emcTrig = (StEmcTriggerMaker*)GetMakerInheritsFrom("StEmcTriggerMaker");
  assert(emcTrig);

  LOG_INFO << "Upsilon(320501) = " << emcTrig->isTrigger(320501) << endm;
  LOG_INFO << "Upsilon(330501) = " << emcTrig->isTrigger(330501) << endm;

  // Get Upsilon trigger
  //mUpsTrig = (StUpsilonTriggerMaker*)GetMakerInheritsFrom("StUpsilonTriggerMaker");
  //assert(mUpsTrig);

  //LOG_INFO << "Upsilon(L2) = " << mUpsTrig->isUpsilonTrigger() << endm;

  // Print particles in event
  LOG_INFO << "Subprocess Id = " << mMcEvent->subProcessId() << endm;
  LOG_INFO << "Ntracks = " << mMcEvent->numberOfPrimaryTracks() << endm;
  cout << "Particles: ";
  for (size_t i = 0; i < mMcEvent->primaryVertex()->numberOfDaughters(); ++i) {
    StMcTrack* mcTrack = mMcEvent->primaryVertex()->daughter(i);
    if (mcTrack->particleDefinition())
      cout << mcTrack->particleDefinition()->name() << " ";
    else
      cout << mcTrack->geantId() << " ";
  }
  cout << endl;

  //Save all the pedestal subtracted ADC values and energies
  // cout<<"loading barrel values"<<endl;
  StEmcDetector* bemcDet = mEvent->emcCollection()->detector(kBarrelEmcTowerId);
  //cout<<"barrel loaded"<<endl;
  for (int i = 0;i<4801;i++){    bemcEnergy[i] = 0;bemcADC[i]=0; }// bemczPos[i]=0; bemcphiPos[i]=0;}
  // cout<<"Set arrays to 0"<<endl;
  // cout<<"# of modules = "<<bemcDet->numberOfModules()<<endl;
  for (unsigned int m = 1; m<=bemcDet->numberOfModules(); ++m){
    //cout<<"module "<<m<<" out of "<<bemcDet->numberOfModules()<<endl;
    StSPtrVecEmcRawHit& hits = bemcDet->module(m)->hits();
    for (StSPtrVecEmcRawHitIterator i = hits.begin(); i != hits.end(); ++i){
      StEmcRawHit* hit = *i;
      if (hit->energy()<=0) continue;
      int softId;
      int modN = hit->module();
      int etaN = hit->eta();
      int subN = hit->sub();
      StEmcGeom *geomBemc = StEmcGeom::getEmcGeom(1);
      geomBemc->getId(modN, etaN, subN, softId);
      int status;
      Float_t pedestal, rms;
      mBemcTables->getStatus(1, softId, status);
      mBemcTables->getPedestal(1,softId,0,pedestal,rms);
      if (hit->adc() < pedestal + 3 * rms) continue;
      bemcADC[softId]= hit->adc();
      //cout << "Id: " << softId << ", ADC: " << hit->adc() << endl;
      bemcEnergy[softId] = hit->energy();

      
      //if (bemcADC[softId]>75)
	//cout<<"bemcADC["<<softId<<"] = "<<bemcADC[softId]<<" energy = "<<bemcEnergy[softId]<<endl;
    }
  }

  // loop over barrelhits to get zpos and phipos
  StSPtrVecEmcPoint& bEmcPoints = mEvent->emcCollection()->barrelPoints();
  int index=0;
  float mindist=1.e9;
  cout << "hmmBLAA" << endl;
  cout << "size is " << bEmcPoints.size() << endl;
  //mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
  for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++, index++) {
    bool associated=false;
    StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
    cout << "neeBLAA" << endl;
    if(bEmcClusters.size()==0 ) continue;
    if(bEmcClusters[0]==NULL) continue;
    cout << "jooBLAA" << endl;
  }


  // Loop over primary tracks
  for (size_t i = 0; i < mMcEvent->primaryVertex()->numberOfDaughters(); ++i) {
    StMcTrack* upsMc = mMcEvent->primaryVertex()->daughter(i);

    // Get Upsilon
    if (upsMc->geantId() >= 160) {
      StMcVertex* UpsStartVertex = upsMc->startVertex();
      if (UpsStartVertex->position() != mMcEvent->primaryVertex()->position()){
	cout<<"Upsilon not at pythia event vertex, bailing out"<<endl;
	return kStOk;}
      
      // Loop over daughters and asign electron and positron
      StMcTrack* eleMc = 0;
      StMcTrack* posMc = 0;
      if (!upsMc->stopVertex())
        {
          cout<<"Found upsilon without stop vertex, bailing out"<<endl;
          return kStOk;
        }

      for (size_t j = 0; j < upsMc->stopVertex()->numberOfDaughters(); ++j) {
	StMcTrack* mcTrack = upsMc->stopVertex()->daughter(j);

	switch (mcTrack->geantId()) {
	case 2:
	  posMc = mcTrack;
	  break;
	case 3:
	  eleMc = mcTrack;
	  break;
	default:
	  LOG_WARN << "Not an electron nor a positron - geantId = " << mcTrack->geantId() << endm;
	  return kStWarn;
	}
      }
      int eleSoftId = 0;int posSoftId = 0;
      // Opening angle between daughters
      float cosThetaMc = eleMc->momentum().unit() * posMc->momentum().unit();
      // cout<<"ele p = "<<eleMc->momentum()<<" ";
      // cout<<"pos p = "<<posMc->momentum()<<" ";
      // cout<<"cosTheta is "<<cosThetaMc<<endl;
      
      int eleAdc = -1000;
      int eleTP = -1000;
      int posAdc = -1000;
      int posTP = -1000;
      bool isL0 = false;
      float eleEcluster = -1000;
      float posEcluster = -1000;
      pair<float,float> eleEtaphi;pair<float,float> posEtaphi;
      eleEtaphi.first = -1000;eleEtaphi.second = -1000;
      posEtaphi.first = -1000;posEtaphi.second = -1000;
      bool isFound = false;
      float L2costheta=-1000;
      float L2invM=-1000;
      int nConeele=0;
      int nConepos=0;

      float eleMcFirstTpcP = -1000;
      float eleMcLastTpcP = -1000;
      float posMcFirstTpcP = -1000;
      float posMcLastTpcP = -1000;
      StThreeVectorF eleMcFirstTpcPos;
      StThreeVectorF eleMcLastTpcPos;
      StThreeVectorF posMcFirstTpcPos;
      StThreeVectorF posMcLastTpcPos;
      float eleEtower2 = -1000;
      float eleEtower3 = -1000;
      float posEtower2 = -1000;
      float posEtower3 = -1000;
      int eleShape = -1000;
      int posShape = -1000;

      if (eleMc->tpcHits().size()){
	eleMcFirstTpcP = eleMc->tpcHits().front()->localMomentum().mag();
	eleMcLastTpcP = eleMc->tpcHits().back()->localMomentum().mag();
	eleMcFirstTpcPos = eleMc->tpcHits().front()->position();
	eleMcLastTpcPos = eleMc->tpcHits().back()->position();
      }
      if (posMc->tpcHits().size()){
	posMcFirstTpcP = posMc->tpcHits().front()->localMomentum().mag();
	posMcLastTpcP = posMc->tpcHits().back()->localMomentum().mag();
	posMcFirstTpcPos = posMc->tpcHits().front()->position();
	posMcLastTpcPos = posMc->tpcHits().back()->position();
      }

      eleSoftId = SetHighTowerVar(eleMc,true);
      posSoftId = SetHighTowerVar(posMc,false);

      if (eleSoftId>0){
	eleAdc = bemcADC[eleSoftId];}
      if (posSoftId>0){
	posAdc = bemcADC[posSoftId];}

      if ( (eleAdc>>4)>L0HTthreshold )
	isL0 = true;
      else if ( (posAdc>>4)>L0HTthreshold )
	isL0 = true;
      if ( eleSoftId>0 ){
	getClusterEtaPhiE(eleSoftId,eleEtaphi,eleEcluster,eleEtower2,eleEtower3,eleShape);}
      if ( posSoftId>0 ){
	getClusterEtaPhiE(posSoftId,posEtaphi,posEcluster,posEtower2,posEtower3,posShape);}
      if ((eleSoftId>0)&&(posSoftId>0)){
	
	StEmcGeom *geomBemc = StEmcGeom::getEmcGeom(1);
	Float_t ex;Float_t ey;Float_t ez;
	Float_t px;Float_t py;Float_t pz;
	geomBemc->getXYZ(eleSoftId,ex,ey,ez);
  cout << "disttt geombemc is " << ez << " and " << eleEtaphi.second << endl;
	geomBemc->getXYZ(posSoftId,px,py,pz);
	StThreeVectorF eleTower3vec(ex,ey,ez);
	StThreeVectorF posTower3vec(px,py,pz);
	L2costheta=eleTower3vec.unit()*posTower3vec.unit();


	
	L2invM=sqrt(2*eleEcluster*posEcluster*(1-L2costheta));}
      

  StEmcGeom *emcGeomBSMDE = StEmcGeom::getEmcGeom("bsmde");
  StEmcGeom *emcGeomBSMDP = StEmcGeom::getEmcGeom("bsmdp");
  StEmcGeom *emcGeomBEMC = StEmcGeom::getEmcGeom(1);

      if ((isL0)&&(((eleEcluster>E1clusterthreshold)&&(posEcluster>E2clusterthreshold))||((eleEcluster>E2clusterthreshold)&&(posEcluster>E1clusterthreshold)))&&(L2costheta<cosThetathreshold)&&(L2invM>invMassLowthreshold)&&(L2invM<invMassHighthreshold)&&(((eleAdc>L2HTthreshold)&&((posAdc>>4)>L0HTthreshold))||((posAdc>L2HTthreshold)&&((eleAdc>>4)>L0HTthreshold))))
	isFound = true;


    
  
      //**********
      
      StTrack* eleRc = findPartner(eleMc);
      StTrack* posRc = findPartner(posMc);
      float cosThetaRc = -1000;
      int eleTPCpos = -1000;
      int posTPCpos = -1000;
      short eleFlag = -1000;
      short posFlag = -1000;
      float eleTPChits = -1000;
      float posTPChits = -1000;
      bool elereco = false;
      bool posreco = false;
      float eleRcDca = -1000, posRcDca = -1000;
      pair<float,float> eleetaphiRcTrack, eleetaphiRcSmd;
      pair<float,float> posetaphiRcTrack, posetaphiRcSmd;
      eleetaphiRcTrack.first = -1000;
      eleetaphiRcTrack.second = -1000;
      posetaphiRcTrack.first = -1000;
      posetaphiRcTrack.second = -1000;
      eleetaphiRcSmd.first = -1000;
      eleetaphiRcSmd.second = -1000;
      posetaphiRcSmd.first = -1000;
      posetaphiRcSmd.second = -1000;
      StThreeVectorF eleRcFirstTpcPos;
      StThreeVectorF eleRcLastTpcPos;
      StThreeVectorF posRcFirstTpcPos;
      StThreeVectorF posRcLastTpcPos;
      float eleR = -1000;
      float posR = -1000;
      float eleRsmd = -1000;
      float posRsmd = -1000;
      int eleRcSoftid = -1000, posRcSoftid = -1000;
      bool eleAccept = false;
      bool posAccept = false;
      int eleRcTowerADC = -1000, posRcTowerADC = -1000;
      float eleRcTowerEnergy = -1000, posRcTowerEnergy = -1000;
      float eleRczDist = -1000, posRczDist = -1000;
      float eleRcphiDist = -1000, posRcphiDist = -1000;
      int eleMaxNeighbor = -1000, posMaxNeighbor = -1000;
      float eleMaxNeighborEnergy = -1000, posMaxNeighborEnergy = -1000;
      int eleMaxNeighborADC = -1000, posMaxNeighborADC = -1000;
      StThreeVectorF eleMom;StThreeVectorF posMom;StThreeVectorF upsMom;

      if (eleRc){
	elereco = true;
	eleMom = eleRc->geometry()->momentum();

  StThreeVectorF primVpos = mMcEvent->primaryVertex()->position();
  //eleRcDca = eleRc->geometry()->helix().distance(primVpos);
  float dca2 = eleRc->outerGeometry()->helix().distance(primVpos);
  eleRcDca = dca2;
  //float dca3 = eleRc->extGeometry()->helix().distance(primVpos);
  //cout << "dca is " << eleRcDca << " vs " << dca2 << " vs " << endl;

	eleFlag = eleRc->flag();
	TrackEtaPhiId(eleRc,eleetaphiRcTrack,eleetaphiRcSmd,eleRcSoftid);


	eleR = sqrt((eleetaphiRcTrack.first-eleEtaphi.first)*(eleetaphiRcTrack.first-eleEtaphi.first)+(eleetaphiRcTrack.second-eleEtaphi.second)*(eleetaphiRcTrack.second-eleEtaphi.second));
	eleRsmd = sqrt((eleetaphiRcSmd.first-eleEtaphi.first)*(eleetaphiRcSmd.first-eleEtaphi.first)+(eleetaphiRcSmd.second-eleEtaphi.second)*(eleetaphiRcSmd.second-eleEtaphi.second));
	eleRcFirstTpcPos = eleRc->detectorInfo()->hits(kTpcId).front()->position();
	eleRcLastTpcPos = eleRc->detectorInfo()->hits(kTpcId).back()->position();
	eleTPCpos = eleRc->numberOfPossiblePoints(kTpcId);	
	eleTPChits = eleRc->fitTraits().numberOfFitPoints(kTpcId);
	
	if(eleRcSoftid>0 && eleRcSoftid<4801){
	  eleMaxNeighbor = getHighestNeighbor(eleRcSoftid);
	  eleRcTowerADC = bemcADC[eleRcSoftid];
	  eleRcTowerEnergy = bemcEnergy[eleRcSoftid];

    float trackPos[2];
    cout << "blaa rc id is " << eleRcSoftid << endl;
    getZPosPhiPos(eleRc,trackPos);
    eleRczDist = trackPos[0];
    eleRcphiDist = trackPos[1];

    //zDist, phiDist
    /*float zPosRaw, xbla, ybla, phiPosRaw;
    float bemczPos, bemcphiPos;
    emcGeomBEMC->getXYZ(eleRcSoftid,xbla,ybla,zPosRaw);
    emcGeomBEMC->getPhi(eleRcSoftid,phiPosRaw);
    bemczPos = zPosRaw ;//hit->position().z();
    bemcphiPos = phiPosRaw;//hit->position().phi();
    eleRczDist = bemczPos - trackPos[0];
    cout << "zdisttt is " << bemczPos << " - " << trackPos[0] << endl;
    eleRcphiDist = bemcphiPos - trackPos[1];
    cout << "phidisttt is " << bemcphiPos << " - " << trackPos[1] << endl;
    cout << "disttt tpc is " << eleRcLastTpcPos.z() << " and " << eleetaphiRcTrack.second << endl; 
    if(eleRcphiDist>=TMath::Pi()) eleRcphiDist=eleRcphiDist-TMath::TwoPi();
    if(eleRcphiDist<-TMath::Pi()) eleRcphiDist=eleRcphiDist+TMath::TwoPi();*/
    //////////
    if(eleMaxNeighbor>0 && eleMaxNeighbor<4801){
	    eleMaxNeighborADC = bemcADC[eleMaxNeighbor];
	    eleMaxNeighborEnergy = bemcEnergy[eleMaxNeighbor];
	  } 
	}
	if (eleFlag>0&&eleFlag!=701&&eleFlag!=801&&eleFlag!=901&&eleFlag<1000&&eleTPChits>20&&eleTPChits/eleTPCpos>0.52&&eleMom.perp()>0.2)
	  eleAccept = true;
      }
      if (posRc){
	posreco = true;
	posMom = posRc->geometry()->momentum();
  StThreeVectorF primVpos = mMcEvent->primaryVertex()->position();
  //eleRcDca = eleRc->geometry()->helix().distance(primVpos);
  float dca2 = posRc->outerGeometry()->helix().distance(primVpos);
  posRcDca = dca2;
	posFlag = posRc->flag();
	TrackEtaPhiId(posRc,posetaphiRcTrack,posetaphiRcSmd,posRcSoftid);
	posR = sqrt((posetaphiRcTrack.first-posEtaphi.first)*(posetaphiRcTrack.first-posEtaphi.first)+(posetaphiRcTrack.second-posEtaphi.second)*(posetaphiRcTrack.second-posEtaphi.second));
	posRsmd = sqrt((posetaphiRcSmd.first-posEtaphi.first)*(posetaphiRcSmd.first-posEtaphi.first)+(posetaphiRcSmd.second-posEtaphi.second)*(posetaphiRcSmd.second-posEtaphi.second));
	posRcFirstTpcPos = posRc->detectorInfo()->hits(kTpcId).front()->position();
	posRcLastTpcPos = posRc->detectorInfo()->hits(kTpcId).back()->position();
	posTPCpos = posRc->numberOfPossiblePoints(kTpcId);
	posTPChits = posRc->fitTraits().numberOfFitPoints(kTpcId);

	if(posRcSoftid>0 && posRcSoftid<4801){
	  posMaxNeighbor = getHighestNeighbor(posRcSoftid);
	  posRcTowerADC = bemcADC[posRcSoftid];
	  posRcTowerEnergy = bemcEnergy[posRcSoftid];
    float trackPos2[2];
    cout << "blaa rc id is " << posRcSoftid << endl;
    getZPosPhiPos(posRc,trackPos2);
    posRczDist = trackPos2[0];
    posRcphiDist = trackPos2[1];
	  if(posMaxNeighbor>0 && posMaxNeighbor<4801){
	    posMaxNeighborADC = bemcADC[posMaxNeighbor];
	    posMaxNeighborEnergy = bemcEnergy[posMaxNeighbor];
	  } 
	}
	if (posFlag>0&&posFlag!=701&&posFlag!=801&&posFlag!=901&&posFlag<1000&&posTPChits>20&&posTPChits/posTPCpos>0.52&&posMom.perp()>0.2)
	  posAccept = true;
      }
      
   
      
      StLorentzVectorF upsRc;
      if (eleRc && posRc) {	
	StLorentzVectorF p1(eleMom,eleMom.massHypothesis(emass));
        StLorentzVectorF p2(posMom,posMom.massHypothesis(emass));
        upsRc = p1 + p2;
	StThreeVectorF upsMom = eleMom+posMom;
	cosThetaRc = eleMom.unit() * posMom.unit();	
      }

      cout<<"Getting vertex daughters"<<endl; 
      StSPtrVecPrimaryTrack& recoPrimTracks = mEvent->primaryVertex()->daughters();
      // cout<<"There are "<<recoPrimTracks.size()<<" reconstructed tracks"<<endl;
      // cout<<"ele eta phi = "<<eleEtaphi.first<<","<<eleEtaphi.second<<endl;
      // cout<<"pos eta phi = "<<posEtaphi.first<<","<<posEtaphi.second<<endl;
      for (int it = 0;it < recoPrimTracks.size();++it){
	StTrack* currTrack = recoPrimTracks[it];
	pair<float,float> currEtaPhi;

	short currFlag = currTrack->flag();
	StThreeVectorF currMom = currTrack->geometry()->momentum();
	float currTPCpos = currTrack->numberOfPossiblePoints(kTpcId);
	float currTPChits = currTrack->fitTraits().numberOfFitPoints(kTpcId);
	if (currFlag>0&&currFlag!=701&&currFlag!=801&&currFlag!=901&&currTPChits>20&&currTPChits/currTPCpos>0.52&&currMom.perp()>0.2){
	  pair<float,float> trash;
	  int trash2;
	  TrackEtaPhiId(currTrack,currEtaPhi,trash,trash2);
	  // cout<<"curr track "<<it<<" eta phi = "<<currEtaPhi.first<<","<<currEtaPhi.second<<endl;
	  float phiele = eleEtaphi.first - currEtaPhi.first;
	  float phipos = posEtaphi.first - currEtaPhi.first;
	  
	  if (phiele > 3.1416)
	    phiele = phiele - 2*3.1416;
	  if (phiele < -3.1416)
	    phiele = phiele + 2*3.1416;
	  if (phipos > 3.1416)
	    phipos = phipos - 2*3.1416;
	  if (phipos < -3.1416)
	    phipos = phipos + 2*3.1416;
	  
	  float Rele = sqrt((eleEtaphi.first-currEtaPhi.first)*(eleEtaphi.first-currEtaPhi.first)+phiele*phiele);
	  float Rpos = sqrt((posEtaphi.first-currEtaPhi.first)*(posEtaphi.first-currEtaPhi.first)+phipos*phipos);
	  if (Rele<0.2)
	    nConeele++;
	  if (Rpos<0.2)
	    nConepos++;
	}
      }
      
      //**********

      // Fill tuple
      float tuple[] = {
  mEvent->id(),
  mEvent->runInfo()->runId(),
	UpsStartVertex->position().x(),
	UpsStartVertex->position().y(),
	UpsStartVertex->position().z(),
	mMcEvent->primaryVertex()->position().x(),
	mMcEvent->primaryVertex()->position().y(),
	mMcEvent->primaryVertex()->position().z(),
	mMcEvent->subProcessId(),
	mMcEvent->numberOfPrimaryTracks(),
	mEvent->primaryVertex()->daughters().size(),
  mgrefmult,

	upsMc->geantId(),
	upsMc->fourMomentum().m(),
	cosThetaMc,
	upsMc->rapidity(),
	upsMc->momentum().mag(),
	upsMc->pt(),
	upsMc->momentum().z(),
	upsMc->pseudoRapidity(),
	upsMc->momentum().phi(),

	eleMc->momentum().mag(),
	eleMc->pt(),
	eleMc->momentum().z(),
	eleMc->pseudoRapidity(),
	eleMc->momentum().phi(),

	posMc->momentum().mag(),
	posMc->pt(),
	posMc->momentum().z(),
	posMc->pseudoRapidity(),
	posMc->momentum().phi(),

	eleEsum,
	eleAdc,
	eleEcluster,
	bemcEnergy[eleSoftId],
	eleEtower2,
	eleEtower3,
	eleShape,
	eleEtaphi.first,
	eleEtaphi.second,
	eleRcSoftid,
	eleRcTowerADC,
	eleRcTowerEnergy,
  eleRczDist,
  eleRcphiDist,
	eleMaxNeighbor,
	eleMaxNeighborADC,
	eleMaxNeighborEnergy,

	posEsum,
	posAdc,
	posEcluster,
	bemcEnergy[posSoftId],
	posEtower2,
	posEtower3,
	posShape,
	posEtaphi.first,
	posEtaphi.second,
	posRcSoftid,
	posRcTowerADC,
	posRcTowerEnergy,
  posRczDist,
  posRcphiDist,
	posMaxNeighbor,
	posMaxNeighborADC,
	posMaxNeighborEnergy,

	L2costheta,
	L2invM,
	isL0,
	isFound,
	
	elereco,
  eleRcDca,
	eleMom.mag(),
	eleMom.perp(),
	eleMom.x(),
	eleMom.y(),
	eleMom.z(),
	eleetaphiRcTrack.first,
	eleetaphiRcTrack.second,
	eleR,
	eleRsmd,
	eleAccept,
	eleMcFirstTpcP,
	eleMcLastTpcP,
	eleMcFirstTpcPos.x(),
	eleMcFirstTpcPos.y(),
	eleMcFirstTpcPos.z(),
	eleMcLastTpcPos.x(),
	eleMcLastTpcPos.y(),
	eleMcLastTpcPos.z(),
	eleRcFirstTpcPos.x(),
	eleRcFirstTpcPos.y(),
	eleRcFirstTpcPos.z(),
	eleRcLastTpcPos.x(),
	eleRcLastTpcPos.y(),
	eleRcLastTpcPos.z(),
	eleTPChits,
	eleTPCpos,

	posreco,
  posRcDca,
	posMom.mag(),
	posMom.perp(),
	posMom.x(),
	posMom.y(),
	posMom.z(),
	posetaphiRcTrack.first,
	posetaphiRcTrack.second,
	posR,
	posRsmd,
	posAccept,
	posMcFirstTpcP,
	posMcLastTpcP,
	posMcFirstTpcPos.x(),
	posMcFirstTpcPos.y(),
	posMcFirstTpcPos.z(),
	posMcLastTpcPos.x(),
	posMcLastTpcPos.y(),
	posMcLastTpcPos.z(),
	posRcFirstTpcPos.x(),
	posRcFirstTpcPos.y(),
	posRcFirstTpcPos.z(),
	posRcLastTpcPos.x(),
	posRcLastTpcPos.y(),
	posRcLastTpcPos.z(),
	posTPChits,
	posTPCpos,

	upsRc.m(),
	upsRc.vect().mag(),
	upsRc.perp(),
	upsRc.pz(),
	cosThetaRc
      };
 

      mTuple->Fill(tuple);
    }
  }

  return kStOk;
}

StTrack* StUpsilonEmbedMaker::findPartner(StMcTrack* mcTrack)
{
  StAssociationMaker* assoc = (StAssociationMaker*)GetMakerInheritsFrom("StAssociationMaker");
  pair<mcTrackMapIter, mcTrackMapIter> p = assoc->mcTrackMap()->equal_range(mcTrack);
  StTrack* maxTrack = 0;
  int maxCommonTpcHits = 0;
  for (mcTrackMapIter k = p.first; k != p.second; ++k) {
    int commonTpcHits = k->second->commonTpcHits();
    StTrack* track = (StTrack*)k->second->partnerTrack()->node()->track(primary);
    if (track && commonTpcHits > maxCommonTpcHits) {
      maxTrack = track;
      maxCommonTpcHits = commonTpcHits;
    }
  }
  return maxTrack;
}

void StUpsilonEmbedMaker::getClusterEtaPhiE(int& id, pair<float,float>& etaphi, float& Ecluster, float& E2, float& E3, int& shape)
{
  if (id == 0){
    cout<<"soft id is zero!!! exit etaphi"<<endl;
    return;}
  float E1 = bemcEnergy[id];

  vector<int> neighbors;
  StEmcPosition* mEmcPosition = new StEmcPosition;
  for (int deta = -1; deta <= 1; ++deta) {
    for (int dphi = -1; dphi <= 1; ++dphi) {
      if (deta || dphi) {
        int id2 = mEmcPosition->getNextTowerId(id, deta, dphi);
        if (id2) {
          neighbors.push_back(id2);
        }
      }
    }
  }
  delete mEmcPosition;

  E2 = 0;E3 = 0;
  int ID2 = 0;
  int ID3 = 0;

  for (int i = 0;i<neighbors.size();i++){
    int idcur = neighbors[i];
    float Ecur = bemcEnergy[idcur];
    if (Ecur > E2){
      E3 = E2;
      ID3 = ID2;
      E2 = Ecur;
      ID2 = idcur;}
    else if (Ecur > E3){
      E3 = Ecur;
      ID3 = idcur;}}
 
  StEmcGeom *geomBemc = StEmcGeom::getEmcGeom(1);
  Ecluster = E1+E2+E3;
  float eta1=0;float eta2=0;float eta3=0;
  float phi1=0;float phi2=0;float phi3=0;
  float xCluster, yCluster;
  geomBemc->getEta(id,eta1);
  geomBemc->getPhi(id,phi1);
  if (ID2 > 0){
    geomBemc->getEta(ID2,eta2);
    geomBemc->getPhi(ID2,phi2);}
  if (ID3 > 0){
    geomBemc->getEta(ID3,eta3);
    geomBemc->getPhi(ID3,phi3);}
  if ((E1+E2+E3)>0){
    etaphi.first = (eta1*E1+eta2*E2+eta3*E3)/Ecluster;
    if(isnan(etaphi.first)){
      cout <<"NAN!NAN!NAN!" << endl;
      cout << "eta1: " << eta1 << ", E1: " << E1 << endl;
      cout << "eta2: " << eta2 << ", E2: " << E2 <<endl;
      cout << "eta3: " << eta3 << ", E3: " << E3 <<endl;
    }
    xCluster = (cos(phi1)*E1+cos(phi2)*E2+cos(phi3)*E3)/Ecluster;
    yCluster = (sin(phi1)*E1+sin(phi2)*E2+sin(phi3)*E3)/Ecluster;
    etaphi.second = atan2(yCluster,xCluster);}
  shape = FindShape(id,ID2,ID3);
  return;
}

void StUpsilonEmbedMaker::getTowerNeighbors(int id, vector<int>& neighbors)
{
  StEmcPosition* mEmcPosition = new StEmcPosition;
  for (int deta = -1; deta <= 1; ++deta) {
    for (int dphi = -1; dphi <= 1; ++dphi) {
      if (deta || dphi) {
        int id2 = mEmcPosition->getNextTowerId(id, deta, dphi);
        if (id2) {
          neighbors.push_back(id2);
        }
      }
    }
  }
  delete mEmcPosition;
}

int StUpsilonEmbedMaker::getHighestNeighbor(int myId){

  vector<int> neighbors;

  getTowerNeighbors(myId,neighbors);
  neighbors.push_back(myId);

  int maxID  = -1;
  int maxAdc = -1;

  for(int i=0; i < neighbors.size(); i++){
    if(bemcADC[neighbors[i]] > maxAdc){
      maxID=neighbors[i];
      maxAdc=bemcADC[neighbors[i]];
    }
  }

  return maxID;
}

int StUpsilonEmbedMaker :: FindShape(int& id1, int& id2, int& id3){
  if ((id2==0)&&(id3==0))
    return 0;
  StEmcPosition* mEmcPosition = new StEmcPosition;
  // for (int deta = -1; deta <= 1; ++deta) {
  //  for (int dphi = -1; dphi <= 1; ++dphi) {
  //    if (deta || dphi) {
  //    int id2 = mEmcPosition->getNextTowerId(id, deta, dphi);
  if (id3==0){
    vector<int> sneighbors;
    sneighbors.push_back(mEmcPosition->getNextTowerId(id1,-1,0));
    sneighbors.push_back(mEmcPosition->getNextTowerId(id1,1,0));
    sneighbors.push_back(mEmcPosition->getNextTowerId(id1,0,-1));
    sneighbors.push_back(mEmcPosition->getNextTowerId(id1,0,1));
    for (int i = 0;i<sneighbors.size();++i){
      if (sneighbors[i] == id2)
        return 1;
    }
  return 2;
  }
  //L shape only happens if id2 and id3 are neighbors
  vector<int> neighbors;
  getTowerNeighbors(id2,neighbors);
  for (int i = 0;i<neighbors.size();++i){
    if (id3 == neighbors[i])
      return 3;//L shape
  }
  if((id3 = mEmcPosition->getNextTowerId(id2,-2,0)&&(id1 = mEmcPosition->getNextTowerId(id2,-1,0)))||(id3 = mEmcPosition->getNextTowerId(id2,2,0)&&(id1 = mEmcPosition->getNextTowerId(id2,1,0)))||(id3 = mEmcPosition->getNextTowerId(id2,0,-2)&&(id1 = mEmcPosition->getNextTowerId(id2,0,-1)))||(id3 = mEmcPosition->getNextTowerId(id2,0,2)&&(id1 = mEmcPosition->getNextTowerId(id2,0,1))))
    return 4; //bar shape
  if ((id3 = mEmcPosition->getNextTowerId(id2,-2,-2))||(id3 = mEmcPosition->getNextTowerId(id2,2,2)))
    return 5; //diagonal
  //only choices left are bar w/ diagonal or two diagonals with same horizontal or vertical location
  if ((id3 = mEmcPosition->getNextTowerId(id2,-2,0))||(id3 = mEmcPosition->getNextTowerId(id2,2,0))||(id3 = mEmcPosition->getNextTowerId(id2,0,-2))||(id3 = mEmcPosition->getNextTowerId(id2,0,-2)))
    return 7; //2 diagonals...
  return 6; //only one left is bar with diagonal!
}

void StUpsilonEmbedMaker :: TrackEtaPhiId(const StTrack* mcTrack, pair<float,float>& etaphi, pair<float,float>& etaphiSmd, int& softid)
{
  //cout<<"in track etaphi StTrack ";
  StEmcPosition* projPos = new StEmcPosition;
  StThreeVectorD* PosVec = new StThreeVectorD;
  StThreeVectorD* MomVec = new StThreeVectorD;
  StEmcGeom *geomBemc = StEmcGeom::getEmcGeom(1);
  Double_t mMagneticField = mEvent->summary()->magneticField()*0.1;//kilogause/tesla
  //cout<<"Defined variables"<<endl;
  projPos->projTrack(PosVec,MomVec,mcTrack,mMagneticField);
  //cout<<"called projection"<<endl;
  etaphi.first = PosVec->pseudoRapidity();
  etaphi.second = PosVec->phi();
  geomBemc->getId(etaphi.second,etaphi.first,softid);

  projPos->projTrack(PosVec,MomVec,mcTrack,mMagneticField,231.72);
  etaphiSmd.first = PosVec->pseudoRapidity();
  etaphiSmd.second = PosVec->phi();

  //cout << "blaa origo id " << softid << " z " << PosVec->z() << " phi " << PosVec->phi() << endl;

  delete projPos;
  delete PosVec;
  delete MomVec;
  //delete geomBemc;
  return;
}

void StUpsilonEmbedMaker :: getZPosPhiPos(const StTrack* mcTrack, float* d)
{
  //cout<<"in track etaphi StTrack ";
  StEmcPosition* projPos = new StEmcPosition;
  StThreeVectorD* PosVecBSMDE = new StThreeVectorD;
  StThreeVectorD* MomVecBSMDE = new StThreeVectorD;
  StEmcGeom *geomBsmde = StEmcGeom::getEmcGeom("bsmde");
  float radiusBSMDE = geomBsmde->Radius(); //1, bemc
  Double_t mMagneticField = mEvent->summary()->magneticField()*0.1;//kilogause/tesla
  //cout<<"Defined variables"<<endl;
  projPos->projTrack(PosVecBSMDE,MomVecBSMDE,mcTrack,mMagneticField,radiusBSMDE); //propagating track to bsmde
  //cout<<"called projection"<<endl;

  pair<float,float> etaphiBSMDE;
  pair<float,float> etaphiSmd;
  int softidBSMDE;
  etaphiBSMDE.first = PosVecBSMDE->pseudoRapidity();
  etaphiBSMDE.second = PosVecBSMDE->phi();
  geomBsmde->getId(etaphiBSMDE.second,etaphiBSMDE.first,softidBSMDE); //getting bsmde id 
  Float_t blax,blay,blaz;
  geomBsmde->getXYZ(softidBSMDE,blax,blay,blaz); //getting z of bsmde cell
  //cout << "blaa 99999999 id " << softidBSMDE << " z " << PosVecBSMDE->z() << " z2 " << blaz << " phi "  << PosVecBSMDE->phi() << endl;
  //d[0] = blaz - PosVecBSMDE->z();
  

  StEmcGeom* geomBemc = StEmcGeom::getEmcGeom("bemc");
  StThreeVectorD* PosVecBEMC = new StThreeVectorD;
  StThreeVectorD* MomVecBEMC = new StThreeVectorD;
  float radiusBEMC = geomBemc->Radius();
  projPos->projTrack(PosVecBEMC,MomVecBEMC,mcTrack,mMagneticField,radiusBEMC);
  pair<float,float> etaphiBEMC;
  int softidBEMC;
  etaphiBEMC.first = PosVecBEMC->pseudoRapidity();
  etaphiBEMC.second = PosVecBEMC->phi();
  geomBemc->getId(etaphiBEMC.second,etaphiBEMC.first,softidBEMC);
  Float_t bemcx,bemcy,bemcz;
  geomBemc->getXYZ(softidBEMC,bemcx,bemcy,bemcz);
  Float_t bemcphi;
  geomBemc->getPhi(softidBEMC,bemcphi);

  d[0] = bemcz - PosVecBSMDE->z();

  StThreeVectorD* PosVecBSMDP = new StThreeVectorD;
  StThreeVectorD* MomVecBSMDP = new StThreeVectorD;
  StEmcGeom *geomBsmdp = StEmcGeom::getEmcGeom("bsmdp");
  float radiusBSMDP = geomBsmdp->Radius(); //1, bemc
  //Double_t mMagneticField = mEvent->summary()->magneticField()*0.1;//kilogause/tesla
  //cout<<"Defined variables"<<endl;
  projPos->projTrack(PosVecBSMDP,MomVecBSMDP,mcTrack,mMagneticField,radiusBSMDP); // propagating track to bsmdp
  //cout<<"called projection"<<endl;

  pair<float,float> etaphiBSMDP;
  int softidBSMDP;
  etaphiBSMDP.first = PosVecBSMDP->pseudoRapidity();
  etaphiBSMDP.second = PosVecBSMDP->phi();
  geomBsmdp->getId(etaphiBSMDP.second,etaphiBSMDP.first,softidBSMDP); // getting id of bsmdp
  Float_t blaphi;
  geomBsmde->getPhi(softidBSMDP,blaphi); //getting phi of bsmdp cell

  //float dphi = blaphi - PosVecBSMDP->phi();
  float dphi = bemcphi - PosVecBSMDP->phi();
  if(dphi>=TMath::Pi()) dphi=dphi-TMath::TwoPi();
    if(dphi<-TMath::Pi()) dphi=dphi+TMath::TwoPi();
  d[1] = dphi;

  //cout << "blaa 99999999 idp " << softidBSMDP << " z " << PosVecBSMDP->z() <<  " phi "  << PosVecBSMDP->phi() << " phi2 " << blaphi << endl;
  

  /*StEmcPosition* projPos = new StEmcPosition;
  Double_t mMagneticField = mEvent->summary()->magneticField()*0.1;//kilogause/tesla
  
  StEmcGeom* emcGeom = StEmcGeom::getEmcGeom("bsmdp");
  StEmcGeom* emcGeom0 = StEmcGeom::getEmcGeom("bsmde");
  StEmcGeom* emcGeomBEMC = StEmcGeom::getEmcGeom("bemc");
  float radius = emcGeom->Radius(); //3, bsmdp
  float radius0 = emcGeom0->Radius(); //3, bsmdp
  float radiusBEMC = emcGeomBEMC->Radius(); //1, bemc

  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;
  StThreeVectorD positionBEMC, momentumBEMC;

  bool okBSMDP = false;
  bool okBSMDE = false;
  bool okBEMC = false;
  okBSMDP = projPos->projTrack(&positionBSMDP, &momentumBSMDP, mcTrack, mMagneticField, radius);
  okBSMDE = projPos->projTrack(&positionBSMDE, &momentumBSMDE, mcTrack, mMagneticField, radius0);
  okBEMC = projPos->projTrack(&positionBEMC, &momentumBEMC, mcTrack, mMagneticField,radiusBEMC);
  cout << "blaa oks are " << okBSMDP << " " << okBSMDE << " " << okBEMC << endl;

  int softidp = 0;
  emcGeom->getId(positionBSMDP.pseudoRapidity(),positionBSMDP.phi(),softidp); //bsmdp plane
  float x = 0, y = 0, z = 0, phi = 0;
  emcGeom->getPhi(softidp,phi);
  int softidz = 0;
  emcGeom0->getId(positionBSMDE.pseudoRapidity(),positionBSMDE.phi(),softidz); //bsmde plane
  emcGeom0->getXYZ(softidz,x,y,z);
  cout << "blaa idz " << softidz << " z " << z << endl;
  cout << "blaa idp " << softidp << " p " << phi << endl;
  int softid = 0;
  emcGeomBEMC->getId(positionBEMC.pseudoRapidity(),positionBEMC.phi(),softid);
  emcGeomBEMC->getPhi(softid,phi);
  emcGeomBEMC->getXYZ(softid,x,y,z);
  cout << "blaa id " << softid << " z " << z << " p " << phi << endl;*/


  /*cout << "disttt id 1 " << softid << endl;
  emcGeom0->getId(positionBSMDE.pseudoRapidity(),positionBSMDE.phi(),softid);
  cout << "disttt id 2 " << softid << endl;
  emcGeom->getId(positionBSMDE.pseudoRapidity(),positionBSMDP.phi(),softid);
  cout << "disttt id 3 " << softid << endl;
  emcGeom0->getId(positionBSMDE.pseudoRapidity(),positionBSMDP.phi(),softid);
  cout << "disttt id 4 " << softid << endl;
  emcGeomBEMC->getId(positionBSMDE.pseudoRapidity(),positionBSMDP.phi(),softid);
  cout << "disttt id 5 " << softid << endl;*/

 /* StEmcCollection* emcColl = mEvent->emcCollection();
  if (!emcColl) LOG_INFO << "EmcCollection doesn't exist for this event." << endm;

  StSPtrVecEmcPoint& bEmcPoints = emcColl->barrelPoints();
  Int_t mod, eta, sub;
  emcGeom0->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub);
  for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++) 
  {
    bool associated=false;

*/
/*
  if (!okBSMDP) LOG_INFO << "Projection for the track failed." << endm;
  d[0] = z - positionBSMDE.z();
  float dphi = phi - positionBSMDP.phi();
    if(dphi>=TMath::Pi()) dphi=dphi-TMath::TwoPi();
    if(dphi<-TMath::Pi()) dphi=dphi+TMath::TwoPi();
  d[1] = dphi;*/
  //d[0] = positionBEMC.z();
  //d[1] = positionBEMC.phi();

  
  delete projPos;
  delete PosVecBSMDP;
  delete MomVecBSMDP;
  delete PosVecBSMDE;
  delete MomVecBSMDE;
  //delete geomBemc;
  return;
}

int StUpsilonEmbedMaker :: getTPADC(int id){
  int TPid;int TPADC = 0;
  StEmcDecoder* EmcDecoder = new StEmcDecoder;
  EmcDecoder->GetTriggerPatchFromTowerId(id,TPid);
  for (int softId = 1;softId<=4800;++softId){
    int TPidcur;
    EmcDecoder->GetTriggerPatchFromTowerId(softId,TPidcur);
    if (TPidcur == TPid) TPADC+=bemcADC[softId];}
  
  return TPADC;
}

Int_t StUpsilonEmbedMaker::getgRefMult(StEvent *event, StThreeVectorF vtxPos)
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

double StUpsilonEmbedMaker::getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos)
{
  if(!globalTrack) return 999;
  StDcaGeometry* trDcaGeom = globalTrack->dcaGeometry();
  if(!trDcaGeom) return 999;
  StPhysicalHelixD dcahh = trDcaGeom->helix();
  double pathlength = dcahh.pathLength(vtxPos, false);
  return (dcahh.at(pathlength)-vtxPos).mag();
}