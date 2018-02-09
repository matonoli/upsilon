// NOTE - chain needs to be declared global so for StHbtEventReader
//==========================================================================================
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class StMuDstMaker;
class StMuDebug;
class StMtdMatchMaker;
class StMtdHitMaker;

void ups_Maker(Int_t nevents=1, 
		const Char_t *fileList, 
		const Char_t *qaflag, 
		const Char_t *outroot="test.root", 
		const Char_t *outtree="tree.root",
		Bool_t debug = 0,
		Bool_t runEmbedding = 0)
{
  load();
  StChain *chain = new StChain; 
  //chain->SetDebug(0);
  StMuDebug::setLevel(-1);  // switch of some debug output
  StMuTimer timer;
  timer.start();

  Int_t dataType = -1;  // 0 - StEvent; 1 - MuDst; 2 - PicoDst
  Bool_t iMC    = 0;
  Bool_t iMB    = 0;

  ifstream infile;
  infile.open(fileList);
  string line;
  getline(infile,line);
  infile.close();

  // determine input data type
  if (line.find("event.root")!=std::string::npos)
    dataType = 0;
  else if (line.find("MuDst.root")!=std::string::npos)
    dataType = 1;
  else if (line.find("picoDst.root")!=std::string::npos)
    dataType = 2;
  else
    {
      printf("Unrecognized file type: %s\n",line.c_str());
      return;
    }

  // check if it is MC
  if (line.find("pythia")!=std::string::npos)
    iMC = 1;

  if (line.find("st_physics")!=std::string::npos)
    iMB = 1;

  // parse to extract run/year
  std::size_t found;
  string name = line;
  found = name.find("st_physics");
  if (found!=std::string::npos)
  name.erase(0,found);

  if(name.find("st_physics_adc") != std::string::npos)
    name.erase(0,15);
  else if(name.find("st_physics") != std::string::npos)
    name.erase(0,11);
  
  name.erase(8,name.length());
  Int_t run = atoi(name.c_str());
  Int_t year = run/1e6 - 1;


  if(iMC==1)    cout << "[i] Simulation" << endl;
  else          cout << "[i] Real data" << endl;
  if(dataType==0)       cout << "[i] Reading StEvent" << endl;
  else if(dataType==1)  cout << "[i] Reading MuDst"   << endl;
  else if(dataType==2)  cout << "[i] Reading PicoDst" << endl;

  cout << "[i] Year " << year << endl;

  if(dataType==0)
    {
      StIOMaker *IOMk = new StIOMaker("IO","r",Form("@%s",fileList),"bfcTree");
      IOMk->SetDebug(0);
      IOMk->SetIOMode("r");
      IOMk->SetBranch("*",0,"0");           //deactivate all branches
      IOMk->SetBranch("eventBranch",0,"r"); //activate event Branch
      if(runEmbedding)      IOMk->SetBranch("geantBranch",0,"r");
      if(runEmbedding||iMC) IOMk->SetBranch("McEventBranch",0,"r"); //activate McEvent Branch
    }
  else if(dataType==1)
    {
    }
  else if(dataType==2)
    {
    }

  St_db_Maker *dbMk;
  StMagFMaker *magfMk;
  if(runEmbedding) 
    {
      dbMk = new St_db_Maker("StarDb", "MySQL:StarDb", "$STAR/StarDb","StarDb");
      if(dataType==2) dbMk->SetDateTime((2000+year)*1e4+101,0);
      // dbMk>SetFlavor("sim","bprsCalib"); // WHAT IS THIS ??? 
      magfMk = new StMagFMaker; 
    }
  
  if(iMC)
    {
      StMcEventMaker *mcEvent = new StMcEventMaker();
    }

  if(dataType==0 && runEmbedding)
    {
      StMcEventMaker *mcEvent = new StMcEventMaker();
      StAssociationMaker *association = new StAssociationMaker();       // TPC association maker
    }

  StRefMultCorr* grefmultCorrUtil = 0x0;
  if(year==14)
    {
      grefmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr();
      grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0); //this one is only for VPDMB5, but for VPD30 and VPDZDCNoVtx, just keep them and use different StRefmultCorr maker
      grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
    }š

  StMtdJpsiUtil *mtdUtil = 0x0;
  if(runJpsiMaker || runEventMixing || runEmbedding || runTrigger || runRunQA || runMtdStudy)
    {
      mtdUtil = new StMtdJpsiUtil();
      mtdUtil->setYear(2000+year);

      if(year==13)
	{
	  /// Run13 pp
	  mtdUtil->setUsePrimVtxClosestToVpd(kFALSE);
	  mtdUtil->setRequireMtdHitForPrimVtx(kFALSE);
	  mtdUtil->setTrigger(1,0,0,0);
	  mtdUtil->setTrackType(0); // 0 - primary tracks; 1 - global tracks
	  mtdUtil->setMaxVtxZ(1e4);
	  mtdUtil->setMaxDiffz(1e4);
	  mtdUtil->setRequireVtxRanking(false);
	  mtdUtil->setMaxDca(3);
	  mtdUtil->setTrackPtLimits(0.2,1e4);
	  mtdUtil->setTrackEtaLimits(-1,1);
	  mtdUtil->setMinNHitsFit(15);
	  mtdUtil->setMinNHitsDedx(10);
	  mtdUtil->setMinFitHitsFaction(0.52);
	  mtdUtil->setMtdHitTrigger(1);
	  mtdUtil->setMuonPt(1.0,1e4);
	  mtdUtil->setNsigmaPiCut(-1,3);
	  mtdUtil->setApplyPtDepCutDz(1);
	  mtdUtil->setApplyPtDepCutDy(1);
	  mtdUtil->setApplyDtofCut(0);
	  mtdUtil->setSigmaDz(2,2);
	  mtdUtil->setSigmaDy(2,2);
	  mtdUtil->setMuonDeltaTof(-1e4,1e4);
	  /*
	 /// Run13 embedding
	 mtdUtil->setMuonDeltaY(-1e4,1e4);
	 mtdUtil->setMuonDeltaTof(-1e4,1e4);
	  */

	  /*
	  // Run13 MB
	  if(iMB) mtdUtil->setTrigger(0,0,0,1);
	  */
	}
      else if(year==14)
	{
	  // Run14 AuAu
	  mtdUtil->setUsePrimVtxClosestToVpd(kTRUE);
	  mtdUtil->setRequireMtdHitForPrimVtx(kFALSE);
	  mtdUtil->setRequireVtxRanking(kFALSE);
	  mtdUtil->setTrigger(1,0,0,0);
	  mtdUtil->setTrackType(0); // 0 - primary tracks; 1 - global tracks
	  mtdUtil->setMaxVtxZ(100);
	  mtdUtil->setMaxVtxR(2);
	  mtdUtil->setMaxDiffz(3);
	  mtdUtil->setMaxDca(1);
	  mtdUtil->setNsigmaPiCut(-1,3);
	  mtdUtil->setMinNHitsFit(15);
	  mtdUtil->setMinNHitsDedx(10);
	  mtdUtil->setMinFitHitsFaction(0.52);
	  mtdUtil->setApplyPtDepCutDz(1);
	  mtdUtil->setApplyPtDepCutDy(1);
	  mtdUtil->setApplyPtDepCutDtof(0);
	  mtdUtil->setSigmaDz(2,2.5);
	  mtdUtil->setSigmaDy(2,2.5);
	  mtdUtil->setMuonDeltaTof(-1e4,0.75);
	  mtdUtil->setMtdHitTrigger(1);
	  mtdUtil->setMuonPt(1.0,1e4);
	  mtdUtil->setRefMultCorr(grefmultCorrUtil);
	  mtdUtil->setVzCorrMethod(1);
	  mtdUtil->setVzCorrFileName("/star/u/marr/mtd/analysis/InputFile/vz_corr.root");
      
	  /*
	  // run MB trigger electronics efficiency
	  mtdUtil->setTrigger(0,0,0,1);
	  mtdUtil->setMtdHitTrigger(0);
	  */

	  // run MB vtx eff
	  //mtdUtil->setTrigger(0,0,0,1);
	}
      else if(year==15)
	{
	  // Run15 pp
	  mtdUtil->setUsePrimVtxClosestToVpd(kTRUE);
	  mtdUtil->setRequireMtdHitForPrimVtx(kFALSE);
	  mtdUtil->setRequireVtxRanking(kTRUE);
	  mtdUtil->setTrigger(1,0,0,0);
	  mtdUtil->setTrackType(0); // 0 - primary tracks; 1 - global tracks
	  mtdUtil->setMaxVtxZ(100);
	  mtdUtil->setMaxVtxR(100);
	  mtdUtil->setMaxDiffz(6);
	  mtdUtil->setMaxDca(3);
	  mtdUtil->setNsigmaPiCut(-2,3);
	  mtdUtil->setMinNHitsFit(15);
	  mtdUtil->setMinNHitsDedx(10);
	  mtdUtil->setMinFitHitsFaction(0.52);
	  mtdUtil->setApplyPtDepCutDz(1);
	  mtdUtil->setApplyPtDepCutDy(1);
	  mtdUtil->setApplyPtDepCutDtof(0);
	  mtdUtil->setSigmaDz(3,3.5);
	  mtdUtil->setSigmaDy(3,3.5);
	  mtdUtil->setMuonDeltaTof(-1.0,1.0);
	  mtdUtil->setMtdHitTrigger(1);
	  mtdUtil->setMuonPt(1.0,1e4);
	}
      else if(year==16)
	{
	  // Run16 AuAu
	  mtdUtil->setUsePrimVtxClosestToVpd(kTRUE);
	  mtdUtil->setRequireMtdHitForPrimVtx(kFALSE);
	  mtdUtil->setRequireVtxRanking(kFALSE);
	  mtdUtil->setTrigger(1,0,0,0);
	  mtdUtil->setTrackType(0); // 0 - primary tracks; 1 - global tracks
	  mtdUtil->setMaxVtxZ(100);
	  mtdUtil->setMaxVtxR(2);
	  mtdUtil->setMaxDiffz(3);
	  mtdUtil->setMaxDca(1);
	  mtdUtil->setNsigmaPiCut(-1,3);
	  mtdUtil->setMinNHitsFit(15);
	  mtdUtil->setMinNHitsDedx(10);
	  mtdUtil->setApplyPtDepCutDz(1);
	  mtdUtil->setApplyPtDepCutDy(1);
	  mtdUtil->setApplyPtDepCutDtof(0);
	  mtdUtil->setSigmaDz(2,2.5);
	  mtdUtil->setSigmaDy(2,2.5);
	  mtdUtil->setMuonDeltaTof(-1e4,1.0);
	  mtdUtil->setMtdHitTrigger(1);
	  mtdUtil->setMuonPt(1.0,1e4);
	  mtdUtil->setRefMultCorr(grefmultCorrUtil);

	  // run MB vtx eff
	  //mtdUtil->setTrigger(0,0,0,1);
	}      

      if(runEmbedding)
	{
	  mtdUtil->setMtdHitTrigger(0);
	}
      
      mtdUtil->Init();

      if(iMC)
	{
	  mtdUtil->setMCflag(1);
	  mtdUtil->setMinNHitsFit(20);
	  mtdUtil->setMinFitHitsFaction(0.52);
	}
      mtdUtil->printConfig();
    }

  if(runTrigger)
    {
      StMtdTrigger *trigger = new StMtdTrigger();
      trigger->setRunYear(year+2000);
      trigger->setDataType(dataType);
      trigger->setMtdUtil(mtdUtil);
    }

  if(runRunQA)
    {
      StMtdRunQa *runQA = new StMtdRunQa();
      runQA->setMtdUtil(mtdUtil);
      runQA->setRunRunQA(true);
      runQA->setRunDetailQA(true);
      runQA->setRunYear(year+2000);
      runQA->setDataType(dataType);
    }

  if(runJpsiMaker)
    {
      if(!iMC && !magfMk && !iMB)    magfMk = new StMagFMaker; 

      StMtdJpsiAnaMaker *ana = new StMtdJpsiAnaMaker();
      ana->setYear(year+2000);
      ana->setDataType(dataType);
      ana->setMtdUtil(mtdUtil);
      ana->setCatchEventsWithGhostVtx(0);
      ana->setCatchJpsiEvent(0);
      ana->setRunMbVtxEff(0);
      ana->setRunMultQA(0);
      ana->setRunVpdQA(0);
      ana->setFindVtxQA(0);
      ana->setRunLumiQA(0);
      ana->setApplyDimuonCut(0);
      ana->setVertexQA(1);
      ana->setRunTrackQA(0);
      ana->setRunProjectionQA(0);
      ana->setStudyProjection(0);
      ana->setRunMatchQA(0);
      ana->setRunTagAndProbe(1);
      ana->setRunJpsiMuonAna(0);
      ana->setRunJpsiAna(1);
      ana->setRunJpsiCutStudy(0);
      ana->setRunPartialTrkEff(0);
      ana->setSaveWeighted(1);
      TString outFileName = Form("jpsi.%s",outroot);
      ana->setOutputFileName(outFileName.Data());
    }

  if(runMtdStudy)
    {
      StMtdStudy *study = new StMtdStudy();
      study->setYear(year+2000);
      study->setDataType(dataType);
      study->setMtdUtil(mtdUtil);
      study->setRunTpcDistortion(0);
      study->setRunSplitTrigUnit(0);
      study->setRunV0(0);
      study->setRunPidQA(0);
      study->setRunRefMult(0);
      study->setRunMtdTrigEff(1);
    }

  
  if(runSys)
    {

      // TPC tracking quality
      // const int nSys = 12;
      // const TString sysname[nSys] = {"default",
      // 				   "dcaUp","dcaDown","NHitsUp","NDedxUp",
      // 				   "dcaUp_NHitsUp", "dcaUp_NDedxUp", "dcaUp_NHitsUp_NDedxUp",
      // 				   "dcaDown_NHitsUp", "dcaDown_NDedxUp", "dcaDown_NHitsUp_NDedxUp",
      // 				   "NHitsUp_NDedxUp"};

      const int nSys = 27;
      const TString sysname[nSys] = {"default",
				     "dzUp","dzDown","dyUp","dyDown","nSigPiUp","nSigPiDown",
				     "dzUp_dyUp","dzUp_dyDown","dzUp_nSigPiUp","dzUp_nSigPiDown","dzUp_dyUp_nSigPiUp","dzUp_dyUp_nSigPiDown","dzUp_dyDown_nSigPiUp","dzUp_dyDown_nSigPiDown",
				     "dzDown_dyUp","dzDown_dyDown","dzDown_nSigPiUp","dzDown_nSigPiDown","dzDown_dyUp_nSigPiUp","dzDown_dyUp_nSigPiDown","dzDown_dyDown_nSigPiUp","dzDown_dyDown_nSigPiDown",
				     "dyUp_nSigPiUp","dyUp_nSigPiDown",
				     "dyDown_nSigPiUp","dyDown_nSigPiDown"};



      for(int isys=0; isys<nSys; isys++)
	{
	  printf("+++++ %s +++++\n",sysname[isys].Data());
	  StMtdJpsiUtil *mtdUtilSys = new StMtdJpsiUtil(*mtdUtil);
	  mtdUtilSys->SetName(Form("MtdUtil_%s",sysname[isys].Data()));
	  if(sysname[isys].Contains("dcaUp"))   mtdUtilSys->setMaxDca(1.5);
	  if(sysname[isys].Contains("dcaDown")) mtdUtilSys->setMaxDca(0.8);
	  if(sysname[isys].Contains("NHitsUp")) mtdUtilSys->setMinNHitsFit(20);
	  if(sysname[isys].Contains("NDedxUp")) mtdUtilSys->setMinNHitsDedx(15);

	  if(sysname[isys].Contains("dzUp"))    mtdUtilSys->setSigmaDz(2.25,2.75);
	  if(sysname[isys].Contains("dzDown"))  mtdUtilSys->setSigmaDz(1.75,2.25);
	  if(sysname[isys].Contains("dyUp"))    mtdUtilSys->setSigmaDy(2.25,2.75);
	  if(sysname[isys].Contains("dyDown"))  mtdUtilSys->setSigmaDy(1.75,2.25);
	  if(sysname[isys].Contains("nSigPiUp"))    mtdUtilSys->setNsigmaPiCut(-1.25, 3.25);
	  if(sysname[isys].Contains("nSigPiDown"))  mtdUtilSys->setNsigmaPiCut(-0.75, 2.75);

	  mtdUtilSys->Init();
	  mtdUtilSys->printConfig();

	  StMtdJpsiAnaMaker *ana = new StMtdJpsiAnaMaker(Form("Sys_JpsiAna_%s",sysname[isys].Data()));
	  ana->setDataType(dataType);
	  ana->setMtdUtil(mtdUtilSys);
	  ana->setRunJpsiAna(1);
	  ana->setRunSystematics(1);
	  TString outFileName = Form("Sys.%s.jpsi.%s",sysname[isys].Data(),outroot);
	  cout << outFileName.Data() << endl;
	  ana->setOutputFileName(outFileName.Data());

	}
    }

  if(runEventMixing)
    {
      StMtdEventMixing *mix = new StMtdEventMixing();
      mix->setMtdUtil(mtdUtil);
      mix->setMixType(0);
      // mix->setNumberOfBufferEvents(1);
      // mix->setCentralityBins(1);
      // mix->setEventPlaneBins(1);
      // mix->setVertexZBins(1);
      mix->setNumberOfBufferEvents(50);
      mix->setCentralityBins(16);
      mix->setEventPlaneBins(12);
      mix->setVertexZBins(20);
      mix->setCentralityRange(0,16);
      mix->setVertexZRange(-100,100); 
      mix->setEventPlaneRange(-0.5*TMath::Pi(),0.5*TMath::Pi());
      TString outFileName1 = Form("jpsi.mix.%s",outroot);
      mix->setOutputFileName(outFileName1.Data());
    }

  if(runEmbedding)
    {
      StMtdEmbedding *embed = new StMtdEmbedding();
      embed->setMtdUtil(mtdUtil);
      //embed->setRunEmbedHadron(1);
    }

  if(runEmulateHLT)
    {
      StMtdEmulateHLT *hlt = new StMtdEmulateHLT();
      hlt->setMtdUtil(mtdUtil);
    }
  
  chain->Init();

  if(nevents<0)
    {
      if(dataType==1)
	nevents = mudstMaker->chain()->GetEntries();
      else if(dataType==2)
	nevents = picoMaker->chain()->GetEntries();
      else
	nevents = 1e6;
    }
  cout<<"total entries: "<<nevents<<endl;

  int istat=0,iev=1;
  EventLoop: if (iev<=nevents && istat!=2) {
    if(iev%1000==0) cout << "Working on eventNumber " << iev << endl;
    chain->Clear();
    istat = chain->Make(iev); // This should call the Make() method in ALL makers
    if (istat == 2) { cout << "Last  Event Processed. Status = " << istat << endl; }
    if (istat == 3) { cout << "Error Event Processed. Status = " << istat << endl; }
    iev++; goto EventLoop;
	  
  } // Event Loop
	
  chain->Finish();

  if(runQAMaker)
    {
      TFile *fout = new TFile(outroot,"RECREATE");
      qaMaker->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  if(runTrigger)
    {
      TFile *fout = new TFile(Form("trigger.%s",outroot),"RECREATE");
      trigger->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  if(runRunQA)
    {
      TFile *fout = new TFile(Form("jpsi.qa.%s",outroot),"RECREATE");
      runQA->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  if(runEmbedding)
    {
      TFile *fout = new TFile(Form("jpsi.embed.%s",outroot),"RECREATE");
      embed->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  if(runMtdStudy)
    {
      TFile *fout = new TFile(Form("jpsi.study.%s",outroot),"RECREATE");
      study->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  if(runEmulateHLT)
    {
      TFile *fout = new TFile(Form("jpsi.hlt.%s",outroot),"RECREATE");
      hlt->GetHistList()->Write();
      fout->Close();
      delete fout;
    }

  printf("========================================\n");
  printf("Input file list: %s\n",fileList);
  printf("Process %d events\n",iev-1);
  printf("Run vpd calib maker: %d\n",runVpdCaliMaker);
  printf("Run hit maker:   %d\n",runHitMaker);
  printf("Run match maker: %d\n",runMatchMaker);
  printf("Run calib maker: %d\n",runMtdCaliMaker);
  printf("Run QA maker:   %d\n",runQAMaker);
  printf("Run trigger maker: %d\n",runTrigger);
  printf("Run jpsi qa maker: %d\n",runRunQA);
  printf("Run jpsi maker: %d\n",runJpsiMaker);
  printf("Run event mixing: %d\n",runEventMixing);
  printf("Run embedding: %d\n",runEmbedding);
  printf("Run HLT emulator: %d\n",runEmulateHLT);
  printf("Run MTD study: %d\n",runMtdStudy);
  printf("Output histo file: %s\n",outroot);
  printf("========================================\n");

  timer.stop();
  cout << "Total time = " << timer.elapsedTime() << endl;
}


void load() {

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  //gROOT->LoadMacro("/star/u/marr/mtd/analysis/StarVMC/Geometry/macros/loadStarGeometry.C");
  //loadStarGeometry("y2016");
  loadSharedLibraries();

  //gSystem->Load("StMagF");
  //gSystem->Load("StTpcDb");
  //gSystem->Load("StDbUtilities");
  gSystem->Load("StDetectorDbMaker.so");
  gSystem->Load("StTpcDb");
  gSystem->Load("StEvent");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StMagF");

  gSystem->Load("StDbUtilities");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("St_g2t.so");
  gSystem->Load("St_geant_Maker.so");
  gSystem->Load("StAssociationMaker");
  gSystem->Load("StMcAnalysisMaker");
  gSystem->Load("libgeometry_Tables");   
  //gSystem->Load("StEmcMixerMaker");
  //gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StTriggerUtilities");


  gSystem->Load("StEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");


  gSystem->Load("StDbLib");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("St_db_Maker");

  gSystem->Load("StEPCalibMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StMtdHitMaker");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StVpdCalibMaker");
  gSystem->Load("StBTofCalibMaker");
  gSystem->Load("StMtdCalibMaker");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StMtdJpsi");
}
