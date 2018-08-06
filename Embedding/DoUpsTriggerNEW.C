void DoUpsTriggerNEW(const char* filelist = "test.list", const char* outpath = ".", const char* outIndex = "00",int nevents = 20)
  
{
  // Load shared libraries
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");

  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEmcSimulatorMaker");
  gSystem->Load("StEmcMixerMaker");
  gSystem->Load("StAssociationMaker");
  //  gSystem->Load("StBbcTriggerMaker");
  gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StUpsilonTriggerMaker");
  gSystem->Load("StUpsilonEmbedMaker");
  gSystem->Load("StUpsilonTreeMaker");

  // Create chain
  StChain* chain = new StChain;
  //cout<<"Chain is made!"<<endl;

  cout<<"Input file is: "<<filelist<<endl;
  // I/O maker

  std::ifstream ifs(filelist);
  StFile *myFiles = new StFile();
  std::string file;
  while(!ifs.eof()){
    getline(ifs,file);
    //cout << file << endl;
    char *a=new char[file.size()+1];
    a[file.size()]=0;
    memcpy(a,file.c_str(),file.size());
    myFiles->AddFile(a);
    delete a;
  }
  //cout << "myFiles made!" << endl;

  
  StIOMaker* ioMaker = new StIOMaker("IO","r",myFiles);
  //ioMaker->SetFile(myFiles);
  ioMaker->SetIOMode("r");
  ioMaker->SetBranch("*",0,"0");             // Deactivate all branches
  ioMaker->SetBranch("geantBranch",0,"r");   // Activate geant branch
  ioMaker->SetBranch("eventBranch",0,"r");   // Activate event branch

  //cout<<"Passed IO maker"<<endl;

  //cout<<"About to init DBmaker."<<endl;

  // STAR database
  St_db_Maker* starDb = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","$PWD/StarDb");
  //St_db_Maker* starDb = new St_db_Maker("StarDb","MySQL:StarDb","$STAR/StarDb","$PWD/StarDb");
  //St_db_Maker* starDb = new St_db_Maker("StarDb","MySQL:StarDb");
  starDb->SetFlavor("sim","bprsCalib");

  // First ADC to energy maker
  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker;
  // this line is important to propagate all the hits into StEvent
  // so, even the pedestals are propagated. In this case
  // the second AdcToEMaker will be responsible for making the
  // cuts (after the simulated hits are embedded)
  adc->saveAllStEvent(kTRUE);

  StEmcPreMixerMaker* preMixer = new StEmcPreMixerMaker("preEmbed"); //this

  // MC event maker
  StMcEventMaker* mcEvent = new StMcEventMaker;
  mcEvent->doPrintEventInfo = true;
  mcEvent->SetDebug(1);
  mcEvent->doPrintMemoryInfo = false;

  // EMC simulator
  StEmcSimulatorMaker* emcSim = new StEmcSimulatorMaker;

  // EMC mixer maker
  // this is where the afterburmer happens.
  //StEmcMixerMaker* emb = new StEmcMixerMaker; //this
  // include the next line if you want to embedd all simuated hits
  // even the ones that do not have a hit in the real data
  //emb->setEmbedAll(kTRUE);

  // Second ADC to energy maker
  StEmcADCtoEMaker* adc2 = new StEmcADCtoEMaker("EReadEmbed"); //this
  adc2->setEmbeddingMode(kTRUE); //this

  // TPC assocation maker
  StAssociationMaker* assoc = new StAssociationMaker;
  assoc->useInTracker();

  // EMC association maker
  StEmcAssociationMaker* emcAssociation = new StEmcAssociationMaker; //this

  // BBC trigger maker
  // StBbcTriggerMaker* bbcTrig = new StBbcTriggerMaker;

  // EMC trigger maker
  StEmcTriggerMaker* emcTrig = new StEmcTriggerMaker("bemctrigger");

  // Upsilon trigger maker
  //StUpsilonTriggerMaker* upsTrig = new StUpsilonTriggerMaker;

  //***New Code from Pibero
  //cout<<"Got this far"<<endl;
  StUpsilonTreeMaker* upsTree = new StUpsilonTreeMaker("upsTree");
  // cout<<"New treemakers"<<endl;
  //upsTree->makeBemcStatusTable("bemcStatus.txt");
  //cout<<"Right before Ups Embedmaker"<<endl;
  // Upsilon embedding maker
  StUpsilonEmbedMaker* upsEmbed = new StUpsilonEmbedMaker;
  //cout<<"Right after Ups EmbedMaker"<<endl;

  int offset=4;
      
  TString IndexCreator(file);

  /*if(IndexCreator.Contains("adc")){
    cout << "found adc" << endl;
    offest=4; //correct for offset from _adc_ files
    }*/

  cout<<"file is "<<IndexCreator.Data()<<endl;
  //IndexCreator.Remove(0,IndexCreator.Index("st_physics")+24+offset);
  //cout << "after remove:" << IndexCreator.Data()<<endl;
  //IndexCreator.Remove(IndexCreator.Index("."));
  //IndexCreator.Remove(IndexCreator.Index("/"));
  //cout<<"after metric indexcreator is "<<IndexCreator.Data()<<endl;
  upsEmbed->SetIndex(outIndex);
 

  //TString RunNumCreator(file);
  //RunNumCreator.Remove(0,RunNumCreator.Index("st_physics")+11+offset);
  //RunNumCreator.Remove(RunNumCreator.Index("_"));
  //upsEmbed->SetRunNumber(RunNumCreator.Data());
  //cout<<"Runnumber is "<<RunNumCreator.Data()<<endl;
  
  TString TSoutfile(outpath);
  //TString TSoutfile("/home/akesich/CurrentUpsTriggerEfficiency/ups_");
  //TString TSoutfile("/eliza13/starhf/rjreed/ups1S_");
  //TSoutfile.Prepend(outpath);
  //TSoutfile.Append(RunNumCreator);
  //TSoutfile.Append("_");
  //TSoutfile.Append(IndexCreator);
  TSoutfile.Append("/ups_");
  TSoutfile.Append(outIndex);
  TSoutfile.Append(".root");
  outpath = TSoutfile.Data();
  cout<<"outpath = "<<outpath<<endl;
  //cout<<"b4 chain init"<<endl;

  chain->lsMakers(chain);

  // Initialize chain
  chain->Init();
  cout<<"Chain init"<<endl;

  // Open output file
  TFile* ofile = TFile::Open(outpath, "recreate");
  assert(ofile);
  cout << "Opened " << outpath << endl;

  // Event loop
  chain->EventLoop(nevents);

  // Save histograms
  ofile->cd();
  int nbytes = upsEmbed->GetTreeList()->Write();
  cout << "Wrote " << nbytes << " bytes" << endl;

  // Close output file
  ofile->Close();
  cout << "Closed " << outpath << endl;
}
