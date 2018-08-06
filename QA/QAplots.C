void QAplots(){

	TString dir("0506qa");

	TCanvas* c1 = new TCanvas("c1","",750,550);
	c1->SetGrid();

	TString name;

	Double_t x1[2], x1b[2], y[2];
	//
	if(1){
	gStyle->SetOptStat(0);
	c1->SetLogz();
	hEventVzvsVzvpdO->SetTitle(";v_{z}^{TPC} [cm];v_{z}^{VPD} [cm]");
	hEventVzvsVzvpdO->Draw("colz");
	TF1* f1 = new TF1("f1","[0]*x+[1]",-100,100);
	TF1* f2 = new TF1("f2","[0]*x+[1]",-100,100);
	f1->SetParameters(1,4);
	f1->Draw("same");
	f2->SetParameters(1,-4);
	f2->Draw("same");

	name = dir;
	name += "/hVzvVz.png";
	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogy();
	TH1D* h1 = hEventVzvsVzvpdO->ProjectionX();
	h1->GetYaxis()->SetTitle("#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();
   	x1[0]=30; x1[1]=x1[0]; x1b[0]=-30; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hVzO.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogy(1);
	TH1D* h1 = hCorrgrefMult->ProjectionY();
	h1->SetTitle(";corrected grefMult; #bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw("colz");

	name = dir;
   	name += "/hMult.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogy(0);
	c1->SetLogz(1);
	hCorrgrefMult->SetTitle(";corrected grefMult; uncorrected grefMult");
	hCorrgrefMult->Draw("colz");

	name = dir;
   	name += "/hMultvMult.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogy(0);
	TH1D* h1 = hCorrCentrality2->ProjectionY();
	h1->SetTitle(";centrality class;#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw("colz");

	name = dir;
   	name += "/hCent.png";
   	c1->SaveAs(name);

	}

	if(1){
	c1->SetLogz(1);
	hTrigAdcId->SetTitle(";ADC value;#bf{BEMC tower ID}");
	hTrigAdcId->Draw("colz");
	name = dir;
   	name += "/hTow.png";
   	c1->SaveAs(name);
	}

	if(1){
	hTrackdEdxvsp->GetXaxis()->SetRangeUser(0.2,10);
	c1->SetLogx(1);
	c1->SetLogz(1);
	hTrackdEdxvsp->SetTitle(";p [GeV]; |dE/dx|");
	hTrackdEdxvsp->Draw("colz");
	name = dir;
   	name += "/hDedx.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(0);
	TH1D* h1 = hTrackEoverpvsp->ProjectionY();
	h1->SetTitle(";E/p;#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=1.5; x1[1]=x1[0]; x1b[0]=0.75; x1b[1]=x1b[0];
   	y[0]=TMath::Power(gPad->GetUymin(),1); y[1]=TMath::Power(gPad->GetUymax(),1);
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hEop.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTrackEtaPhiPtP->ProjectionY();
	h1->SetTitle(";#eta;#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=1.0; x1[1]=x1[0]; x1b[0]=-1.0; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hEta.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTrackEtaPhiPtP->ProjectionZ();
	h1->SetTitle(";#varphi;#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=1.0; x1[1]=x1[0]; x1b[0]=-1.0; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	//g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	//g1b->Draw("SAME");

   	name = dir;
   	name += "/hPhi.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTrackpMomgMom->ProjectionX();
	h1->SetTitle(";p [GeV];#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=4.5; x1[1]=x1[0]; x1b[0]=3.5; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(-202);
   	g1->SetLineColor(kBlack);
   	g1->SetFillColor(kBlack);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hP.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTracknSigmaElectronvsp->ProjectionY();
	h1->SetTitle(";n#sigma_{e};#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=3.0; x1[1]=x1[0]; x1b[0]=-1.5; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hnS.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTrackzDistphiDist->ProjectionX();
	h1->SetTitle(";zDist [cm];#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Rebin();
	h1->Draw();
	c1->Update();

	x1[0]=5.0; x1[1]=x1[0]; x1b[0]=-5.0; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hZdist.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hTrackzDistphiDist->ProjectionY();
	h1->SetTitle(";phiDist [cm];#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Rebin();
	h1->Draw();
	c1->Update();

	x1[0]=0.050; x1[1]=x1[0]; x1b[0]=-0.050; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hPhidist.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hMotherPtEtaPhi->ProjectionX();
	h1->SetTitle(";pair p_{T} [GeV];#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Rebin();
	h1->Draw();
	c1->Update();

	x1[0]=1.0; x1[1]=x1[0]; x1b[0]=10.0; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	//g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hPairPt.png";
   	c1->SaveAs(name);
	}

	if(1){
	c1->SetLogx(0);
	c1->SetLogy(1);
	TH1D* h1 = hMotherYvEta->ProjectionX();
	h1->SetTitle(";y;#bf{counts}");
	h1->SetLineColor(kAzure);
	h1->SetFillColor(kBlue-10);
	h1->Draw();
	c1->Update();

	x1[0]=0.5; x1[1]=x1[0]; x1b[0]=-0.5; x1b[1]=x1b[0];
   	y[0]=TMath::Power(10,gPad->GetUymin()); y[1]=TMath::Power(10,gPad->GetUymax());
   	cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(202);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-202);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");

   	name = dir;
   	name += "/hPairY.png";
   	c1->SaveAs(name);
	}

  if(1){
  c1->SetLogy(0);
  hTrackR->SetTitle(";TPC-BEMC matching distance R; counts");
  hTrackR->SetLineColor(kAzure);
  hTrackR->SetFillColor(kBlue-10);
  hTrackR->GetXaxis()->SetRangeUser(0.0,0.1);
  hTrackR->Draw("");
  c1->Update();

  x1[0]=0.025; x1[1]=x1[0]; 
    y[0]=gPad->GetUymin(); y[1]=gPad->GetUymax();
    cout << "y0 " << y[0] << " y1 " << y[1] <<endl;
   
    TGraph* g1 = new TGraph(2, x1, y);
    g1->SetLineWidth(202);
    g1->SetLineColor(kRed);
    g1->SetFillColor(kRed);
    g1->SetFillStyle(3003);
    g1->Draw("SAME");

  name = dir;
    name += "/R.png";
    c1->SaveAs(name);

  }
}
/*
  KEY: TH1F     hEventSelection;1
  KEY: TH1F     hEventCuts;1
  KEY: TH1F     hEventTrigger;1
  KEY: TH2F     hEventVzvsVzvpdO;1      hEventVzvsVzvpdO
  KEY: TH1F     hEventdVzO;1    hEventdVzO
  KEY: TH1F     hEventVrO;1     hEventVrO
  KEY: TH1F     hEventnEmcPidsO;1       hEventnEmcPidsO
  KEY: TH1F     hEventrunId;1   hEventrunId
  KEY: TH1F     hEventeventId;1 hEventeventId
  KEY: TH1F     hEventfillId;1  hEventfillId
  KEY: TH1F     hEventrefMult;1 hEventrefMult
  KEY: TH2F     hEventVzvsVzvpd;1       hEventVzvsVzvpd
  KEY: TH1F     hEventdVz;1     hEventdVz
  KEY: TH1F     hEventVr;1      hEventVr
  KEY: TH2F     hEventnTrigTowers;1     hEventnTrigTowers
  KEY: TH2F     hEventnBtowEmc;1        hEventnBtowEmc
  KEY: TH2F     hEventnEmcTracks;1      hEventnEmcTracks
  KEY: TH2F     hEventnBtowTracks;1     hEventnBtowTracks
  KEY: TH1F     hEventAllTracks;1       hEventAllTracks
  KEY: TH1F     hEventTriggerTracks;1
  KEY: TH2F     hCorrgrefMult;1
  KEY: TH2F     hCorrgrefMultWT;1
  KEY: TH2F     hCorrCentrality;1
  KEY: TH2F     hCorrCentrality2;1
  KEY: TH1F     hCorrWeight;1
  KEY: TH1F     hTrigFlag;1     hTrigFlag
  KEY: TH2F     hTrigAdcId;1    hTrigAdcId
  KEY: TH2F     hBtowAdcId;1    hBtowAdcId
  KEY: TH2F     hEmcAdcId;1     hEmcAdcId
  KEY: TH1F     hEmcE;1 hEmcE
  KEY: TH1F     hEmcE0vE;1      hEmcE0vE
  KEY: TH2F     hEmczDistphiDist;1      hEmczDistphiDist
  KEY: TH1F     hEmcTrackIndex;1        hEmcTrackIndex
  KEY: TH1F     hTracksSelection;1      hTracksSelection
  KEY: TH2F     hTrackEoverpvsp;1       hTrackEoverpvsp
  KEY: TH2F     hTrackdEdxvsp;1 hTrackdEdxvsp
  KEY: TH2F     hTracknSigmaElectronvsp;1       hTracknSigmaElectronvsp
  KEY: TH2F     hTracknSigmaPionvsp;1   hTracknSigmaPionvsp
  KEY: TH3F     hTrackEtaPhiPtG;1       hTrackEtaPhiPtG
  KEY: TH3F     hTrackEtaPhiPtP;1       hTrackEtaPhiPtP
  KEY: TH1F     hTrackDca;1     hTrackDca
  KEY: TH1F     hTracknHitsRatio;1      hTracknHitsRatio
  KEY: TH2F     hTrackAdc0vsTower;1     hTrackAdc0vsTower
  KEY: TH2F     hTrackzDistphiDist;1    hTrackzDistphiDist
  KEY: TH1F     hTrackE0vE;1    hTrackE0vE
  KEY: TH1F     hTracknHitsFit;1        hTracknHitsFit
  KEY: TH2F     hTrackphiDist2;1        hTrackphiDist2
  KEY: TH2F     hTrackSigmavsEoverp;1
  KEY: TH2F     hTrackSigmavsE0vE;1
  KEY: TH2F     hTrackpMomgMom;1
  KEY: TH2F     hElectronPhivdEdx;1
  KEY: TH2F     hElectronTowervp;1
  KEY: TH1F     hElectronPt;1
  KEY: TH1F     hElectronP;1
  KEY: TH1F     hPairCuts;1
  KEY: TH3F     hMotherPtEtaPhi;1
  KEY: TH2F     hMotherYvEta;1
  KEY: TH3F     hUpsPtYPhi;1
  KEY: TH3F     hIMpp;1
  KEY: TH3F     hIMmm;1
  KEY: TH3F     hIMpm;1
  KEY: TH1F     hIMun;1
  KEY: TH1F     hIMli;1
  KEY: TH1F     hIMmix;1
*/