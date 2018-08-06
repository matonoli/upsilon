#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void calcRaa() {

  TFile *fin1 = new TFile("combine_11_14_16.root","READ");

  TGraphAsymmErrors* comb1Ssys = (TGraphAsymmErrors*)fin1->Get("gr_1S_Sys");
  TGraphAsymmErrors* comb2S3Ssys = (TGraphAsymmErrors*)fin1->Get("gr_2S3S_Sys");
  TGraphErrors* comb1Sstat = (TGraphErrors*)fin1->Get("gr_1S_Stat");
  TGraphErrors* comb2S3Sstat = (TGraphErrors*)fin1->Get("gr_2S3S_Stat");

	// constants for raa calculation
	float lumi = 5201.188E3;		// mb-1
	float nEvQuoted = 149.417;
	float nEvSeen = 114.209;//e+08,115.364;//e+08;
	float sigmaPP = 42;
	float sigmaAA = 6000;
	float dy = 1.0;

	//float dsigPP = 72.7E-9;//81.0E-9;//72.1E-9;
	//float dsigPPerr = 11E-9;

  //anthony
  //float dsigPP = 64.4E-9;//81.0E-9;//72.1E-9;
  //float dsigPPerr = 15.8E-9;

  //zaochen
  float dsigPP = 81.0E-9;//81.0E-9;//72.1E-9;
  float dsigPPerr = 9.2E-9;

	float dsigPP1S = 0.704*dsigPP;
	float dsigPP23S = 0.296*dsigPP;

	// GLAUBER-----------------
	float nbin09 = 301.44186; // +-        23.16514
	float nbin29 = 392.28700; // +-        27.62289
	float nbin05 = (222.68354+124.779+64.2989+30.3083+13.30831)/5;
	float nbin25 = (222.68354+124.779+64.2989)/3;
	float nbin25err = sqrt(18.03825*18.03825+24.57679*24.57679+30.38000*30.38000)/sqrt(3.);
	float nbin57 = (607.739+375.511)*0.5;
	float nbin57err = sqrt(33.49483*33.49483+30.46155*30.46155)/sqrt(2.);
	float nbin79 = 960.99459;
	float nbin[] = {nbin29, nbin25, nbin57, nbin79};
	float nbinerr[] = {27.62289, nbin25err, nbin57err, 27.58884};
	cout << "nbin05 " << nbin05 << endl;
	cout << "nbin25 " << nbin25 << endl;

	float npart09 = 126.568;
	float npart29 = 160.961;
	float npart05 = (115.34+76.499+47.489+27.338+14.269)/5;
	float npart25 = (115.34+76.499+47.489)/3;
	float npart57 = (235.473+167.04)/2;
	float npart79 = 324.662;
  double shift = 10;
	Double_t npart[] = {npart29+shift, npart25+shift, npart57+shift, npart79+shift};
  cout << "npart is " << npart29 << " " << npart25 << " " << npart57 << " " << npart79 << endl;

	// LOAD EFFICIENCIES FROM A FILE---------------
	/*double effi1S[] = {0.0451608, 0.0534371, 0.0440506, 0.0388117};//0.0274 ;
	double effi2S[] = {0.0521197, 0.0625701, 0.0511656, 0.0468241};//0.0316339;
	double effi3S[] = {1.2162*effi1S[0],1.2162*effi1S[1],1.2162*effi1S[2],1.2162*effi1S[3]};//0.0333246;
	double effi23S[] = {(689.*effi2S[0]+311.*effi3S[0])/(689.+311.),(689.*effi2S[1]+311.*effi3S[1])/(689.+311.),(689.*effi2S[2]+311.*effi3S[2])/(689.+311.),(689.*effi2S[3]+311.*effi3S[3])/(689.+311.)};//(689*effi2S+311*effi3S)/(689.+311.);
	cout << "effi 2s3s is " << effi23S[0] << endl;*/
	const double coffLM = 0.457395797;
	const double coffH  = 0.542604203;
	fstream feffiLM1S("data/effiLM1S.txt");
	fstream feffiLM2S("data/effiLM2S.txt");
	fstream feffiLM3S("data/effiLM3S.txt");
	fstream feffiH1S("data/effiH1S.txt");
	fstream feffiH2S("data/effiH2S.txt");
	fstream feffiH3S("data/effiH3S.txt");
    double effiLM1S[4], effiLM1SL[4], effiLM1ST[4];
    double effiLM2S[4], effiLM2SL[4], effiLM2ST[4];
    double effiLM3S[4], effiLM3SL[4], effiLM3ST[4];
    double effiH1S[4], effiH1SL[4], effiH1ST[4];
    double effiH2S[4], effiH2SL[4], effiH2ST[4];
    double effiH3S[4], effiH3SL[4], effiH3ST[4];
    double effi1S[4], effi1SL[4], effi1ST[4];
    double effi2S[4], effi2SL[4], effi2ST[4];
    double effi3S[4], effi3SL[4], effi3ST[4];
    double effi23S[4], effi23SL[4], effi23ST[4];
	string t; int lc = 0;
	while ( !feffiLM1S.eof() && lc < 12 )	{
   		getline(feffiLM1S,t);
   		if (lc%3==0) effiLM1S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiLM1SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiLM1ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	lc = 0;
   	while ( !feffiLM2S.eof() && lc < 12 )	{
   		getline(feffiLM2S,t);
   		if (lc%3==0) effiLM2S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiLM2SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiLM2ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	lc = 0;
   	while ( !feffiLM3S.eof() && lc < 12 )	{
   		getline(feffiLM3S,t);
   		if (lc%3==0) effiLM3S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiLM3SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiLM3ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	lc = 0;
	while ( !feffiH1S.eof() && lc < 12 )	{
   		getline(feffiH1S,t);
   		if (lc%3==0) effiH1S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiH1SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiH1ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	lc = 0;
   	while ( !feffiH2S.eof() && lc < 12 )	{
   		getline(feffiH2S,t);
   		if (lc%3==0) effiH2S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiH2SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiH2ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	lc = 0;
   	while ( !feffiH3S.eof() && lc < 12 )	{
   		getline(feffiH3S,t);
   		if (lc%3==0) effiH3S[lc/3] = atof(t.c_str());
   		if (lc%3==1) effiH3SL[lc/3] = atof(t.c_str());
   		if (lc%3==2) effiH3ST[lc/3] = atof(t.c_str());
   		lc++;	}
   	cout << "blaa " << effiLM1S[0] << endl;
   	for (int i = 0; i < 4; ++i)
   	{
   		effi1S[i] 	= coffLM*effiLM1S[i] + coffH*effiH1S[i];
   		effi1SL[i] 	= coffLM*effiLM1SL[i] + coffH*effiH1SL[i];
   		effi1ST[i] 	= coffLM*effiLM1ST[i] + coffH*effiH1ST[i];
   		effi2S[i] 	= coffLM*effiLM2S[i] + coffH*effiH2S[i];
   		effi2SL[i] 	= coffLM*effiLM2SL[i] + coffH*effiH2SL[i];
   		effi2ST[i] 	= coffLM*effiLM2ST[i] + coffH*effiH2ST[i];
   		effi3S[i] 	= coffLM*effiLM3S[i] + coffH*effiH3S[i];
   		effi3SL[i] 	= coffLM*effiLM3SL[i] + coffH*effiH3SL[i];
   		effi3ST[i] 	= coffLM*effiLM3ST[i] + coffH*effiH3ST[i];
   		effi23S[i] 	= 0.689*effi2S[i]+0.311*effi3S[i];
   		effi23SL[i] = 0.689*effi2SL[i]+0.311*effi3SL[i];
   		effi23ST[i] = 0.689*effi2ST[i]+0.311*effi3ST[i];
   	}

   	printf("Working with following efficiencies:\n");
   	for (int iter = 0; iter < 4; ++iter)
   	{
   		printf("cent %i | eff_Y1S : %.4f + %.4f - %.4f ,  eff_Y2S : %.4f + %.4f - %.4f, eff_Y3S : %.4f + %.4f - %.4f, eff_Y23S : %.4f + %.4f - %.4f \n",
   			iter,effi1S[iter],effi1SL[iter],effi1ST[iter],
   			effi2S[iter],effi2SL[iter],effi2ST[iter],
   			effi3S[iter],effi3SL[iter],effi3ST[iter],
   			effi23S[iter],effi23SL[iter],effi23ST[iter]);
   	}

    //zao crosscheck
    //effi1S[0] =  0.0326725;
    //lumi = 4.69979e+09/6000.;
    //nEvSeen = 1;
    //nEvQuoted = 1;


	// LOAD YIELDS FROM A FILE----------------
    fstream fyields("data/yields.txt");
    double yield1S[4];
    double err1S[4];
    double yield23S[4];
    double err23S[4];
	lc = 0;
	while ( !fyields.eof() && lc < 16 )	{
   		getline(fyields,t);
   		if (lc%4==0) yield1S[lc/4] = atof(t.c_str());
   		if (lc%4==1) err1S[lc/4] = atof(t.c_str());
   		if (lc%4==2) yield23S[lc/4] = atof(t.c_str());
   		if (lc%4==3) err23S[lc/4] = atof(t.c_str());
   		lc++;	}
   	printf("------------------------------------------- \n");
   	printf("Working with following yields:\n");
   	for (int iter = 0; iter < 4; ++iter)
   	{
   		printf("cent %i | Y1S : %.1f +- %.1f   ,  Y23S : %.1f +- %.1f \n",iter,yield1S[iter],err1S[iter],yield23S[iter],err23S[iter]);
   	}

   	double raa1S[4];
   	double raa1Sstat[4];
   	double raa1SsysL[4];
   	double raa1SsysT[4];
   	double raa1Sncol[4];

    double ny1S[4];
    double ny1SsysL[4];
    double ny1SsysT[4];
    double ny1Sstat[4];
    double ny23S[4];
    double ny23SsysL[4];
    double ny23SsysT[4];
    double ny23Sstat[4];


   	double raa23S[4];
   	double raa23Sstat[4];
   	double raa23SsysL[4];
   	double raa23SsysT[4];
   	double raa23Sncol[4];

   	double dcent[] = {0.6, 0.3, 0.2, 0.1};
   	double furtherL[4] = {sqrt(0.08*0.08 + 0.0488*0.0488), sqrt(0.08*0.08 + 0.015*0.015), sqrt(0.08*0.08 + 0.061*0.061), sqrt(0.08*0.08 + 0.082*0.082)};
   	double furtherT[4] = {sqrt(0.08*0.08 + 0.0348*0.0348), sqrt(0.08*0.08 + 0.025*0.025), sqrt(0.08*0.08 + 0.061*0.061), sqrt(0.08*0.08 + 0.036*0.036)};

   	for (int iC = 0; iC < 4; ++iC)
   	{
   		raa1S[iC] = yield1S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[iC]*sigmaAA*dcent[iC]*dsigPP1S*nbin[iC]*dy);
   		raa1Sstat[iC] = err1S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[iC]*sigmaAA*dcent[iC]*dsigPP1S*nbin[iC]*dy);

      ny1S[iC] = yield1S[iC]*nEvQuoted/(lumi*nEvSeen*effi1S[iC]*sigmaAA*dcent[iC]*dy);
      ny1Sstat[iC] = err1S[iC]*nEvQuoted/(lumi*nEvSeen*effi1S[iC]*sigmaAA*dcent[iC]*dy);

   		double yieldE = raa1S[iC]*effi1S[iC]*(1./(effi1S[iC]+effi1SL[iC]));
   		raa1SsysL[iC] = sqrt( (raa1S[iC]-yieldE)*(raa1S[iC]-yieldE) + raa1S[iC]*furtherL[iC]*raa1S[iC]*furtherL[iC]);
   		yieldE = raa1S[iC]*effi1S[iC]*(1./(effi1S[iC]-effi1ST[iC]));
   		raa1SsysT[iC] = sqrt( (raa1S[iC]-yieldE)*(raa1S[iC]-yieldE) + raa1S[iC]*furtherT[iC]*raa1S[iC]*furtherT[iC]);
   		yieldE = raa1S[iC]*nbin[iC]*(1./(nbin[iC]-nbinerr[iC]));
   		raa1Sncol[iC] = fabs(yieldE-raa1S[iC]);

   		raa23S[iC] = yield23S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[iC]*sigmaAA*dcent[iC]*dsigPP23S*nbin[iC]*dy);
   		raa23Sstat[iC] = err23S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[iC]*sigmaAA*dcent[iC]*dsigPP23S*nbin[iC]*dy);

   		yieldE = raa23S[iC]*effi23S[iC]*(1./(effi23S[iC]+effi23SL[iC]));
   		raa23SsysL[iC] = sqrt( (raa23S[iC]-yieldE)*(raa23S[iC]-yieldE) + 4*raa23S[iC]*furtherL[iC]*raa23S[iC]*furtherL[iC]);
   		yieldE = raa23S[iC]*effi23S[iC]*(1./(effi23S[iC]-effi23ST[iC]));
   		raa23SsysT[iC] = sqrt( (raa23S[iC]-yieldE)*(raa23S[iC]-yieldE) + 4*raa23S[iC]*furtherT[iC]*raa23S[iC]*furtherT[iC]);
   		yieldE = raa23S[iC]*nbin[iC]*(1./(nbin[iC]-nbinerr[iC]));
   		raa23Sncol[iC] = fabs(yieldE-raa23S[iC]);

   		raa23S[iC] = yield23S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[iC]*sigmaAA*dcent[iC]*dsigPP23S*nbin[iC]*dy);
   		raa23Sstat[iC] = err23S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[iC]*sigmaAA*dcent[iC]*dsigPP23S*nbin[iC]*dy);
   		//raa23Sstat[iC] = err23S[iC]*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[iC]*sigmaAA*dcent[iC]*dsigPP23S*nbin[iC]*dy);


   		printf("raa 1s is : %f +- %f \n", raa1S[iC], raa1Sstat[iC]);
   		printf("raa 23s is : %f +- %f \n", raa23S[iC], raa23Sstat[iC]);
   	}

   	double normX[] = {398.};
   	double normY[] = {1.0 };
   	double normEL[] = {0.704*dsigPPerr/dsigPP1S};
   	double normET[] = {0.704*dsigPPerr/dsigPP1S};

   	double normE2L[] = {0.296*dsigPPerr/dsigPP1S};
   	double normE2T[] = {0.296*dsigPPerr/dsigPP1S};


   	// DRAW EFFICIENCIES
   	{
   	double effX[] = {-0.1, 0.9, 1.9, 2.9};
   	double effX2[] = {0, 1, 2, 3};
   	double effX3[] = {0.1, 1.1, 2.1, 3.1};
	double effWid[] = {0.05, 0.05, 0.05, 0.05};
   	TGraphAsymmErrors* pEff = new TGraphAsymmErrors(4,effX,effi1S,effWid,effWid,effi1ST,effi1SL);
   	TGraphAsymmErrors* pEff2 = new TGraphAsymmErrors(4,effX2,effi2S,effWid,effWid,effi2ST,effi2SL);
   	TGraphAsymmErrors* pEff3 = new TGraphAsymmErrors(4,effX3,effi3S,effWid,effWid,effi3ST,effi3SL);
	pEff->SetTitle(";;efficiency");
	TAxis * xax = pEff->GetXaxis();
  	TAxis * yax = pEff->GetYaxis();
  xax->SetBinLabel(xax->FindBin(0.0),"0-60%");
  xax->SetBinLabel(xax->FindBin(1.0),"30-60%");
  xax->SetBinLabel(xax->FindBin(2.0),"10-30%");
  xax->SetBinLabel(xax->FindBin(3.0),"0-10%");
  xax->LabelsOption("h");
  xax->CenterTitle();
  xax->SetLabelSize(0.075);
  xax->SetTitleSize(0.055);
  yax->CenterTitle();
  yax->SetLabelSize(0.045);
  yax->SetTitleSize(0.055);
  yax->SetTitleOffset(+0.88);
  pEff->SetMinimum(0.0);
  pEff->SetMaximum(0.15);
  pEff->SetMarkerSize(0.7);
  pEff->SetMarkerStyle(21);
  pEff->SetMarkerColor(kBlack);
  pEff->SetLineColor(kPink-5);
  pEff->SetLineWidth(2);
  pEff->SetFillColor(kPink-5);
  pEff->SetFillStyle(1001);
  pEff2->SetMarkerSize(0.7);
  pEff2->SetMarkerStyle(21);
  pEff2->SetMarkerColor(kBlack);
  pEff2->SetLineColor(kSpring-5);
  pEff2->SetLineWidth(2);
  pEff2->SetFillColor(kSpring-5);
  pEff2->SetFillStyle(1001);
  pEff3->SetMarkerSize(0.7);
  pEff3->SetMarkerStyle(21);
  pEff3->SetMarkerColor(kBlack);
  pEff3->SetLineColor(kAzure+1);
  pEff3->SetLineWidth(2);
  pEff3->SetFillColor(kAzure+1);
  pEff3->SetFillStyle(1001);
  pEff->Draw("A2P");
  pEff->Draw("2P");
  pEff2->Draw("2P");
  pEff3->Draw("2P");
  TLegend* leff = new TLegend(0.40,0.64,0.81,0.85);
  leff->SetNColumns(1);
  leff->SetFillStyle(0);
  leff->SetBorderSize(0);
  leff->SetTextSize(0.042);
  leff->AddEntry((TObject*)0,"#bf{2014 Au+Au #sqrt{s_{NN}} = 200 GeV}"," ");
   	leff->AddEntry(pEff,"#Upsilon(1S) efficiency #times acceptance","plf");
   	leff->AddEntry(pEff2,"#Upsilon(2S) efficiency #times acceptance","plf");
   	leff->AddEntry(pEff3,"#Upsilon(3S) efficiency #times acceptance","plf");
   	leff->Draw();
   	}

   	
   	// DRAW RAA

   	printf("Raa is %f  +- %f (stat) + %f - %f (syst) \n", raa1S[0], raa1Sstat[0], raa1SsysL[0], raa1SsysT[0]);
   	printf("Raa 23 is %f  +- %f (stat) + %f - %f (syst) \n", raa23S[0], raa23Sstat[0], raa23SsysL[0], raa23SsysT[0]);
   	
   	new TCanvas;
   	double raaWid[] = {5., 5., 5., 5.};
   	TGraphAsymmErrors* pRaa1S = new TGraphAsymmErrors(1,npart,raa1S,raaWid,raaWid,raa1Sstat,raa1Sstat);
   	TGraphAsymmErrors* pRaa1S_2 = new TGraphAsymmErrors(3,&npart[1],&raa1S[1],raaWid,raaWid,&raa1Sstat[1],&raa1Sstat[1]);
   	TGraphAsymmErrors* pRaa1S_sys = new TGraphAsymmErrors(4,npart,raa1S,raaWid,raaWid,raa1SsysL,raa1SsysT);
   	TGraphAsymmErrors* pRaa1S_ncol = new TGraphAsymmErrors(4,npart,raa1S,raaWid,raaWid,raa1Sncol,raa1Sncol);
   	TGraphAsymmErrors* pRaa1S_norm = new TGraphAsymmErrors(1,normX,normY,raaWid,raaWid,normEL,normET);
  	pRaa1S->SetTitle(";<N_{part}>;R_{AA}");
  	xax = pRaa1S->GetXaxis();
  	yax = pRaa1S->GetYaxis();
  	//xax->CenterTitle();
  	xax->SetLabelSize(0.04);
  	xax->SetTitleSize(0.05);
  	xax->SetTitleOffset(+0.9);
  	yax->CenterTitle();
  	yax->SetLabelSize(0.04);
  	yax->SetTitleSize(0.05);
  	yax->SetTitleOffset(+0.8);
  	xax->SetLimits(0.,400.);
  	pRaa1S->SetMinimum(0.0);
  	pRaa1S->SetMaximum(1.4);
  	pRaa1S->SetMarkerSize(2.2);
  	pRaa1S->SetMarkerStyle(30);
  	pRaa1S->SetMarkerColor(2);
  	pRaa1S->SetLineColor(2);
  	pRaa1S->SetLineWidth(1);
   	pRaa1S->Draw("AP");
  	pRaa1S_ncol->SetLineColor(kGray);
  	pRaa1S_ncol->SetLineWidth(1);
  	pRaa1S_ncol->SetFillColor(kGray);
  	pRaa1S_ncol->SetFillStyle(1001);
  	pRaa1S_ncol->Draw("2");
   	pRaa1S->Draw("P");
  	pRaa1S_2->SetMarkerSize(2.2);
  	pRaa1S_2->SetMarkerStyle(29);
  	pRaa1S_2->SetMarkerColor(2);
  	pRaa1S_2->SetLineColor(2);
  	pRaa1S_2->SetLineWidth(1);
  	pRaa1S_2->Draw("P");
  	pRaa1S_sys->SetLineColor(2);
  	pRaa1S_sys->SetLineWidth(1);
  	pRaa1S_sys->SetFillColor(2);
  	pRaa1S_sys->SetFillStyle(2001);
  	pRaa1S_sys->Draw("2");
  	pRaa1S_norm->SetLineColor(2);
  	pRaa1S_norm->SetLineWidth(1);
  	pRaa1S_norm->SetFillColor(2);
  	pRaa1S_norm->SetFillStyle(3001);
  	pRaa1S_norm->Draw("2");

    comb1Ssys->Draw("P2");
    comb1Sstat->Draw("P");

  	TF1* unity = new TF1("unity","pol0",0,400.);
  	unity->SetParameter(0,1.0);
  	unity->SetLineColor(kBlack);
  	unity->SetLineStyle(2);
  	unity->Draw("same");

  	TLegend* lraa1S = new TLegend(0.13,0.69,0.41,0.88);
   	lraa1S->SetNColumns(1);
   	lraa1S->SetFillStyle(0);
   	lraa1S->SetBorderSize(0);
   	lraa1S->SetTextSize(0.032);
   	lraa1S->AddEntry((TObject*)0,"#bf{Au+Au 2014 #sqrt{s_{NN}} = 200 GeV}"," ");
   	lraa1S->AddEntry(pRaa1S_2,"#Upsilon(1S)->e^{+}e^{-}, |y|<0.5, 60-30%, 30-10%, 10-0%","p");
   	lraa1S->AddEntry(pRaa1S,"0-60%","p");
   	//lraa1S->AddEntry(pRaa1S_sys,"systematic uncertainty","f");
   	lraa1S->AddEntry(pRaa1S_ncol,"STAR N_{coll} uncertainty","f");
   	lraa1S->Draw();


   	// 23S

   	new TCanvas;
   	TGraphAsymmErrors* pRaa23S = new TGraphAsymmErrors(1,npart,raa23S,raaWid,raaWid,raa23Sstat,raa23Sstat);
   	TGraphAsymmErrors* pRaa23S_2 = new TGraphAsymmErrors(3,&npart[1],&raa23S[1],raaWid,raaWid,&raa23Sstat[1],&raa23Sstat[1]);
   	TGraphAsymmErrors* pRaa23S_sys = new TGraphAsymmErrors(4,npart,raa23S,raaWid,raaWid,raa23SsysL,raa23SsysT);
   	TGraphAsymmErrors* pRaa23S_ncol = new TGraphAsymmErrors(4,npart,raa23S,raaWid,raaWid,raa23Sncol,raa23Sncol);
   	TGraphAsymmErrors* pRaa23S_norm = new TGraphAsymmErrors(1,normX,normY,raaWid,raaWid,normEL,normET);
  	pRaa23S->SetTitle(";<N_{part}>;R_{AA}");
  	xax = pRaa23S->GetXaxis();
  	yax = pRaa23S->GetYaxis();
  	//xax->CenterTitle();
  	xax->SetLabelSize(0.04);
  	xax->SetTitleSize(0.05);
  	xax->SetTitleOffset(+0.9);
  	yax->CenterTitle();
  	yax->SetLabelSize(0.04);
  	yax->SetTitleSize(0.05);
  	yax->SetTitleOffset(+0.8);
  	xax->SetLimits(0.,400.);
  	pRaa23S->SetMinimum(0.0);
  	pRaa23S->SetMaximum(1.4);
  	pRaa23S->SetMarkerSize(2.2);
  	pRaa23S->SetMarkerStyle(30);
  	pRaa23S->SetMarkerColor(2);
  	pRaa23S->SetLineColor(2);
  	pRaa23S->SetLineWidth(1);
   	pRaa23S->Draw("AP");
  	pRaa23S_ncol->SetLineColor(kGray);
  	pRaa23S_ncol->SetLineWidth(1);
  	pRaa23S_ncol->SetFillColor(kGray);
  	pRaa23S_ncol->SetFillStyle(1001);
  	pRaa23S_ncol->Draw("2");
   	pRaa23S->Draw("P");
  	pRaa23S_2->SetMarkerSize(2.2);
  	pRaa23S_2->SetMarkerStyle(29);
  	pRaa23S_2->SetMarkerColor(2);
  	pRaa23S_2->SetLineColor(2);
  	pRaa23S_2->SetLineWidth(1);
  	pRaa23S_2->Draw("P");
  	pRaa23S_sys->SetLineColor(2);
  	pRaa23S_sys->SetLineWidth(1);
  	pRaa23S_sys->SetFillColor(2);
  	pRaa23S_sys->SetFillStyle(2001);
  	pRaa23S_sys->Draw("2");
  	pRaa23S_norm->SetLineColor(2);
  	pRaa23S_norm->SetLineWidth(1);
  	pRaa23S_norm->SetFillColor(2);
  	pRaa23S_norm->SetFillStyle(3001);
  	pRaa23S_norm->Draw("2");

  	unity->Draw("same");

  	TLegend* lraa23S = new TLegend(0.13,0.69,0.41,0.88);
   	lraa23S->SetNColumns(1);
   	lraa23S->SetFillStyle(0);
   	lraa23S->SetBorderSize(0);
   	lraa23S->SetTextSize(0.032);
   	lraa23S->AddEntry((TObject*)0,"#bf{Au+Au 2014 #sqrt{s_{NN}} = 200 GeV}"," ");
   	lraa23S->AddEntry(pRaa23S_2,"#Upsilon(2S+3S)->e^{+}e^{-}, |y|<0.5, 60-30%, 30-10%, 10-0%","p");
   	lraa23S->AddEntry(pRaa23S,"0-60%","p");
   	//lraa23S->AddEntry(pRaa23S_sys,"systematic uncertainty","f");
   	lraa23S->AddEntry(pRaa23S_ncol,"STAR N_{coll} uncertainty","f");
   	lraa23S->Draw();
  	


	//29
	/*float yield1S = 180.0;
	float yield23S = 74.4;
	float err1S = 21.3;
	float err23S = 17.8;*/
}
	#if 0
	float raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[0]*sigmaAA*0.6*dsigPP1S*nbin29*dy);
	float raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[0]*sigmaAA*0.6*dsigPP23S*nbin29*dy);
	float check = 200./(5201.188E3*(nEvSeen/nEvQuoted)*effi1S[0]);
	check = check*42./(6000*0.8);
	check = check/nbin09;
	check = check/dsigPP1S;//(0.0000000521);
	printf("raa 1s is %f \n", raa1S);
	//printf("raa 23s is %f \n", raa23S);
	//printf("check 1s is %f \n", check);

	float raapoint1S_0[1];
	float raapoint1SL_0[1];
	float raapoint1ST_0[1];

	float raapoint23S_0[1];
	float raapoint23SL_0[1];
	float raapoint23ST_0[1];

	raapoint1S_0[0] = raa1S;
	raapoint1SL_0[0] = raa1S*(yield1S+err1S)/yield1S - raa1S;
	raapoint1ST_0[0] = raa1S - raa1S*(yield1S-err1S)/yield1S;

	raapoint23S_0[0] = raa23S;
	raapoint23SL_0[0] = raa23S*(yield23S+err23S)/yield23S - raa23S;
	raapoint23ST_0[0] = raa23S - raa23S*(yield23S-err23S)/yield23S;



	int npoints = 3;
	float raapoint1S[npoints];
	float raapoint1SL[npoints];
	float raapoint1ST[npoints];
	float raapoint23S[npoints];
	float raapoint23SL[npoints];
	float raapoint23ST[npoints];
	float npart0[] = {npart25,npart57,npart79};
	float npart2[] = {npart29};

	//25
	yield1S = 37.9;
	err1S 	=  11.5;
	yield23S = 19.4;
	err23S = 9.3;

	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[1]*sigmaAA*0.3*dsigPP1S*nbin25*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[1]*sigmaAA*0.3*dsigPP23S*nbin25*dy);
	raapoint1S[0] = raa1S;
	raapoint23S[0] = raa23S;
	printf("05 raa 1s is %f \n", raa1S);
	//printf("raa 23s is %f \n", raa23S);

	raapoint1SL[0] = raa1S*(yield1S+err1S)/yield1S - raa1S;
	raapoint1ST[0] = raa1S - raa1S*(yield1S-err1S)/yield1S;

	raapoint23SL[0] = raa23S*(yield23S+err23S)/yield23S - raa23S;
	raapoint23ST[0] = raa23S - raa23S*(yield23S-err23S)/yield23S;

	/*yield1S = yield1S + err1S;
	yield23S = yield23S + err23S;
	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.5*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.5*dsigPP23S*nbin05*dy);
	//raapoint1SL[0] = raa1S - raapoint1S[0];
	//raapoint23SL[0] = raa23S  - raapoint23S[0];

	yield1S = yield1S - 2.*err1S;
	yield23S = yield23S - 2.*err23S;
	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.5*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.5*dsigPP23S*nbin05*dy);
	//raapoint1ST[0] = (-1.)*(raa1S - raapoint1S[0]);
	//raapoint23ST[0] = (-1.)*(raa23S  - raapoint23S[0]);*/
	
	//57
	yield1S = 98.5;
	err1S = 14.8;
	yield23S = 26.8;
	err23S = 12.7;

	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[2]*sigmaAA*0.2*dsigPP1S*nbin57*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[2]*sigmaAA*0.2*dsigPP23S*nbin57*dy);
	raapoint1S[1] = raa1S;
	raapoint23S[1] = raa23S;
	printf("57 raa 1s is %f \n", raa1S);
	//printf("raa 23s is %f \n", raa23S);

	raapoint1SL[1] = raa1S*(yield1S+err1S)/yield1S - raa1S;
	raapoint1ST[1] = raa1S - raa1S*(yield1S-err1S)/yield1S;


	raapoint23SL[1] = raa23S*(yield23S+err23S)/yield23S - raa23S;
	raapoint23ST[1] = raa23S - raa23S*(yield23S-err23S)/yield23S;

	/*cout << "yield " << yield1S << endl;
	//yield1S = yield1S + err1S;

	cout << "yield " << yield1S << endl;
	yield23S = yield23S + err23S;
	raa1S = (yield1S+err1S)*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.2*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.2*dsigPP23S*nbin05*dy);
	printf("57 raa 1s is %f \n", raa1S);
	//printf("raa 23s is %f \n", raa23S);
	raapoint1SL[1] = raa1S - raapoint1S[1];
	raapoint23SL[1] = raa23S  - raapoint23S[1];

	yield1S = yield1S - 2.*err1S;
	yield23S = yield23S - 2.*err23S;
	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.2*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.2*dsigPP23S*nbin05*dy);
	printf("57 raa 1s is %f \n", raa1S);
	//printf("raa 23s is %f \n", raa23S);
	raapoint1ST[1] = (-1.)*(raa1S - raapoint1S[1]);
	raapoint23ST[1] = (-1.)*(raa23S  - raapoint23S[1]);*/

	//79
	yield1S = 44.9;
	err1S = 10.3;
	yield23S = 28.1;
	err23S = 8.4;

	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S[3]*sigmaAA*0.1*dsigPP1S*nbin79*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S[3]*sigmaAA*0.1*dsigPP23S*nbin79*dy);
	raapoint1S[2] = raa1S;
	raapoint23S[2] = raa23S;
	printf("raa 1s is %f \n", raa1S);
	printf("raa 23s is %f \n", raa23S);

	raapoint1SL[2] = raa1S*(yield1S+err1S)/yield1S - raa1S;
	raapoint1ST[2] = raa1S - raa1S*(yield1S-err1S)/yield1S;



	raapoint23SL[2] = raa23S*(yield23S+err23S)/yield23S - raa23S;
	raapoint23ST[2] = raa23S - raa23S*(yield23S-err23S)/yield23S;

/*
	yield1S = yield1S + err1S;
	yield23S = yield23S + err23S;
	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.1*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.1*dsigPP23S*nbin05*dy);
	printf("raa 1s is %f \n", raa1S);
	printf("raa 23s is %f \n", raa23S);
	raapoint1SL[2] = raa1S - raapoint1S[2];
	raapoint23SL[2] = raa23S  - raapoint23S[2];

	yield1S = yield1S - 2*err1S;
	yield23S = yield23S - 2*err23S;
	raa1S = yield1S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi1S*sigmaAA*0.1*dsigPP1S*nbin05*dy);
	raa23S = yield23S*sigmaPP*nEvQuoted/(lumi*nEvSeen*effi23S*sigmaAA*0.1*dsigPP23S*nbin05*dy);
	printf("raa 1s is %f \n", raa1S);
	printf("raa 23s is %f \n", raa23S);
	raapoint1ST[2] = (-1.)*(raa1S - raapoint1S[2]);
	raapoint23ST[2] = (-1.)*(raa23S  - raapoint23S[2]);*/


	float wid[] = {5., 5., 5.};

	TCanvas* ceff = new TCanvas("ceff","ceff",800,600);
	ceff->SetGrid();
	TGraphAsymmErrors* plot = new TGraphAsymmErrors(npoints,npart0,raapoint1S,wid,wid,raapoint1SL,raapoint1ST);
	plot->SetTitle(";N_{part};R_{AA}");
	TAxis * xax = plot->GetXaxis();
  	TAxis * yax = plot->GetYaxis();
  	xax->CenterTitle();
  	xax->SetLabelSize(0.04);
  	xax->SetTitleSize(0.04);
  	yax->CenterTitle();
  	yax->SetLabelSize(0.04);
  	yax->SetTitleSize(0.05);
  	yax->SetTitleOffset(+0.9);
  	plot->SetMinimum(0.0);
  	plot->SetMaximum(1.4);
  	plot->SetMarkerSize(2.2);
  	plot->SetMarkerStyle(29);
  	plot->SetMarkerColor(2);
  	plot->SetLineColor(2);
  	plot->SetLineWidth(2);
  	plot->SetFillColor(20);
    plot->SetFillStyle(1001);
   	plot->Draw("AP");
  	xax->SetLimits(0.,400.);
   	plot->Draw("AP");
   	ceff->SaveAs("1s.png");


	TGraphAsymmErrors* plot2 = new TGraphAsymmErrors(1,npart2,raapoint1S_0,wid,wid,raapoint1SL_0,raapoint1ST_0);
	plot2->SetTitle(";N_{part};R_{AA}");
	TAxis * xax2 = plot2->GetXaxis();
  	TAxis * yax2 = plot2->GetYaxis();
  	xax2->CenterTitle();
  	xax2->SetLabelSize(0.04);
  	xax2->SetTitleSize(0.04);
  	yax2->CenterTitle();
  	yax2->SetLabelSize(0.04);
  	yax2->SetTitleSize(0.05);
  	yax2->SetTitleOffset(+0.9);
  	plot2->SetMinimum(0.0);
  	plot2->SetMaximum(1.4);
  	plot2->SetMarkerSize(2.2);
  	plot2->SetMarkerStyle(30);
  	plot2->SetMarkerColor(2);
  	plot2->SetLineColor(2);
  	plot2->SetLineWidth(2);
  	plot2->SetFillColor(20);
    plot2->SetFillStyle(1001);
  	xax2->SetLimits(0.,400.);
   	plot2->Draw("P");
   	ceff->SaveAs("1s_0.png");

   	TGraphAsymmErrors* plot3 = new TGraphAsymmErrors(npoints,npart0,raapoint23S,wid,wid,raapoint23SL,raapoint23ST);
	plot3->SetTitle(";N_{part};R_{AA}");
	TAxis * xax3 = plot3->GetXaxis();
  	TAxis * yax3 = plot3->GetYaxis();
  	xax3->CenterTitle();
  	xax3->SetLabelSize(0.04);
  	xax3->SetTitleSize(0.04);
  	yax3->CenterTitle();
  	yax3->SetLabelSize(0.04);
  	yax3->SetTitleSize(0.05);
  	yax3->SetTitleOffset(+0.9);
  	plot3->SetMinimum(0.0);
  	plot3->SetMaximum(1.4);
  	plot3->SetMarkerSize(2.2);
  	plot3->SetMarkerStyle(29);
  	plot3->SetMarkerColor(2);
  	plot3->SetLineColor(2);
  	plot3->SetLineWidth(2);
  	plot3->SetFillColor(20);
    plot3->SetFillStyle(1001);
   	plot3->Draw("AP");
  	xax3->SetLimits(0.,400.);
   	plot3->Draw("AP");
   	ceff->SaveAs("23s.png");


	TGraphAsymmErrors* plot4 = new TGraphAsymmErrors(1,npart2,raapoint23S_0,wid,wid,raapoint23SL_0,raapoint23ST_0);
	plot4->SetTitle(";N_{part};R_{AA}");
	TAxis * xax4 = plot4->GetXaxis();
  	TAxis * yax4 = plot4->GetYaxis();
  	xax4->CenterTitle();
  	xax4->SetLabelSize(0.04);
  	xax4->SetTitleSize(0.04);
  	yax4->CenterTitle();
  	yax4->SetLabelSize(0.04);
  	yax4->SetTitleSize(0.05);
  	yax4->SetTitleOffset(+0.9);
  	plot4->SetMinimum(0.0);
  	plot4->SetMaximum(1.4);
  	plot4->SetMarkerSize(2.2);
  	plot4->SetMarkerStyle(30);
  	plot4->SetMarkerColor(2);
  	plot4->SetLineColor(2);
  	plot4->SetLineWidth(2);
  	plot4->SetFillColor(20);
    plot4->SetFillStyle(1001);
  	xax4->SetLimits(0.,400.);
   	plot4->Draw("P");
   	ceff->SaveAs("23s_0.png");






}
#endif