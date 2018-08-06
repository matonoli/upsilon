#include <iostream>
#include <cmath>

void nSigEff() {

	// LOAD FILES
	TFile* fin = new TFile("hNSigEff.root","READ");
	TFile* fhadrons = new TFile("hadronMeans.root","READ");

	// OUTPUT FILE
	TFile* fout = new TFile("nSigEff.root","RECREATE");

	// SELECT BINNING
	int nbins = 12;
	//float xbins[] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.5, 10.0}; //12
	//float xbins[] = {2.4, 2.8, 3.2, 3.6, 4.0, 4.75, 5.5, 7.0, 8.5, 10.0};//10
	//3.5, 4.25, 5.0, 6.0, 8.0, 10.0};
	//float xbins[] = {2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.25,6.0,7.0,8.5,11}; //13
	//float xbins[] = {2.0,2.5,3.0,3.5,4.0,4.5,5.25,6.0,7.0,8.5,11.0};   //11
	float xbins[] = {2.0,2.5,3.0,3.5,4.0,4.5,5.25,6.0,6.75,8.0,9.5,11.5};   //12
	//float xbins[] = {2.0,2.3,2.8,3.2,3.6,4.0,4.5,5.25,6.0,6.75,8.0,9.5,11.5};   //13
	//float xbins[] = {};
	//float pimu[9] = {-4.5, -4, -3.5, -3, -2.7, -2.5, -2.3, -2.1, -2.0};
	//float kpmu[9] = {}
	
	// OTHER VARIABLES
	const float left = -1.5;
	const float right = 3;
	int nRepeats = 3; // refit how many times?
	int KPflag = 1; // show: 0 - kaons and protons, 1 - kaons + protons, 2 - both

	// LOAD HISTOS
	fhadrons->cd();
	TH2F* hPion = (TH2F*)gDirectory->Get("hPionnSigmaEvsP");
	TH2F* hKaon = (TH2F*)gDirectory->Get("hKaonnSigmaEvsP");
	TH2F* hProton = (TH2F*)gDirectory->Get("hProtonnSigmaEvsP");

	// CONSTRUCT HADRON MEANS SPECTRUM
	int graphBins = hPion->GetNbinsX();
	float piMeans[graphBins], kaMeans[graphBins], prMeans[graphBins], kpMeans[graphBins];
	float piMeansErr[graphBins], kaMeansErr[graphBins], prMeansErr[graphBins], kpMeansErr[graphBins];
	float ptpos0[graphBins];
	for (int i = 0; i < graphBins-1; ++i)
	{
		piMeans[i] = hPion->ProjectionY("",i,i)->GetMean();
		piMeansErr[i] = hPion->ProjectionY("",i,i)->GetRMS();
		kaMeans[i] = hKaon->ProjectionY("",i,i)->GetMean();
		kaMeansErr[i] = hKaon->ProjectionY("",i,i)->GetRMS();
		prMeans[i] = hProton->ProjectionY("",i,i)->GetMean();
		prMeansErr[i] = hProton->ProjectionY("",i,i)->GetRMS();

		TH1D* hKP1 = (TH1D*)hKaon->ProjectionY("",i,i)->Clone("hKP1");
		TH1D* hKP2 = (TH1D*)hProton->ProjectionY("",i,i)->Clone("hKP2");
		hKP1->Add(hKP2);
		kpMeans[i] = hKP1->GetMean();
		kpMeansErr[i] = hKP1->GetRMS();

		ptpos0[i] = hPion->ProjectionX()->GetBinCenter(i);
	}

	TCanvas* cmeans = new TCanvas("cmeans","cmeans",800,600);
	TGraphErrors* piMean = new TGraphErrors(graphBins-1,ptpos0,piMeans,0,piMeansErr);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	piMean->SetTitle("hadron means of n#sigma_e;p (GeV/c);#mu_{n#sigma}");
	piMean->SetMarkerColor(kGreen+2);
	piMean->SetMarkerStyle(20);
	piMean->SetMarkerSize(1.2);
	piMean->SetLineColor(kGreen+2);
	piMean->SetLineWidth(2);
	piMean->Fit("pol3","","",2,15);
	piMean->Draw("");
	TGraphErrors* kaMean = new TGraphErrors(graphBins-1,ptpos0,kaMeans,0,kaMeansErr);
	kaMean->SetMarkerColor(kMagenta+2);
	kaMean->SetMarkerStyle(20);
	kaMean->SetMarkerSize(1.2);
	kaMean->SetLineColor(kMagenta+2);
	kaMean->SetLineWidth(2);
	kaMean->Fit("pol3","","",2,15);
	if (KPflag%2 == 0) kaMean->Draw("same");
	TGraphErrors* prMean = new TGraphErrors(graphBins-1,ptpos0,prMeans,0,prMeansErr);
	prMean->SetMarkerColor(kRed+2);
	prMean->SetMarkerStyle(20);
	prMean->SetMarkerSize(1.2);
	prMean->SetLineColor(kRed+2);
	prMean->SetLineWidth(2);
	prMean->Fit("pol3","","",2,15);
	if (KPflag%2 == 0) prMean->Draw("same");
	TGraphErrors* kpMean = new TGraphErrors(graphBins-1,ptpos0,kpMeans,0,kpMeansErr);
	kpMean->SetMarkerColor(kCyan+2);
	kpMean->SetMarkerStyle(20);
	kpMean->SetMarkerSize(1.2);
	kpMean->SetLineColor(kCyan+2);
	kpMean->SetLineWidth(2);
	kpMean->Fit("pol3","","",2,15);
	if (KPflag>0) kpMean->Draw("same");

	fin->cd();
	TH3F* hSigP0 = (TH3F*)gDirectory->Get("hElectronnSigmavsP");
	TH3F* hSigP = (TH3F*)hSigP0->Clone("hSigP");
	hSigP->Sumw2();
	TH1D* hSigZero[nbins-1];	
	TH1D* hSig[nbins-1];
	TCanvas* can[nbins-1];

	// DIVIDE HISTOGRAMS INTO EQUAL BINS
	/*TH1D* hP = hSigP->ProjectionX("",0,-1);
	int nentries = hP->Integral(0,-1);
	float step = (float)nentries/(nbins-1);

	float xbins[nbins];
	int ibins[nbins];
	xbins[0] = 2.0;
	int binj = 1;
	ibins[0] = hP->FindBin(xbins[0]);

	int sum = 0;
	for (int i = 0; i < hP->GetNbinsX(); ++i)
	{
		sum += hP->GetBinContent(i);
		if (sum >= step) {
			xbins[binj] = hP->GetBinCenter(i+1);
			ibins[binj] = i+1;
			sum = sum-step;
			binj++;
		}
	}

	for (int i = binj+1; i < nbins+1; ++i)
	{
		xbins[i] = 12.0;
		ibins[i] = 120;
	}

	// CHECK
	cout << "total entries: " << nentries << endl;
	cout << "step: " << step << endl;
	for (int i = 0; i < nbins; ++i)
	{
		cout << " interval: " << i << endl;
		cout << " first bin: " << ibins[i] << " end bin: " << ibins[i+1] << endl;
		cout << " coordinates: " << xbins[i] << " and " << xbins[i+1] << endl;
		cout << " integral1: " << hP->Integral(ibins[i],ibins[i+1]-1) << endl;
		int a = hSigP->ProjectionX()->FindBin(xbins[i]);
		int b = hSigP->ProjectionX()->FindBin(xbins[i+1]);
		cout << " integral2: " << hSigP->ProjectionY("",a,b-1); 
	}*/


	//float xbins[] = GetEqualBins(hSigP->ProjectionX(),12);
	/*TH1D* hP = hSigP->ProjectionX();
	double integral = hP->Integral();
	cout << "integral " << integral <<endl;
	xbins[0] = 2.0;
	double step = integral/nbins;
	cout << "step " << step << endl;
	int sum = 0;
	for (int i = 1; i < nbins; ++i)
	{
		int bin = hP->FindBin(xbins[i-1]);
		for (int j = bin; j < hP->GetNbinsX(); ++j)
		{

				cout << "sum " << sum << " bin " << j << " x " << hP->GetBinCenter(j) << endl;
			if (sum >= i*step) {
				cout << "br sum " << sum << " bin " << j << " x " << hP->GetBinCenter(j) << endl;
				xbins[i] = hP->GetBinCenter(j);
				break; 
			}
			sum += hP->GetBinContent(j);
		}
	}*/

	// FIT VARIABLES
	float emu0 = -0.3;
	TF1* fAll[nbins-1];

	// RESULTS
	float emu[nbins];
	float emuerr[nbins];
	float esig[nbins];
	float esigerr[nbins];
	float ptpos[nbins];
	float ptl[nbins];
	float ptr[nbins];
	float eff[nbins];
	float pur[nbins];
	float effT[nbins];
	float effL[nbins];

	// STYLING
	TLatex latex;
	latex.SetTextSize(0.045);
   	latex.SetNDC(kTRUE);
	gStyle->SetOptStat(0);/*
	gStyle->SetOptFit(10);
	gStyle->SetStatW(0.3);
	gStyle->SetFitFormat("5.3g");
	gStyle->SetStatFontSize(0.13);*/
	TLegend* l1 = new TLegend(0.74,0.58,0.85,0.75);
   	l1->SetNColumns(1);
   	l1->SetFillStyle(0);
   	l1->SetBorderSize(0);
   	l1->SetTextSize(0.042);

	// LOOPING PT BINS
	for (int i = 0; i < nbins-1; ++i)
	{
		int a = hSigP->ProjectionX()->FindBin(xbins[i]);
		int b = hSigP->ProjectionX()->FindBin(xbins[i+1]);

		TString name = Form("hSig%i",i);
		hSigZero[i] = hSigP->ProjectionY("",a,b-1,0,-1);
		//hSigZero[i] = hSigP->ProjectionY("",a+1,b,0,-1);
		hSigZero[i]->Rebin(6);
		hSigZero[i]->SetMarkerStyle(20);
		hSigZero[i]->SetMarkerSize(0.8);

		hSig[i] = (TH1D*)hSigZero[i]->Clone(name);

		cout << xbins[i] << " - " << xbins[i+1] << ", step " << i << ", " << a << " : " << b << endl;
		cout << "entries: " << hSig[i]->Integral() << endl;

		TString cname = Form("can%i",i);
		can[i] = new TCanvas(cname,cname,500,400);
		can[i]->SetLogy();
		can[i]->SetGridx();
		hSig[i]->Draw();
		can[i]->Update();
		//cout << pow(10,can[i]->GetFrame()->GetY1()) << " and " << pow(10,can[i]->GetFrame()->GetY2()) << endl;
		
		float pimu = hPion->ProjectionY("",a,b-1)->GetMean();
		TH1D* hKP1 = (TH1D*)hKaon->ProjectionY("",a,b-1)->Clone("hKP1");
		TH1D* hKP2 = (TH1D*)hProton->ProjectionY("",a,b-1)->Clone("hKP2");
		hKP1->Add(hKP2);
		float kpmu = hKP1->GetMean();
		float emuerr0 = 0;
		double pars[9];
		TString faname = Form("fAll%i",i);
		float hmax = hSig[i]->GetMaximum();
		cout << "max is " << hmax << endl;

		fAll[i] = new TF1(faname,"gaus(0)+gaus(3)+gaus(6)",-10,5);
		fAll[i]->SetParameters(hmax,emu0,1.00,hmax,pimu+emu0,1.2,hmax,kpmu+emu0,1.25);
		fAll[i]->SetParLimits(0,0,100*hmax);
		fAll[i]->SetParLimits(3,0,100*hmax);
		fAll[i]->SetParLimits(6,0,100*hmax);
		fAll[i]->SetParLimits(1,-0.6,0.1);
		fAll[i]->SetParLimits(2,0.80,1.20);
		fAll[i]->SetParLimits(4,-6,-1.5);
		fAll[i]->SetParLimits(5,1.1,1.3);
		fAll[i]->SetParLimits(7,-7,-4);
		fAll[i]->SetParLimits(8,1.15,1.25);
		TFitResultPtr fr = hSig[i]->Fit(faname,"BIMSQ","",-10,2.5);
		fAll[i]->GetParameters(pars);
		for (int iR = 0; iR < nRepeats; ++iR)
		{
			cout << "pimu " << pimu << " emu " << pars[1] << " fixing to " << pimu+pars[1] << endl;
			emuerr0 = fAll[i]->GetParError(1);
			fAll[i]->SetParameter(4,pimu+pars[1]);
			fAll[i]->SetParLimits(4,pimu+pars[1]-3*emuerr0,pimu+pars[1]+3*emuerr0);
			fAll[i]->SetParameter(7,kpmu+pars[1]);
			fAll[i]->SetParLimits(7,kpmu+pars[1]-3*emuerr0,kpmu+pars[1]+3*emuerr0);
			fr = hSig[i]->Fit(faname,"BIMSQ","",-10,2.5);
			fAll[i]->GetParameters(pars);
		}

		//fr->PrintCovMatrix(cout);
		TMatrixDSym cov = fr->GetCovarianceMatrix();
		TMatrixDSym newCov(2);
		newCov(0,0) = cov(1,1);		newCov(0,1) = cov(1,2);
		newCov(1,0) = cov(2,1);		newCov(1,1) = cov(2,2);
		//newCov.Print();
		double params[2];
		params[0] = fAll[i]->GetParameter(1);
		params[1] = fAll[i]->GetParameter(2);
		double cov2[4];
		cov2[0] = cov(1,1); cov2[1] = cov(1,2);
		cov2[2] = cov(2,1); cov2[3] = cov(2,2);

		//for (int j = 3; j < 9; ++j) fAll[i]->FixParameter(j,pars[j]);
		//hSig[i]->Fit(faname,"BIQM","",pars[1]-6*pars[2],3);

		// zoom out
		hSig[i]->GetYaxis()->SetRangeUser(pow(10,gPad->GetUymin()),5*pow(10,gPad->GetUymax()));

		// TAKE CONTRIBUTIONS
		TF1* fE = new TF1("fE","gaus(0)",-10,5);
		fE->SetLineColor(kBlue);
		TF1* fPi = new TF1("fPi","gaus(0)",-10,5);
		fPi->SetLineColor(kGreen+1);
		TF1* fKP = new TF1("fKP","gaus(0)",-10,5);
		fKP->SetLineColor(kMagenta+2);
		fE->SetParameters(fAll[i]->GetParameters());
		fPi->SetParameters(&fAll[i]->GetParameters()[3]);
		fKP->SetParameters(&fAll[i]->GetParameters()[6]);
		fE->Draw("same");
		fPi->Draw("same");
		fKP->Draw("same");
		fAll[i]->Draw("same");

		float sc = (float)hSig[i]->GetNbinsX()/15;
		cout << sc << endl;
		float iE = fE->Integral(left,right);
		float iE0 = fE->Integral(-10,5);
		float iA = fAll[i]->Integral(left,right);

		cout << "int then " << fAll[i]->Integral(left,right) << endl;

		TF1* fGaus = new TF1("fG","exp(-0.5*((x-[0])/[1])**2)",-10,5);
		fGaus->SetParameters(params);
		float efferr = (iE/iE0)*fGaus->IntegralError(left,right,params,cov2)/fGaus->Integral(-10,5);
		cout << "err ERROROROOROR is " << efferr << endl;
		cout << "erooro " << fGaus->IntegralError(-10,5,params,cov2)/fGaus->Integral(-10,5) << endl; // can i neglect error of denom

		//fG->SetParameters(1,)

		/*fAll[i]->FixParameter(0,pars[0]);
		for (int j = 3; j < 9; ++j)
		{
			fAll[i]->FixParameter(j,pars[j]);
		}
		hSig[i]->Fit(faname,"BIM","",-10,4);
		cout << "int now " << fAll[i]->Integral(left,right) << endl;
		float efferr = fAll[i]->IntegralError(left,right)*sc/iE0;
		cout << "err is " << efferr << endl;*/

		//cout << "iE = " << iE << " +- " << fE->IntegralError(left,right)*sc << endl;

		eff[i] = iE/iE0;
		pur[i] = iE/iA;

		cout << "eff: " << eff[i] << " , pur: " << pur[i] << endl;

		emu[i] = fAll[i]->GetParameter(1);
		emuerr[i] = fAll[i]->GetParError(1);
		cout << emuerr[i] << endl;
		esig[i] = fabs(fAll[i]->GetParameter(2));
		esigerr[i] = fAll[i]->GetParError(2);
		ptpos[i] = 0.5*(xbins[i]+xbins[i+1]);
		ptl[i] = ptpos[i] - xbins[i];
		ptr[i] = xbins[i+1] - ptpos[i];
		effT[i] = 0;
		effL[i] = 0;

		TF1* fG = new TF1("fG","gaus(0)",-10,5);


		for (int j=-1;j<=1; j++)
		{
			for (int k=-1;k<=1;k++) 
			{
				if (j==0 || k==0) continue;
				fG->SetParameters(1,emu[i]+j*emuerr[i],esig[i]+k*esigerr[i]);
				float eff0 = fG->Integral(left,right) / fG->Integral(-10,5);
				//cout << "eff0 is " << eff0 << endl;
				if ( eff[i]-eff0 > effT[i] ) effT[i] = eff[i] - eff0;
				if ( eff0-eff[i] > effL[i] ) effL[i] = eff0 - eff[i];
			}
		}

		cout << "errs old: " << effT[i] << " and " << effL[i] << endl;
		cout << " vs " << efferr/2 << endl;

		//or from integralError()
		bool interr = 1;
		if (interr){
		effT[i] = efferr/2;
		effL[i] = efferr/2; }

		TString texstr(Form("%.1f < p < %.1f",xbins[i],xbins[i+1]));
		//texstr += xbins[i];
		//texstr += " < p_{T} < ";
		//texstr = xbins[i] + texstr;
		//texstr += xbins[i+1];

  		latex.DrawLatex(0.15, 0.80, texstr);
  		if (i == 0)
  		{
  			l1->AddEntry(fE, " e", "l");
  			l1->AddEntry(fPi, " #pi", "l");
  			l1->AddEntry(fKP, " K+p", "l");
  		}
  		l1->Draw();

  		can[i]->SaveAs(Form("nsigeff_phe_%.1f_%.1f.png",xbins[i],xbins[i+1]));

		//cout << " + " << eff[i] - effT[i] << ", - " << eff[i] + effL[i] << endl;
	}

	// PLOT EFF VS PT
	TCanvas* ceff = new TCanvas("ceff","ceff",800,600);
	ceff->SetGrid();
	TGraphAsymmErrors* plot = new TGraphAsymmErrors(nbins-1,ptpos,eff,ptl,ptr,effL,effT);
	plot->SetTitle(";p (GeV/c);n#sigma_{e} efficiency");
	TAxis * xax = plot->GetXaxis();
  	TAxis * yax = plot->GetYaxis();
  	xax->CenterTitle();
  	xax->SetLabelSize(0.04);
  	xax->SetTitleSize(0.04);
  	yax->CenterTitle();
  	yax->SetLabelSize(0.04);
  	yax->SetTitleSize(0.05);
  	yax->SetTitleOffset(+0.9);
  	plot->SetMinimum(0.75);
  	plot->SetMaximum(1.05);
  	plot->SetMarkerSize(0.8);
  	plot->SetMarkerStyle(20);
  	plot->SetMarkerColor(2);
  	plot->SetLineColor(2);
  	plot->SetLineWidth(2);
  	plot->SetFillColor(20);
    plot->SetFillStyle(1001);
   	plot->Draw("A2");
   	plot->Draw("P SAME");

   	TLatex tex;
  	tex.SetTextSize(0.035);
   	tex.SetNDC(kTRUE);
  	tex.DrawLatex(0.55, 0.81, "m<100 MeV electron efficiency");
  	tex.DrawLatex(0.55, 0.76, "-1.5 < n#sigma_{e}< 3.0 cut");
  	tex.DrawLatex(0.55, 0.71, "|y|<1  0-100%  AuAu14 BHT2");

  	ceff->SaveAs("nsigeff_phe.png");

  	TCanvas* cmu = new TCanvas("cmu","cmu",800,600);
	cmu->SetGrid();
	TGraphAsymmErrors* plotmu = new TGraphAsymmErrors(nbins-1,ptpos,emu,ptl,ptr,emuerr,emuerr);
	plotmu->SetTitle(";p (GeV/c);n#sigma_{e} mean");
	plotmu->SetMarkerSize(0.8);
  	plotmu->SetMarkerStyle(20);
  	plotmu->SetMarkerColor(2);
  	plotmu->SetLineColor(2);
  	plotmu->SetLineWidth(2);
  	plotmu->SetMinimum(-0.8);
  	plotmu->SetMaximum(0.2);
	plotmu->Draw("AP");
  	cmu->SaveAs("nsigmu_phe.png");

	TCanvas* csig = new TCanvas("csig","csig",800,600);
	csig->SetGrid();
	TGraphAsymmErrors* plotsig = new TGraphAsymmErrors(nbins-1,ptpos,esig,ptl,ptr,esigerr,esigerr);
	plotsig->SetTitle(";p (GeV/c);n#sigma_{e} sigma");
	plotsig->SetMarkerSize(0.8);
  	plotsig->SetMarkerStyle(20);
  	plotsig->SetMarkerColor(2);
  	plotsig->SetLineColor(2);
  	plotsig->SetLineWidth(2);
  	plotsig->SetMinimum(0.75);
  	plotsig->SetMaximum(1.5);
	plotsig->Draw("AP");
  	csig->SaveAs("nsigsig_phe.png");

  	TCanvas* cfit = new TCanvas("cfit","cfit",800,600);
	cfit->SetGrid();

	fout->cd();

   	TF1* efit = new TF1("efit","pol1",2,12);
   	plot->SetMarkerColor(1);
  	plot->SetLineColor(15);
  	plot->SetLineWidth(1);
   	plot->Fit(efit,"M0");
   	plot->Fit(efit,"M0");
   	plot->Draw("AP");
   	efit->SetLineWidth(2);
   	efit->Draw("same");
   	efit->Write();

   	TF1* efitT = (TF1*)efit->Clone("efitT");
   	efitT->SetParameter(0,efit->GetParameter(0)+efit->GetParError(0));
   	efitT->SetParameter(1,efit->GetParameter(1)+efit->GetParError(1));
   	efitT->SetLineWidth(2);
   	efitT->SetLineStyle(2);
   	efitT->Draw("same");
   	efitT->Write();

   	TF1* efitL = (TF1*)efit->Clone("efitL");
   	efitL->SetParameter(0,efit->GetParameter(0)-efit->GetParError(0));
   	efitL->SetParameter(1,efit->GetParameter(1)-efit->GetParError(1));
   	efitL->SetLineWidth(2);
   	efitL->SetLineStyle(2);
   	efitL->Draw("same");
   	efitL->Write();

   	//efit->Draw("same");
   	//efit->Draw();
   	tex.DrawLatex(0.55, 0.81, "m<100 MeV electron efficiency");
  	tex.DrawLatex(0.55, 0.76, "-1.5 < n#sigma_{e}< 3.0 cut");
  	tex.DrawLatex(0.55, 0.71, "|y|<1  0-100%  AuAu14 BHT2");

  	fout->Write();

  }