#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void drawHFTeffect() {

	TFile* fin = new TFile("hist_tmpH.root","READ");
	//TFile* fin = new TFile("hist_tmpE99lm.root","READ");
	//TFile* fin = new TFile("hist_tmpnsiglm_1.root","READ");
	TH3F* hVvP = (TH3F*)fin->Get("hEventVzvsNPrimHard");

	/*TFile* fin = new TFile("rootOutputs/15lm_0403_tree.root","READ");
	TTree* upsTree = (TTree*)fin->Get("upsTree");
	TH3F* hVvP = new TH3F("hVvP","",400,-50,50,50,0,50,10,0,10);
	upsTree->Draw("tEventCent9:tNElectrons:tEventTpcVz>>hVvP");*/


	const int nbins = 60;
	float vz[nbins];
	float errvz[nbins];
	float nPrim29[nbins];
	float err2nP29[nbins];
	float nPrim25[nbins];
	float nPrim57[nbins];
	float nPrim79[nbins];

	
	for (int i = 0; i < nbins; i++)
	{
		vz[i] = -30. + (float)i*(60.)/(nbins-1);
		errvz[i] = 0.01;
	}

	for (int i = 0; i < nbins-1; ++i)
	{
		int a = hVvP->ProjectionX()->FindBin(vz[i]);
		int b = hVvP->ProjectionX()->FindBin(vz[i+1]);
		nPrim29[i] = hVvP->ProjectionY("",a,b,3,9)->GetMean();
		err2nP29[i] = hVvP->ProjectionY("",a,b,3,9)->GetMeanError();
		//nPrim25[i] = hVvP->ProjectionY("",a,b,3,5)->GetMean();
		//nPrim57[i] = hVvP->ProjectionY("",a,b,6,7)->GetMean();
		//nPrim79[i] = hVvP->ProjectionY("",a,b,8,9)->GetMean();
	}

	TGraphErrors* plot = new TGraphErrors(nbins-1,vz,nPrim29,errvz,err2nP29);
	plot->SetTitle("; v_{z} [cm]; <primary tracks> (norm.)");
	plot->SetMarkerStyle(20);
	plot->SetMarkerSize(1.0);
	plot->SetMinimum(0.6);
	plot->SetMaximum(1.2);

	TF1* fpol0 = new TF1("fpol0","pol0",-30,30);
	TF1* fpol0B = new TF1("fpol0B","pol0",-30,30);
	fpol0B->SetLineColor(kBlue);
	plot->Fit("fpol0","0","",-30,-16);
	float scale = fpol0->GetParameter(0);
	plot->Fit("fpol0","0","",16,30);
	scale = 0.5*( scale + fpol0->GetParameter(0));
	for (int i = 0; i < plot->GetN(); i++) plot->GetY()[i] *= 1./scale;
	for (int i = 0; i < plot->GetN(); i++) plot->GetEY()[i] *= 1./scale;

	plot->Fit("fpol0B","0","",-10,10);
	
	fpol0->SetParameter(0,1);
	//fpol0B->SetParameter(0,fpol0B->GetParameter(0)*1./scale);
	plot->Draw("AP");
	fpol0->Draw("same");
	fpol0B->Draw("same");


}