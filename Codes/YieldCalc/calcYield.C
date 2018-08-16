#include <iostream>
#include <string>
#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

#include "style.C"
//#include "LechLabels.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "TLatex.h"

//#include "../pp2011constants.h"
#include "LoadXYscanned.h"
#include "GraphOperations.h"
#include "PropagateErrors.h"

#include "TMatrixDSym.h"

using namespace std;
using namespace RooFit;

Int_t padw = 800, padh = 600;

const Int_t nPadX = 1;
const Int_t nPadY = 1;

const Int_t nHist = 1;

//TGraphAsymmErrors *LoadRawYield(TString name);
TGraphAsymmErrors *LoadEfficiency(TString nameL, TString nameH, double cL = 0.457395797, double cH = 0.542604203);
TGraphAsymmErrors *CalculateYield(TGraphAsymmErrors *g, TGraphAsymmErrors *eff,
		double lum, double nEvSeen, double nEvAll, double *dcent, double dy = 1.0);
/////////////////////////////////////////////////////////////
void calcYield() {

	style();
	gStyle->SetOptStat(0);
	gStyle->SetOptDate(0);
	gStyle->SetLineWidth(1);
	gROOT->ForceStyle();
	//gStyle->SetHistLineColor(kBlue);

	// Constants for corrected yield calculation
	double lumi = 5.201188;		// [nb^{-1}]
	double nEvQuoted = 149.417;
	double nEvSeen = 114.209;//e+08,115.364;//e+08;
	double dcent[] = {0.3,0.2,0.1}; // centrality bin width 30-60%, 10-30%, 0-10%
	double dcentInt[] = {0.6}; // centrality integrated width 0-60%
	double dy = 1.0; // |y|<0.5

	double frac2Sover2S3S = 0.689;
	double frac3Sover2S3S = 0.311;

	gSystem->cd("data");


	//TFile *file = new TFile("combine_11_14_16.root","READ");

	TGraphAsymmErrors* rawYield_1S_2S_3S = LoadXYscanned("yield_1S.txt");
	TGraphAsymmErrors* rawYield_1S = LoadXYscanned("yield_1S.txt");
	TGraphAsymmErrors* rawYield_2S_3S = LoadXYscanned("yield_2S3S.txt");

	TGraphAsymmErrors* rawYield_1S_2S_3S_int = LoadXYscanned("yield_1S_int.txt");
	TGraphAsymmErrors* rawYield_1S_int = LoadXYscanned("yield_1S_int.txt");
	TGraphAsymmErrors* rawYield_2S_3S_int = LoadXYscanned("yield_2S3S_int.txt");

	gSystem->cd("../");


	gSystem->cd("efficiency");

	// Efficiency
/*
	TGraphAsymmErrors* eff_1S = (TGraphAsymmErrors*)file->Get("gr_1S_Sys");
	TGraphAsymmErrors* eff_2S = (TGraphAsymmErrors*)file->Get("gr_1S_Stat");
	TGraphAsymmErrors* eff_3S = (TGraphAsymmErrors*)file->Get("gr_2S3S_Stat");
*/
	TGraphAsymmErrors* eff_1S = LoadEfficiency("eff_1S_L.txt", "eff_1S_H.txt");
	TGraphAsymmErrors* eff_2S = LoadEfficiency("eff_2S_L.txt", "eff_2S_H.txt");
	TGraphAsymmErrors* eff_3S = LoadEfficiency("eff_3S_L.txt", "eff_3S_H.txt");

	TGraphAsymmErrors* eff_1S_int = LoadEfficiency("eff_1S_L_int.txt", "eff_1S_H_int.txt");
	TGraphAsymmErrors* eff_2S_int = LoadEfficiency("eff_2S_L_int.txt", "eff_2S_H_int.txt");
	TGraphAsymmErrors* eff_3S_int = LoadEfficiency("eff_3S_L_int.txt", "eff_3S_H_int.txt");

	TGraphAsymmErrors *eff_2S_scaled = (TGraphAsymmErrors *)eff_2S->Clone(Form("Raa_%s", eff_2S->GetName()));
	TGraphAsymmErrors *eff_3S_scaled = (TGraphAsymmErrors *)eff_3S->Clone(Form("Raa_%s", eff_3S->GetName()));

	TGraphAsymmErrors *eff_2S_int_scaled = (TGraphAsymmErrors *)eff_2S_int->Clone(Form("Raa_%s", eff_2S_int->GetName()));
	TGraphAsymmErrors *eff_3S_int_scaled = (TGraphAsymmErrors *)eff_3S_int->Clone(Form("Raa_%s", eff_3S_int->GetName()));

	MultiplyGraph(eff_2S_scaled, frac2Sover2S3S);
	MultiplyGraph(eff_3S_scaled, frac3Sover2S3S);

	MultiplyGraph(eff_2S_int_scaled, frac2Sover2S3S);
	MultiplyGraph(eff_3S_int_scaled, frac3Sover2S3S);

	TGraphAsymmErrors *eff_2S_3S = AddGraphs(eff_2S_scaled, eff_3S_scaled);
	TGraphAsymmErrors *eff_2S_3S_int = AddGraphs(eff_2S_int_scaled, eff_3S_int_scaled);

	gSystem->cd("../");

	//-----------------------

	//TGraphAsymmErrors *yield_STAR_1S_2S_3S = CalculateYield(rawYield_1S_2S_3S, eff_1S_2S_3S, lumi, nEvQuoted, nEvSeen, dcent, dy);
	TGraphAsymmErrors *yield_STAR_1S = CalculateYield(rawYield_1S, eff_1S, lumi, nEvSeen, nEvQuoted, dcent, dy);
	TGraphAsymmErrors *yield_STAR_2S_3S = CalculateYield(rawYield_2S_3S, eff_2S_3S, lumi, nEvSeen, nEvQuoted, dcent, dy);

	TGraphAsymmErrors *yield_STAR_1S_int = CalculateYield(rawYield_1S_int, eff_1S_int, lumi, nEvSeen, nEvQuoted, dcentInt, dy);
	TGraphAsymmErrors *yield_STAR_2S_3S_int = CalculateYield(rawYield_2S_3S_int, eff_2S_3S_int, lumi, nEvSeen, nEvQuoted, dcentInt, dy);


	MultiplyGraph(yield_STAR_1S, 1.0);
	MultiplyGraph(yield_STAR_2S_3S, 1.0);

	//-----------------------

	//yield_STAR_1S_2S_3S->SetLineColor(kRed+1);
	yield_STAR_1S->SetLineColor(kGreen+1);
	yield_STAR_2S_3S->SetLineColor(kBlue+1);

	//yield_STAR_1S_2S_3S->SetMarkerColor(kRed+1);
	yield_STAR_1S->SetMarkerColor(kGreen+1);
	yield_STAR_2S_3S->SetMarkerColor(kBlue+1);

	//yield_STAR_1S_2S_3S->SetLineWidth(2);
	yield_STAR_1S->SetLineWidth(2);
	yield_STAR_2S_3S->SetLineWidth(2);


	//yield_STAR_1S_2S_3S_int->SetLineColor(kRed+1);
	yield_STAR_1S_int->SetLineColor(kGreen+1);
	yield_STAR_2S_3S_int->SetLineColor(kBlue+1);

	//yield_STAR_1S_2S_3S_int->SetMarkerColor(kRed+1);
	yield_STAR_1S_int->SetMarkerColor(kGreen+1);
	yield_STAR_2S_3S_int->SetMarkerColor(kBlue+1);

	//yield_STAR_1S_2S_3S_int->SetLineWidth(2);
	yield_STAR_1S_int->SetLineWidth(2);
	yield_STAR_2S_3S_int->SetLineWidth(2);

/*
	syst_xsec_STAR_1S_2S_3S->SetLineColor(kRed+2);
	syst_xsec_STAR_1S->SetLineColor(kGreen+2);
	syst_xsec_STAR_2S_3S->SetLineColor(kBlue+3);

	syst_xsec_STAR_1S_2S_3S->SetFillColor(kRed);
	syst_xsec_STAR_1S->SetFillColor(kGreen);
	syst_xsec_STAR_2S_3S->SetFillColor(kBlue);

	syst_xsec_STAR_1S_2S_3S->SetFillStyle(0);
	syst_xsec_STAR_1S->SetFillStyle(0);
	syst_xsec_STAR_2S_3S->SetFillStyle(0);

	syst_xsec_STAR_1S_2S_3S->SetLineWidth(2);
	syst_xsec_STAR_1S->SetLineWidth(2);
	syst_xsec_STAR_2S_3S->SetLineWidth(2);

	glob_xsec_STAR_1S_2S_3S->SetLineColor(kRed-10);
	glob_xsec_STAR_1S->SetLineColor(kGreen-10);
	glob_xsec_STAR_2S_3S->SetLineColor(kBlue-9);

	glob_xsec_STAR_1S_2S_3S->SetFillColor(kRed-10);
	glob_xsec_STAR_1S->SetFillColor(kGreen-10);
	glob_xsec_STAR_2S_3S->SetFillColor(kBlue-9);

	glob_xsec_STAR_1S_2S_3S->SetFillStyle(1001);
	glob_xsec_STAR_1S->SetFillStyle(1001);
	glob_xsec_STAR_2S_3S->SetFillStyle(1001);

	glob_xsec_STAR_1S_2S_3S->SetLineWidth(2);
	glob_xsec_STAR_1S->SetLineWidth(2);
	glob_xsec_STAR_2S_3S->SetLineWidth(2);*/

	//yield_STAR_1S_2S_3S->SetMarkerStyle(20);
	yield_STAR_1S->SetMarkerStyle(33);
	yield_STAR_2S_3S->SetMarkerStyle(21);

	//yield_STAR_1S_2S_3S->SetMarkerSize(1.4);
	yield_STAR_1S->SetMarkerSize(2);
	yield_STAR_2S_3S->SetMarkerSize(1.2);


	//yield_STAR_1S_2S_3S_int->SetMarkerStyle(24);
	yield_STAR_1S_int->SetMarkerStyle(27);
	yield_STAR_2S_3S_int->SetMarkerStyle(25);

	//yield_STAR_1S_2S_3S_int->SetMarkerSize(1.4);
	yield_STAR_1S_int->SetMarkerSize(2);
	yield_STAR_2S_3S_int->SetMarkerSize(1.2);

/*	//----------------------------------------------

	gSystem->cd("worldData");


	gSystem->cd("../");

	//----------------------------------------------


	//----------------------------------------------

	gSystem->cd("models");


	gSystem->cd("../");

*/
	//----------------------------------------------

	TCanvas *cnv = new TCanvas("cnv", "cnv", nPadX*padw, nPadY*padh);

	//TH1D *hBack = new TH1D("hBack", "Corrected #varUpsilon yield vs. N_{part}", 400, 0.0, 400.0);
	TH1D *hBack = new TH1D("hBack", "#varUpsilon cross section vs. N_{part}", 400, 0.0, 400.0);

	hBack->SetTitle("");
	//hBack->GetYaxis()->SetTitle("B_{ee} #frac{dN}{dy} [1]");
	hBack->GetYaxis()->SetTitle("B_{ee} #frac{d#sigma}{dy} [nb]");
	hBack->GetXaxis()->SetTitle("N_{part} [1]");

	hBack->GetYaxis()->SetRangeUser( 5e1, 1e6);
	hBack->GetXaxis()->SetRangeUser( 0.0, 400.0);


	cnv->cd()->SetLogy();	//----------------

	hBack->Draw();

/*
	glob_xsec_STAR_1S_2S_3S->Draw("samee2");
	glob_xsec_STAR_1S->Draw("samee2");
	glob_xsec_STAR_2S_3S->Draw("samee2");

	syst_xsec_STAR_1S_2S_3S->Draw("samee2");
	syst_xsec_STAR_1S->Draw("samee2");
	syst_xsec_STAR_2S_3S->Draw("samee2");
*/
	//yield_STAR_1S_2S_3S->Draw("psame1");
	yield_STAR_1S->Draw("psame1");
	yield_STAR_2S_3S->Draw("psame1");

	//yield_STAR_1S_2S_3S_int->Draw("psame1");
	yield_STAR_1S_int->Draw("psame1");
	yield_STAR_2S_3S_int->Draw("psame1");


	//---------------


	//---------------

	TLatex tl;
	tl.SetTextFont(42);
	tl.SetTextSize(0.042);
	tl.SetNDC();
	//if(Preliminary) tl.DrawLatex(0.22, 0.2, "#color[2]{STAR Preliminary}");
	tl.DrawLatex(0.19, 0.8, "STAR Au+Au @ 200 GeV");
	tl.DrawLatex(0.19, 0.74, " Run14 |y|<0.5");

	TLegend *leg = new TLegend(0.55, 0.5, 0.86, 0.86);
	TLegend *leg2 = new TLegend(0.3, 0.66, 0.52, 0.78);
	//leg->AddEntry(yield_STAR_1S_2S_3S, "STAR #varUpsilon(1S+2S+3S)", "p");
	leg->AddEntry(yield_STAR_1S, "STAR #varUpsilon(1S)", "p");
	leg->AddEntry(yield_STAR_2S_3S, "STAR #varUpsilon(2S+3S)", "p");
	leg->AddEntry(yield_STAR_1S_int, "STAR #varUpsilon(1S) 0-60%", "p");
	leg->AddEntry(yield_STAR_2S_3S_int, "STAR #varUpsilon(2S+3S) 0-60%", "p");
/*
	leg2->AddEntry(syst_yield_STAR_1S_2S_3S, "Uncorrelated syst.", "f");
	leg2->AddEntry(glob_yield_STAR_1S_2S_3S, "Correlated syst.", "f");
*/
	leg->SetFillColor(kWhite);
	leg->SetFillStyle(1000);
	leg2->SetFillColor(kWhite);
	leg2->SetFillStyle(1000);

	leg->Draw("same");
	//leg2->Draw("same");

	//------------------------------------

	TString mkdir = "output";
	gSystem->MakeDirectory(mkdir.Data());
	gSystem->cd(mkdir.Data());

	cnv->cd()->SaveAs("UpsCorrYield.png");
	cnv->cd()->SaveAs("UpsCorrYield.pdf");
	cnv->cd()->SaveAs("UpsCorrYield.eps");

	yield_STAR_1S->SetName("UpsCorrYield_STAR_1S");
	yield_STAR_2S_3S->SetName("UpsCorrYield_STAR_2S3S");

	yield_STAR_1S_int->SetName("UpsCorrYield_STAR_1S_int");
	yield_STAR_2S_3S_int->SetName("UpsCorrYield_STAR_2S3S_int");

	WriteGraph(yield_STAR_1S);
	WriteGraph(yield_STAR_2S_3S);

	WriteGraph(yield_STAR_1S_int);
	WriteGraph(yield_STAR_2S_3S_int);

	gSystem->cd("../");


}


TGraphAsymmErrors *LoadEfficiency(TString nameL, TString nameH, double cL, double cH)
{

	TGraphAsymmErrors *effL = LoadXYscanned(nameL);
	TGraphAsymmErrors *effH = LoadXYscanned(nameH);

	int nL = effL->GetN();
	int nH = effH->GetN();

	if(nL != nH)
	{
		cout<<"WARNING! Efficiency graphs have different number of points. Exiting."<<endl;
		return NULL;
	}

	MultiplyGraph(effL, cL);
	MultiplyGraph(effH, cH);

	TGraphAsymmErrors* eff = AddGraphs(effL, effH);

	eff->SetName(Form("combEff_%s", effL->GetName()));

	return eff;
}

TGraphAsymmErrors *CalculateYield(TGraphAsymmErrors *g, TGraphAsymmErrors *eff,
		double lum, double nEvSeen, double nEvAll, double *dcent, double dy)
{

	int n = g->GetN();
	int neff = eff->GetN();

	TGraphAsymmErrors *yield = (TGraphAsymmErrors *)g->Clone(Form("corrYield_%s", g->GetName()));

	if(n != neff)
	{
		cout<<"WARNING! Yield and efficiency graphs have different number of points. Exiting."<<endl;
		return yield;
	}

	for (int i = 0; i < n; ++i) {

		// yield
		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);
		double ey = (eh+el)/2.0;

		// eff
		double eff_x, eff_y, eff_eh, eff_el;

		eff->GetPoint(i, eff_x, eff_y);
		eff_eh = eff->GetErrorYhigh(i);
		eff_el = eff->GetErrorYlow(i);
		double eff_ey = (eff_eh+eff_el)/2.0;


		RooRealVar rawYield("rawYield","Raw yield", 1.0, 0.0, 1e6);
		RooRealVar eff("eff","Efficiency", 1.0, 0.0, 1.0);
		RooRealVar effLum("effLum","Effective integrated luminosity", 1.0, 0.0, 1e12);
		RooRealVar dCent("dCent","Centrality bin width", 1.0, 0.0, 1e6);
		RooRealVar rapWidth("rapWidth","Rapidity bin width", 1.0, 0.0, 1e6);

		effLum.setUnit("nb^{-1}");

		rawYield.setVal(y);
		eff.setVal(eff_y);
		effLum.setVal(lum*(nEvSeen/nEvAll)); // effective integrated luminosity
		dCent.setVal(dcent[i]);
		rapWidth.setVal(dy);

		// set errors
		rawYield.setError(ey);
		eff.setError(eff_ey);
		eff.setError(0.2*eff_y);
		effLum.setError(0.0);
		dCent.setError(0.0);
		rapWidth.setError(0.0);

		// constants
		effLum.setConstant(kTRUE);
		dCent.setConstant(kTRUE);
		rapWidth.setConstant(kTRUE);

		cout<<endl;
		cout<<"rawYield"<<"\t";
		cout<<"eff"<<"\t";
		cout<<"effLum"<<"\t";
		cout<<"dCent"<<"\t";
		cout<<"rapWidth"<<"\t";
		cout<<endl;

		cout<<rawYield.getVal()<<"\t";
		cout<<eff.getVal()<<"\t";
		cout<<effLum.getVal()<<"\t";
		cout<<dCent.getVal()<<"\t";
		cout<<rapWidth.getVal()<<"\t";
		cout<<endl;

		cout<<"rawYieldErr"<<"\t";
		cout<<"effErr"<<"\t";
		cout<<endl;

		cout<<rawYield.getError()<<"\t";
		cout<<eff.getError()<<"\t";
		cout<<endl;

		RooGenericPdf corrYield("corrYield", "Corrected yield",
				"rawYield/(eff*effLum*dCent*rapWidth)",RooArgSet(rawYield,eff,effLum,dCent,rapWidth));

		RooArgSet argSet(rawYield,eff);
		TMatrixDSym covMat(2);

		GetCovMatUncorr(argSet, covMat);

		covMat.Print();

		double  corrYield_val = corrYield.getVal();
		double  corrYield_err = getPropagatedError(corrYield, argSet, covMat);

		cout<<"corrYield = "<<corrYield_val<<endl;
		cout<<"corrYield_err = "<<corrYield_err<<endl;

		yield->SetPoint(i, x, corrYield_val);
		yield->SetPointEYhigh(i, corrYield_err);
		yield->SetPointEYlow(i, corrYield_err);

	}

	return yield;
}
