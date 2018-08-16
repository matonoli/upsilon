//#include <exception>
//#include <assert.h>
#include <iostream>
//#include <vector>
//#include <algorithm>
//#include <math.h>
#include <string>
#include <fstream>

#include "TString.h"
//#include "TClonesArray.h"
//#include "TRefArray.h"
//#include "TRef.h"
#include "TFile.h"
//#include "TArrayI.h"
//#include "TTree.h"
#include "TH1.h"
//#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
//#include "TBranch.h"
//#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
//#include "TPaveLabel.h"
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

TGraphAsymmErrors *CalculateRaa(TGraphAsymmErrors *g, TGraphAsymmErrors *gncoll,
		double ups_xsec_pp, double ups_xsec_pp_err, double sig_pp, double sig_AA);
/////////////////////////////////////////////////////////////
void calcRaa() {

	style();
	gStyle->SetOptStat(0);
	gStyle->SetOptDate(0);
	gStyle->SetLineWidth(1);
	gROOT->ForceStyle();
	//gStyle->SetHistLineColor(kBlue);

	// Constants for Raa calculation
/*	double lumi = 5201.188E3;		// [mb^{-1}]
	double nEvQuoted = 149.417;
	double nEvSeen = 114.209;//e+08,115.364;//e+08;*/
	double UpsXsec_pp = 81.0e-3; // [nb]
	double UpsXsec_pp_err = 9.2e-3; // [nb]
	double sigmaPP = 42; // [mb]
	double sigmaAA = 6000; // [mb]
	//double dy = 1.0; // |y|<0.5

	//sigmaPP = 1.0; sigmaAA = 1.0;

	double Ups1SXsec_pp = 0.704*UpsXsec_pp;
	double Ups2S3SXsec_pp = 0.296*UpsXsec_pp;

	// propagate uncertainty on the ratio
	double Ups1SXsec_pp_err = 0.704*UpsXsec_pp_err;
	double Ups2S3SXsec_pp_err = 0.296*UpsXsec_pp_err;

	gSystem->cd("data");


	//TGraphAsymmErrors* rawYield_1S_2S_3S = LoadXYscanned("UpsCorrYield_STAR_1S2S3S.txt");
	TGraphAsymmErrors* rawYield_1S = LoadXYscanned("UpsCorrYield_STAR_1S.txt");
	TGraphAsymmErrors* rawYield_2S_3S = LoadXYscanned("UpsCorrYield_STAR_2S3S.txt");

	//TGraphAsymmErrors* rawYield_1S_2S_3S_int = LoadXYscanned("UpsCorrYield_STAR_1S2S3S_int.txt");
	TGraphAsymmErrors* rawYield_1S_int = LoadXYscanned("UpsCorrYield_STAR_1S_int.txt");
	TGraphAsymmErrors* rawYield_2S_3S_int = LoadXYscanned("UpsCorrYield_STAR_2S3S_int.txt");

	TGraphAsymmErrors* Ncoll = LoadXYscanned("Ncoll.txt");
	TGraphAsymmErrors* NcollInt = LoadXYscanned("NcollInt.txt"); // 392.287


	gSystem->cd("../");

	//-----------------------

	//TGraphAsymmErrors *yield_STAR_1S_2S_3S = CalculateRaa(rawYield_1S_2S_3S, Ncoll, UpsXsec_pp, UpsXsec_pp_err, sigmaPP, sigmaAA);
	TGraphAsymmErrors *yield_STAR_1S = CalculateRaa(rawYield_1S, Ncoll, Ups1SXsec_pp, Ups1SXsec_pp_err, sigmaPP, sigmaAA);
	TGraphAsymmErrors *yield_STAR_2S_3S = CalculateRaa(rawYield_2S_3S, Ncoll, Ups2S3SXsec_pp, Ups2S3SXsec_pp_err, sigmaPP, sigmaAA);

	//TGraphAsymmErrors *yield_STAR_1S_2S_3S_int = CalculateRaa(rawYield_1S_2S_3S_int, NcollInt, UpsXsec_pp, UpsXsec_pp_err, sigmaPP, sigmaAA);
	TGraphAsymmErrors *yield_STAR_1S_int = CalculateRaa(rawYield_1S_int, NcollInt, Ups1SXsec_pp, Ups1SXsec_pp_err, sigmaPP, sigmaAA);
	TGraphAsymmErrors *yield_STAR_2S_3S_int = CalculateRaa(rawYield_2S_3S_int, NcollInt, Ups2S3SXsec_pp, Ups2S3SXsec_pp_err, sigmaPP, sigmaAA);


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
	glob_xsec_STAR_2S_3S->SetLineWidth(2);
*/
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

	//----------------------------------------------

	gSystem->cd("worldData");


	gSystem->cd("../");

	//----------------------------------------------


	//----------------------------------------------

	gSystem->cd("models");


	gSystem->cd("../");


	//----------------------------------------------

	TCanvas *cnv = new TCanvas("cnv", "cnv", nPadX*padw, nPadY*padh);

	TH1D *hBack = new TH1D("hBack", "Corrected #varUpsilon yield vs. N_{part}", 400, 0.0, 400.0);

	hBack->SetTitle("");
	hBack->GetYaxis()->SetTitle("B_{ee} #frac{dN}{dy} [1]");
	hBack->GetXaxis()->SetTitle("N_{part} [1]");

	hBack->GetYaxis()->SetRangeUser( 0.0, 1.8);
	hBack->GetXaxis()->SetRangeUser( 0.0, 400.0);


	//cnv->cd()->SetLogy();	//----------------

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
	tl.SetTextSize(0.045);
	tl.SetNDC();
	//if(Preliminary) tl.DrawLatex(0.22, 0.2, "#color[2]{STAR Preliminary}");
	tl.DrawLatex(0.19, 0.8, "STAR Au+Au @ 200 GeV");
	tl.DrawLatex(0.19, 0.74, "Run14 |y|<0.5");

	TLegend *leg = new TLegend(0.52, 0.5, 0.86, 0.86);
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

	cnv->cd()->SaveAs("UpsRaaVsNpart.png");
	cnv->cd()->SaveAs("UpsRaaVsNpart.pdf");
	cnv->cd()->SaveAs("UpsRaaVsNpart.eps");

	yield_STAR_1S->SetName("UpsRaa_STAR_1S");
	yield_STAR_2S_3S->SetName("UpsRaa_STAR_2S3S");

	yield_STAR_1S_int->SetName("UpsRaa_STAR_1S");
	yield_STAR_2S_3S_int->SetName("UpsRaa_STAR_2S3S");

	WriteGraph(yield_STAR_1S);
	WriteGraph(yield_STAR_2S_3S);

	WriteGraph(yield_STAR_1S_int);
	WriteGraph(yield_STAR_2S_3S_int);

	gSystem->cd("../");


}

TGraphAsymmErrors *CalculateRaa(TGraphAsymmErrors *g, TGraphAsymmErrors *gncoll,
		double ups_xsec_pp, double ups_xsec_pp_err, double sig_pp, double sig_AA)
{

	int n = g->GetN();
	int n_ncoll = gncoll->GetN();

	TGraphAsymmErrors *raa = (TGraphAsymmErrors *)g->Clone(Form("Raa_%s", g->GetName()));

	if(n != n_ncoll)
	{
		cout<<"WARNING! Yield and <N_{coll}> graphs have different number of points. Exiting."<<endl;
		return raa;
	}

	for (int i = 0; i < n; ++i) {

		// yield
		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);
		double ey = (eh+el)/2.0;

		// N_{coll}
		double nc_x, nc_y, nc_eh, nc_el;

		gncoll->GetPoint(i, nc_x, nc_y);
		nc_eh = gncoll->GetErrorYhigh(i);
		nc_el = gncoll->GetErrorYlow(i);
		double nc_ey = (nc_eh+nc_el)/2.0;

		RooRealVar yield("yield","Corrected yield", 1.0, 0.0, 1e6);
		RooRealVar upsXsecpp("upsXsecpp","Upsilon cross section in p+p", 1.0, 0.0, 1e6);
		RooRealVar sigpp("sigpp","Inelastic p+p cross section", 1.0, 0.0, 1e6);
		RooRealVar sigAA("sigAA","Inelastic A+A cross section", 1.0, 0.0, 1e6);
		RooRealVar Ncoll("Ncoll","Number of binary collisions", 1.0, 0.0, 1e3);

		upsXsecpp.setUnit("pb");
		sigpp.setUnit("mb");
		sigAA.setUnit("mb");

		yield.setVal(y);
		upsXsecpp.setVal(ups_xsec_pp);
		sigpp.setVal(sig_pp);
		sigAA.setVal(sig_AA);
		Ncoll.setVal(nc_y);

		// set errors
		yield.setError(ey);
		upsXsecpp.setError(ups_xsec_pp_err);
		Ncoll.setError(nc_ey);
		Ncoll.setError(0.0);

		// constants
		sigpp.setConstant(kTRUE);
		sigAA.setConstant(kTRUE);
		Ncoll.setConstant(kTRUE);


		cout<<endl;
		cout<<"yield"<<"\t";
		cout<<"upsXsecpp"<<"\t";
		cout<<"sigpp"<<"\t";
		cout<<"sigAA"<<"\t";
		cout<<"Ncoll"<<"\t";
		cout<<endl;

		cout<<yield.getVal()<<"\t";
		cout<<upsXsecpp.getVal()<<"\t";
		cout<<sigpp.getVal()<<"\t";
		cout<<sigAA.getVal()<<"\t";
		cout<<Ncoll.getVal()<<"\t";
		cout<<endl;

		cout<<"yield"<<"\t";
		cout<<"upsXsecpp"<<"\t";
		cout<<"Ncoll"<<"\t";
		cout<<endl;

		cout<<yield.getError()<<"\t";
		cout<<upsXsecpp.getError()<<"\t";
		cout<<Ncoll.getError()<<"\t";
		cout<<endl;

		RooGenericPdf Raa("Raa", "Nuclear modification factor",
				"yield/upsXsecpp*(sigpp/sigAA)*(1.0/Ncoll)",RooArgSet(yield,upsXsecpp,sigpp,sigAA,Ncoll));


		//RooArgSet argSet(yield,upsXsecpp,Ncoll);
		RooArgSet argSet(yield,upsXsecpp);
		TMatrixDSym covMat(2); //3

		GetCovMatUncorr(argSet, covMat);

		double  raa_val = Raa.getVal();
		double  raa_err = getPropagatedError(Raa, argSet, covMat);

		cout<<"raa_val = "<<raa_val<<endl;
		cout<<"raa_err = "<<raa_err<<endl;

		raa->SetPoint(i, x, raa_val);
		raa->SetPointEYhigh(i, raa_err);
		raa->SetPointEYlow(i, raa_err);

	}

	return raa;
}
