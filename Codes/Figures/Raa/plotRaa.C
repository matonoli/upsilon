#include <exception>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string>
#include <fstream>

#include "TString.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveLabel.h"
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

using namespace std;
using namespace RooFit;

Int_t padw = 800, padh = 600;

const Int_t nPadX = 1;
const Int_t nPadY = 1;

const Int_t nHist = 1;

const Int_t nWorld = 13;

void CalculateInvariantYield(TGraphAsymmErrors *g);
void CalculateInvariantYieldCEM(TGraphAsymmErrors *g);
TGraphAsymmErrors* LoadCEM(const char *name = "", double scale = 1.0, int npoints = 1, int line_offset = 4);
TGraphAsymmErrors* LoadCGCNRQCD(const char *name = "", double scale = 1.0, int npoints = 1, int line_offset = 6);

/////////////////////////////////////////////////////////////
void plotPtXsec() {

	style();
	gStyle->SetOptStat(0);
	gStyle->SetOptDate(0);
	gStyle->SetLineWidth(1);
	gROOT->ForceStyle();
	//gStyle->SetHistLineColor(kBlue);

	Bool_t drawCEM = 1;
	Bool_t drawCGCNRQCD = 1;
	Bool_t Preliminary = 1;

	if(!Preliminary) gSystem->cd("data");
	if(Preliminary) gSystem->cd("dataPreliminary");

	TGraphAsymmErrors* xsec_STAR_1S_2S_3S = LoadXYscanned("Xsec_all_allRap_0mult100_cent.txt");
	TGraphAsymmErrors* xsec_STAR_1S = LoadXYscanned("Xsec_1S_allRap_0mult100_cent.txt");
	TGraphAsymmErrors* xsec_STAR_2S_3S = LoadXYscanned("Xsec_2S_3S_allRap_0mult100_cent.txt");

	TGraphAsymmErrors* syst_xsec_STAR_1S_2S_3S = LoadXYscanned("gXsec_all_allRap_0mult100_syst.txt");
	TGraphAsymmErrors* syst_xsec_STAR_1S = LoadXYscanned("gXsec_1S_allRap_0mult100_syst.txt");
	TGraphAsymmErrors* syst_xsec_STAR_2S_3S = LoadXYscanned("gXsec_2S_3S_allRap_0mult100_syst.txt");

	TGraphAsymmErrors* glob_xsec_STAR_1S_2S_3S = LoadXYscanned("gXsec_all_allRap_0mult100_glob.txt");
	TGraphAsymmErrors* glob_xsec_STAR_1S = LoadXYscanned("gXsec_1S_allRap_0mult100_glob.txt");
	TGraphAsymmErrors* glob_xsec_STAR_2S_3S = LoadXYscanned("gXsec_2S_3S_allRap_0mult100_glob.txt");

	TFile *file = new TFile("PtSpectrumFits.root", "read");

	TF1 *fit_STAR_1S_2S_3S = (TF1 *)file->Get("fitInv_Ups_1S+2S+3S");
	TF1 *fit_STAR_1S = (TF1 *)file->Get("fitInv_Ups_1S");
	TF1 *fit_STAR_2S_3S = (TF1 *)file->Get("fitInv_Ups_2S+3S");

	gSystem->cd("../");

	//-----------------------

	CalculateInvariantYield(xsec_STAR_1S_2S_3S);
	CalculateInvariantYield(xsec_STAR_1S);
	CalculateInvariantYield(xsec_STAR_2S_3S);

	CalculateInvariantYield(syst_xsec_STAR_1S_2S_3S);
	CalculateInvariantYield(syst_xsec_STAR_1S);
	CalculateInvariantYield(syst_xsec_STAR_2S_3S);

	CalculateInvariantYield(glob_xsec_STAR_1S_2S_3S);
	CalculateInvariantYield(glob_xsec_STAR_1S);
	CalculateInvariantYield(glob_xsec_STAR_2S_3S);

	//-----------------------

	xsec_STAR_1S_2S_3S->SetLineColor(kRed+1);
	xsec_STAR_1S->SetLineColor(kGreen+1);
	xsec_STAR_2S_3S->SetLineColor(kBlue+1);

	xsec_STAR_1S_2S_3S->SetMarkerColor(kRed+1);
	xsec_STAR_1S->SetMarkerColor(kGreen+1);
	xsec_STAR_2S_3S->SetMarkerColor(kBlue+1);

	xsec_STAR_1S_2S_3S->SetLineWidth(2);
	xsec_STAR_1S->SetLineWidth(2);
	xsec_STAR_2S_3S->SetLineWidth(2);


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

	xsec_STAR_1S_2S_3S->SetMarkerStyle(20);
	xsec_STAR_1S->SetMarkerStyle(33);
	xsec_STAR_2S_3S->SetMarkerStyle(21);

	xsec_STAR_1S_2S_3S->SetMarkerSize(1.4);
	xsec_STAR_1S->SetMarkerSize(2);
	xsec_STAR_2S_3S->SetMarkerSize(1.2);


	fit_STAR_1S_2S_3S->SetLineColor(kRed-7);
	fit_STAR_1S->SetLineColor(kGreen-7);
	fit_STAR_2S_3S->SetLineColor(kBlue-7);

	fit_STAR_1S_2S_3S->SetLineWidth(1);
	fit_STAR_1S->SetLineWidth(1);
	fit_STAR_2S_3S->SetLineWidth(1);

	//----------------------------------------------

	gSystem->cd("worldData");


	gSystem->cd("../");

	//----------------------------------------------


	//----------------------------------------------

	gSystem->cd("models");

	// Warning the 2S and 3S ratios are for direct Ups - not correct!
	TGraphAsymmErrors* xsec_CEM_1S = LoadCEM("ups_pp_CEM.dat", 1000*0.0238, 49);
	TGraphAsymmErrors* xsec_CEM_2S = LoadCEM("ups_pp_CEM.dat", 1000*0.0238*0.51, 49);
	TGraphAsymmErrors* xsec_CEM_3S = LoadCEM("ups_pp_CEM.dat", 1000*0.0238*0.35, 49);
	// nb->pb, B_ee, deltay?

	CalculateInvariantYieldCEM(xsec_CEM_1S);
	CalculateInvariantYieldCEM(xsec_CEM_2S);
	CalculateInvariantYieldCEM(xsec_CEM_3S);

	TGraphAsymmErrors* xsec_CEM_2S_3S = AddGraphs(xsec_CEM_2S, xsec_CEM_3S);
	TGraphAsymmErrors* xsec_CEM_1S_2S_3S = AddGraphs(xsec_CEM_1S, xsec_CEM_2S_3S);

	//CGC+NRQCD

	TGraphAsymmErrors* xsec_CGC_NRQCD_1S = LoadCGCNRQCD("ups_pp_CGC+NRQCD_1S.dat", 1.0, 98);
	TGraphAsymmErrors* xsec_CGC_NRQCD_2S = LoadCGCNRQCD("ups_pp_CGC+NRQCD_2S.dat", 1.0, 98);
	TGraphAsymmErrors* xsec_CGC_NRQCD_3S = LoadCGCNRQCD("ups_pp_CGC+NRQCD_3S.dat", 1.0, 98);

	TGraphAsymmErrors* xsec_CGC_NRQCD_2S_3S = AddGraphs(xsec_CGC_NRQCD_2S, xsec_CGC_NRQCD_3S);
	TGraphAsymmErrors* xsec_CGC_NRQCD_1S_2S_3S = AddGraphs(xsec_CGC_NRQCD_1S, xsec_CGC_NRQCD_2S_3S);

	gSystem->cd("../");

	xsec_CEM_1S_2S_3S->SetFillStyle(3013);	// 1001
	xsec_CEM_1S_2S_3S->SetLineColor(kRed+2); // kGray
	xsec_CEM_1S_2S_3S->SetFillColor(kRed+2);

	xsec_CEM_1S->SetFillStyle(1001);	// 1001
	xsec_CEM_1S->SetLineColor(kGray+1); // kGray, kGreen+2
	xsec_CEM_1S->SetFillColor(kGray+1); // kGreen+2

	xsec_CEM_2S_3S->SetFillStyle(3013);	// 1001
	xsec_CEM_2S_3S->SetLineColor(kBlue+3); // kGray
	xsec_CEM_2S_3S->SetFillColor(kBlue+3);


	xsec_CGC_NRQCD_1S->SetFillStyle(3001);
	xsec_CGC_NRQCD_1S->SetLineColor(kViolet-3); // kGreen-5
	xsec_CGC_NRQCD_1S->SetFillColor(kViolet-3);

	xsec_CGC_NRQCD_2S->SetFillStyle(3001);
	xsec_CGC_NRQCD_2S->SetLineColor(kRed-5);
	xsec_CGC_NRQCD_2S->SetFillColor(kRed-5);

	xsec_CGC_NRQCD_3S->SetFillStyle(3001);
	xsec_CGC_NRQCD_3S->SetLineColor(kBlue-5);
	xsec_CGC_NRQCD_3S->SetFillColor(kBlue-5);

	xsec_CGC_NRQCD_1S_2S_3S->SetFillStyle(3001);
	xsec_CGC_NRQCD_1S_2S_3S->SetLineColor(kRed-5);
	xsec_CGC_NRQCD_1S_2S_3S->SetFillColor(kRed-5);

	xsec_CGC_NRQCD_2S_3S->SetFillStyle(3001);
	xsec_CGC_NRQCD_2S_3S->SetLineColor(kBlue-5);
	xsec_CGC_NRQCD_2S_3S->SetFillColor(kBlue-5);



	//----------------------------------------------

	TCanvas *cnv = new TCanvas("cnv", "cnv", nPadX*padw, nPadY*padh);

	TH1D *hBack = new TH1D("hBack", "#varUpsilon invariant cross section vs. p_{T}", 1000, 0.0, 10.0);

	hBack->SetTitle("");
	hBack->GetYaxis()->SetTitle("B_{ee} #frac{1}{2#pi p_{T}} #frac{d^{2}#sigma}{dp_{T}dy} [#frac{pb}{(GeV/c)^{2}}]");
	hBack->GetXaxis()->SetTitle("p_{T} [GeV/c]");

	hBack->GetYaxis()->SetRangeUser( 5e-3, 500.0);
	hBack->GetXaxis()->SetRangeUser( 0.0, 10.0);


	cnv->cd()->SetLogy();	//----------------

	hBack->Draw();

	if (drawCEM){
		//xsec_CEM_1S_2S_3S->Draw("samee3");
		xsec_CEM_1S->Draw("samee3");
		//xsec_CEM_2S_3S->Draw("samee3");
	}

	if (drawCGCNRQCD){
		//xsec_CGC_NRQCD_1S_2S_3S->Draw("samee3");
		xsec_CGC_NRQCD_1S->Draw("samee3");
		//xsec_CGC_NRQCD_2S->Draw("samee3");
		//xsec_CGC_NRQCD_3S->Draw("samee3");
		//xsec_CGC_NRQCD_2S_3S->Draw("samee3");
	}


	//fit_STAR_1S_2S_3S->Draw("same");
	//fit_STAR_1S->Draw("same");
	//fit_STAR_2S_3S->Draw("same");

	glob_xsec_STAR_1S_2S_3S->Draw("samee2");
	glob_xsec_STAR_1S->Draw("samee2");
	glob_xsec_STAR_2S_3S->Draw("samee2");

	syst_xsec_STAR_1S_2S_3S->Draw("samee2");
	syst_xsec_STAR_1S->Draw("samee2");
	syst_xsec_STAR_2S_3S->Draw("samee2");

	xsec_STAR_1S_2S_3S->Draw("psame1");
	xsec_STAR_1S->Draw("psame1");
	xsec_STAR_2S_3S->Draw("psame1");


	//---------------


	//---------------

	TLatex tl;
	tl.SetTextFont(42);
	tl.SetTextSize(0.045);
	tl.SetNDC();
	if(Preliminary) tl.DrawLatex(0.22, 0.2, "#color[2]{STAR Preliminary}");
	tl.DrawLatex(0.19, 0.8, "STAR p+p @ 500 GeV");
	tl.DrawLatex(0.19, 0.74, "|y|<1");

	TLegend *leg = new TLegend(0.52, 0.5, 0.86, 0.86);
	TLegend *leg2 = new TLegend(0.3, 0.66, 0.52, 0.78);
	leg->AddEntry(xsec_STAR_1S_2S_3S, "STAR #varUpsilon(1S+2S+3S)", "p");
	leg->AddEntry(xsec_STAR_1S, "STAR #varUpsilon(1S)", "p");
	leg->AddEntry(xsec_STAR_2S_3S, "STAR #varUpsilon(2S+3S)", "p");

	if (drawCEM){
		//leg->AddEntry(xsec_CEM_1S_2S_3S, "CEM R. Vogt #varUpsilon(1S+2S+3S)", "f");
		leg->AddEntry(xsec_CEM_1S, "CEM R. Vogt #varUpsilon(1S)", "f");
		//leg->AddEntry(xsec_CEM_2S_3S, "CEM R. Vogt #varUpsilon(2S+3S)", "f");
	}
	if (drawCGCNRQCD){
		//leg->AddEntry(xsec_CGC_NRQCD_1S_2S_3S, "CGC+NRQCD #varUpsilon(1S+2S+3S)", "f");
		leg->AddEntry(xsec_CGC_NRQCD_1S, "CGC+NRQCD #varUpsilon(1S)", "f");
		//leg->AddEntry(xsec_CGC_NRQCD_2S_3S, "CGC+NRQCD #varUpsilon(2S+3S)", "f");
	}
	//leg->AddEntry(xsec_CGC_NRQCD_2S, "CGC+NRQCD #varUpsilon(2S)", "f");
	//leg->AddEntry(xsec_CGC_NRQCD_3S, "CGC+NRQCD #varUpsilon(3S)", "f");
	leg2->AddEntry(syst_xsec_STAR_1S_2S_3S, "Uncorrelated syst.", "f");
	leg2->AddEntry(glob_xsec_STAR_1S_2S_3S, "Correlated syst.", "f");

	leg->SetFillColor(kWhite);
	leg->SetFillStyle(1000);
	leg2->SetFillColor(kWhite);
	leg2->SetFillStyle(1000);

	leg->Draw("same");
	leg2->Draw("same");

	//------------------------------------

	TString mkdir = "output";
	if(Preliminary) mkdir = "outputPreliminary";
	gSystem->MakeDirectory(mkdir.Data());
	gSystem->cd(mkdir.Data());

	if (!drawCEM && !drawCGCNRQCD){
		cnv->cd()->SaveAs("UpsInvXsecPt.png");
		cnv->cd()->SaveAs("UpsInvXsecPt.pdf");
		cnv->cd()->SaveAs("UpsInvXsecPt.eps");
	}

	if (drawCEM && drawCGCNRQCD){
		cnv->cd()->SaveAs("UpsInvXsecPt_models.png");
		cnv->cd()->SaveAs("UpsInvXsecPt_models.pdf");
		cnv->cd()->SaveAs("UpsInvXsecPt_models.eps");
	}

	if (drawCEM && !drawCGCNRQCD){
		cnv->cd()->SaveAs("UpsInvXsecPt_CEM.png");
		cnv->cd()->SaveAs("UpsInvXsecPt_CEM.pdf");
		cnv->cd()->SaveAs("UpsInvXsecPt_CEM.eps");
	}

	if (drawCGCNRQCD && !drawCEM){
		cnv->cd()->SaveAs("UpsInvXsecPt_CGCNRQCD.png");
		cnv->cd()->SaveAs("UpsInvXsecPt_CGCNRQCD.pdf");
		cnv->cd()->SaveAs("UpsInvXsecPt_CGCNRQCD.eps");
	}

	gSystem->cd("../");


}

void CalculateInvariantYield(TGraphAsymmErrors *g)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		y /= 2.0*TMath::Pi()*x;
		y /= 2.0; // rapidity bin width
		y /= 2.0; // pt bin width

		eh /= 2.0*TMath::Pi()*x;
		eh /= 2.0; // rapidity bin width
		eh /= 2.0; // pt bin width

		el /= 2.0*TMath::Pi()*x;
		el /= 2.0; // rapidity bin width
		el /= 2.0; // pt bin width

		g->SetPoint(i, x, y);
		g->SetPointEYhigh(i, eh);
		g->SetPointEYlow(i, el);

	}

}

void CalculateInvariantYieldCEM(TGraphAsymmErrors *g)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		y /= 2.0*TMath::Pi()*x;
		//y /= 2.0; // rapidity bin width

		eh /= 2.0*TMath::Pi()*x;
		//eh /= 2.0; // rapidity bin width

		el /= 2.0*TMath::Pi()*x;
		//el /= 2.0; // rapidity bin width

		g->SetPoint(i, x, y);
		g->SetPointEYhigh(i, eh);
		g->SetPointEYlow(i, el);

	}

}


TGraphAsymmErrors* LoadCEM(const char *name, double scale, int npoints, int line_offset)
{

	std::ifstream file;
	string string;

	int n = 0;
	const int nmax = 1000;

	double array_x[nmax];
	double array_y[nmax];
	double array_dx_minus[nmax];
	double array_dx_plus[nmax];
	double array_dy_minus[nmax];
	double array_dy_plus[nmax];

	file.open(name, std::ifstream::in);

	if(!file.is_open()) {
		cout<<"ERROR : No input file!"<<endl;
		cout<<name<<endl;
		return 0;
	}

	cout<<endl;

	for (int i = 0; i < line_offset; ++i) {
		std::getline(file, string);
		cout<<string;
	}

	cout<<endl;

	while(!file.eof())	{

		if(n>=nmax)	{
			cout<<"WARNING!! Maximum number of points reached!"<<endl;
			break;
		}

		std::getline(file, string);
		//cout<<string<<endl;

		if (string=="# EoF") break;
		if (string=="\043 EoF") break;
		if (n == npoints) break;

		//cout<<endl;

		//cout<<string.data()<<endl;

		double x = 0.0;
		double y = 0.0;
		double dx_hi = 0.0;
		double dx_lo = 0.0;
		double dy_hi = 0.0;
		double dy_lo = 0.0;

		//sscanf(string.data(), "%f%f%f%f%f%f", array_x[x], array_y[x], array_dx_minus[x], array_dx_plus[x], array_dy_minus[x], array_dy_plus[x]);
		//sscanf(string.data(), "%f\t%f\t%f\t%f\t%f\t%f", &array_x[x], &array_y[x], &array_dx_minus[x], &array_dx_plus[x], &array_dy_minus[x], &array_dy_plus[x]);
		sscanf(string.data(), "%lf\t%lf\t%lf\t%lf\t", &x, &y, &dy_hi, &dy_lo);

/*
		cout<<x<<"\t";
		cout<<y<<"\t";
		cout<<dx_lo<<"\t";
		cout<<dx_hi<<"\t";
		cout<<dy_lo<<"\t";
		cout<<dy_hi;

		cout<<endl;
*/
		dy_lo = y-dy_lo;
		dy_hi = dy_hi-y;

		y *= scale;
		dy_lo *= scale;
		dy_hi *= scale;

		array_x[n] = x;
		array_y[n] = y;
		array_dy_minus[n] = dy_lo;
		array_dy_plus[n] = dy_hi;

		cout<<array_x[n]<<"\t";
		cout<<array_y[n]<<"\t";
		cout<<array_dy_minus[n]<<"\t";
		cout<<array_dy_plus[n];

		cout<<endl;
		++n;
	}

	TGraphAsymmErrors *graph = new TGraphAsymmErrors(n);

	for (int i = 0; i < n; ++i) {


		graph->SetPoint(i, array_x[i], array_y[i]);
		graph->SetPointEYlow(i, array_dy_minus[i]);
		graph->SetPointEYhigh(i, array_dy_plus[i]);

	}

	return graph;
}


TGraphAsymmErrors* LoadCGCNRQCD(const char *name, double scale, int npoints, int line_offset)
{

	std::ifstream file;
	string string;

	int n = 0;
	const int nmax = 1000;

	double array_x[nmax];
	double array_y[nmax];
	double array_dx_minus[nmax];
	double array_dx_plus[nmax];
	double array_dy_minus[nmax];
	double array_dy_plus[nmax];

	file.open(name, std::ifstream::in);

	if(!file.is_open()) {
		cout<<"ERROR : No input file!"<<endl;
		cout<<name<<endl;
		return 0;
	}

	cout<<endl;

	for (int i = 0; i < line_offset; ++i) {
		std::getline(file, string);
		cout<<string;
	}

	cout<<endl;

	while(!file.eof())	{

		if(n>=nmax)	{
			cout<<"WARNING!! Maximum number of points reached!"<<endl;
			break;
		}

		std::getline(file, string);
		//cout<<string<<endl;

		if (string=="# EoF") break;
		if (string=="\043 EoF") break;
		if (n == npoints) break;

		//cout<<endl;

		//cout<<string.data()<<endl;

		double x = 0.0;
		double y = 0.0;
		double dx_hi = 0.0;
		double dx_lo = 0.0;
		double dy_hi = 0.0;
		double dy_lo = 0.0;
		double y_hi = 0.0;
		double y_lo = 0.0;

		//sscanf(string.data(), "%f%f%f%f%f%f", array_x[x], array_y[x], array_dx_minus[x], array_dx_plus[x], array_dy_minus[x], array_dy_plus[x]);
		//sscanf(string.data(), "%f\t%f\t%f\t%f\t%f\t%f", &array_x[x], &array_y[x], &array_dx_minus[x], &array_dx_plus[x], &array_dy_minus[x], &array_dy_plus[x]);
		sscanf(string.data(), "%lf\t%lf\t%lf\t", &x, &y_lo, &y_hi);

/*
		cout<<x<<"\t";
		cout<<y<<"\t";
		cout<<dx_lo<<"\t";
		cout<<dx_hi<<"\t";
		cout<<dy_lo<<"\t";
		cout<<dy_hi;

		cout<<endl;
*/
		y = 0.5*(y_lo+y_hi);

		dy_lo = y-y_lo;
		dy_hi = y_hi-y;

		y *= scale;
		dy_lo *= scale;
		dy_hi *= scale;

		array_x[n] = x;
		array_y[n] = y;
		array_dy_minus[n] = dy_lo;
		array_dy_plus[n] = dy_hi;

		cout<<array_x[n]<<"\t";
		cout<<array_y[n]<<"\t";
		cout<<array_dy_minus[n]<<"\t";
		cout<<array_dy_plus[n];

		cout<<endl;
		++n;
	}

	TGraphAsymmErrors *graph = new TGraphAsymmErrors(n);

	for (int i = 0; i < n; ++i) {


		graph->SetPoint(i, array_x[i], array_y[i]);
		graph->SetPointEYlow(i, array_dy_minus[i]);
		graph->SetPointEYhigh(i, array_dy_plus[i]);

	}

	return graph;
}



