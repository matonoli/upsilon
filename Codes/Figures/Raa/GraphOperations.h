/*
 * GraphOperations.h
 *
 *  Created on: 3 lis 2017
 *      Author: Khaless
 */

#ifndef GRAPHOPERATIONS_H_
#define GRAPHOPERATIONS_H_

#include <iostream>
#include <fstream>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>


void MultiplyGraph(TGraphAsymmErrors *g, double scale = 1.0);
void MultiplyGraph(TGraphErrors *g, double scale = 1.0);
void MultiplyGraphX(TGraphAsymmErrors *g, double scale = 1.0);
void MultiplyGraphX(TGraphErrors *g, double scale = 1.0);
TGraphAsymmErrors* AddGraphs(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2);
TGraphErrors* AddGraphs(TGraphErrors* gr1, TGraphErrors* gr2);
TGraphAsymmErrors* SubtractGraphs(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2);
TGraphErrors* SubtractGraphs(TGraphErrors* gr1, TGraphErrors* gr2);
TGraphAsymmErrors* AddErrors(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2);
void SetEx(TGraphAsymmErrors *g, double ex = 0.0);

void printGraph(TGraphAsymmErrors *gr);
void printGraphLatex(TGraphAsymmErrors *gr);

void writeGraph(TGraphAsymmErrors *gr);

void MultiplyGraph(TGraphAsymmErrors *g, double scale)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		y *= scale;
		eh *= scale;
		el *= scale;

		g->SetPoint(i, x, y);
		g->SetPointEYhigh(i, eh);
		g->SetPointEYlow(i, el);

	}

}

void MultiplyGraph(TGraphErrors *g, double scale)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		y *= scale;
		eh *= scale;
		el *= scale;

		g->SetPoint(i, x, y);
		g->SetPointError(i, 0.0, eh);
		//g->SetPointEYlow(i, el);

	}

}

void MultiplyGraphX(TGraphAsymmErrors *g, double scale)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		x *= scale;

		g->SetPoint(i, x, y);
		g->SetPointEYhigh(i, eh);
		g->SetPointEYlow(i, el);

	}

}

void MultiplyGraphX(TGraphErrors *g, double scale)
{

	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, eh, el;

		g->GetPoint(i, x, y);
		eh = g->GetErrorYhigh(i);
		el = g->GetErrorYlow(i);

		x *= scale;

		g->SetPoint(i, x, y);
		g->SetPointError(i, 0.0, eh);
		//g->SetPointEYlow(i, el);

	}

}

TGraphAsymmErrors* AddGraphs(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2)
{

	double x = 0.0;
	double y = 0.0;
	double dx_hi = 0.0;
	double dx_lo = 0.0;
	double dy_hi = 0.0;
	double dy_lo = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;
	double dx_hi2 = 0.0;
	double dx_lo2 = 0.0;
	double dy_hi2 = 0.0;
	double dy_lo2 = 0.0;

	int n = gr1->GetN();
	//TGraphAsymmErrors* gr = new TGraphAsymmErrors(*gr1);
	TGraphAsymmErrors* gr = (TGraphAsymmErrors*)gr1->Clone();

	for (int i = 0; i < n; ++i) {

		gr1->GetPoint(i, x, y);
		gr2->GetPoint(i, x2, y2);

		dy_hi = gr1->GetErrorYhigh(i);
		dy_lo = gr1->GetErrorYlow(i);
		dy_hi2 = gr2->GetErrorYhigh(i);
		dy_lo2 = gr2->GetErrorYlow(i);

		gr->SetPoint(i, x, y+y2);
		gr->SetPointEYhigh(i, sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2));
		gr->SetPointEYlow(i, sqrt(dy_lo*dy_lo+dy_lo2*dy_lo2));

	}


	return gr;
}

TGraphErrors* AddGraphs(TGraphErrors* gr1, TGraphErrors* gr2)
{

	double x = 0.0;
	double y = 0.0;
	double dx_hi = 0.0;
	double dx_lo = 0.0;
	double dy_hi = 0.0;
	double dy_lo = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;
	double dx_hi2 = 0.0;
	double dx_lo2 = 0.0;
	double dy_hi2 = 0.0;
	double dy_lo2 = 0.0;

	int n = gr1->GetN();
	TGraphErrors* gr = (TGraphErrors*)gr1->Clone();

	for (int i = 0; i < n; ++i) {

		gr1->GetPoint(i, x, y);
		gr2->GetPoint(i, x2, y2);

		dy_hi = gr1->GetErrorYhigh(i);
		dy_lo = gr1->GetErrorYlow(i);
		dy_hi2 = gr2->GetErrorYhigh(i);
		dy_lo2 = gr2->GetErrorYlow(i);

		gr->SetPoint(i, x, y+y2);
		gr->SetPointError(i, 0.0, sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2));

	}


	return gr;
}

TGraphAsymmErrors* SubtractGraphs(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2)
{

	double x = 0.0;
	double y = 0.0;
	double dx_hi = 0.0;
	double dx_lo = 0.0;
	double dy_hi = 0.0;
	double dy_lo = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;
	double dx_hi2 = 0.0;
	double dx_lo2 = 0.0;
	double dy_hi2 = 0.0;
	double dy_lo2 = 0.0;

	int n = gr1->GetN();
	//TGraphAsymmErrors* gr = new TGraphAsymmErrors(*gr1);
	TGraphAsymmErrors* gr = (TGraphAsymmErrors*)gr1->Clone();

	for (int i = 0; i < n; ++i) {

		gr1->GetPoint(i, x, y);
		gr2->GetPoint(i, x2, y2);

		dy_hi = gr1->GetErrorYhigh(i);
		dy_lo = gr1->GetErrorYlow(i);
		dy_hi2 = gr2->GetErrorYhigh(i);
		dy_lo2 = gr2->GetErrorYlow(i);

		gr->SetPoint(i, x, y-y2);
		gr->SetPointEYhigh(i, sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2));
		gr->SetPointEYlow(i, sqrt(dy_lo*dy_lo+dy_lo2*dy_lo2));

	}


	return gr;
}


TGraphErrors* SubtractGraphs(TGraphErrors* gr1, TGraphErrors* gr2)
{

	double x = 0.0;
	double y = 0.0;
	double dx_hi = 0.0;
	double dx_lo = 0.0;
	double dy_hi = 0.0;
	double dy_lo = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;
	double dx_hi2 = 0.0;
	double dx_lo2 = 0.0;
	double dy_hi2 = 0.0;
	double dy_lo2 = 0.0;

	int n = gr1->GetN();
	//TGraphAsymmErrors* gr = new TGraphAsymmErrors(*gr1);
	TGraphErrors* gr = (TGraphErrors*)gr1->Clone();

	for (int i = 0; i < n; ++i) {

		gr1->GetPoint(i, x, y);
		gr2->GetPoint(i, x2, y2);

		dy_hi = gr1->GetErrorYhigh(i);
		dy_lo = gr1->GetErrorYlow(i);
		dy_hi2 = gr2->GetErrorYhigh(i);
		dy_lo2 = gr2->GetErrorYlow(i);

		gr->SetPoint(i, x, y-y2);
		gr->SetPointError(i, 0.0, sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2));

	}


	return gr;
}


TGraphAsymmErrors* AddErrors(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2)
{

	double x = 0.0;
	double y = 0.0;
	double dx_hi = 0.0;
	double dx_lo = 0.0;
	double dy_hi = 0.0;
	double dy_lo = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;
	double dx_hi2 = 0.0;
	double dx_lo2 = 0.0;
	double dy_hi2 = 0.0;
	double dy_lo2 = 0.0;

	int n = gr1->GetN();
	//TGraphAsymmErrors* gr = new TGraphAsymmErrors(*gr1);
	TGraphAsymmErrors* gr = (TGraphAsymmErrors*)gr1->Clone();

	for (int i = 0; i < n; ++i) {

		gr1->GetPoint(i, x, y);
		gr2->GetPoint(i, x2, y2);

		dy_hi = gr1->GetErrorYhigh(i);
		dy_lo = gr1->GetErrorYlow(i);
		dy_hi2 = gr2->GetErrorYhigh(i);
		dy_lo2 = gr2->GetErrorYlow(i);

		//cout<<"dy_hi = "<<dy_hi<<"\t dy_hi2 = "<<dy_hi2<<"\t dy_hi sum = "<<sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2)<<endl;

		gr->SetPoint(i, x, y);
		gr->SetPointEYhigh(i, sqrt(dy_hi*dy_hi+dy_hi2*dy_hi2));
		gr->SetPointEYlow(i, sqrt(dy_lo*dy_lo+dy_lo2*dy_lo2));

	}


	return gr;
}

void SetEx(TGraphAsymmErrors *g, double ex)
{
	int n = g->GetN();

	for (int i = 0; i < n; ++i) {

		double x, y, exh, exl;

		g->GetPoint(i, x, y);
		//exh = g->GetErrorXhigh(i);
		//exl = g->GetErrorXlow(i);


		g->SetPoint(i, x, y);
		g->SetPointEXhigh(i, ex);
		g->SetPointEXlow(i, ex);

	}

}


void printGraph(TGraphAsymmErrors *graph)
{

	int np = graph->GetN();

	//cout<<"np = "<<np<<endl;
	cout<<endl;
	cout<<graph->GetName()<<endl;

	for (int i = 0; i < np; ++i) {

		double x, y, exh, exl, eyh, eyl;

		graph->GetPoint(i, x, y);

		exh = graph->GetErrorXhigh(i);
		exl = graph->GetErrorXlow(i);

		eyh = graph->GetErrorYhigh(i);
		eyl = graph->GetErrorYlow(i);

		cout<<x<<" "<<y<<" "<<" "<<eyl<<" "<<eyh<<endl;
	}

	cout<<endl;

}


void printGraphLatex(TGraphAsymmErrors *graph)
{

	int np = graph->GetN();

	//cout<<"np = "<<np<<endl;
	cout<<endl;
	cout<<graph->GetName()<<endl;

	for (int i = 0; i < np; ++i) {

		double x, y, exh, exl, eyh, eyl;

		graph->GetPoint(i, x, y);

		exh = graph->GetErrorXhigh(i);
		exl = graph->GetErrorXlow(i);

		eyh = graph->GetErrorYhigh(i);
		eyl = graph->GetErrorYlow(i);

		//cout<<x<<" "<<y<<" "<<" "<<eyl<<" "<<eyh<<endl;
		cout<<" & $"<<Form("%.1f",eyh*100)<<"\\%$";
	}

	cout<<" \\\\";
	cout<<endl;

	for (int i = 0; i < np; ++i) {

		double x, y, exh, exl, eyh, eyl;

		graph->GetPoint(i, x, y);

		exh = graph->GetErrorXhigh(i);
		exl = graph->GetErrorXlow(i);

		eyh = graph->GetErrorYhigh(i);
		eyl = -graph->GetErrorYlow(i);

		//cout<<x<<" "<<y<<" "<<" "<<eyl<<" "<<eyh<<endl;
		cout<<" & $"<<Form("%.1f",eyl*100)<<"\\%$";
	}

	cout<<" \\\\";
	cout<<endl;

}


void writeGraph(TGraphAsymmErrors *graph)
{

	int np = graph->GetN();

	//cout<<"np = "<<np<<endl;

	std::ofstream file;
	file.open(Form("./output/%s_syst.txt", graph->GetName()));

	file<<"# xyscan format data"<<endl;
	file<<"# Date:"<<endl;
	file<<"# Scanned by:"<<endl;
	file<<"# Source:"<<endl;
	file<<"# Comment:"<<endl;
	file<<"# Format: x y -dx +dx -dy +dy"<<endl;

	for (int i = 0; i < np; ++i) {

		double x, y, exh, exl, eyh, eyl;

		graph->GetPoint(i, x, y);

		exh = graph->GetErrorXhigh(i);
		exl = graph->GetErrorXlow(i);

		eyh = graph->GetErrorYhigh(i);
		eyl = graph->GetErrorYlow(i);

		file<<x<<" "<<y<<" "<<exl<<" "<<exh<<" "<<eyl<<" "<<eyh<<endl;

	}
	file<<"# EoF"<<endl;
	file.close();

}

#endif /* GRAPHOPERATIONS_H_ */
