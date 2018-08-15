/*
 * LoadXYscanned.h
 *
 *  Created on: 13 pa� 2014
 *      Author: Khaless
 */

#ifndef LOADXYSCANNED_H_
#define LOADXYSCANNED_H_

#include <iostream>
//#include <sstream>
#include <fstream>
//#include <stdio.h>

#include <TGraphAsymmErrors.h>

//TGraphAsymmErrors* LoadXYscanned(const char *name, int line_offset);

TGraphAsymmErrors* LoadXYscanned(const char *name = "", int line_offset = 6)	{

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
		sscanf(string.data(), "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", &x, &y, &dx_lo, &dx_hi, &dy_lo, &dy_hi);

/*
		cout<<x<<"\t";
		cout<<y<<"\t";
		cout<<dx_lo<<"\t";
		cout<<dx_hi<<"\t";
		cout<<dy_lo<<"\t";
		cout<<dy_hi;

		cout<<endl;
*/
		array_x[n] = x;
		array_y[n] = y;
		array_dx_minus[n] = dx_lo;
		array_dx_plus[n] = dx_hi;
		array_dy_minus[n] = dy_lo;
		array_dy_plus[n] = dy_hi;

		cout<<array_x[n]<<"\t";
		cout<<array_y[n]<<"\t";
		cout<<array_dx_minus[n]<<"\t";
		cout<<array_dx_plus[n]<<"\t";
		cout<<array_dy_minus[n]<<"\t";
		cout<<array_dy_plus[n];

		cout<<endl;
		++n;
	}

	TGraphAsymmErrors *graph = new TGraphAsymmErrors(n);

	for (int i = 0; i < n; ++i) {


		graph->SetPoint(i, array_x[i], array_y[i]);
		graph->SetPointEXlow(i, array_dx_minus[i]);
		graph->SetPointEXhigh(i, array_dx_plus[i]);
		graph->SetPointEYlow(i, array_dy_minus[i]);
		graph->SetPointEYhigh(i, array_dy_plus[i]);

	}

	return graph;
}


#endif /* LOADXYSCANNED_H_ */
