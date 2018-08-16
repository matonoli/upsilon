#include <iostream>
#include <fstream>
#include <cmath>
#include "../Codes/RaaCalc/GraphOperations.h"

using namespace std;

void makeXYfiles(const Char_t *inputFile="yields.txt", Int_t option = 0) {

	fstream fyields(inputFile);
	if (option == 0 || 1) {
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


	float npart09 = 126.568;
	float npart29 = 160.961;
	float npart05 = (115.34+76.499+47.489+27.338+14.269)/5;
	float npart25 = (115.34+76.499+47.489)/3;
	float npart57 = (235.473+167.04)/2;
	float npart79 = 324.662;
	float shift = 0;
	Double_t npart[] = {npart29+shift, npart25+shift, npart57+shift, npart79+shift};
	Double_t zero[] = {0, 0, 0, 0};

   	//producing graph now
   	TGraphAsymmErrors* p1Sint 		= new TGraphAsymmErrors(1,npart,yield1S,zero,zero,err1S,err1S);
   	TGraphAsymmErrors* p1S 	  		= new TGraphAsymmErrors(3,&npart[1],&yield1S[1],zero,zero,&err1S[1],&err1S[1]);
   	TGraphAsymmErrors* p23Sint 		= new TGraphAsymmErrors(1,npart,yield23S,zero,zero,err23S,err23S);
   	TGraphAsymmErrors* p23S 	  	= new TGraphAsymmErrors(3,&npart[1],&yield23S[1],zero,zero,&err23S[1],&err23S[1]);

   	//produing XY files now
   	p1Sint->SetName("yield_1S_int");
   	p1Sint->SetName("yield_1S");
   	p1Sint->SetName("yield_2S3S_int");
   	p1Sint->SetName("yield_2S3S");
   	
   	}


}