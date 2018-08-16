#include <iostream>
#include <fstream>
#include <cmath>

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

   		//producing graph now
   		TGraph* y1S = new TGraph()
   	}


}