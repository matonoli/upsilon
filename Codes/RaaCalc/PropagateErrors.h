/*
 * PropagateErrors.h
 *
 *  Created on: 19 sty 2017
 *      Author: Khaless
 */

#ifndef PROPAGATEERRORS_H_
#define PROPAGATEERRORS_H_

//#include <vector>

#include "TList.h"
#include "TVectorD.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGlobalFunc.h"
#include "TMatrixDSym.h"
#include "RooGenericPdf.h"
#include "RooArgList.h"


using namespace RooFit;
using namespace std;

Double_t getPropagatedError(RooGenericPdf &f, RooArgSet parSet, TMatrixDSym covMat);
Double_t getPropagatedErrorUncorr(RooGenericPdf &f, RooArgSet parSet);

int GetCovMat(RooArgSet parSet, TMatrixDSym &covMat);
int GetCovMatUncorr(RooArgSet parSet, TMatrixDSym &covMat);

// Warning!! Use covariance matrix with parameters in the same order as in parSet!
Double_t getPropagatedError(RooGenericPdf &f, RooArgSet parSet, TMatrixDSym covMat)
{

	RooAbsReal* cloneFunc = (RooAbsReal*) f.cloneTree() ;
	//RooArgSet* errorParams = cloneFunc->getObservables(parSet);
	RooArgSet* nset = cloneFunc->getObservables(parSet);
	//RooArgSet* nset = cloneFunc->getParameters(*errorParams);
	RooArgList* paramList = new RooArgList(*nset);

	Double_t *plusVar = new Double_t[paramList->getSize()];
	Double_t *minusVar = new Double_t[paramList->getSize()];
/*
    cout<<endl;

    //errorParams->Print();
	nset->Print();
	paramList->Print();

    cout<<endl;
*/
	TMatrixDSym V(covMat);
	//V.Print();
	for (Int_t ivar=0 ; ivar<paramList->getSize() ; ivar++) {

	    RooRealVar *rrv = (RooRealVar*)paramList->at(ivar);

	    Double_t cenVal = rrv->getVal();
	    Double_t errVal = sqrt(V(ivar,ivar));
/*
	    cout<<"ivar = "<<ivar<<"\t";
	    cout<<"cenVal = "<<cenVal<<"\t";
	    cout<<"errVal = "<<errVal;
	    cout<<endl;
*/
	    // Make Plus variation
	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal+errVal);
	    //plusVar[ivar] = ((RooAbsReal*)cloneFunc->getVal(nset));
	    plusVar[ivar] = cloneFunc->getVal();

	    // Make Minus variation
	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal-errVal);
	    //minusVar[ivar] =  ((RooAbsReal*)cloneFunc->getVal(nset));
	    minusVar[ivar] =  cloneFunc->getVal();

	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal);

	    //cout<<"plusVar = "<<plusVar[ivar]<<endl;
	    //cout<<"minusVar = "<<minusVar[ivar]<<endl;
	}	// paramList loop


	TMatrixDSym C(paramList->getSize());
	//TList errVec() ;
	for (int i=0 ; i<paramList->getSize(); i++) {
		//errVec[i] = sqrt(V(i,i)) ;
	    for (int j=i ; j<paramList->getSize(); j++) {
	      C(i,j) = V(i,j)/sqrt(V(i,i)*V(j,j));
	      C(j,i) = C(i,j);
	    }
	}

	// Make vector of variations
	TVectorD F(paramList->getSize());
	for (unsigned int j=0 ; j<paramList->getSize(); j++) {
		F[j] = (plusVar[j]-minusVar[j])/2;
	}

	//C.Print();
	//F.Print();

	// Calculate error in linear approximation from variations and correlation coefficient
	Double_t sum = F*(C*F) ;

	delete cloneFunc ;
	//delete errorParams ;
	delete nset ;

	delete [] plusVar;
	delete [] minusVar;

	return sqrt(sum) ;

	return 0.0;

}

// Warning!! Use covariance matrix with parameters in the same order as in parSet!
Double_t getPropagatedErrorUncorr(RooGenericPdf &f, RooArgSet parSet)
{

	RooAbsReal* cloneFunc = (RooAbsReal*) f.cloneTree() ;
	//RooArgSet* errorParams = cloneFunc->getObservables(parSet);
	RooArgSet* nset = cloneFunc->getObservables(parSet);
	//RooArgSet* nset = cloneFunc->getParameters(*errorParams);
	RooArgList* paramList = new RooArgList(*nset);

	Double_t *plusVar = new Double_t[paramList->getSize()];
	Double_t *minusVar = new Double_t[paramList->getSize()];
/*
    cout<<endl;

    //errorParams->Print();
	nset->Print();
	paramList->Print();

    cout<<endl;
*/
	TMatrixDSym covMat(paramList->getSize());
	for (Int_t i=0 ; i<paramList->getSize() ; i++) {

	    RooRealVar *rrv = (RooRealVar*)paramList->at(i);

	    covMat(i,i) = rrv->getError()*rrv->getError();
	}

	TMatrixDSym V(covMat);
	//V.Print();
	for (Int_t ivar=0 ; ivar<paramList->getSize() ; ivar++) {

	    RooRealVar *rrv = (RooRealVar*)paramList->at(ivar);

	    Double_t cenVal = rrv->getVal();
	    Double_t errVal = sqrt(V(ivar,ivar));
/*
	    cout<<"ivar = "<<ivar<<"\t";
	    cout<<"cenVal = "<<cenVal<<"\t";
	    cout<<"errVal = "<<errVal;
	    cout<<endl;
*/
	    // Make Plus variation
	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal+errVal);
	    //plusVar[ivar] = ((RooAbsReal*)cloneFunc->getVal(nset));
	    plusVar[ivar] = cloneFunc->getVal();

	    // Make Minus variation
	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal-errVal);
	    //minusVar[ivar] =  ((RooAbsReal*)cloneFunc->getVal(nset));
	    minusVar[ivar] =  cloneFunc->getVal();

	    ((RooRealVar*)paramList->at(ivar))->setVal(cenVal);

	    //cout<<"plusVar = "<<plusVar[ivar]<<endl;
	    //cout<<"minusVar = "<<minusVar[ivar]<<endl;
	}	// paramList loop


	TMatrixDSym C(paramList->getSize());
	//TList errVec() ;
	for (int i=0 ; i<paramList->getSize(); i++) {
		//errVec[i] = sqrt(V(i,i)) ;
	    for (int k=i ; k<paramList->getSize(); k++) {
	      C(i,k) = V(i,k)/sqrt(V(i,i)*V(k,k));
	      C(k,i) = C(i,k);
	    }
	}

	// Make vector of variations
	TVectorD F(paramList->getSize());
	for (unsigned int j=0 ; j<paramList->getSize(); j++) {
		F[j] = (plusVar[j]-minusVar[j])/2;
	}

	//C.Print();
	//F.Print();

	// Calculate error in linear approximation from variations and correlation coefficient
	Double_t sum = F*(C*F) ;

	delete cloneFunc ;
	//delete errorParams ;
	delete nset ;

	delete [] plusVar;
	delete [] minusVar;

	return sqrt(sum) ;

	return 0.0;

}


int GetCovMat(RooArgSet parSet, TMatrixDSym &covMat)
{

	RooArgList* paramList = new RooArgList(parSet);

	if (covMat.GetNrows() != covMat.GetNcols())
		{

			cout<<"ERROR! Matrix not square!"<<endl;
			return 0;
		}
	if (paramList->getSize() != covMat.GetNcols())
		{

			cout<<"ERROR! Matrix dimension differs from parSet!"<<endl;
			return 0;
		}

	for (int i=0 ; i<paramList->getSize() ; i++) {
	    RooRealVar *rrv_i = (RooRealVar*)paramList->at(i);

		for (int j=0 ; j<paramList->getSize() ; j++) {
		    RooRealVar *rrv_j = (RooRealVar*)paramList->at(j);

		    covMat(i,j) = rrv_i->getError()*rrv_j->getError();
		}
	}

	//covMat.Print();

	return 1;

}


int GetCovMatUncorr(RooArgSet parSet, TMatrixDSym &covMat)
{

	RooArgList* paramList = new RooArgList(parSet);

	if (covMat.GetNrows() != covMat.GetNcols())
		{

			cout<<"ERROR! Matrix not square!"<<endl;
			return 0;
		}
	if (paramList->getSize() != covMat.GetNcols())
		{

			cout<<"ERROR! Matrix dimension differs from parSet!"<<endl;
			return 0;
		}

	for (int i=0 ; i<paramList->getSize() ; i++) {
	    RooRealVar *rrv_i = (RooRealVar*)paramList->at(i);

		for (int j=0 ; j<paramList->getSize() ; j++) {
		    RooRealVar *rrv_j = (RooRealVar*)paramList->at(j);

		    if(i!=j) continue;

		    covMat(i,j) = rrv_i->getError()*rrv_j->getError();
		}
	}

	//covMat.Print();

	return 1;

}

#endif /* PROPAGATEERRORS_H_ */
