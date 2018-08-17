#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"

#include <vector>

using namespace std;

const float elmass = 0.000510998910;

void readUpsTree(const Char_t *inputFile="test.root", Int_t cutSet = 0) {
	
	TFile* fin = new TFile(inputFile,"READ");
	//TFile *fsim = new TFile("simBB_0419_2eq6.root","READ");
	//TFile *fsim = new TFile("simBB_0405_2e.root","READ");
	TFile *fsim = new TFile("../upstree/simBB_0505_2e_m.root","READ");
	//TFile *fsim = new TFile("simBB_0405_2eq20.root","READ");
	//TFile *fsim = new TFile("simBB_0420_DYq6.root","READ");
	//TFile *fsim = new TFile("simBB_0409_2eq15.root","READ");


	TTree* upsTree = (TTree*)fin->Get("upsTree");
	TTree* bbbar = (TTree*)fsim->Get("bbbar");

	Float_t fMass = 0;
	Int_t fCentrality = 0;
	Int_t fSign = 0;
	TTree* massTree = new TTree("massTree", "massTree");
	massTree->Branch("fMass", &fMass, "fMass/F");
	massTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	massTree->Branch("fSign", &fSign, "fSign/I");

	fin->cd();
	int nEntries = upsTree->GetEntries();
	cout << "upsTree has " << nEntries << " entries " << endl;

	const int maxTracks = 100;
    Int_t tEventId;
    Int_t tEventRunId;
    Float_t tEventTpcVx;
    Float_t tEventTpcVy;
    Float_t tEventTpcVz;
    Float_t tEventVpdVz;
    Int_t tEventCent9;
    Int_t tEventCent16;
    Int_t tEventRefMult;
    Int_t tEventRefMultC;
    Int_t tNElectrons;
    Float_t tElePt[maxTracks];
    Float_t tEleP[maxTracks];
    Float_t tElePx[maxTracks];
    Float_t tElePy[maxTracks];
    Float_t tElePz[maxTracks];
    Float_t tEleEta[maxTracks];
    Int_t tEleCharge[maxTracks];
    Int_t tEleNHits[maxTracks];
    Float_t tEleNRat[maxTracks];
    Int_t tEleNDedx[maxTracks];
    Float_t tEleDca[maxTracks];
    Float_t tEleDedx[maxTracks];
    Float_t tEleNSigE[maxTracks];
    Float_t tEleE[maxTracks];
    Float_t tEleEcl[maxTracks];
    Float_t tEleZDist[maxTracks];
    Float_t tElePhiDist[maxTracks];
    Float_t tEleR[maxTracks];
    Int_t tEleTrig[maxTracks];

	upsTree->SetBranchAddress("tEventId",&tEventId);
	upsTree->SetBranchAddress("tEventRunId",&tEventRunId);
	upsTree->SetBranchAddress("tEventTpcVx",&tEventTpcVx);
	upsTree->SetBranchAddress("tEventTpcVy",&tEventTpcVy);
	upsTree->SetBranchAddress("tEventTpcVz",&tEventTpcVz);
	upsTree->SetBranchAddress("tEventVpdVz",&tEventVpdVz);
	upsTree->SetBranchAddress("tEventCent9",&tEventCent9);
	upsTree->SetBranchAddress("tEventCent16",&tEventCent16);
	upsTree->SetBranchAddress("tEventRefMult",&tEventRefMult);
	upsTree->SetBranchAddress("tEventRefMultC",&tEventRefMultC);
	upsTree->SetBranchAddress("tNElectrons",&tNElectrons);
	upsTree->SetBranchAddress("tElePt",tElePt);
	upsTree->SetBranchAddress("tEleP",tEleP);
	upsTree->SetBranchAddress("tElePx",tElePx);
	upsTree->SetBranchAddress("tElePy",tElePy);
	upsTree->SetBranchAddress("tElePz",tElePz);
	upsTree->SetBranchAddress("tEleEta",tEleEta);
	upsTree->SetBranchAddress("tEleCharge",tEleCharge);
	upsTree->SetBranchAddress("tEleNHits",tEleNHits);
	upsTree->SetBranchAddress("tEleNRat",tEleNRat);
	upsTree->SetBranchAddress("tEleNDedx",tEleNDedx);
	upsTree->SetBranchAddress("tEleDca",tEleDca);
	upsTree->SetBranchAddress("tEleDedx",tEleDedx);
	upsTree->SetBranchAddress("tEleNSigE",tEleNSigE);
	upsTree->SetBranchAddress("tEleE",tEleE);
	upsTree->SetBranchAddress("tEleEcl",tEleEcl);
	upsTree->SetBranchAddress("tEleZDist",tEleZDist);
	upsTree->SetBranchAddress("tElePhiDist",tElePhiDist);
	upsTree->SetBranchAddress("tEleR",tEleR);
	upsTree->SetBranchAddress("tEleTrig",tEleTrig);
	//Float_t px[kMaxTrack];

	TH1F* hEventVz               	= new TH1F("hEventVz","hEventVz",100,-30,30);
    TH1F* hEleP              		= new TH1F("hEleP","",200,0,15);
    TH1F* hIMun   		        	= new TH1F("hIMun","",200,0,20);
    hIMun->Sumw2();
    TH1F* hIMli    		        	= new TH1F("hIMli","",200,0,20);
    hIMli->Sumw2();

    TH1F* hRbg   		        	= new TH1F("hRbg","",400,0,20);
    hRbg->Sumw2();

    hIMmix                  		= new TH1F("hIMmix","",200,0,20);
    hIMmix->Sumw2();
    hIMmixC                 		= new TH2F("hIMmixC","",200,0,20,10,0,10);
    hIMmixC->Sumw2();

    TH1F* hUpsPt 	             	= new TH1F("hUpsPt","",200,0,20);
    hUpsPt->Sumw2();

    int ncuts = 8;
    if (cutSet+1>ncuts)	{
    	cout << "Invalid set of cuts was chosen! Setting to default." << endl;
    	cutSet = 0;	}

    	//6 - optimised for p17lm y<1 (sig)
    	//7 - optimised for p17lm y<.5 (sig)

    const float cutVzTpc[] 		= {30., 	30.,	30.,	30.,	30.,	30.,	30.,	30,};	//vertex
    const float cutVzDif[] 		= {4., 		4., 	4.,		4.,		4.,		4.,		4.,		4.};

    const int cutNHits[] 		= {20, 		30, 	25,		20,		20,		20,		20,		20};	//nhits
    const float cutNRat[] 		= {0.52, 	0.52, 	0.52,	0.52,	0.52,	0.52,	0.52,	0.52};
    const int cutNDedx[] 		= {10, 		10, 	10,		10,		10,		10,		10,		10};
    const float cutNsigL[] 		= {-1.5,	-1.5,	-1.5,	-1.5,	-1.5,	-1.5,	-1.4,	-1.3};
    const float cutNsigT[] 		= {3.0, 	3.0, 	3.0,	3.0,	3.0,	3.0,	3.0,	3.0};
    const float cutEta[] 		= {1.0, 	1.0, 	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};
    const float cutEpL[] 		= {0.75, 	0.75, 	0.3,	0.3,	0.75,	0.75,	0.65,	0.70};
    const float cutEpT[] 		= {1.5, 	1.5, 	1.8,	1.8,	1.5,	1.5,	1.6,	1.4};
    const float cutR[] 			= {0.025, 	0.025, 	0.025,	0.03,	0.025,	0.025,	0.026,	0.026};
    const float cutPlow[] 		= {3.25, 	3.5, 	3.5,	3.0,	3.5,	3.5,	3.5,	3.5};
    const float cutDca[] 		= {0.75, 	0.75, 	1.5,	1.5,	0.75,	0.75,	3.0,	3.0};
    const float cutPlead[] 		= {4.5, 	4.5, 	4.5,	4.5,	4.5,	4.5,	4.5,	4.5};
    const float cutPairY[] 		= {1.0, 	1.0, 	1.0,	1.0,	1.0,	0.5,	1.0,	0.5};
    const float cutPairPt[] 	= {10.0, 	10.0, 	10.0,	10.0,	10.0,	10.0,	10.,	10.};

    // make rbg histo
    TString bbcuts = Form("eleGlDca<%f&&posGlDca<%f&&",cutDca[cutSet],cutDca[cutSet]);
    bbcuts += Form("elePrP>%f&&posPrP>%f&&",cutPlow[cutSet],cutPlow[cutSet]);
    bbcuts += Form("elePrEta>-%f&&posPrEta>-%f&&",cutEta[cutSet],cutEta[cutSet]);
    bbcuts += Form("elePrEta<%f&&posPrEta<%f&&",cutEta[cutSet],cutEta[cutSet]);
    bbcuts += Form("(elePrP>%f||posPrP>%f)&&",cutPlead[cutSet],cutPlead[cutSet]);
    bbcuts += Form("motherPt<%f&&",cutPairPt[cutSet]);
    bbcuts += Form("motherY>-%f&&",cutPairY[cutSet]);
    bbcuts += Form("motherY<%f",cutPairY[cutSet]);
    //bbcuts += Form("&&motherY<%f",-10.);
    bbbar->Draw("motherM>>hRbg",bbcuts.Data(),"");

    	//"eleGlDca<0.75&&posGlDca<0.75&&elePrP>3.25&&posPrP>3.25&&elePrEta<1.0&&elePrEta>-1.0&&posPrEta<1.0&&posPrEta>-1.0&&(elePrP>4.5||posPrP>4.5)&&motherPt<10&&motherY<1.0&&motherY>-1.0")")

    cout << "done with rbg" << endl;

    // mixing events
    int mixCounter = 0;
    int mixSaved = 5000;
    vector<TLorentzVector> trigElectrons;
    vector<Float_t> trigElectronsVz;
    vector<Int_t> trigElectronsC;
    vector<TLorentzVector> allElectrons;
    vector<Float_t> allElectronsVz;
    vector<Int_t> allElectronsC;

    //nEntries = 50000;
	for (int iEv = 0; iEv < nEntries; ++iEv)
	{
		upsTree->GetEntry(iEv);
		hEventVz->Fill(tEventTpcVz);

		//event cuts
		if (fabs(tEventTpcVz) > cutVzTpc[cutSet]) continue;
		if (fabs(tEventTpcVz-tEventVpdVz) > cutVzDif[cutSet]) continue; 
		//if (tEventCent9 > 4) continue;

		for (int iTr1 = 0; iTr1 < tNElectrons; ++iTr1)
		{
			//track cuts
			if (tEleNHits[iTr1] < cutNHits[cutSet]) continue;
			if (tEleNRat[iTr1] < cutNRat[cutSet]) continue;
			//if (tEleNDedx[iTr1] < cutNdDedx[cutSet]) continue;
			if (tEleNSigE[iTr1] < cutNsigL[cutSet]) continue;
			if (tEleNSigE[iTr1] > cutNsigT[cutSet]) continue;
			if (fabs(tEleEta[iTr1]) > cutEta[cutSet]) continue;
			if (tEleEcl[iTr1]/tEleP[iTr1] < cutEpL[cutSet]) continue;
			if (tEleEcl[iTr1]/tEleP[iTr1] > cutEpT[cutSet]) continue;
			if (tEleR[iTr1] > cutR[cutSet]) continue;
			if (tEleP[iTr1] < cutPlow[cutSet]) continue;
			if (tEleDca[iTr1] > cutDca[cutSet]) continue;
			hEleP->Fill(tEleP[iTr1]);

			for (int iTr2 = iTr1+1; iTr2 < tNElectrons; ++iTr2)
			{
				//partner cuts
				if (tEleNHits[iTr2] < cutNHits[cutSet]) continue;
				if (tEleNRat[iTr2] < cutNRat[cutSet]) continue;
				//if (tEleNDedx[iTr2] < cutNDedx[cutSet]) continue;
				if (tEleNSigE[iTr2] < cutNsigL[cutSet]) continue;
				if (tEleNSigE[iTr2] > cutNsigT[cutSet]) continue;
				if (fabs(tEleEta[iTr2]) > cutEta[cutSet]) continue;
				if (tEleEcl[iTr2]/tEleP[iTr2] < cutEpL[cutSet]) continue;
				if (tEleEcl[iTr2]/tEleP[iTr2] > cutEpT[cutSet]) continue;
				if (tEleR[iTr2] > cutR[cutSet]) continue;
				if (tEleP[iTr2] < cutPlow[cutSet]) continue;
				if (tEleDca[iTr2] > cutDca[cutSet]) continue;

				//pair cuts
				if (!(tEleTrig[iTr1]%2)&&!(tEleTrig[iTr2]%2)) continue;

				TLorentzVector* tl1 = new TLorentzVector(tElePx[iTr1],tElePy[iTr1],tElePz[iTr1],sqrt(tEleP[iTr1]*tEleP[iTr1]+elmass*elmass));
				TLorentzVector* tl2 = new TLorentzVector(tElePx[iTr2],tElePy[iTr2],tElePz[iTr2],sqrt(tEleP[iTr2]*tEleP[iTr2]+elmass*elmass));
				TLorentzVector* q = new TLorentzVector(*tl1+*tl2);

				if (mixCounter < 2*mixSaved) {
					if ((mixCounter+1)%2) {
						if (tEleTrig[iTr1]%2) trigElectrons.push_back(*tl1);
						else if (tEleTrig[iTr2]%2) trigElectrons.push_back(*tl2);
						trigElectronsVz.push_back(tEventTpcVz);
						trigElectronsC.push_back(tEventCent9);	}
					if (mixCounter%2) {
						if (!tEleTrig[iTr1]%2) allElectrons.push_back(*tl1);
						else if (!tEleTrig[iTr2]%2) allElectrons.push_back(*tl2);
						allElectronsVz.push_back(tEventTpcVz);
						allElectronsC.push_back(tEventCent9);	}
					mixCounter++;
				}

				if (tEleP[iTr1]<cutPlead[cutSet] &&tEleP[iTr2]<cutPlead[cutSet]) continue;

				if (q->M() < 0.2) continue;
				if (q->Perp() > cutPairPt[cutSet]) continue;
				if (fabs(q->Rapidity()) > cutPairY[cutSet]) continue;

				fMass = q->M();
				fCentrality = tEventCent9;
				fSign = tEleCharge[iTr1]*tEleCharge[iTr2];
				massTree->Fill();

				if (tEleCharge[iTr1]*tEleCharge[iTr2] < 0) {
					hIMun->Fill(q->M()); 
					if (q->M() < 9.9 && q->M() > 8.9) hUpsPt->Fill(q->Perp());}
				else {
					hIMli->Fill(q->M()); }
			}
		}
	}

	cout << "done reading tree" << endl;

	// do mixed events mass here
	for (int iET = 0; iET < trigElectrons.size(); ++iET)
	{
		TLorentzVector* p1 = new TLorentzVector(trigElectrons[iET]);
		for (int iEA = 0; iEA < allElectrons.size(); ++iEA)
		{
			TLorentzVector* p2 = new TLorentzVector(allElectrons[iEA]);
			if (fabs(trigElectronsVz[iET]-allElectronsVz[iEA])>7) continue;
			if (trigElectronsC[iET] != allElectronsC[iEA]) continue;

			TLorentzVector* q = new TLorentzVector(*p1 + *p2);
			if (p1->P() < cutPlead[cutSet] && p2->P() < cutPlead[cutSet]) continue;
			if (q->M() < 0.2) continue;
			if (q->Perp() > cutPairPt[cutSet]) continue;
			if (fabs(q->Rapidity()) > cutPairY[cutSet]) continue;

			hIMmix->Fill(q->M());
			hIMmixC->Fill(q->M(),trigElectronsC[iET]);
		}
	}

	cout << "done mixing" << endl;

	int abin = hIMli->FindBin(5.0);
	int bbin = hIMli->FindBin(15.0);
	hIMmix->Scale((float)hIMli->Integral(abin,bbin)/hIMmix->Integral(abin,bbin));
	
	TFile *frbg = new TFile("hRbg.root","RECREATE");
	frbg->cd();
	hRbg->Write();
	frbg->Close();


	TFile* fmass = new TFile("masstree.root","RECREATE");
	fmass->cd();
	hIMmixC->Write();
	massTree->Write();
	fmass->Close();

	hIMun->Draw();
	hIMli->SetLineColor(kRed);
	hIMli->Draw("same");

	hIMmix->SetLineColor(kGreen+2);
	hIMmix->Draw("same");

	new TCanvas;
	hUpsPt->Draw();

	int a = hIMun->FindBin(8.9);
	int b = hIMun->FindBin(9.9);
	float un = hIMun->Integral(a,b);
	float li = hIMli->Integral(a,b);
	cout << "un " << un << " li " << li << " ups " << un-li << " s/b " << (un-li)/li << " sig " << (un-li)/sqrt(un) << endl;

	new TCanvas;
	hIMun->Add(hIMli,-1);
	//hIMun->Add(hIMmix,-1);
	hIMun->Rebin();
	hIMun->GetXaxis()->SetRangeUser(5.,15.);
	hIMun->Draw();

	

}
