// apologies to anyone who ever has to go through this

#include <iostream>
#include <fstream>
#include <cmath>

// defining smearing functions globally
TF1* fMuPt;
TF1* fSigPt;
TF1* fMuPhi;
TF1* fSigPhi;
TF1* fMuEta;
TF1* fSigEta;

float Eerr(float E) {
	float Eres = sqrt(1.5*1.5+14.0*14.0/E)*0.01*E;
	return Eres;
}

TLorentzVector SmearTrack(TLorentzVector* t){
    

    double ptT = t->Pt();
    double phiT = t->Phi();
    double etaT = t->Eta();

    double muPt = fMuPt->Eval(ptT);
    double sigPt = fSigPt->Eval(ptT);
    double muPhi = fMuPhi->Eval(ptT);
    double sigPhi = fSigPhi->Eval(ptT);
    double muEta = fMuEta->Eval(ptT);
    double sigEta = fSigEta->Eval(ptT);

    double ptM = ptT*(1+gRandom->Gaus(muPt,sigPt));
    double phiM = phiT + gRandom->Gaus(muPhi,sigPhi);
    double etaM = etaT + gRandom->Gaus(muEta,sigEta);

    TLorentzVector ts;// = new TLorentzVector;
    ts.SetPtEtaPhiM(ptM,etaM,phiM,0.00051);

    return ts;
}

void nTupEff_sym(const Char_t *inputFile="upsE_0209.root", Int_t cutSet = 0, Int_t dataSet = 0, Int_t cMin = 0, Int_t cMax = 9, const Char_t *outFile="effiLM1S.txt") {

	// LOAD FILES
	TFile* fin = new TFile(inputFile,"READ");
	//TFile* fin = new TFile("upsE2S_0322.root","READ");
	//TFile* fin = new TFile("upsE3S_0322.root","READ");
	TFile* fnsig = new TFile("nSigEff.root","READ");
	TFile* fdca = new TFile("hDca.root","READ");
	//TFile* fdca = new TFile("dcaEffNPE.root","READ");
	TFile* fSmear = new TFile("hSmearFunc.root","READ");
  
	// OUTPUT FILE
	TFile* fout = new TFile("Y1S_eff.root","RECREATE");
	TString dir = "ntupeff/";

	// SELECT BINNING
	int nbinspt = 21;
	float hleft = 0.0;
	float hright = 10.5;

	int nbinspt1E = 40;
	float hleft1E = 1.0;
	float hright1E = 13.0;

	// CUT VALUES

	//const float centrality = 400;

	// event
	int ncuts = 6;
	const float cutvzl[] 			= {30.0, 	30., 	30.,	30.,	30.,	30.};
	const float cutvdiff[] 			= {4.0, 	4., 	4.,		4.,		4.,		4. };			// vpd not in ntuple :/

	// tracking
	const int cutnhitsf[] 			= {20, 		30, 	25,		20,		30,		20};
	const float cutnratio[] 		= {0.52, 	0.52, 	0.52,	0.52,	0.52,	0.52};
	const int cutndedx[] 			= {10,	 	10, 	10,		10,		10,		10};		// nHitsdEdx not in ntuple :/
	const float cutdca[] 			= {0.75, 	0.75, 	1.5,	1.5,	0.75,	0.75};			// dca not in ntuple :/
	/*const float cutdcaprob = 0.9689;	// !!!!!!
	const float cutdcaprobL = 0.9608;
	const float cutdcaprobT = 0.9725;*/

	// kinematics
	const float cutetal[] 			= {1.0, 	1.0, 	1.0,	1.0,	1.0,	1.0};
	const float cutplow[] 			= {3.25,	3.5,	3.5,	3.0,	3.25,	3.5};		//sigma for 3.5 = 0.15		// for 3 it's 0.13
	const float cutplowL[] 			= {3.11,	3.35,	3.35,	2.87,	3.11,	3.35};		//2.87;
	const float cutplowT[] 			= {3.39,	3.65,	3.65,	3.13,	3.39,	3.65};
	const float cutplead[] 			= {4.5,		4.5,	4.5,	4.5,	4.5,	4.5};		//sigma = 0.25
	const float cutpleadL[] 		= {4.25,	4.25,	4.25,	4.25,	4.25,	4.25};
	const float cutpleadT[] 		= {4.75,	4.75,	4.75,	4.75,	4.75,	4.75};

	// emc
	const float cutadc = 297;//19<<4;
		const float cutadcL = cutadc-9;
		const float cutadcT = cutadc+9;
	
	const float cutepl[] 			= {0.75,	0.75,	0.3,	0.3,	0.75,	0.75};		//sigma for 0.3 = 0.04     // divide by sqrt(12) ?
		const float cuteplL[] 		= {0.69,	0.69,	0.26,	0.26,	0.69,	0.69};		//for 0.75 = 0.06, for 1.5 = 0.1
		const float cuteplT[] 		= {0.81,	0.81,	0.34,	0.34,	0.81,	0.81};
	const float cutepr[] 			= {1.5,		1.5,	1.8,	1.8,	1.5,	1.5};		//sigma for 1.8 = 0.08
		const float cuteprL[] 		= {1.6,		1.6,	1.88,	1.88,	1.6,	1.6};
		const float cuteprT[] 		= {1.4,		1.4,	1.72,	1.72,	1.4,	1.4};
	const float cutzdl = 5.0;			// zDist not in ntuple :/
	const float cutpdl = 0.05;			// phiDist not in ntuple :/
	const float cutrd[]  			= {0.025,	0.025,	0.03,	0.03,	0.025,	0.025};		// sigma = 0.0005, anthony has 0.008, but rc vs mc lasttpc point sigma= max 0.8cm->0.002
		const float cutrdL[]  		= {0.027,	0.027,	0.032,	0.032,	0.027,	0.027}; // divide by sqrt(12) ?
		const float cutrdT[]  		= {0.023,	0.023,	0.028,	0.028,	0.023,	0.023};

	// pair
	const float cuturap[]			= {1.0,		1.0,	1.0,	1.0,	1.0,	0.5};
	const float cutupt[] 			= {10.0,	10.0,	10.0,	10.0,	10.0,	10.0};
	const float elmass = 0.000510998910;

	// GREFMULT BINS
	int grefmultLM[] 	= {10, 22, 42, 74, 120, 183, 268, 379, 447, 9999};
	int grefmultH[] 	= {10, 22, 43, 75, 123, 188, 275, 393, 467, 9999};

	// LOAD SHIT IN
	TNtuple* ups = (TNtuple*)fin->Get("ups");
	TF1* efit = (TF1*)fnsig->Get("efit");
	TF1* efitT = (TF1*)fnsig->Get("efitT");
	TF1* efitL = (TF1*)fnsig->Get("efitL");

    fMuPt = (TF1*)fSmear->Get("fMuPt");
    fSigPt = (TF1*)fSmear->Get("fSigPt");
    fMuPhi = (TF1*)fSmear->Get("fMuPhi");
    fSigPhi = (TF1*)fSmear->Get("fSigPhi");
    fMuEta = (TF1*)fSmear->Get("fMuEta");
    fSigEta = (TF1*)fSmear->Get("fSigEta");

    TH1F* hDcaGl = (TH1F*)fdca->Get("hDcaGl");
    int dcabin = hDcaGl->FindBin(cutdca[cutSet]);
    float cutdcaprob = (float)hDcaGl->Integral(1,dcabin)/hDcaGl->Integral(1,999);
    float cutdcaprobL = (float)hDcaGl->Integral(1,dcabin+1)/hDcaGl->Integral(1,999);
    float cutdcaprobT = (float)hDcaGl->Integral(1,dcabin-1)/hDcaGl->Integral(1,999);

    // GET RUN NUMBERS OF LOWMID/HIGH EVENTS
    fstream fileLowmid("RunListLM.txt");
    fstream fileHigh("RunListH.txt");
    vector<int> runIds[2];
	string t;
	while ( !fileLowmid.eof() )	{
   		getline(fileLowmid,t);
   		runIds[0].push_back(atoi(t.c_str()));	}
	while ( !fileHigh.eof() )	{
   		getline(fileHigh,t);
   		runIds[1].push_back(atoi(t.c_str()));	}

	// SET UP TREE
	int nentries = ups->GetEntries();
	cout << "nTuple has " << nentries << " entries." << endl;
	float vx, vy, vz, vzRc;
	float runId, ntracksRec, grefmult;
	float eleReco, eleTpcHits, eleTpcPos;
	float eleRcFirstTpcX, eleRcFirstTpcY, eleRcFirstTpcZ;
	float posRcFirstTpcX, posRcFirstTpcY, posRcFirstTpcZ; 
	float posReco, posTpcHits, posTpcPos;
	float eleP, elePt, elePx, elePy, elePz, elePhi, eleEta, elePRc, elePtRc, elePxRc, elePyRc, elePzRc, elePhiRc, eleEtaRc;
	float posP, posPt, posPx, posPy, posPz, posPhi, posEta, posPRc, posPtRc, posPxRc, posPyRc, posPzRc, posPhiRc, posEtaRc;
	float eleAdc, eleEcluster, eleEtower, eleRcAdc, eleRcEnergy, eleR;
	float posAdc, posEcluster, posEtower, posRcAdc, posRcEnergy, posR;
	float isL0;
	float upsPt, upsRap, upsPtRc, upsMRc;

	ups->SetBranchAddress("vz",&vz);
	ups->SetBranchAddress("vx",&vz);
	ups->SetBranchAddress("vy",&vz);
	ups->SetBranchAddress("vzRc",&vzRc);
	ups->SetBranchAddress("ntracksRec",&ntracksRec);
	ups->SetBranchAddress("runId",&runId);
	ups->SetBranchAddress("grefmult",&grefmult);
	ups->SetBranchAddress("isL0",&isL0);

	ups->SetBranchAddress("eleReco",&eleReco);
	ups->SetBranchAddress("eleTpcHits",&eleTpcHits);
	ups->SetBranchAddress("eleTpcPos",&eleTpcPos);
	ups->SetBranchAddress("eleP",&eleP);
	ups->SetBranchAddress("elePt",&elePt);
	ups->SetBranchAddress("elePz",&elePz);
	ups->SetBranchAddress("elePhi",&elePhi);
	ups->SetBranchAddress("eleEta",&eleEta);
	ups->SetBranchAddress("elePRc",&elePRc);
	ups->SetBranchAddress("elePtRc",&elePtRc);
	ups->SetBranchAddress("elePxRc",&elePxRc);
	ups->SetBranchAddress("elePyRc",&elePyRc);
	ups->SetBranchAddress("elePzRc",&elePzRc);
	ups->SetBranchAddress("elePhiRc",&elePhiRc);
	ups->SetBranchAddress("eleEtaRc",&eleEtaRc);
	ups->SetBranchAddress("eleAdc",&eleAdc);
	ups->SetBranchAddress("eleEcluster",&eleEcluster);
	ups->SetBranchAddress("eleEtower",&eleEtower);
	ups->SetBranchAddress("eleRcAdc",&eleRcAdc);
	ups->SetBranchAddress("eleRcEnergy",&eleRcEnergy);
	ups->SetBranchAddress("eleR",&eleR);

	ups->SetBranchAddress("posReco",&posReco);
	ups->SetBranchAddress("posTpcHits",&posTpcHits);
	ups->SetBranchAddress("posTpcPos",&posTpcPos);
	ups->SetBranchAddress("posP",&posP);
	ups->SetBranchAddress("posPt",&posPt);
	ups->SetBranchAddress("posPz",&posPz);
	ups->SetBranchAddress("posPhi",&posPhi);
	ups->SetBranchAddress("posEta",&posEta);
	ups->SetBranchAddress("posPRc",&posPRc);
	ups->SetBranchAddress("posPtRc",&posPtRc);
	ups->SetBranchAddress("posPxRc",&posPxRc);
	ups->SetBranchAddress("posPyRc",&posPyRc);
	ups->SetBranchAddress("posPzRc",&posPzRc);
	ups->SetBranchAddress("posPhiRc",&posPhiRc);
	ups->SetBranchAddress("posEtaRc",&posEtaRc);
	ups->SetBranchAddress("posAdc",&posAdc);
	ups->SetBranchAddress("posEcluster",&posEcluster);
	ups->SetBranchAddress("posEtower",&posEtower);
	ups->SetBranchAddress("posRcAdc",&posRcAdc);
	ups->SetBranchAddress("posRcEnergy",&posRcEnergy);
	ups->SetBranchAddress("posR",&posR);

	ups->SetBranchAddress("upsPt",&upsPt);
	ups->SetBranchAddress("upsRap",&upsRap);
	ups->SetBranchAddress("upsPtRc",&upsPtRc);
	ups->SetBranchAddress("upsMRc",&upsMRc);

	// do histograms here
	TH1F::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetHistLineColor(kAzure);
	gStyle->SetHistFillColor(kBlue-10);
	gStyle->SetTitleOffset(1.4,"Y");

	TH1F* hNhits 	= new TH1F("hNhits",";nHitsFit;counts",60,0,60);
	TH1F* hNrat 	= new TH1F("hNrat",";nHitsRatio;counts",50,0,1);
	TH1F* hEta 		= new TH1F("hEta",";eta;counts",100,-2,2);
	TH1F* hEp 		= new TH1F("hEp",";E/p;counts",300,0,3);
	TH1F* hR 		= new TH1F("hR",";R;counts",100,0,0.2);

	TH1F* hUpt		= new TH1F("hUpt",";upsilon p_{T};counts",50,0,12);
	TH1F* hUrap		= new TH1F("hUrap",";upsilon rapidity;counts",60,-1.5,1.5);
	TH1F* hUphi		= new TH1F("hUphi",";upsilon phi;counts",60,-3.2,3.2);



	//gStyle->SetHistLineColor(kRed);
	gStyle->SetMarkerSize(0.8);
  	gStyle->SetMarkerStyle(20);
  	gStyle->SetMarkerColor(2);
  	gStyle->SetHistLineColor(2);
  	gStyle->SetHistLineWidth(2);

  	int nhistos = 14;
  	int nsyst = 18;

	TH1F* hNom[nhistos][nsyst];
	TH1F* hDen[nsyst];
	TH1F* hTotCnt[nsyst];
 	TH1F* hAccE[nsyst];
 	TH1F* hTrigE[nsyst];
	TH1F* hRecoE[nsyst];
	TH1F* hNhitsE[nsyst];
	TH1F* hNratE[nsyst];
	TH1F* hEtaE[nsyst];
	TH1F* hPlowE[nsyst];
	TH1F* hPleadE[nsyst];
	TH1F* hEpE[nsyst];
	TH1F* hRE[nsyst];
	TH1F* hDcaE[nsyst];
	TH1F* hUrapE[nsyst];
	TH1F* hNsigE[nsyst];
	TH1F* hTotE[nsyst];
	TH1F* hDca[nsyst];

	TH1F* h1ENom[nhistos][nsyst];

	float err[nbinspt][nsyst];

	for (int j = 0; j < nsyst; ++j)
	{
		for (int i = 0; i < nhistos; ++i)
		{
			hNom[i][j] 	= new TH1F(Form("hNom%i_%i",i,j),"",nbinspt,hleft,hright);
			h1ENom[i][j]		= new TH1F(Form("h1ENom%i_%i",i,j),"",nbinspt1E,hleft1E,hright1E);	
		}

		hDen[j]	 		= new TH1F(Form("hDen_%i",j),"",nbinspt,hleft,hright);
		hTotCnt[j]	 	= new TH1F(Form("hTotCnt_%i",j),"",10,0,10);
 		hAccE[j]	 	= new TH1F(Form("hAccE_%i",j),"",nbinspt,hleft,hright);
 		hTrigE[j]	 	= new TH1F(Form("hTrigE_%i",j),"",nbinspt,hleft,hright);
		hRecoE[j]	 	= new TH1F(Form("hRecoE_%i",j),"",nbinspt,hleft,hright);


		hNhitsE[j]	 	= new TH1F(Form("hNhitsE_%i",j),"",nbinspt,hleft,hright);
		hNratE[j]	 	= new TH1F(Form("hNratE_%i",j),"",nbinspt,hleft,hright);

		hEtaE[j]	 	= new TH1F(Form("hEtaE_%i",j),"",nbinspt,hleft,hright);
		hPlowE[j]	 	= new TH1F(Form("hPlowE_%i",j),"",nbinspt,hleft,hright);
		hPleadE[j]	 	= new TH1F(Form("hPleadE_%i",j),"",nbinspt,hleft,hright);

		hEpE[j]	 		= new TH1F(Form("hEpE_%i",j),"",nbinspt,hleft,hright);
		hRE[j]	 		= new TH1F(Form("hRE_%i",j),"",nbinspt,hleft,hright);

		hDcaE[j]	    = new TH1F(Form("hDcaE_%i",j),"",nbinspt,hleft,hright);

		hUrapE[j]		= new TH1F(Form("hUrapE_%i",j),"",nbinspt,hleft,hright);

		hNsigE[j]		= new TH1F(Form("hNsigE_%i",j),"",nbinspt,hleft,hright);

		hTotE[j]		= new TH1F(Form("hTotE_%i",j),"",nbinspt,hleft,hright);
		hDca[j]	 		= new TH1F(Form("hDca_%i",j),"",60,0,15);
	
	}


	TH1F*	hTotEL		= new TH1F(Form("hTotEL"),"",nbinspt,hleft,hright);
	TH1F*	hTotET		= new TH1F(Form("hTotET"),"",nbinspt,hleft,hright);

	TH1F* h1ETrigPhi1	= new TH1F("h1ETrigPhi1","",200,-3.4,3.4);
	TH1F* h1ETrigPhi2	= new TH1F("h1ETrigPhi2","",200,-3.4,3.4);
	TH1F* h1ERecoPhi1	= new TH1F("h1ERecoPhi1","",50,-3.4,3.4);
	TH1F* h1ERecoPhi2	= new TH1F("h1ERecoPhi2","",50,-3.4,3.4);
	TH1F* h1ETrigPhi10	= new TH1F("h1ETrigPhi10","",200,-3.4,3.4);
	TH1F* h1ETrigPhi20	= new TH1F("h1ETrigPhi20","",200,-3.4,3.4);
	TH1F* h1ERecoPhi10	= new TH1F("h1ERecoPhi10","",50,-3.4,3.4);
	TH1F* h1ERecoPhi20	= new TH1F("h1ERecoPhi20","",50,-3.4,3.4);

	TH2F* hEpErr 		  
	= new TH2F("hEpErr","",100,0,3,100,0,3);
  	TH2F* hRcminMcPt      = new TH2F("hRcminMcPt",";p_{T}^{MC};(1/p_{T}^{MC})(p_{T}^{RC}-p_{T}^{MC})",200,0,20,200,-1.0,1.0);
  	TH2F* hRcminMcPhi     = new TH2F("hRcminMcPhi",";p_{T}^{MC};(#phi_{T}^{RC}-#phi_{T}^{MC}",200,0,20,200,-1.0,1.0);
  	TH2F* hRcminMcEta     = new TH2F("hRcminMcEta",";p_{T}^{MC};(#eta_{T}^{RC}-#eta_{T}^{MC}",200,0,20,200,-1.0,1.0);


	TH1F* h1EDen;
	TH1F* h1ETrig[nsyst];
	TH1F* h1EEp[nsyst];

	h1EDen		= new TH1F(Form("h1EDen"),"",nbinspt1E,hleft1E,hright1E);
	for (int j = 0; j < nsyst; ++j)
	{
		h1ETrig[j]		= new TH1F(Form("h1ETrig_%i",j),"",nbinspt1E,hleft1E,hright1E);
		h1EEp[j]		= new TH1F(Form("h1EEp_%i",j),"",nbinspt1E,hleft1E,hright1E);
	}

	TH3F* hUpsPepP 	= new TH3F("hUpsPepP",";Upsilon p_{T};electron p [GeV/c];positron p",nbinspt,hleft,hright,150,0,15,150,0,15);	
    hEMCEFFepUS     = new TH3F("hEMCEFFepUS","",150,0,15,200,0,4,10,0,10);

    // MC random uniform
    TRandom3* r3 = new TRandom3();

	//nentries = 50000;
	for (int i = 0; i < nentries; ++i)
	{
		ups->GetEntry(i);
		if (i%50000==0) cout << "Working on entry: " << i << endl;

		if (dataSet==0||dataSet==1) {
			bool isInDataset = false;
			for (int iRun = 0; iRun < runIds[dataSet].size(); ++iRun)		
			{
				if (runId==runIds[dataSet][iRun]) isInDataset = true;		
			}
			if (!isInDataset) continue;	
		}

		if (dataSet!=1) {
			if (grefmult <= grefmultLM[cMin] || grefmult > grefmultLM[cMax]) continue;		}
		else {
			if (grefmult <= grefmultH[cMin] || grefmult > grefmultH[cMax]) continue;		}

		//if (grefmult < 379) continue;

		if (fabs(upsRap)>cuturap[cutSet]) continue;

		//if (ntracksRec<1000) continue;

		float dcar1 = r3->Uniform(0,1);		// use MC for external efficiencies
		float dcar2 = r3->Uniform(0,1);
		float nsigr1 = r3->Uniform(0,1);
		float nsigr2 = r3->Uniform(0,1);
		for (int j = 0; j < nsyst; ++j)
    	{

			hDen[j]->Fill(upsPt);
			hNom[0][j]->Fill(upsPt);
			hTotCnt[j]->Fill(1);	
			if (!j) h1EDen->Fill(eleP);
			if (!j) h1EDen->Fill(posP);
			h1ENom[0][j]->Fill(eleP);
			h1ENom[0][j]->Fill(posP);		

			// event cuts
			TVector3 primV(vx,vy,vz);	

			//tpc acceptance?
			if (fabs(eleEta)<cutetal[cutSet]) h1ENom[1][j]->Fill(eleP);
			if (fabs(posEta)<cutetal[cutSet]) h1ENom[1][j]->Fill(posP);
			if (fabs(eleEta)>cutetal[cutSet]||fabs(posEta)>cutetal[cutSet]) continue;
			h1ENom[1][j]->SetTitle("single e TPC acceptance;electron p [GeV/c];efficiency");
			hAccE[j]->Fill(upsPt);
			hNom[1][j]->Fill(upsPt);
			hNom[1][j]->SetTitle("TPC acceptance;upsilon p_{T};efficiency");

			// could the trigger be fired?

			h1ENom[2][j]->Fill(eleP);
			h1ENom[2][j]->Fill(posP);	

			//if (!j && elePt > 10) cout << "pt " << elePt << " adc " << eleAdc << endl;

			if (!j && eleAdc >= (int)cutadcL) h1ENom[3][1]->Fill(eleP);
			if (!j && eleAdc >= (int)cutadc) h1ENom[3][0]->Fill(eleP);

			if (!j && eleAdc >= (int)cutadc && eleP < 9 && eleP > 6) h1ETrigPhi1->Fill(eleEta);
			if (!j && eleP < 9 && eleP > 6) h1ETrigPhi10->Fill(eleEta);
			if (!j && eleAdc >= (int)cutadc && eleP > 11) h1ETrigPhi2->Fill(eleEta);
			if (!j && eleP > 11) h1ETrigPhi20->Fill(eleEta);
			if (!j && posAdc >= (int)cutadc && posP < 9 && posP > 6) h1ETrigPhi1->Fill(posEta);
			if (!j && posP < 9 && posP > 6) h1ETrigPhi10->Fill(posEta);
			if (!j && posAdc >= (int)cutadc && posP > 11) h1ETrigPhi2->Fill(posEta);
			if (!j && posP > 11) h1ETrigPhi20->Fill(posEta);

			if (!j && eleAdc >= (int)cutadcT) h1ENom[3][2]->Fill(eleP);
			if (!j && posAdc >= (int)cutadcL) h1ENom[3][1]->Fill(posP);
			if (!j && posAdc >= (int)cutadc) h1ENom[3][0]->Fill(posP);
			if (!j && posAdc >= (int)cutadcT) h1ENom[3][2]->Fill(posP);
			if (j==1) if (eleAdc < (int)cutadcL && posAdc < (int)cutadcL) continue;
			if (j==2) if (eleAdc < (int)cutadcT && posAdc < (int)cutadcT) continue;
			if (j!=1&&j!=2) if (eleAdc < (int)cutadc && posAdc < (int)cutadc) continue;

			hTrigE[j]->Fill(upsPt);
			hNom[2][j]->Fill(upsPt);
			hNom[2][j]->SetTitle("trigger efficiency;upsilon p_{T}; efficiency");
			h1ENom[4][j]->Fill(eleP);
			h1ENom[4][j]->Fill(posP);	

			// were the daughters reconstructed?
			if (eleReco==1) h1ENom[5][j]->Fill(eleP);

			if (!j && 1 && eleReco==1 && eleP < 6) h1ERecoPhi1->Fill(elePhi);
			if (!j && 1 && eleP < 6) h1ERecoPhi10->Fill(elePhi);
			if (!j && 1 && eleReco==1 && eleP > 11) h1ERecoPhi2->Fill(elePhi);
			if (!j && 1 && eleP > 11) h1ERecoPhi20->Fill(elePhi);

			if (posReco==1) h1ENom[5][j]->Fill(posP);
			if (eleReco!=1||posReco!=1) continue;
			h1ENom[5][j]->SetTitle("single e reconstruction;electron p [GeV/c];efficiency");
			hRecoE[j]->Fill(upsPt);
			hNom[3][j]->Fill(upsPt);
			hNom[3][j]->SetTitle(" daughter reconstruction;upsilon p_{T};efficiency");

			// good tracking?
			h1ENom[6][j]->Fill(eleP);
			h1ENom[6][j]->Fill(posP);
			if (!j) hNhits->Fill(eleTpcHits); 
			if (!j) hNhits->Fill(posTpcHits);
			if (eleTpcHits > cutnhitsf[cutSet]) h1ENom[7][j]->Fill(eleP);
			if (posTpcHits > cutnhitsf[cutSet]) h1ENom[7][j]->Fill(posP);
			h1ENom[7][j]->SetTitle("single e track quality;electron p [GeV/c];efficiency");
			if (eleTpcHits < cutnhitsf[cutSet] || posTpcHits < cutnhitsf[cutSet]) continue;
			if (!j) hNrat->Fill(eleTpcHits/eleTpcPos); 
			if (!j) hNrat->Fill(posTpcHits/posTpcPos);
			hNhitsE[j]->Fill(upsPt);
			hNom[4][j]->Fill(upsPt);
			hNom[4][j]->SetTitle("nHitsFit cut;upsilon p_{T}; efficiency");
			if (eleTpcHits/eleTpcPos < cutnratio[cutSet] || posTpcHits/posTpcPos < cutnratio[cutSet]) continue;
			hNratE[j]->Fill(upsPt);
			hNom[5][j]->Fill(upsPt);
			hNom[5][j]->SetTitle("nHitsRatio cut;upsilon p_{T};efficiency ");	

			//-----------------------------------------------------
			// +dca
			TVector3 eleRcFirstTpc(eleRcFirstTpcX,eleRcFirstTpcY,eleRcFirstTpcZ);
			TVector3 posRcFirstTpc(posRcFirstTpcX,posRcFirstTpcY,posRcFirstTpcZ);
			float eleDca = (primV - eleRcFirstTpc).Mag();
			float posDca = (primV - posRcFirstTpc).Mag();
			hDca[j]->Fill(eleDca);
			hDca[j]->Fill(posDca);		// ???	
			// momentum resolution study
			if (!j) {
				if (elePt != 0) hRcminMcPt->Fill(elePt,(elePtRc-elePt)/elePt);
      			hRcminMcPhi->Fill(elePt,elePhiRc-elePhi);
      			hRcminMcEta->Fill(elePt,eleEtaRc-eleEta); 
				if (posPt != 0) hRcminMcPt->Fill(posPt,(posPtRc-posPt)/posPt);
      			hRcminMcPhi->Fill(posPt,posPhiRc-posPhi);
      			hRcminMcEta->Fill(posPt,posEtaRc-posEta); 
      			/*cout << "elept " << elePt << " " << elePtRc << endl;
      			cout << "elephi " << elePhi << " " << elePhiRc << endl;
      			cout << "eleeta " << eleEta << " " << eleEtaRc << endl;*/
      		}
			// ep err study
			if (!j) {
			float eperr = (eleEcluster+gRandom->Gaus(0,Eerr(eleEcluster)))/(elePRc);
			hEpErr->Fill(eleEcluster/eleP,eperr);
			}
			//-----------------------------------------------------
					

			// emc
			h1ENom[8][j]->Fill(eleP);
			h1ENom[8][j]->Fill(posP);
			if (!j) hEp->Fill(eleEcluster/elePRc); 
			if (!j) hEp->Fill(posEcluster/posPRc);
			if (!j) hEMCEFFepUS->Fill(elePRc,eleEcluster/elePRc,1);
			if (!j) hEMCEFFepUS->Fill(posPRc,posEcluster/posPRc,1);
			if (j==3 && eleEcluster/elePRc > cutepl[cutSet] && eleEcluster/elePRc < cutepr[cutSet]) h1ENom[9][3]->Fill(eleP);
			if (j==4 && eleEcluster/elePRc > cuteplL[cutSet] && eleEcluster/elePRc < cuteprL[cutSet]) h1ENom[9][4]->Fill(eleP);
			if (j==5 && eleEcluster/elePRc > cuteplT[cutSet] && eleEcluster/elePRc < cuteprT[cutSet]) h1ENom[9][5]->Fill(eleP);
			if (j==3 && posEcluster/posPRc > cutepl[cutSet] && posEcluster/posPRc < cutepr[cutSet]) h1ENom[9][3]->Fill(posP);
			if (j==4 && posEcluster/posPRc > cuteplL[cutSet] && posEcluster/posPRc < cuteprL[cutSet]) h1ENom[9][4]->Fill(posP);
			if (j==5 && posEcluster/posPRc > cuteplT[cutSet] && posEcluster/posPRc < cuteprT[cutSet]) h1ENom[9][5]->Fill(posP);
			if (j!=4&&j!=5) if (eleEcluster/elePRc < cutepl[cutSet] || posEcluster/posPRc < cutepl[cutSet]) continue;
			if (j!=4&&j!=5) if (eleEcluster/elePRc > cutepr[cutSet] || posEcluster/posPRc > cutepr[cutSet]) continue;
			if (j==4) if (eleEcluster/elePRc < cuteplL[cutSet] || posEcluster/posPRc < cuteplL[cutSet]) continue;
			if (j==4) if (eleEcluster/elePRc > cuteprL[cutSet] || posEcluster/posPRc > cuteprL[cutSet]) continue;
			if (j==5) if (eleEcluster/elePRc < cuteplT[cutSet] || posEcluster/posPRc < cuteplT[cutSet]) continue;
			if (j==5) if (eleEcluster/elePRc > cuteprT[cutSet] || posEcluster/posPRc > cuteprT[cutSet]) continue;
			h1ENom[9][j]->SetTitle("single e E/p;electron p [GeV/c];efficiency");
			hEpE[j]->Fill(upsPt);
			hNom[6][j]->Fill(upsPt);
			hNom[6][j]->SetTitle("E/p cut;upsilon p_{T};efficiency");

			//r
			h1ENom[10][j]->Fill(eleP);
			h1ENom[10][j]->Fill(posP);
			if (!j) hR->Fill(eleR);
			if (!j) hR->Fill(posR);
			if (j==6 && eleR<cutrd[cutSet]) h1ENom[11][j]->Fill(eleP);
			if (j==6 && posR<cutrd[cutSet]) h1ENom[11][j]->Fill(posP);
			if (j==7 && eleR<cutrdL[cutSet]) h1ENom[11][j]->Fill(eleP);
			if (j==7 && posR<cutrdL[cutSet]) h1ENom[11][j]->Fill(posP);
			if (j==8 && eleR<cutrdT[cutSet]) h1ENom[11][j]->Fill(eleP);
			if (j==8 && posR<cutrdT[cutSet]) h1ENom[11][j]->Fill(posP);
			if (j!=7&&j!=8) if (eleR>cutrd[cutSet] || posR>cutrd[cutSet]) continue;
			if (j==7) if (eleR>cutrdL[cutSet] || posR>cutrdL[cutSet]) continue;
			if (j==8) if (eleR>cutrdT[cutSet] || posR>cutrdT[cutSet]) continue;
			hRE[j]->Fill(upsPt);
			hNom[7][j]->Fill(upsPt);
			hNom[7][j]->SetTitle("R cut;upsilon p_{T};efficiency");


			// good kinematics?
			h1ENom[12][j]->Fill(eleP);
			h1ENom[12][j]->Fill(posP);
			if (!j) hUpsPepP->Fill(upsPt,eleP,posP);
			if (!j) hEta->Fill(eleEtaRc); 
			if (!j) hEta->Fill(posEtaRc);
			if (fabs(eleEtaRc) > cutetal[cutSet] || fabs(posEtaRc) > cutetal[cutSet]) continue;
			hEtaE[j]->Fill(upsPt);
			hNom[8][j]->Fill(upsPt);
			hNom[8][j]->SetTitle("eta cut;upsilon p_{T}; efficiency");
			if (j!=10&&j!=11) if (elePRc < cutplow[cutSet] || posPRc < cutplow[cutSet]) continue;
			if (j==10) if (elePRc < cutplowL[cutSet] || posPRc < cutplowL[cutSet]) continue;
			if (j==11) if (elePRc < cutplowT[cutSet] || posPRc < cutplowT[cutSet]) continue;
			hPlowE[j]->Fill(upsPt);
			hNom[9][j]->Fill(upsPt);
			hNom[9][j]->SetTitle("low momentum cut;upsilon p_{T}; efficiency");
			if (j!=10&&j!=11) if (elePRc < cutplead[cutSet] && posPRc < cutplead[cutSet]) continue;
			if (j==10) if (elePRc < cutpleadL[cutSet] && posPRc < cutpleadT[cutSet]) continue;
			if (j==11) if (elePRc < cutpleadL[cutSet] && posPRc < cutpleadT[cutSet]) continue;
			hPleadE[j]->Fill(upsPt);
			hNom[10][j]->Fill(upsPt);
			hNom[10][j]->SetTitle("high momentum cut;upsilon p_{T};efficiency ");

			// dca
			if (j!=13&&j!=14) if (dcar1 > cutdcaprob || dcar2 > cutdcaprob) continue; 	// value taken from minimc
			if (j==13) if (dcar1 > cutdcaprobL || dcar2 > cutdcaprobL) continue;
			if (j==14) if (dcar1 > cutdcaprobT || dcar2 > cutdcaprobT) continue;
			hDcaE[j]->Fill(upsPt);
			hNom[11][j]->Fill(upsPt);
			hNom[11][j]->SetTitle("DCA cut;upsilon p_{T};efficiency ");	
	
			// pair
			TLorentzVector* emom = new TLorentzVector(elePxRc,elePyRc,elePzRc,elePRc*elePRc+elmass*elmass);
			TLorentzVector* pmom = new TLorentzVector(posPxRc,posPyRc,posPzRc,posPRc*posPRc+elmass*elmass);
			TLorentzVector* pair = new TLorentzVector(*emom+*pmom);
			if (pair->Rapidity() > cuturap[cutSet]) continue;
			hUrapE[j]->Fill(upsPt);
			hNom[12][j]->Fill(upsPt);
			hNom[12][j]->SetTitle("pair rapidity cut;upsilon p_{T}; efficiency");	

			// nsigeff
			if (j!=16&&j!=17) if (nsigr1 > efit->Eval(elePRc) || nsigr2 > efit->Eval(posPRc)) continue;
			if (j==16) if (nsigr1 > efitT->Eval(elePRc) || nsigr2 > efitT->Eval(posPRc)) continue;
			if (j==17) if (nsigr1 > efitL->Eval(elePRc) || nsigr2 > efitL->Eval(posPRc)) continue;
			
			hNsigE[j]->Fill(upsPt);
			hNom[13][j]->Fill(upsPt);
			hNom[13][j]->SetTitle("nSigmaE cut;upsilon p_{T};efficiency");		
				

			//printf("vz is %f ,", vz);
			//printf("vzRc is %f \n", vzRc);
			//if (eleReco==0) continue;
			//printf("elePt is %f ,", elePt);
			//printf("elePtRc is %f \n", elePtRc);
			
			hTotCnt[j]->Fill(2);
			hTotE[j]->Fill(upsPt);
		}

	}
	
	cout << "htot1 " << hTotE[1]->GetBinContent(hTotE[1]->FindBin(9.8)) << endl;
	cout << "htot0 " << hTotE[0]->GetBinContent(hTotE[1]->FindBin(9.8)) << endl;
	cout << "htot2 " << hTotE[2]->GetBinContent(hTotE[1]->FindBin(9.8)) << endl;


	for (int j = 0; j < nsyst; ++j)
	{
		hAccE[j]->Divide(hDen[j]);
		hTrigE[j]->Divide(hDen[j]);
		hRecoE[j]->Divide(hDen[j]);
		hNhitsE[j]->Divide(hDen[j]);
		hNratE[j]->Divide(hDen[j]);
		hEtaE[j]->Divide(hDen[j]);
		hPlowE[j]->Divide(hDen[j]);
		hPleadE[j]->Divide(hDen[j]);
		hEpE[j]->Divide(hDen[j]);
		hRE[j]->Divide(hDen[j]);
		hDcaE[j]->Divide(hDen[j]);
		hUrapE[j]->Divide(hDen[j]);
		hNsigE[j]->Divide(hDen[j]);
		hTotE[j]->Divide(hDen[j]);
	
		for (int i = nhistos-1; i > 0; --i) 
		{
  			hNom[i][j]->Divide(hNom[i-1][j]);
  		}
  	}

  	// calculate syst. error
  	int nerrs = 7;
  	float errT[nbinspt][nerrs];
  	float errL[nbinspt][nerrs];
  	float errEmc[nbinspt];

  	// fill flats
  	const float errAccT = 0.;		// assuming symmetric this time (orig 0.017, 0.03)
  	const float errAccL = 0.;		// lets add them in calcRaa

  	const float errTrackL = 0.0655;		// gives 11.6 for yield 
  	const float errTrackT = 0.0655; 	// gives 15.1 for yield (+HFT effect)

  	// fill default
  	for (int i = 1; i < nbinspt; ++i)
  	{
  		errL[i][0] = hTotE[0]->GetBinContent(i);
  		
  		errL[i][1] = (!errL[i][0]) ? 0 : (float)(hTotE[1]->GetBinContent(i) - errL[i][0])/errL[i][0];		//trig
  		errL[i][2] = (!errL[i][0]) ? 0 : (float)(hTotE[4]->GetBinContent(i) - errL[i][0])/errL[i][0];		//ep
  		errL[i][3] = (!errL[i][0]) ? 0 : (float)(hTotE[7]->GetBinContent(i) - errL[i][0])/errL[i][0];		//r
  		errL[i][4] = (!errL[i][0]) ? 0 : (float)(hTotE[10]->GetBinContent(i) - errL[i][0])/errL[i][0];		//kine
  		errL[i][5] = (!errL[i][0]) ? 0 : (float)(hTotE[13]->GetBinContent(i) - errL[i][0])/errL[i][0];		//dca
  		errL[i][6] = (!errL[i][0]) ? 0 : (float)(hTotE[16]->GetBinContent(i) - errL[i][0])/errL[i][0];		//nsig

  		errT[i][1] = (!errL[i][0]) ? 0 : -(float)(hTotE[2]->GetBinContent(i) - errL[i][0])/errL[i][0];		//trig
  		errT[i][2] = (!errL[i][0]) ? 0 : -(float)(hTotE[5]->GetBinContent(i) - errL[i][0])/errL[i][0];		//ep
  		errT[i][3] = (!errL[i][0]) ? 0 : -(float)(hTotE[8]->GetBinContent(i) - errL[i][0])/errL[i][0];		//r
  		errT[i][4] = (!errL[i][0]) ? 0 : -(float)(hTotE[11]->GetBinContent(i) - errL[i][0])/errL[i][0];		//kine
  		errT[i][5] = (!errL[i][0]) ? 0 : -(float)(hTotE[14]->GetBinContent(i) - errL[i][0])/errL[i][0];		//dca
  		errT[i][6] = (!errL[i][0]) ? 0 : -(float)(hTotE[17]->GetBinContent(i) - errL[i][0])/errL[i][0];		//nsig
  	}

  	//lets symmetrise here, actually no, cant symmetrise here yet (it's bin by bin)
  	/*for (int i = 0; i < nbinspt; ++i)
  	{
  		for (int j = 1; j < 7; ++j)
  		{
  			if (errL[i][j] > errT[i][j]) errT[i][j] = errL[i][j];
  			else errL[i][j] = errT[i][j];
  		}
  	}*/

  	for (int i = 1; i < nbinspt; ++i)
  	{
  		errL[i][0] = errAccL*errAccL + 4*errTrackL*errTrackL;
  		errT[i][0] = errAccT*errAccT + 4*errTrackT*errTrackT;
  		//cout << "--0 err L is " << sqrt(errL[i][0]) << " err T is " << sqrt(errT[i][0]) << endl;
  		

  		// for simplicity lets symmetrise here based on posteriori knowledge
  		errL[i][0] += errT[i][1]*errT[i][1];// + errT[i][2] + errT[i][3] + errL[i][4] + errL[i][5] + errL[i][6];
  		errL[i][0] += errT[i][2]*errT[i][2];
  		errL[i][0] += errT[i][3]*errT[i][3];
  		errL[i][0] += errL[i][4]*errL[i][4];
  		errL[i][0] += errL[i][5]*errL[i][5];
  		errL[i][0] += errL[i][6]*errL[i][6];

  		errT[i][0] += errT[i][1]*errT[i][1];// + errT[i][2] + errT[i][3] + errL[i][4] + errL[i][5] + errL[i][6];
  		errT[i][0] += errT[i][2]*errT[i][2];
  		errT[i][0] += errT[i][3]*errT[i][3];
  		errT[i][0] += errL[i][4]*errL[i][4];
  		errT[i][0] += errL[i][5]*errL[i][5];
  		errT[i][0] += errL[i][6]*errL[i][6];

  		/*for (int k = 1; k < 7; ++k)
  		{
  			//cout << "k " << k << " err L is " << errL[i][k] << " err T is " << errT[i][k] << endl;
  			//cout << "sq is " << errL[i][k]*errL[i][k] << endl;
  			//cout << errL[i][0] << " + " << errL[i][k]*errL[i][k] << " = " << errL[i][0]+errL[i][k]*errL[i][k] << endl;

  			errL[i][0] += errL[i][k]*errL[i][k];
  			errT[i][0] += errT[i][k]*errT[i][k];
  			//cout << "tot err L is " << sqrt(errL[i][0]) << " err T is " << sqrt(errT[i][0]) << endl;
  		}*/

  		errL[i][0] = sqrt(errL[i][0]);
  		errT[i][0] = sqrt(errT[i][0]);
  		//cout << "--FINAL err L is " << errL[i][0] << " err T is " << errT[i][0] << endl;
  	}  

  	cout << "SYMMETRISED" << endl;

  	// drawing here
  	

  	//distributions
  	TCanvas* c1 = new TCanvas("c1","",750,600);
  	c1->SetGridx();
  	hNhits->Draw("H");
  	hNhits->Draw("SAME");

  	Double_t x1[2], x1b[2], y[2];
	x1[0]=0; x1[1]=x1[0]; x1b[0]=cutnhitsf[cutSet]; x1b[1]=x1b[0];
   	c1->Update();
   	y[0]=gPad->GetUymin(); y[1]=gPad->GetUymax();
   	cout << "y0 " << y[0] << " y1 " << y[1] << endl;
   	
   
   	TGraph* g1 = new TGraph(2, x1, y);
   	g1->SetLineWidth(802);
   	g1->SetLineColor(kRed);
   	g1->SetFillColor(kRed);
   	g1->SetFillStyle(3003);
   	//g1->Draw("SAME");

   	TGraph* g1b = new TGraph(2, x1b, y);
   	g1b->SetLineWidth(-802);
   	g1b->SetLineColor(kRed);
   	g1b->SetFillColor(kRed);
   	g1b->SetFillStyle(3003);
   	g1b->Draw("SAME");
   	
   	c1->SaveAs(dir+"nhitsf.png");

   	hEta->Draw("H");
   	hEta->Draw("SAME");
   	c1->Update();
   	g1->SetPoint(0,cutetal[cutSet],gPad->GetUymin());
   	g1->SetPoint(1,cutetal[cutSet],gPad->GetUymax());
   	g1b->SetPoint(0,-cutetal[cutSet],gPad->GetUymin());
   	g1b->SetPoint(1,-cutetal[cutSet],gPad->GetUymax());
   	g1->Draw("SAME");
   	g1b->Draw("SAME");
   	c1->SaveAs(dir+"eta.png");

   	hEp->Draw("H");
   	hEp->Draw("SAME");
   	c1->Update();
   	g1->SetPoint(0,cutepr[cutSet],gPad->GetUymin());
   	g1->SetPoint(1,cutepr[cutSet],gPad->GetUymax());
   	g1b->SetPoint(0,cutepl[cutSet],gPad->GetUymin());
   	g1b->SetPoint(1,cutepl[cutSet],gPad->GetUymax());
   	g1->Draw("SAME");
   	g1b->Draw("SAME");
   	c1->SaveAs(dir+"ep.png");

   	hR->Draw("H");
   	hR->Draw("SAME");
   	c1->Update();
   	g1->SetPoint(0,cutrd[cutSet],gPad->GetUymin());
   	g1->SetPoint(1,cutrd[cutSet],gPad->GetUymax());
   	g1->Draw("SAME");
   	c1->SaveAs(dir+"r.png");

   	hUpsPepP->GetXaxis()->SetRange(1,hUpsPepP->GetXaxis()->FindBin(4.0));
   	hUpsPepP->Project3D("zy")->Draw("colz");
   	g1->SetPoint(0,cutplow[cutSet],15);
   	g1->SetPoint(1,cutplow[cutSet],0);
   	g1b->SetPoint(0,15,cutplow[cutSet]);
   	g1b->SetPoint(1,0,cutplow[cutSet]);
   	g1->Draw("SAME");
   	g1b->Draw("SAME");
   	c1->SaveAs(dir+"plowpt.png");

   	hUpsPepP->GetXaxis()->SetRange(hUpsPepP->GetXaxis()->FindBin(6.0),150);
   	hUpsPepP->Project3D("zy")->Draw("colz");
   	g1->SetPoint(0,cutplow[cutSet],15);
   	g1->SetPoint(1,cutplow[cutSet],0);
   	g1b->SetPoint(0,15,cutplow[cutSet]);
   	g1b->SetPoint(1,0,cutplow[cutSet]);
   	g1->Draw("SAME");
   	g1b->Draw("SAME");
   	c1->SaveAs(dir+"phighpt.png");

   	// different systematics
  	TH1F* hSysL[nerrs+1];
  	TH1F* hSysT[nerrs+1];
  	for (int k = 0; k < nerrs+1; ++k)
  	{
  		hSysL[k] = new TH1F(Form("hSysL_%i",k),"",nbinspt,hleft,hright);
  		hSysT[k] = new TH1F(Form("hSysT_%i",k),"",nbinspt,hleft,hright);
  		//hSysL[k]->Sumw2();
  	}

  	hSysL[0]->SetLineColor(kBlack);
  	hSysL[0]->SetMarkerColor(kBlack);
  	hSysL[1]->SetLineColor(kRed);
  	hSysL[1]->SetMarkerColor(kRed);
  	hSysL[2]->SetLineColor(kGreen+2);
  	hSysL[2]->SetMarkerColor(kGreen+2);
  	hSysL[3]->SetLineColor(kGreen+4);
  	hSysL[3]->SetMarkerColor(kGreen+4);
  	hSysL[4]->SetLineColor(kCyan+1);
  	hSysL[4]->SetMarkerColor(kCyan+1);
  	//hSysL[5]->SetLineColor(kBlack);
  	//hSysL[5]->SetMarkerColor(kBlack);
  	//hSysL[6]->SetLineColor(kBlue+2);
  	//hSysL[6]->SetMarkerColor(kBlue+2);
  	hSysL[7]->SetLineColor(kBlue+2);
  	hSysL[7]->SetMarkerColor(kBlue+2);

  	for (int k = 0; k < nerrs; ++k)
  	{
  		for (int i = 1; i < nbinspt; ++i)
  		{
  			hSysL[k]->SetBinContent(i,errT[i][k]*100);
  			hSysL[k]->SetBinError(i,0.001);
  			hSysT[k]->SetBinContent(i,errT[i][k]*100);
  			hSysT[k]->SetBinError(i,0.001);
  		}
  	}

	for (int i = 1; i < nbinspt; ++i) // fill tpc flat
  		{
  			hSysL[7]->SetBinContent(i,2*errTrackL*100);
  			hSysL[7]->SetBinError(i,0.001);
  			hSysT[7]->SetBinContent(i,2*errTrackT*100);
  			hSysT[7]->SetBinError(i,0.001);
  		}


  	ofstream outfile_sys;
  	TString sysfile = "sys_";
  	sysfile += outFile;
  	outfile_sys.open(sysfile);
  	double sysErrL[8];
  	TF1* fp0 = new TF1("fp0","pol0",hleft,hright);
  	for (int iF = 1; iF < 7; ++iF)
  	{
  		hSysL[iF]->Fit("fp0");
  		double result = fp0->GetParameter(0);
  		result = 0.01*result;
  		double result2 = 1./(1-result)-1;
  		result = 1-1./(1+result);
  		result = 0.5*result + 0.5*result2;
  		outfile_sys << result2 << endl;
  	}
  	outfile_sys.close();

   	hSysL[0]->GetYaxis()->SetRangeUser(0.0,20.0);
   	hSysL[0]->GetXaxis()->SetRangeUser(0.0,10.0);
  	hSysL[0]->SetFillColorAlpha(kBlack,0.05);	// total
  	hSysL[0]->SetTitle("systematic uncertainty contributions;upsilon p_{T};percentage");
  	hSysL[0]->Draw("][ hist");
  	hSysL[0]->Draw("same");
  	hSysL[7]->SetFillColorAlpha(kBlue+2,0.05);	// track
  	hSysL[7]->Draw("][ hist same");
  	hSysL[7]->Draw("same");
  	hSysL[1]->SetFillColorAlpha(kRed,0.05); // trigger
  	hSysL[1]->Draw("][ hist same");
  	hSysL[1]->Draw("same");
  	cout << "Trigger L" << endl;
  	hSysL[1]->Fit("pol0");
  	hSysL[2]->SetFillColorAlpha(kGreen+2,0.05); //ep
  	hSysL[2]->Draw("][ hist same");
  	hSysL[2]->Draw("same");
  	cout << "Ep L" << endl;
  	hSysL[2]->Fit("pol0");
  	hSysL[3]->SetFillColorAlpha(kGreen+4,0.05); //matching
  	hSysL[3]->Draw("][ hist same");
  	hSysL[3]->Draw("same");
  	cout << "R L" << endl;
  	hSysL[3]->Fit("pol0");
  	cout << "drawing" << endl;
  	hSysL[4]->SetFillColorAlpha(kCyan+1,0.05); // kine
  	hSysL[4]->Draw("][ hist same");
  	hSysL[4]->Draw("same");
  	cout << "Kine L" << endl;
  	hSysL[4]->Fit("pol0");


  	cout << "Trigger T" << endl;
  	hSysT[1]->Fit("pol0");

  	cout << "Ep T" << endl;
  	hSysT[2]->Fit("pol0");

  	cout << "R T" << endl;
  	hSysT[3]->Fit("pol0");

  	cout << "Kine T" << endl;
  	hSysT[4]->Fit("pol0");

  	TLegend* l3 = new TLegend(0.45,0.74,0.96,0.89);
   	l3->SetNColumns(2);
   	l3->SetFillStyle(0);
   	l3->SetBorderSize(0);
   	l3->SetTextSize(0.032);
   	l3->AddEntry(hSysL[0],"total","pl");
   	l3->AddEntry(hSysL[1],"trigger","pl");
   	l3->AddEntry(hSysL[2],"E/p","pl");
   	l3->AddEntry(hSysL[3],"BEMC matching","pl");
   	l3->AddEntry(hSysL[4],"kinematics","pl");
   	l3->AddEntry(hSysL[7],"TPC tracking","pl");
  	l3->Draw();
  	c1->SaveAs(dir+"systs.png");

  	hSysL[0]->Fit("pol0");


   	//efficiencies
   	TLegend* l1 = new TLegend(0.14,0.65,0.86,0.85);
   	l1->SetNColumns(2);
   	l1->SetFillStyle(0);
   	l1->SetBorderSize(0);
   	l1->SetTextSize(0.032);

   	l1->AddEntry((TObject*)0,"#bf{Au+Au 2014 #sqrt{s_{NN}} = 200 GeV}"," ");
   	l1->AddEntry((TObject*)0," "," ");

   	//c1->SetLogy();
   	c1->SetGridx();
   	hAccE[0]->GetYaxis()->SetRangeUser(0.0,0.75);
	hAccE[0]->SetTitle(";upsilon p_{T} [GeV/c];efficiency");
   	hAccE[0]->SetLineColor(kOrange+2);
   	hAccE[0]->SetMarkerColor(kOrange+2);
   	l1->AddEntry(hAccE[0],"acceptance","pl");
   	hAccE[0]->Draw();
   	//l1->Draw();
   	c1->SaveAs(dir+"eff0.png");

   	//hTrigE->SetTitle(";upsilon p_{T};efficiency");
   	l1->AddEntry(hTrigE[0],"trigger","pl");
   	hTrigE[0]->Draw("same");
   	//l1->Draw();
   	c1->SaveAs(dir+"eff1.png");

   	hRecoE[0]->SetLineColor(kMagenta+2);
   	hRecoE[0]->SetMarkerColor(kMagenta+2);
   	l1->AddEntry(hRecoE[0],"daughters reconstructed","pl");
   	hRecoE[0]->Draw("same");
   	//l1->Draw();
   	c1->SaveAs(dir+"eff2.png");

   	hNratE[0]->SetLineColor(kBlue+2);
   	hNratE[0]->SetMarkerColor(kBlue+2);
   	l1->AddEntry(hNratE[0],"track quality","pl");
   	hNratE[0]->Draw("same");
   	//l1->Draw();
   	c1->SaveAs(dir+"eff3.png");

   	hPleadE[0]->SetLineColor(kCyan+1);
   	hPleadE[0]->SetMarkerColor(kCyan+1);
   	l1->AddEntry(hPleadE[0],"kinematics","pl");
   	hPleadE[0]->Draw("same");
   	//l1->Draw();
   	c1->SaveAs(dir+"eff4.png");

   	hRE[0]->SetLineColor(kGreen+2);
   	hRE[0]->SetMarkerColor(kGreen+2);
   	l1->AddEntry(hRE[0],"EMC pid","pl");
   	hRE[0]->Draw("same");
   	//l1->Draw();
   	c1->SaveAs(dir+"eff5.png");


   	hTotE[0]->SetLineColor(kBlack);
   	hTotE[0]->SetMarkerColor(kBlack);
   	l1->AddEntry(hTotE[0],"total reconstruction","pl");
   	hTotE[0]->Draw("same");
   	l1->Draw();
   	c1->SaveAs(dir+"eff6.png");

   	//separate efficiencies
   	//apply summed errors
   	for (int i = 1; i < nbinspt; ++i)
   	{
   		float bc = hTotE[0]->GetBinContent(i);
   		float err = hTotE[0]->GetBinError(i);
   		hTotEL->SetBinContent(i,bc + bc*errL[i][0]);
   		hTotEL->SetBinError(i,err);
   		hTotET->SetBinContent(i,bc - bc*errT[i][0]);
   		hTotET->SetBinError(i,err);
   	}
   	hTotE[0]->SetTitle(";upsilon p_{T};reconstruction efficiency");
   	hTotE[0]->GetYaxis()->SetRangeUser(0.0,0.1);
   	hTotEL->SetFillColorAlpha(kBlack, 0.20);
   	hTotET->SetFillColorAlpha(kWhite, 1.0);
   	hTotEL->SetLineColorAlpha(kWhite, 0.0);
   	hTotET->SetLineColorAlpha(kWhite, 0.0);
   	hTotE[0]->Draw();
   	hTotEL->Draw("hist same");
   	hTotE[0]->Draw("same");
   	hTotET->Draw("hist same");
   	hTotE[0]->Draw("axis axig same");
   	gPad->RedrawAxis();
   	c1->SaveAs(dir+"effTOT.png");

   	// SAVE EFFI IN FILE
   	ofstream outfile;
  	outfile.open(outFile);
   	TF1* fpol0 = new TF1("fpol0","pol0",hleft,hright);
   	float center, allErr;

   	hTotE[0]->Fit("fpol0");
   	center = fpol0->GetParameter(0);
  	outfile << center << endl;
   	hTotEL->Fit("fpol0");
   	allErr = (fpol0->GetParameter(0)-center)*(fpol0->GetParameter(0)-center)+(fpol0->GetParError(0)*fpol0->GetParError(0));
   	outfile << sqrt(allErr) << endl;
   	hTotET->Fit("fpol0");
   	allErr = (fpol0->GetParameter(0)-center)*(fpol0->GetParameter(0)-center)+(fpol0->GetParError(0)*fpol0->GetParError(0));
	outfile << sqrt(allErr) << endl;
	outfile.close();
   	

   	//reco
   	for (int i = 0; i < nhistos; ++i)
   	{
   		hNom[i][0]->Draw();
   		c1->SaveAs(dir+Form("sepeff%i.png",i));
   	}

   	//qa
   	ups->Draw("upsPt>>hUpt");
   	ups->Draw("upsPhi>>hUphi");
   	ups->Draw("upsRap>>hUrap");

   	hUpt->Draw("H");
   	hUpt->Draw("SAME");
   	c1->SaveAs(dir+"upt.png");
   	hUphi->Draw("H");
   	hUphi->Draw("SAME");
   	c1->SaveAs(dir+"uphi.png");
   	hUrap->Draw("H");
   	hUrap->Draw("SAME");
   	c1->SaveAs(dir+"urap.png");

   	// single electron efficiencies

   	for (int j = 0; j < nsyst; ++j)
	{	
		for (int i = 12; i > 0; --i) 
		{
  			h1ENom[i][j]->Divide(h1ENom[i-1][j]);
  			h1ENom[i][j]->GetYaxis()->SetRangeUser(0.0,1.0);
  		}
  	}

   	for (int j = 0; j < nsyst; ++j)
   	{
   		h1ETrig[j]->Divide(h1ENom[1][j]);
   		h1EEp[j]->Divide(h1ENom[5][3+j]);
   	}


   	hNom[2][1]->SetFillColorAlpha(kRed,0.15);
   	hNom[2][2]->SetFillColorAlpha(kWhite,1.0);
   	hNom[2][1]->SetLineColorAlpha(kWhite,0.0);
   	hNom[2][2]->SetLineColorAlpha(kWhite,0.0);
   	hNom[2][0]->Draw();
	hNom[2][1]->Draw("hist same");
	hNom[2][0]->Draw("same");
   	hNom[2][2]->Draw("hist same");
   	hNom[2][0]->Draw("axis axig same");
   	c1->SaveAs(dir+"trigsys.png");


   	//fill 1e errors

   	h1ENom[1][1]->SetFillColorAlpha(kOrange+2,0.15);		//sys 012
   	h1ENom[1][2]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[1][1]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[1][2]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[1][0]->SetLineColor(kOrange+2);
   	h1ENom[1][0]->SetMarkerColor(kOrange+2);
   	h1ENom[1][1]->SetTitle(";electron p [GeV/c];acceptance");
   	h1ENom[1][1]->Draw("hist");
   	h1ENom[1][0]->Draw("same");
   	h1ENom[1][2]->Draw("hist same");
   	h1ENom[1][0]->Draw("same");
   	h1ENom[1][0]->Draw("axis same");
   	h1ENom[1][0]->Draw("axig same");
   	c1->SaveAs(dir+"1Eacc.png");

   	h1ENom[3][1]->SetFillColorAlpha(kRed,0.15);		//sys 012
   	h1ENom[3][2]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[3][1]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[3][2]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[3][0]->SetLineColor(kRed);
   	h1ENom[3][0]->SetMarkerColor(kRed);
   	h1ENom[3][1]->SetTitle(";electron p [GeV/c];trigger efficiency");
   	h1ENom[3][1]->Draw("hist");
   	h1ENom[3][0]->Draw("same");
   	h1ENom[3][2]->Draw("hist same");
   	h1ENom[3][0]->Draw("same");
   	h1ENom[3][0]->Draw("axis same");
   	h1ENom[3][0]->Draw("axig same");
   	c1->SaveAs(dir+"1Etrig.png");

   	h1ENom[5][1]->SetFillColorAlpha(kMagenta+2,0.15);		//sys 012
   	h1ENom[5][2]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[5][1]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[5][2]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[5][0]->SetLineColor(kMagenta+2);
   	h1ENom[5][0]->SetMarkerColor(kMagenta+2);
   	h1ENom[5][1]->SetTitle(";electron p [GeV/c];track reconstruction");
   	h1ENom[5][1]->Draw("hist");
   	h1ENom[5][0]->Draw("same");
   	h1ENom[5][2]->Draw("hist same");
   	h1ENom[5][0]->Draw("same");
   	h1ENom[5][0]->Draw("axis same");
   	h1ENom[5][0]->Draw("axig same");
   	c1->SaveAs(dir+"1Ereco.png");

   	h1ENom[7][1]->Scale(1+errTrackL); 	// = 0.058
   	h1ENom[7][2]->Scale(1-errTrackL);
   	TLegend* l2 = new TLegend(0.14,0.14,0.46,0.30);
   	l2->SetNColumns(1);
   	l2->SetFillStyle(0);
   	l2->SetBorderSize(0);
   	l2->SetTextSize(0.032);
   	l2->AddEntry(h1ENom[5][0],"track reconstruction","pl");

   	h1ENom[7][1]->SetFillColorAlpha(kBlue+2,0.15);		//sys 012
   	h1ENom[7][2]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[7][1]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[7][2]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[7][0]->SetLineColor(kBlue+2);
   	h1ENom[7][0]->SetMarkerColor(kBlue+2);
   	h1ENom[7][1]->SetTitle(";electron p [GeV/c];track quality");
   	l2->AddEntry(h1ENom[7][0],"track quality","pl");
   	h1ENom[7][1]->Multiply(h1ENom[5][1]);
   	h1ENom[7][0]->Multiply(h1ENom[5][0]);
   	h1ENom[7][2]->Multiply(h1ENom[5][2]);
   	h1ENom[7][1]->Draw("hist same");
   	h1ENom[7][0]->Draw("same");
   	h1ENom[7][2]->Draw("hist same");
   	h1ENom[7][0]->Draw("same");
   	h1ENom[7][0]->Draw("axis same");
   	h1ENom[7][0]->Draw("axig same");
   	l2->Draw();
   	c1->SaveAs(dir+"1Etrac.png");


   	h1ENom[9][4]->SetFillColorAlpha(kGreen+2,0.25);		//sys 
   	h1ENom[9][5]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[9][4]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[9][5]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[9][3]->SetLineColor(kGreen+2);
   	h1ENom[9][3]->SetMarkerColor(kGreen+2);
   	h1ENom[9][4]->SetTitle(";electron p [GeV/c];E/p cut efficiency");
   	h1ENom[9][4]->Draw("hist");
   	h1ENom[9][3]->Draw("same");
   	h1ENom[9][5]->Draw("hist same");
   	h1ENom[9][3]->Draw("same");
   	h1ENom[9][3]->Draw("axis same");
   	h1ENom[9][3]->Draw("axig same");
   	c1->SaveAs(dir+"1Eep.png");

   	h1ENom[11][7]->SetFillColorAlpha(kGreen+4,0.15);		//sys 012
   	h1ENom[11][8]->SetFillColorAlpha(kWhite,1.0);
   	h1ENom[11][7]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[11][8]->SetLineColorAlpha(kWhite,0.0);
   	h1ENom[11][6]->SetLineColor(kGreen+4);
   	h1ENom[11][6]->SetMarkerColor(kGreen+4);
   	h1ENom[11][7]->SetTitle(";electron p [GeV/c];BEMC matching efficiency");
   	h1ENom[11][7]->Draw("hist");
   	h1ENom[11][6]->Draw("same");
   	h1ENom[11][8]->Draw("hist same");
   	h1ENom[11][6]->Draw("same");
   	h1ENom[11][6]->Draw("axis same");
   	h1ENom[11][6]->Draw("axig same");
   	c1->SaveAs(dir+"1Er.png");

   	h1ETrigPhi1->SetTitle("e pt < 6; electron phi; trigger efficiency");
   	//h1ETrigPhi1->Divide(h1ETrigPhi10);
   	//h1ETrigPhi1->GetYaxis()->SetRangeUser(0.0,1.0);
   	h1ETrigPhi1->Draw("hist");
   	c1->SaveAs(dir+"1eTrPhi1.png");
   	h1ETrigPhi2->SetTitle("e pt > 11; electron phi; trigger efficiency");
   	//h1ETrigPhi2->Divide(h1ETrigPhi20);
   	//h1ETrigPhi2->GetYaxis()->SetRangeUser(0.0,1.0);
   	h1ETrigPhi2->Draw("hist");
   	c1->SaveAs(dir+"1eTrPhi2.png");
   	h1ERecoPhi1->SetTitle("e pt < 6; electron phi; reco efficiency");
   	h1ERecoPhi1->Divide(h1ERecoPhi10);
   	h1ERecoPhi1->GetYaxis()->SetRangeUser(0.0,1.0);
   	h1ERecoPhi1->Draw("hist");
   	c1->SaveAs(dir+"1eRePhi1.png");
   	h1ERecoPhi2->SetTitle("e pt > 11; electron phi; reco efficiency");
   	h1ERecoPhi2->Divide(h1ERecoPhi20);
   	h1ERecoPhi2->GetYaxis()->SetRangeUser(0.0,1.0);
   	h1ERecoPhi2->Draw("hist");
   	c1->SaveAs(dir+"1eRePhi2.png");

   	hRcminMcPt->Draw("colz");
   	c1->SaveAs(dir+"rcmcpt.png");
   	hRcminMcPhi->Draw("colz");
   	c1->SaveAs(dir+"rcmcphi.png");
   	hRcminMcEta->Draw("colz");
   	c1->SaveAs(dir+"rcmceta.png");

}