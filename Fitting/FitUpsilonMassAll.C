static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

using namespace RooFit;

Int_t nBins = 40;
Float_t fitMin = 6.0;
Float_t fitMax = 14.0;

double yields[4];
double yields2[4];



void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}

void FitUpsilonMass(Int_t centralityMin = 7, Int_t centralityMax = 9) {

TFile *DT_file = TFile::Open("masstree.root", "read");

TFile *fTemplates1S = TFile::Open("hEmb1S.root", "read");
TFile *fTemplates2S = TFile::Open("hEmb2S.root", "read");
TFile *fTemplates3S = TFile::Open("hEmb3S.root", "read");

TFile *fRbg = TFile::Open("hRbg.root","read");

TH2F *hTemplate1S_2D = (TH2F*) fTemplates1S->Get("hUnIM");
TH1D *hTemplate1S = hTemplate1S_2D->ProjectionY("hTemplate1S");
hTemplate1S->Rebin();

TH2F *hTemplate2S_2D = (TH2F*) fTemplates2S->Get("hUnIM");
TH1D *hTemplate2S = hTemplate2S_2D->ProjectionY("hTemplate2S");
hTemplate2S->Rebin();

TH2F *hTemplate3S_2D = (TH2F*) fTemplates3S->Get("hUnIM");
TH1D *hTemplate3S = hTemplate3S_2D->ProjectionY("hTemplate3S");
hTemplate3S->Rebin();

TH1D* hRbg = (TH1D*)fRbg->Get("hRbg");
hRbg->Rebin();
hRbg->Smooth(1000);
//hRbg->Rebin();


//________________________________________________DATA_________________________________________________________
Int_t fCentrality;
Float_t fMass;
Int_t fSign;
Float_t fM;

Int_t centValues[10] = {80,70,60,50,40,30,20,10,5,0};

TTree *massTree = (TTree*)DT_file->Get("massTree");
massTree->SetBranchAddress("fMass",&fMass);
massTree->SetBranchAddress("fCentrality",&fCentrality);
massTree->SetBranchAddress("fSign",&fSign);

UInt_t nDT = massTree->GetEntries();

TH1D *hDiLeptonMassDT_Hist = new TH1D("hDiLeptonMassDT_Hist","Invariant mass of dileptons candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_Hist->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}) (GeV/c)");
hDiLeptonMassDT_Hist->Sumw2();

TH1D *hDiLeptonMassDT_LS = new TH1D("hDiLeptonMassDT_LS","Invariant mass like-sign candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_LS->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{+} l^{-}l^{-}) (GeV/c)");
hDiLeptonMassDT_LS->Sumw2();

TTree *fTreeDileptonMass = new TTree("fTreeDileptonMass", "fTreeDileptonMass");
fTreeDileptonMass->Branch("MassDT", &fM, "MassDT/F");
TTree *fTreeDileptonMassLS = new TTree("fTreeDileptonMassLS", "fTreeDileptonMassLS");
fTreeDileptonMassLS->Branch("MassDT", &fM, "MassDT/F");

Bool_t ZNA0n, ZNC0n;

    //event loop
int howmany = 0;
for(Int_t i=0; i<nDT; i++){
	massTree->GetEntry(i);
	if (fCentrality<centralityMin || fCentrality>=centralityMax) continue;
	fM = fMass;
	if (fSign < 0) {
		fTreeDileptonMass->Fill();
		hDiLeptonMassDT_Hist->Fill(fM); }
	else {
    fTreeDileptonMassLS->Fill();
    hDiLeptonMassDT_LS->Fill(fM);
  }
	howmany++;
}
cout << "howmany is " << howmany << endl;

// Declare observable - Mass for data
RooRealVar w("w","w",-2,2);
w=1;	
RooRealVar MassDT("MassDT","M_{ ee} (GeV/#it{c}^{2})",fitMin,fitMax);
RooDataSet DT_UnBin("DT_UnBin","DT_UnBin",RooArgSet(MassDT,w),Import(*fTreeDileptonMass),WeightVar(w));
//RooRealVar MassDTLS("MassDTLS","M_{ l^{+-}l^{+-}} (GeV/#it{c}^{2})",fitMin,fitMax);
//w=-1;
RooDataSet DT_UnBinLS("DT_UnBinLS","DT_UnBinLS",RooArgSet(MassDT,w),Import(*fTreeDileptonMassLS),WeightVar("w"));
RooDataHist DT_histLS("DT_histLS","DT_histLS",MassDT,Import(*hDiLeptonMassDT_LS));
RooDataHist DT_hist("DT_hist","DT_hist",MassDT,Import(*hDiLeptonMassDT_Hist));

//US-LS
/*RooRealVar w("w","w",-1,1) ;
w.setVal(1); 
DT_UnBin.addColumn(w);
DT_UnBin.setWeightVar(w);
w.setVal(-1); 
DT_UnBinLS.addColumn(w);
DT_UnBinLS.setWeightVar(w);

// Append datasets and interpret w as event weight*/
//DT_UnBin.append(DT_UnBinLS);

// Declare observable - Masses for MC
RooRealVar MassMC_Signal_1S("MassMC_Signal_1S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",fitMin,fitMax);
RooRealVar MassMC_Signal_2S("MassMC_Signal_2S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",fitMin,fitMax);
RooRealVar MassMC_Signal_3S("MassMC_Signal_3S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",fitMin,fitMax);

RooDataHist MC_hist_Signal_1S("MC_hist_Signal_1S","MC_hist_Signal_1S",MassMC_Signal_1S,Import(*hTemplate1S));
RooDataHist MC_hist_Signal_2S("MC_hist_Signal_2S","MC_hist_Signal_2S",MassMC_Signal_2S,Import(*hTemplate2S));
RooDataHist MC_hist_Signal_3S("MC_hist_Signal_3S","MC_hist_Signal_3S",MassMC_Signal_3S,Import(*hTemplate3S));

// Histogram for rbg
RooDataHist DT_Background("DT_Background","DT_Background",MassDT,Import(*hRbg));

// Crystal Ball P.D.F. variables for 1S,2S,3S
RooRealVar Mean_1S("Mean_1S","m_{1S}",9.388,8,12);
Mean_1S.setError(0.011);
//RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.1938,0.05,0.4);
RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.2289,0.05,0.4);
Sigma_1S.setError(0.009);
RooRealVar Alpha_1S("Alpha_1S","Alpha_1S",1.331,0.0,2);
Alpha_1S.setError(0.1015);
RooRealVar N_1S("N_1S","N_1S",0.9619,0,4);
N_1S.setError(0.224);

RooRealVar Mean_2S("Mean_2S","m_{2S}",9.934,8,12);
Mean_2S.setError(0.010);
//RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2041,0.05,0.4);
RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2540,0.05,0.4);
Sigma_2S.setError(0.0078);
RooRealVar Alpha_2S("Alpha_2S","Alpha_2S",1.285,0.0,2);
Alpha_2S.setError(0.0848);
RooRealVar N_2S("N_2S","N_2S",1.279,0,4);
N_2S.setError(0.160);

RooRealVar Mean_3S("Mean_3S","m_{3S}",10.27,8,12);
Mean_3S.setError(0.010);
//RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2303,0.05,0.4);
RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2655,0.05,0.4);
Sigma_3S.setError(0.0088);
RooRealVar Alpha_3S("Alpha_3S","Alpha_3S",1.073,0.0,2);
Alpha_3S.setError(0.090);
RooRealVar N_3S("N_3S","N_3S",1.035,0,4);
N_3S.setError(0.163);


//Fit Crystal Ball to 1S MC signal
RooCBShape CrystalBall_1S_MC("CrystalBall_1S_MC","CrystalBall_1S_MC",MassMC_Signal_1S,Mean_1S,Sigma_1S,Alpha_1S,N_1S);
//CrystalBall_1S_MC.fitTo(MC_hist_Signal_1S,Range(5,15));
TCanvas* cMC_1S = new TCanvas("cMC_1S","cMC_1S",700,700) ;
cMC_1S->cd();
RooPlot* plot_MC_1S = MassMC_Signal_1S.frame(Title(" "));
MC_hist_Signal_1S.plotOn(plot_MC_1S);
CrystalBall_1S_MC.plotOn(plot_MC_1S,Name("CrystalBall"),LineColor(myDarkBlue));
plot_MC_1S->Draw();

TLegend *myLegend_1S = new TLegend(0.41,0.55,0.893,0.85);
myLegendSetUp(myLegend_1S,0.03,1);
myLegend_1S->AddEntry((TObject*)0,TString::Format("m_{1S} = %5.3f #pm %5.3f GeV/#it{c}^{2}",Mean_1S.getVal(),Mean_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("#sigma_{1S} = %5.3f #pm %5.3f GeV/#it{c}^{2}",Sigma_1S.getVal(),Sigma_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("#alpha_{1S} = %5.3f #pm %5.3f",Alpha_1S.getVal(),Alpha_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("n_{1S} = %5.3f #pm %5.3f",N_1S.getVal(),N_1S.getError())," ");
myLegend_1S->Draw();
cMC_1S->SaveAs("cMC_1S.png");

//Fit Crystal Ball to 2S MC signal
RooCBShape CrystalBall_2S_MC("CrystalBall_2S_MC","CrystalBall_2S_MC",MassMC_Signal_2S,Mean_2S,Sigma_2S,Alpha_2S,N_2S);
//CrystalBall_2S_MC.fitTo(MC_hist_Signal_2S,Range(5.0,15));

//Fit Crystal Ball to 3S MC signal
RooCBShape CrystalBall_3S_MC("CrystalBall_3S_MC","CrystalBall_3S_MC",MassMC_Signal_3S,Mean_3S,Sigma_3S,Alpha_3S,N_3S);
//CrystalBall_3S_MC.fitTo(MC_hist_Signal_3S,Range(5,15));

// Fit a sum of the functions to the data
RooRealVar nSignal_1S("nSignal_1S","N_{1S}",1,0,1e06);
RooRealVar nSignal_2S("nSignal_2S","N_{2S}",1,0,1e06);
RooRealVar nSignal_3S("nSignal_3S","N_{3S}",1,0,1e06);
RooRealVar nSignal_2S3S("nSignal_2S3S","N_{2S3S}",1,0,1e06);
RooRealVar nBackground("nBackground","nBackground",1,0,1e06);
RooRealVar nCB2("nCB2","nCB2",1,0,1e06);

RooRealVar Sigma_1S_DT("Sigma_1S_DT","Sigma_1S_DT",Sigma_1S.getVal(),Sigma_1S.getVal()-0*Sigma_1S.getError(),Sigma_1S.getVal()+0*Sigma_1S.getError());
RooRealVar Sigma_2S_DT("Sigma_2S_DT","Sigma_2S_DT",Sigma_2S.getVal(),Sigma_2S.getVal()-0*Sigma_2S.getError(),Sigma_2S.getVal()+0*Sigma_2S.getError());
RooRealVar Sigma_3S_DT("Sigma_3S_DT","Sigma_3S_DT",Sigma_3S.getVal(),Sigma_3S.getVal()-0*Sigma_3S.getError(),Sigma_3S.getVal()+0*Sigma_3S.getError());

RooCBShape CrystalBall_1S_DT("CrystalBall_1S_DT","CrystalBall_1S_DT",MassDT,RooConst(Mean_1S.getVal()),Sigma_1S_DT,RooConst(Alpha_1S.getVal()),RooConst(N_1S.getVal()));
RooCBShape CrystalBall_2S_DT("CrystalBall_2S_DT","CrystalBall_2S_DT",MassDT,RooConst(Mean_2S.getVal()),Sigma_2S_DT,RooConst(Alpha_2S.getVal()),RooConst(N_2S.getVal()));
RooCBShape CrystalBall_3S_DT("CrystalBall_3S_DT","CrystalBall_3S_DT",MassDT,RooConst(Mean_3S.getVal()),Sigma_3S_DT,RooConst(Alpha_3S.getVal()),RooConst(N_3S.getVal()));
RooRealVar BLAOne("BLAOne","BLAOne",689,689,689);
RooRealVar BLATwo("BLATwo","BLATwo",311,311,311);
//RooRealVar BLAOne("BLAOne","BLAOne",0.689,0.689,0.689);
//RooRealVar BLATwo("BLATwo","BLATwo",0.311,0.311,0.311);
RooAddPdf CrystalBall_2S3S_DT("CrystalBall_2S3S_DT","CrystalBall_2S3S_DT",RooArgList(CrystalBall_2S_DT,CrystalBall_3S_DT),RooArgList(BLAOne,BLATwo));

//Combinatorial background
RooRealVar CB_DT_ParA("CB_DT_ParA","CB_DT_ParA",1,-100,100);
RooRealVar CB_DT_ParB("CB_DT_ParB","CB_DT_ParB",0.1,-100,100);
RooRealVar CB_DT_ParC("CB_DT_ParC","CB_DT_ParC",-0.1,-100,100);
RooRealVar CB_DT_ParD("CB_DT_ParD","CB_DT_ParD",-0.1,-100,100);
RooRealVar CB_DT_ParE("CB_DT_ParE","CB_DT_ParE",-0.1,-100,100);
RooChebychev Pol3("Pol3","Pol3",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC));
RooChebychev Pol5("Pol5","Pol5",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
RooChebychev CBShapeLS(Pol3,"CBShapeLS");
//RooChebychev CBShapeLS("CBShapeLS","CBShapeLS",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
if (centralityMin==0&&centralityMax==9)		// for central use pol3
{
	//CBShapeLS = RooChebychev("CBShapeLS","CBShapeLS",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
}
//RooGaussian CBShape_DT("CBShape_DT","CBShape_DT",MassDT,CB_DT_ParA,CB_DT_ParB);
RooAddPdf CBShapeLS_2("CBShapeLS_2","CBShapeLS_2",RooArgList(CBShapeLS),RooArgList(nCB2));
CBShapeLS_2.fitTo(DT_UnBinLS);
cout << "ncb2 is " << nCB2.getVal() << endl;
RooRealVar nCB("nCB","nCB",nCB2.getVal(),nCB2.getVal(),nCB2.getVal());

RooChebychev CBShapeUS("CBShapeUS","CBShapeUS",MassDT,RooArgSet(RooConst(CB_DT_ParA.getVal()),RooConst(CB_DT_ParB.getVal()),RooConst(CB_DT_ParC.getVal()),RooConst(CB_DT_ParD.getVal()),RooConst(CB_DT_ParE.getVal())));

//RooHistPdf CBShapeUS("CBShapeUS","CBShapeUS",MassDT,DT_histLS);

//Fit Polynom with Root because RooFit is useless
//TF1* fPol5 = new TF1("fPol5","[0]*([1]*x**4 +[2]*x**3 +[3]*x**2 +[4]*x +[5])",fitMin,fitMax);
/*TF1* fPol5 = new TF1("fPol5","ROOT::Math::Chebyshev5(x,[0],[1],[2],[3],[4],[5])",fitMin,fitMax);
fPol5->SetParameters(1,0.1,0.1,0.1,0.1,0.1);
fPol5->FixParameter(0,1);
hDiLeptonMassDT_LS->Fit("fPol5","B","",fitMin,fitMax);
CB_DT_ParA.setVal(fPol5->GetParameter(1));
CB_DT_ParB.setVal(fPol5->GetParameter(2));
CB_DT_ParC.setVal(fPol5->GetParameter(3));
CB_DT_ParD.setVal(fPol5->GetParameter(4));
CB_DT_ParE.setVal(fPol5->GetParameter(5));*/
//nCB2.setVal(fPol5->GetParameter(0));


//RooChebychev CBShapeLS("CBShapeLS","CBShapeLS",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
//CBShapeLS.fitTo(DT_UnBinLS,Range(fitMin,fitMax));
//RooRealVar nCB2("nCB2","nCB2",1,0,1e06);
//RooAddPdf CBShapeLS_2("CBShapeLS_2","CBShapeLS_2",RooArgSet(CBShapeLS),RooArgList(nCB2));

//CBShapeLS_2.fitTo(DT_histLS);
//CBShapeLS_2.fitTo(DT_UnBinLS,Range(fitMin,fitMax));
//RooRealVar nCB("nCB","nCB",nCB2.getVal(),nCB2.getVal(),nCB2.getVal());
//RooChebychev CBShapeUS("CBShapeUS","CBShapeUS",MassDT,RooArgSet(RooConst(CB_DT_ParA.getVal()),RooConst(CB_DT_ParB.getVal()),RooConst(CB_DT_ParC.getVal()),RooConst(CB_DT_ParD.getVal()),RooConst(CB_DT_ParE.getVal())));

TCanvas* cCB_DT = new TCanvas("cCB_DT","cCB_DT",700,700) ;
cCB_DT->cd();
RooPlot* plot_CB_DT = MassDT.frame(Bins(nBins),Title(" "));
DT_UnBinLS.plotOn(plot_CB_DT);
CBShapeLS.plotOn(plot_CB_DT,Name("pol5"),LineColor(myDarkRed));
plot_CB_DT->Draw();

//Background function
//RooRealVar BG_DT_ParA("BG_DT_ParA","BG_DT_ParA",1,1,1);
RooRealVar BG_DT_ParB("BG_DT_ParB","BG_DT_ParB",50,0,300);
RooRealVar BG_DT_ParC("BG_DT_ParC","BG_DT_ParC",10,0,300);
RooRealVar BG_DT_ParD("BG_DT_ParD","BG_DT_ParD",70,0,300);

RooHistPdf BackgroundShape("BackgroundShape","BackgroundShape",MassDT,DT_Background);
//RooGenericPdf BackgroundShape("BackgroundShape","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
//							RooArgSet(MassDT,BG_DT_ParB,BG_DT_ParC,BG_DT_ParD));
//BackgroundShape.fitTo(DT_Background,Range(5,15));

TCanvas* cRbg = new TCanvas("cRbg","cRbg",750,550) ;
cRbg->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_BG_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_Background.plotOn(plot_BG_DT); 
BackgroundShape.plotOn(plot_BG_DT,LineStyle(2),LineWidth(3),LineColor(1));
DT_Background.plotOn(plot_BG_DT,MarkerSize(1.1),MarkerColor(kGreen+2),LineColor(kGreen+2));
plot_BG_DT->GetYaxis()->SetTitleOffset(1.4);
plot_BG_DT->GetYaxis()->SetTitleSize(0.04);
plot_BG_DT->GetXaxis()->SetTitleSize(0.04);
plot_BG_DT->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_BG_DT->Draw();


RooRealVar BG_DT_ParBfix("BG_DT_ParBfix","BG_DT_ParBfix",BG_DT_ParB.getVal(),BG_DT_ParB.getVal(),BG_DT_ParB.getVal());
RooRealVar BG_DT_ParCfix("BG_DT_ParCfix","BG_DT_ParCfix",BG_DT_ParC.getVal(),BG_DT_ParC.getVal(),BG_DT_ParC.getVal());
RooRealVar BG_DT_ParDfix("BG_DT_ParDfix","BG_DT_ParDfix",BG_DT_ParD.getVal(),BG_DT_ParD.getVal(),BG_DT_ParD.getVal());

//RooGenericPdf BackgroundShape_DT("BackgroundShape","TMath::Power(MassDT,BG_DT_ParBfix)/TMath::Power(1+MassDT/BG_DT_ParCfix,BG_DT_ParDfix)",
//							RooArgSet(MassDT,BG_DT_ParBfix,BG_DT_ParCfix,BG_DT_ParDfix));

//RooGenericPdf BackgroundShape_DT("BackgroundShape_DT","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
//							MassDT,RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal()));
//							RooArgSet(MassDT,RooConst(BG_DT_ParA.getVal()),RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal())));


//RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S_DT,CrystalBall_3S_DT,BackgroundShape_DT,CBShape_DT),
//							   RooArgList(nSignal_1S,nSignal_2S,nSignal_3S,nBackground,nCB));
RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S3S_DT,BackgroundShape,CBShapeUS),
							   RooArgList(nSignal_1S,nSignal_2S3S,nBackground,nCB));

//RooRealVar nBackground("nBackground2","nBackground2",1,0,1e06);
DT_FitFunction.fitTo(DT_UnBin);
if (centralityMax==5) nBackground.setMax(1.5*nSignal_1S.getVal());
DT_FitFunction.fitTo(DT_UnBin);
//if (centralityMax==5) nBackground.setMax(1.7*nSignal_1S.getVal());
DT_FitFunction.fitTo(DT_UnBin);

TCanvas* cDT = new TCanvas("cDT","cDT",750,550) ;
cDT->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
//DT_UnBin.plotOn(plot_DT); 
DT_UnBin.plotOn(plot_DT,MarkerSize(1.3),MarkerColor(kRed),LineColor(kRed),LineWidth(1.5),MarkerStyle(24),RooFit::Name("usp"));

DT_FitFunction.plotOn(plot_DT,Components(BackgroundShape),LineStyle(2),LineWidth(3),LineColor(1),RooFit::Name("bg3"));
DT_FitFunction.plotOn(plot_DT,Components(CBShapeUS),LineStyle(2),LineWidth(3),LineColor(kBlue),RooFit::Name("bg2"));
DT_FitFunction.plotOn(plot_DT,LineWidth(2),LineColor(1),RooFit::Name("fiit"));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_1S_DT),LineStyle(1),LineColor(kPink-5),LineWidth(3),RooFit::Name("u1s"));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_2S3S_DT),LineStyle(1),LineColor(kSpring-5),LineWidth(3),RooFit::Name("u23s"));
//DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_2S_DT),LineStyle(1),LineColor(myDarkGreen),LineWidth(3));
//DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_3S_DT),LineStyle(1),LineColor(myDarkGreen),LineWidth(3));

DT_UnBin.plotOn(plot_DT,MarkerSize(1.25),MarkerColor(kRed),LineColor(kRed),LineWidth(1.5),MarkerStyle(24),RooFit::Name("usp"));
DT_UnBinLS.plotOn(plot_DT,MarkerSize(1.25),MarkerColor(kBlue),LineColor(kBlue),LineWidth(1.5),MarkerStyle(24),RooFit::Name("lsp"));

plot_DT->GetYaxis()->SetTitleOffset(1.2);
plot_DT->GetYaxis()->SetTitleSize(0.04);
plot_DT->GetXaxis()->SetTitleSize(0.04);
plot_DT->GetXaxis()->CenterTitle();
plot_DT->GetYaxis()->CenterTitle();
plot_DT->GetYaxis()->SetTitle(TString::Format("Counts per %.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_DT->Draw();
RooChi2Var chi2ndf("chi2ndf","chi2ndf",DT_FitFunction,DT_hist);

TLegend *myLegend_DT = new TLegend(0.064,0.608,0.547,0.87);
myLegendSetUp(myLegend_DT,0.035,1);
//myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{combBg} = %.1f #pm %.1f",nCB2.getVal(),nCB2.getError())," ");
//myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{corrBg} = %.1f #pm %.1f",nBackground.getVal(),nBackground.getError())," ");
//myLegend_DT->AddEntry((TObject*)0,TString::Format("#chi^{2}/ndf = %.1f / %i",chi2ndf.getVal(),nBins)," ");
TString textstr = Form("#bf{%i-%i % |y|<0.5}",centValues[centralityMax],centValues[centralityMin]);
cout << "string is " << textstr << endl;
myLegend_DT->AddEntry((TObject*)0,"#bf{2014 Au+Au @200 GeV}"," ");
myLegend_DT->AddEntry((TObject*)0,textstr," ");
myLegend_DT->AddEntry((TObject*)0,"#bf{p_{T}^{#Upsilon}<10 GeV/c}"," ");
myLegend_DT->AddEntry((TObject*)0,""," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{1S} = %.1f #pm %.1f",nSignal_1S.getVal(),nSignal_1S.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{2S+3S} = %.1f #pm %.1f",nSignal_2S3S.getVal(),nSignal_2S3S.getError())," ");

myLegend_DT->Draw();

TLegend *myLegend_BG = new TLegend(0.55,0.60,0.78,0.87);
myLegendSetUp(myLegend_BG,0.033,1);
myLegend_BG->AddEntry(plot_DT->findObject("usp"), "unlike sign", "pl");
myLegend_BG->AddEntry(plot_DT->findObject("lsp"), "like sign", "pl");
myLegend_BG->AddEntry((TObject*)0,""," ");
myLegend_BG->AddEntry(plot_DT->findObject("fiit"), "total fit", "l");
myLegend_BG->AddEntry(plot_DT->findObject("u1s"), "#Upsilon(1S)", "l");
myLegend_BG->AddEntry(plot_DT->findObject("u23s"), "#Upsilon(2S+3S)", "l");
myLegend_BG->AddEntry(plot_DT->findObject("bg2"), "combinatorial background", "l");
myLegend_BG->AddEntry(plot_DT->findObject("bg3"), "correlated background", "l");
/*TLegend *myLegend_BG = new TLegend(0.46,0.55,0.94,0.87);
myLegendSetUp(myLegend_BG,0.03,1);
myLegend_BG->AddEntry((TObject*)0,"BG = A*x^B/(1+x/C)^D"," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("A = %.1f #pm %.1f",BG_DT_ParA.getVal(),BG_DT_ParA.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("B = %.1f #pm %.1f",BG_DT_ParB.getVal(),BG_DT_ParB.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("C = %.1f #pm %.1f",BG_DT_ParC.getVal(),BG_DT_ParC.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("D = %.1f #pm %.1f",BG_DT_ParD.getVal(),BG_DT_ParD.getError())," ");
*/myLegend_BG->Draw();

cDT->SaveAs(TString::Format("new_DataFit_Cent%d_%d.png",centralityMin,centralityMax));
cDT->SaveAs(TString::Format("new_DataFit_Cent%d_%d.eps",centralityMin,centralityMax));

yields[0] = nSignal_1S.getVal();
yields[1] = nSignal_1S.getError();
yields[2] = nSignal_2S3S.getVal();
yields[3] = nSignal_2S3S.getError();

}



///////////////////////////////////////////////////////
void FitUpsilonMassUSLS(Int_t centralityMin = 7, Int_t centralityMax = 9) {

TFile *DT_file = TFile::Open("masstree.root", "read");

TFile *fTemplates1S = TFile::Open("hEmb1S.root", "read");
TFile *fTemplates2S = TFile::Open("hEmb2S.root", "read");
TFile *fTemplates3S = TFile::Open("hEmb3S.root", "read");

TFile *fRbg = TFile::Open("hRbg.root","read");

TH2F *hTemplate1S_2D = (TH2F*) fTemplates1S->Get("hUnIM");
TH1D *hTemplate1S = hTemplate1S_2D->ProjectionY("hTemplate1S");
hTemplate1S->Rebin();

TH2F *hTemplate2S_2D = (TH2F*) fTemplates2S->Get("hUnIM");
TH1D *hTemplate2S = hTemplate2S_2D->ProjectionY("hTemplate2S");
hTemplate2S->Rebin();

TH2F *hTemplate3S_2D = (TH2F*) fTemplates3S->Get("hUnIM");
TH1D *hTemplate3S = hTemplate3S_2D->ProjectionY("hTemplate3S");
hTemplate3S->Rebin();

TH1D* hRbg = (TH1D*)fRbg->Get("hRbg");
hRbg->Smooth();


//________________________________________________DATA_________________________________________________________
Int_t fCentrality;
Float_t fMass;
Int_t fSign;
Float_t fM;

Int_t centValues[10] = {80,70,60,50,40,30,20,10,5,0};

TTree *massTree = (TTree*)DT_file->Get("massTree");
massTree->SetBranchAddress("fMass",&fMass);
massTree->SetBranchAddress("fCentrality",&fCentrality);
massTree->SetBranchAddress("fSign",&fSign);

UInt_t nDT = massTree->GetEntries();

TH1D *hDiLeptonMassDT_Hist = new TH1D("hDiLeptonMassDT_Hist","Invariant mass of dileptons candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_Hist->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}) (GeV/c)");
hDiLeptonMassDT_Hist->Sumw2();

TH1D *hDiLeptonMassDT_LS = new TH1D("hDiLeptonMassDT_LS","Invariant mass like-sign candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_LS->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{+} l^{-}l^{-}) (GeV/c)");
hDiLeptonMassDT_LS->Sumw2();

TTree *fTreeDileptonMass = new TTree("fTreeDileptonMass", "fTreeDileptonMass");
fTreeDileptonMass->Branch("MassDT", &fM, "MassDT/F");
TTree *fTreeDileptonMassLS = new TTree("fTreeDileptonMassLS", "fTreeDileptonMassLS");
fTreeDileptonMassLS->Branch("MassDT", &fM, "MassDT/F");

Bool_t ZNA0n, ZNC0n;

    //event loop
int howmany = 0;
for(Int_t i=0; i<nDT; i++){
	massTree->GetEntry(i);
	if (fCentrality<centralityMin || fCentrality>=centralityMax) continue;
	fM = fMass;
	if (fSign < 0) {
		fTreeDileptonMass->Fill();
		hDiLeptonMassDT_Hist->Fill(fM); }
	else {
    fTreeDileptonMassLS->Fill();
    hDiLeptonMassDT_LS->Fill(fM);
  }
	howmany++;
}
cout << "howmany is " << howmany << endl;

// Declare observable - Mass for data
RooRealVar w("w","w",-2,2);
w=1;	
RooRealVar MassDT("MassDT","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",fitMin,fitMax);
RooDataSet DT_UnBin("DT_UnBin","DT_UnBin",RooArgSet(MassDT,w),Import(*fTreeDileptonMass),WeightVar(w));
//RooRealVar MassDTLS("MassDTLS","M_{ l^{+-}l^{+-}} (GeV/#it{c}^{2})",fitMin,fitMax);
w=-1;
RooDataSet DT_UnBinLS("DT_UnBinLS","DT_UnBinLS",RooArgSet(MassDT,w),Import(*fTreeDileptonMassLS),WeightVar("w"));

hDiLeptonMassDT_Hist->Add(hDiLeptonMassDT_LS,-1);
RooDataHist DT_hist("DT_hist","DT_hist",MassDT,Import(*hDiLeptonMassDT_Hist));

//US-LS
/*RooRealVar w("w","w",-1,1) ;
w.setVal(1); 
DT_UnBin.addColumn(w);
DT_UnBin.setWeightVar(w);
w.setVal(-1); 
DT_UnBinLS.addColumn(w);
DT_UnBinLS.setWeightVar(w);

// Append datasets and interpret w as event weight*/
DT_UnBin.append(DT_UnBinLS);

// Declare observable - Masses for MC
RooRealVar MassMC_Signal_1S("MassMC_Signal_1S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);
RooRealVar MassMC_Signal_2S("MassMC_Signal_2S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);
RooRealVar MassMC_Signal_3S("MassMC_Signal_3S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);

RooDataHist MC_hist_Signal_1S("MC_hist_Signal_1S","MC_hist_Signal_1S",MassMC_Signal_1S,Import(*hTemplate1S));
RooDataHist MC_hist_Signal_2S("MC_hist_Signal_2S","MC_hist_Signal_2S",MassMC_Signal_2S,Import(*hTemplate2S));
RooDataHist MC_hist_Signal_3S("MC_hist_Signal_3S","MC_hist_Signal_3S",MassMC_Signal_3S,Import(*hTemplate3S));

// Histogram for rbg
RooDataHist DT_Background("DT_Background","DT_Background",MassDT,Import(*hRbg));

// Crystal Ball P.D.F. variables for 1S,2S,3S
RooRealVar Mean_1S("Mean_1S","m_{1S}",9.388,8,12);
Mean_1S.setError(0.011);
//RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.1938,0.05,0.4);
RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.2289,0.05,0.4);
Sigma_1S.setError(0.009);
RooRealVar Alpha_1S("Alpha_1S","Alpha_1S",1.331,0.0,2);
Alpha_1S.setError(0.1015);
RooRealVar N_1S("N_1S","N_1S",0.9619,0,4);
N_1S.setError(0.224);

RooRealVar Mean_2S("Mean_2S","m_{2S}",9.934,8,12);
Mean_2S.setError(0.010);
//RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2041,0.05,0.4);
RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2540,0.05,0.4);
Sigma_2S.setError(0.0078);
RooRealVar Alpha_2S("Alpha_2S","Alpha_2S",1.285,0.0,2);
Alpha_2S.setError(0.0848);
RooRealVar N_2S("N_2S","N_2S",1.279,0,4);
N_2S.setError(0.160);

RooRealVar Mean_3S("Mean_3S","m_{3S}",10.27,8,12);
Mean_3S.setError(0.010);
//RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2303,0.05,0.4);
RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2655,0.05,0.4);
Sigma_3S.setError(0.0088);
RooRealVar Alpha_3S("Alpha_3S","Alpha_3S",1.073,0.0,2);
Alpha_3S.setError(0.090);
RooRealVar N_3S("N_3S","N_3S",1.035,0,4);
N_3S.setError(0.163);


//Fit Crystal Ball to 1S MC signal
RooCBShape CrystalBall_1S_MC("CrystalBall_1S_MC","CrystalBall_1S_MC",MassMC_Signal_1S,Mean_1S,Sigma_1S,Alpha_1S,N_1S);
//CrystalBall_1S_MC.fitTo(MC_hist_Signal_1S,Range(5,15));

//Fit Crystal Ball to 2S MC signal
RooCBShape CrystalBall_2S_MC("CrystalBall_2S_MC","CrystalBall_2S_MC",MassMC_Signal_2S,Mean_2S,Sigma_2S,Alpha_2S,N_2S);
//CrystalBall_2S_MC.fitTo(MC_hist_Signal_2S,Range(5.0,15));

//Fit Crystal Ball to 3S MC signal
RooCBShape CrystalBall_3S_MC("CrystalBall_3S_MC","CrystalBall_3S_MC",MassMC_Signal_3S,Mean_3S,Sigma_3S,Alpha_3S,N_3S);
//CrystalBall_3S_MC.fitTo(MC_hist_Signal_3S,Range(5,15));

// Fit a sum of the functions to the data
RooRealVar nSignal_1S("nSignal_1S","N_{1S}",1,0,1e06);
RooRealVar nSignal_2S("nSignal_2S","N_{2S}",1,0,1e06);
RooRealVar nSignal_3S("nSignal_3S","N_{3S}",1,0,1e06);
RooRealVar nSignal_2S3S("nSignal_2S3S","N_{2S3S}",1,0,1e06);
RooRealVar nBackground("nBackground","nBackground",1,0,1e06);

RooRealVar Sigma_1S_DT("Sigma_1S_DT","Sigma_1S_DT",Sigma_1S.getVal(),Sigma_1S.getVal()-0*Sigma_1S.getError(),Sigma_1S.getVal()+0*Sigma_1S.getError());
RooRealVar Sigma_2S_DT("Sigma_2S_DT","Sigma_2S_DT",Sigma_2S.getVal(),Sigma_2S.getVal()-0*Sigma_2S.getError(),Sigma_2S.getVal()+0*Sigma_2S.getError());
RooRealVar Sigma_3S_DT("Sigma_3S_DT","Sigma_3S_DT",Sigma_3S.getVal(),Sigma_3S.getVal()-0*Sigma_3S.getError(),Sigma_3S.getVal()+0*Sigma_3S.getError());

RooCBShape CrystalBall_1S_DT("CrystalBall_1S_DT","CrystalBall_1S_DT",MassDT,RooConst(Mean_1S.getVal()),Sigma_1S_DT,RooConst(Alpha_1S.getVal()),RooConst(N_1S.getVal()));
RooCBShape CrystalBall_2S_DT("CrystalBall_2S_DT","CrystalBall_2S_DT",MassDT,RooConst(Mean_2S.getVal()),Sigma_2S_DT,RooConst(Alpha_2S.getVal()),RooConst(N_2S.getVal()));
RooCBShape CrystalBall_3S_DT("CrystalBall_3S_DT","CrystalBall_3S_DT",MassDT,RooConst(Mean_3S.getVal()),Sigma_3S_DT,RooConst(Alpha_3S.getVal()),RooConst(N_3S.getVal()));
RooRealVar BLAOne("BLAOne","BLAOne",689,689,689);
RooRealVar BLATwo("BLATwo","BLATwo",311,311,311);
//RooRealVar BLAOne("BLAOne","BLAOne",0.689,0.689,0.689);
//RooRealVar BLATwo("BLATwo","BLATwo",0.311,0.311,0.311);
RooAddPdf CrystalBall_2S3S_DT("CrystalBall_2S3S_DT","CrystalBall_2S3S_DT",RooArgList(CrystalBall_2S_DT,CrystalBall_3S_DT),RooArgList(BLAOne,BLATwo));


//Background function
//RooRealVar BG_DT_ParA("BG_DT_ParA","BG_DT_ParA",1,1,1);
RooRealVar BG_DT_ParB("BG_DT_ParB","BG_DT_ParB",50,0,300);
RooRealVar BG_DT_ParC("BG_DT_ParC","BG_DT_ParC",10,0,300);
RooRealVar BG_DT_ParD("BG_DT_ParD","BG_DT_ParD",70,0,300);

RooHistPdf BackgroundShape("BackgroundShape","BackgroundShape",MassDT,DT_Background);

//RooGenericPdf BackgroundShape("BackgroundShape","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
	//						RooArgSet(MassDT,BG_DT_ParB,BG_DT_ParC,BG_DT_ParD));
//BackgroundShape.fitTo(DT_Background,Range(5,15));

TCanvas* cRbg = new TCanvas("cRbg","cRbg",750,550) ;
cRbg->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_BG_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_Background.plotOn(plot_BG_DT); 
BackgroundShape.plotOn(plot_BG_DT,LineStyle(2),LineWidth(3),LineColor(1));
DT_Background.plotOn(plot_BG_DT,MarkerSize(1.1),MarkerColor(kGreen+2),LineColor(kGreen+2));
plot_BG_DT->GetYaxis()->SetTitleOffset(1.4);
plot_BG_DT->GetYaxis()->SetTitleSize(0.04);
plot_BG_DT->GetXaxis()->SetTitleSize(0.04);
plot_BG_DT->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_BG_DT->Draw();


RooRealVar BG_DT_ParBfix("BG_DT_ParBfix","BG_DT_ParBfix",BG_DT_ParB.getVal(),BG_DT_ParB.getVal(),BG_DT_ParB.getVal());
RooRealVar BG_DT_ParCfix("BG_DT_ParCfix","BG_DT_ParCfix",BG_DT_ParC.getVal(),BG_DT_ParC.getVal(),BG_DT_ParC.getVal());
RooRealVar BG_DT_ParDfix("BG_DT_ParDfix","BG_DT_ParDfix",BG_DT_ParD.getVal(),BG_DT_ParD.getVal(),BG_DT_ParD.getVal());


//RooGenericPdf BackgroundShape_DT("BackgroundShape","TMath::Power(MassDT,BG_DT_ParBfix)/TMath::Power(1+MassDT/BG_DT_ParCfix,BG_DT_ParDfix)",
//							RooArgSet(MassDT,BG_DT_ParBfix,BG_DT_ParCfix,BG_DT_ParDfix));

//RooGenericPdf BackgroundShape_DT("BackgroundShape_DT","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
//							MassDT,RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal()));
//							RooArgSet(MassDT,RooConst(BG_DT_ParA.getVal()),RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal())));


//RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S_DT,CrystalBall_3S_DT,BackgroundShape_DT,CBShape_DT),
//							   RooArgList(nSignal_1S,nSignal_2S,nSignal_3S,nBackground,nCB));
RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S3S_DT,BackgroundShape),
							   RooArgList(nSignal_1S,nSignal_2S3S,nBackground));
DT_FitFunction.fitTo(DT_UnBin);

if (centralityMax==5) nBackground.setMax(1.5*nSignal_1S.getVal());
DT_FitFunction.fitTo(DT_UnBin);

TCanvas* cDT = new TCanvas("cDT2","cDT2",750,550) ;
cDT->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_UnBin.plotOn(plot_DT); 
DT_FitFunction.plotOn(plot_DT,Components(BackgroundShape),LineStyle(2),LineWidth(3),LineColor(1),RooFit::Name("bg3"));
DT_FitFunction.plotOn(plot_DT,LineWidth(2),LineColor(1),RooFit::Name("fiit"));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_1S_DT),LineStyle(1),LineColor(myDarkRed),LineWidth(3));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_2S3S_DT),LineStyle(1),LineColor(myDarkGreen),LineWidth(3));
DT_UnBin.plotOn(plot_DT,MarkerSize(1.1),MarkerColor(kRed),LineColor(kRed));
plot_DT->GetYaxis()->SetTitleOffset(1.4);
plot_DT->GetYaxis()->SetTitleSize(0.04);
plot_DT->GetXaxis()->SetTitleSize(0.04);
plot_DT->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_DT->Draw();
RooChi2Var chi2ndf("chi2ndf","chi2ndf",DT_FitFunction,DT_hist);

TLegend *myLegend_DT = new TLegend(0.064,0.588,0.547,0.87);
myLegendSetUp(myLegend_DT,0.03,1);
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{1S} = %.1f #pm %.1f",nSignal_1S.getVal(),nSignal_1S.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{2S+3S} = %.1f #pm %.1f",nSignal_2S3S.getVal(),nSignal_2S3S.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{corrbg} = %.1f #pm %.1f",nBackground.getVal(),nBackground.getError())," ");
//myLegend_DT->AddEntry((TObject*)0,TString::Format("#chi^{2}/ndf = %.1f / %i",chi2ndf.getVal(),nBins)," ");
myLegend_DT->Draw();

TLegend *myLegend_BG = new TLegend(0.57,0.60,0.94,0.87);
myLegendSetUp(myLegend_BG,0.03,1);
myLegend_BG->AddEntry(plot_DT->findObject("fiit"), "sig. + corr. bg", "l");
myLegend_BG->AddEntry(plot_DT->findObject("bg3"), "corr. bg", "l");
TString textstr = Form("#bf{%i-%i %}",centValues[centralityMax],centValues[centralityMin]);
cout << "string is " << textstr << endl;
myLegend_BG->AddEntry((TObject*)0,textstr," ");
/*TLegend *myLegend_BG = new TLegend(0.46,0.55,0.94,0.87);
myLegendSetUp(myLegend_BG,0.03,1);
myLegend_BG->AddEntry((TObject*)0,"BG = A*x^B/(1+x/C)^D"," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("A = %.1f #pm %.1f",BG_DT_ParA.getVal(),BG_DT_ParA.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("B = %.1f #pm %.1f",BG_DT_ParB.getVal(),BG_DT_ParB.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("C = %.1f #pm %.1f",BG_DT_ParC.getVal(),BG_DT_ParC.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("D = %.1f #pm %.1f",BG_DT_ParD.getVal(),BG_DT_ParD.getError())," ");
*/myLegend_BG->Draw();

cDT->SaveAs(TString::Format("USLSDataFit_Cent%d_%d.png",centralityMin,centralityMax));
cout << "n2s vs n3s " << BLAOne.getVal() << " " << BLATwo.getVal() << endl;

yields2[0] = nSignal_1S.getVal();
yields2[1] = nSignal_1S.getError();
yields2[2] = nSignal_2S3S.getVal();
yields2[3] = nSignal_2S3S.getError();

}


/////////////////////////////////////////////////////////////////////////////////////////
void FitUpsilonMassME(Int_t centralityMin = 7, Int_t centralityMax = 9) {

TFile *DT_file = TFile::Open("masstree.root", "read");

TFile *fTemplates1S = TFile::Open("hEmb1S.root", "read");
TFile *fTemplates2S = TFile::Open("hEmb2S.root", "read");
TFile *fTemplates3S = TFile::Open("hEmb3S.root", "read");

TFile *fRbg = TFile::Open("hRbg.root","read");

TH2F* hIMmixC	= (TH2F*)DT_file->Get("hIMmixC");
TH1F* hIMmix 	= (TH1F*)hIMmixC->ProjectionX("mixed events; m; counts", centralityMin, centralityMax)->Clone("hIMmix");
hIMmix->Sumw2();
hIMmix->Smooth();
//hIMmix->Rebin();

TH2F *hTemplate1S_2D = (TH2F*) fTemplates1S->Get("hUnIM");
TH1D *hTemplate1S = hTemplate1S_2D->ProjectionY("hTemplate1S");
hTemplate1S->Rebin();

TH2F *hTemplate2S_2D = (TH2F*) fTemplates2S->Get("hUnIM");
TH1D *hTemplate2S = hTemplate2S_2D->ProjectionY("hTemplate2S");
hTemplate2S->Rebin();

TH2F *hTemplate3S_2D = (TH2F*) fTemplates3S->Get("hUnIM");
TH1D *hTemplate3S = hTemplate3S_2D->ProjectionY("hTemplate3S");
hTemplate3S->Rebin();

TH1D* hRbg = (TH1D*)fRbg->Get("hRbg");
hRbg->Smooth();
//hRbg->Rebin();


//________________________________________________DATA_________________________________________________________
Int_t fCentrality;
Float_t fMass;
Int_t fSign;
Float_t fM;

Int_t centValues[10] = {80,70,60,50,40,30,20,10,5,0};

TTree *massTree = (TTree*)DT_file->Get("massTree");
massTree->SetBranchAddress("fMass",&fMass);
massTree->SetBranchAddress("fCentrality",&fCentrality);
massTree->SetBranchAddress("fSign",&fSign);

UInt_t nDT = massTree->GetEntries();

TH1D *hDiLeptonMassDT_Hist = new TH1D("hDiLeptonMassDT_Hist","Invariant mass of dileptons candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_Hist->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{-}) (GeV/c)");
hDiLeptonMassDT_Hist->Sumw2();

TH1D *hDiLeptonMassDT_LS = new TH1D("hDiLeptonMassDT_LS","Invariant mass like-sign candidates DT",nBins,fitMin,fitMax);
hDiLeptonMassDT_LS->GetXaxis()->SetTitle("Invariant mass(l^{+}l^{+} l^{-}l^{-}) (GeV/c)");
hDiLeptonMassDT_LS->Sumw2();

TTree *fTreeDileptonMass = new TTree("fTreeDileptonMass", "fTreeDileptonMass");
fTreeDileptonMass->Branch("MassDT", &fM, "MassDT/F");
TTree *fTreeDileptonMassLS = new TTree("fTreeDileptonMassLS", "fTreeDileptonMassLS");
fTreeDileptonMassLS->Branch("MassDT", &fM, "MassDT/F");

Bool_t ZNA0n, ZNC0n;

    //event loop
int howmany = 0;
for(Int_t i=0; i<nDT; i++){
	massTree->GetEntry(i);
	if (fCentrality<centralityMin || fCentrality>=centralityMax) continue;
	fM = fMass;
	if (fSign < 0) {
		fTreeDileptonMass->Fill();
		hDiLeptonMassDT_Hist->Fill(fM); }
	else {
    fTreeDileptonMassLS->Fill();
    hDiLeptonMassDT_LS->Fill(fM);
  }
	howmany++;
}
cout << "howmany is " << howmany << endl;

// Declare observable - Mass for data
RooRealVar w("w","w",-2,2);
w=1;	
RooRealVar MassDT("MassDT","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",fitMin,fitMax);
RooDataSet DT_UnBin("DT_UnBin","DT_UnBin",RooArgSet(MassDT,w),Import(*fTreeDileptonMass),WeightVar(w));
//RooRealVar MassDTLS("MassDTLS","M_{ l^{+-}l^{+-}} (GeV/#it{c}^{2})",fitMin,fitMax);
//w=-1;
RooDataSet DT_UnBinLS("DT_UnBinLS","DT_UnBinLS",RooArgSet(MassDT,w),Import(*fTreeDileptonMassLS),WeightVar("w"));
RooDataHist DT_histLS("DT_histLS","DT_histLS",MassDT,Import(*hDiLeptonMassDT_LS));
RooDataHist DT_hist("DT_hist","DT_hist",MassDT,Import(*hDiLeptonMassDT_Hist));

//US-LS
/*RooRealVar w("w","w",-1,1) ;
w.setVal(1); 
DT_UnBin.addColumn(w);
DT_UnBin.setWeightVar(w);
w.setVal(-1); 
DT_UnBinLS.addColumn(w);
DT_UnBinLS.setWeightVar(w);

// Append datasets and interpret w as event weight*/
//DT_UnBin.append(DT_UnBinLS);

// Declare observable - Masses for MC
RooRealVar MassMC_Signal_1S("MassMC_Signal_1S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);
RooRealVar MassMC_Signal_2S("MassMC_Signal_2S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);
RooRealVar MassMC_Signal_3S("MassMC_Signal_3S","M_{ l^{+}l^{-}} (GeV/#it{c}^{2})",5.0,15.0);

RooDataHist MC_hist_Signal_1S("MC_hist_Signal_1S","MC_hist_Signal_1S",MassMC_Signal_1S,Import(*hTemplate1S));
RooDataHist MC_hist_Signal_2S("MC_hist_Signal_2S","MC_hist_Signal_2S",MassMC_Signal_2S,Import(*hTemplate2S));
RooDataHist MC_hist_Signal_3S("MC_hist_Signal_3S","MC_hist_Signal_3S",MassMC_Signal_3S,Import(*hTemplate3S));

// Histogram for mixed events
RooDataHist ME_Background("ME_Background","ME_Background",MassDT,Import(*hIMmix));

// Histogram for rbg
RooDataHist DT_Background("DT_Background","DT_Background",MassDT,Import(*hRbg));

// Crystal Ball P.D.F. variables for 1S,2S,3S
RooRealVar Mean_1S("Mean_1S","m_{1S}",9.388,8,12);
Mean_1S.setError(0.011);
//RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.1938,0.05,0.4);
RooRealVar Sigma_1S("Sigma_1S","#sigma_{1S}",0.2289,0.05,0.4);
Sigma_1S.setError(0.009);
RooRealVar Alpha_1S("Alpha_1S","Alpha_1S",1.331,0.0,2);
Alpha_1S.setError(0.1015);
RooRealVar N_1S("N_1S","N_1S",0.9619,0,4);
N_1S.setError(0.224);

RooRealVar Mean_2S("Mean_2S","m_{2S}",9.934,8,12);
Mean_2S.setError(0.010);
//RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2041,0.05,0.4);
RooRealVar Sigma_2S("Sigma_2S","#sigma_{2S}",0.2540,0.05,0.4);
Sigma_2S.setError(0.0078);
RooRealVar Alpha_2S("Alpha_2S","Alpha_2S",1.285,0.0,2);
Alpha_2S.setError(0.0848);
RooRealVar N_2S("N_2S","N_2S",1.279,0,4);
N_2S.setError(0.160);

RooRealVar Mean_3S("Mean_3S","m_{3S}",10.27,8,12);
Mean_3S.setError(0.010);
//RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2303,0.05,0.4);
RooRealVar Sigma_3S("Sigma_3S","#sigma_{3S}",0.2655,0.05,0.4);
Sigma_3S.setError(0.0088);
RooRealVar Alpha_3S("Alpha_3S","Alpha_3S",1.073,0.0,2);
Alpha_3S.setError(0.090);
RooRealVar N_3S("N_3S","N_3S",1.035,0,4);
N_3S.setError(0.163);


//Fit Crystal Ball to 1S MC signal
RooCBShape CrystalBall_1S_MC("CrystalBall_1S_MC","CrystalBall_1S_MC",MassMC_Signal_1S,Mean_1S,Sigma_1S,Alpha_1S,N_1S);
//CrystalBall_1S_MC.fitTo(MC_hist_Signal_1S,Range(5,15));
TCanvas* cMC_1S = new TCanvas("cMC_1S","cMC_1S",700,700) ;
cMC_1S->cd();
RooPlot* plot_MC_1S = MassMC_Signal_1S.frame(Title(" "));
MC_hist_Signal_1S.plotOn(plot_MC_1S);
CrystalBall_1S_MC.plotOn(plot_MC_1S,Name("CrystalBall"),LineColor(myDarkBlue));
plot_MC_1S->Draw();

TLegend *myLegend_1S = new TLegend(0.41,0.55,0.893,0.85);
myLegendSetUp(myLegend_1S,0.03,1);
myLegend_1S->AddEntry((TObject*)0,TString::Format("m_{1S} = %5.3f #pm %5.3f GeV/#it{c}^{2}",Mean_1S.getVal(),Mean_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("#sigma_{1S} = %5.3f #pm %5.3f GeV/#it{c}^{2}",Sigma_1S.getVal(),Sigma_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("#alpha_{1S} = %5.3f #pm %5.3f",Alpha_1S.getVal(),Alpha_1S.getError())," ");
myLegend_1S->AddEntry((TObject*)0,TString::Format("n_{1S} = %5.3f #pm %5.3f",N_1S.getVal(),N_1S.getError())," ");
myLegend_1S->Draw();
cMC_1S->SaveAs("cMC_1S.png");

//Fit Crystal Ball to 2S MC signal
RooCBShape CrystalBall_2S_MC("CrystalBall_2S_MC","CrystalBall_2S_MC",MassMC_Signal_2S,Mean_2S,Sigma_2S,Alpha_2S,N_2S);
//CrystalBall_2S_MC.fitTo(MC_hist_Signal_2S,Range(5.0,15));

//Fit Crystal Ball to 3S MC signal
RooCBShape CrystalBall_3S_MC("CrystalBall_3S_MC","CrystalBall_3S_MC",MassMC_Signal_3S,Mean_3S,Sigma_3S,Alpha_3S,N_3S);
//CrystalBall_3S_MC.fitTo(MC_hist_Signal_3S,Range(5,15));

// Fit a sum of the functions to the data
RooRealVar nSignal_1S("nSignal_1S","N_{1S}",1,0,1e06);
RooRealVar nSignal_2S("nSignal_2S","N_{2S}",1,0,1e06);
RooRealVar nSignal_3S("nSignal_3S","N_{3S}",1,0,1e06);
RooRealVar nSignal_2S3S("nSignal_2S3S","N_{2S3S}",1,0,1e06);
RooRealVar nBackground("nBackground","nBackground",1,0,1e06);
RooRealVar nCB2("nCB2","nCB2",1,0,1e06);

RooRealVar Sigma_1S_DT("Sigma_1S_DT","Sigma_1S_DT",Sigma_1S.getVal(),Sigma_1S.getVal()-0*Sigma_1S.getError(),Sigma_1S.getVal()+0*Sigma_1S.getError());
RooRealVar Sigma_2S_DT("Sigma_2S_DT","Sigma_2S_DT",Sigma_2S.getVal(),Sigma_2S.getVal()-0*Sigma_2S.getError(),Sigma_2S.getVal()+0*Sigma_2S.getError());
RooRealVar Sigma_3S_DT("Sigma_3S_DT","Sigma_3S_DT",Sigma_3S.getVal(),Sigma_3S.getVal()-0*Sigma_3S.getError(),Sigma_3S.getVal()+0*Sigma_3S.getError());

RooCBShape CrystalBall_1S_DT("CrystalBall_1S_DT","CrystalBall_1S_DT",MassDT,RooConst(Mean_1S.getVal()),Sigma_1S_DT,RooConst(Alpha_1S.getVal()),RooConst(N_1S.getVal()));
RooCBShape CrystalBall_2S_DT("CrystalBall_2S_DT","CrystalBall_2S_DT",MassDT,RooConst(Mean_2S.getVal()),Sigma_2S_DT,RooConst(Alpha_2S.getVal()),RooConst(N_2S.getVal()));
RooCBShape CrystalBall_3S_DT("CrystalBall_3S_DT","CrystalBall_3S_DT",MassDT,RooConst(Mean_3S.getVal()),Sigma_3S_DT,RooConst(Alpha_3S.getVal()),RooConst(N_3S.getVal()));
RooRealVar BLAOne("BLAOne","BLAOne",689,689,689);
RooRealVar BLATwo("BLATwo","BLATwo",311,311,311);
//RooRealVar BLAOne("BLAOne","BLAOne",0.689,0.689,0.689);
//RooRealVar BLATwo("BLATwo","BLATwo",0.311,0.311,0.311);
RooAddPdf CrystalBall_2S3S_DT("CrystalBall_2S3S_DT","CrystalBall_2S3S_DT",RooArgList(CrystalBall_2S_DT,CrystalBall_3S_DT),RooArgList(BLAOne,BLATwo));

//Combinatorial background
/*RooRealVar CB_DT_ParA("CB_DT_ParA","CB_DT_ParA",1,-100,100);
RooRealVar CB_DT_ParB("CB_DT_ParB","CB_DT_ParB",0.1,-100,100);
RooRealVar CB_DT_ParC("CB_DT_ParC","CB_DT_ParC",-0.1,-100,100);
RooRealVar CB_DT_ParD("CB_DT_ParD","CB_DT_ParD",-0.1,-100,100);
RooRealVar CB_DT_ParE("CB_DT_ParE","CB_DT_ParE",-0.1,-100,100);*/
//RooChebychev CBShapeLS("CBShapeLS","CBShapeLS",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
//RooGaussian CBShape_DT("CBShape_DT","CBShape_DT",MassDT,CB_DT_ParA,CB_DT_ParB);
RooHistPdf CBShapeLS("CBShapeLS","CBShapeLS",MassDT,ME_Background);
RooAddPdf CBShapeLS_2("CBShapeLS_2","CBShapeLS_2",RooArgList(CBShapeLS),RooArgList(nCB2));
CBShapeLS_2.fitTo(DT_UnBinLS,Range(fitMin,fitMax));
cout << "ncb2 is " << nCB2.getVal() << endl;
RooRealVar nCB("nCB","nCB",nCB2.getVal(),nCB2.getVal(),nCB2.getVal());

//RooChebychev CBShapeUS("CBShapeUS","CBShapeUS",MassDT,RooArgSet(RooConst(CB_DT_ParA.getVal()),RooConst(CB_DT_ParB.getVal()),RooConst(CB_DT_ParC.getVal()),RooConst(CB_DT_ParD.getVal()),RooConst(CB_DT_ParE.getVal())));

//RooHistPdf CBShapeUS("CBShapeUS","CBShapeUS",MassDT,DT_histLS);

//Fit Polynom with Root because RooFit is useless
//TF1* fPol5 = new TF1("fPol5","[0]*([1]*x**4 +[2]*x**3 +[3]*x**2 +[4]*x +[5])",fitMin,fitMax);
/*TF1* fPol5 = new TF1("fPol5","ROOT::Math::Chebyshev5(x,[0],[1],[2],[3],[4],[5])",fitMin,fitMax);
fPol5->SetParameters(1,0.1,0.1,0.1,0.1,0.1);
fPol5->FixParameter(0,1);
hDiLeptonMassDT_LS->Fit("fPol5","B","",fitMin,fitMax);
CB_DT_ParA.setVal(fPol5->GetParameter(1));
CB_DT_ParB.setVal(fPol5->GetParameter(2));
CB_DT_ParC.setVal(fPol5->GetParameter(3));
CB_DT_ParD.setVal(fPol5->GetParameter(4));
CB_DT_ParE.setVal(fPol5->GetParameter(5));*/
//nCB2.setVal(fPol5->GetParameter(0));


//RooChebychev CBShapeLS("CBShapeLS","CBShapeLS",MassDT,RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC,CB_DT_ParD,CB_DT_ParE));
//CBShapeLS.fitTo(DT_UnBinLS,Range(fitMin,fitMax));
//RooRealVar nCB2("nCB2","nCB2",1,0,1e06);
//RooAddPdf CBShapeLS_2("CBShapeLS_2","CBShapeLS_2",RooArgSet(CBShapeLS),RooArgList(nCB2));

//CBShapeLS_2.fitTo(DT_histLS);
//CBShapeLS_2.fitTo(DT_UnBinLS,Range(fitMin,fitMax));
//RooRealVar nCB("nCB","nCB",nCB2.getVal(),nCB2.getVal(),nCB2.getVal());
//RooChebychev CBShapeUS("CBShapeUS","CBShapeUS",MassDT,RooArgSet(RooConst(CB_DT_ParA.getVal()),RooConst(CB_DT_ParB.getVal()),RooConst(CB_DT_ParC.getVal()),RooConst(CB_DT_ParD.getVal()),RooConst(CB_DT_ParE.getVal())));

TCanvas* cCB_DT = new TCanvas("cCB_DT3","cCB_DT3",700,700) ;
cCB_DT->cd();
RooPlot* plot_CB_DT = MassDT.frame(Bins(nBins),Title(" "));
DT_UnBinLS.plotOn(plot_CB_DT);
CBShapeLS.plotOn(plot_CB_DT,Name("pol5"),LineColor(myDarkRed));
plot_CB_DT->Draw();

#if 1

//Background function
//RooRealVar BG_DT_ParA("BG_DT_ParA","BG_DT_ParA",1,1,1);
RooRealVar BG_DT_ParB("BG_DT_ParB","BG_DT_ParB",50,0,300);
RooRealVar BG_DT_ParC("BG_DT_ParC","BG_DT_ParC",10,0,300);
RooRealVar BG_DT_ParD("BG_DT_ParD","BG_DT_ParD",70,0,300);

RooHistPdf BackgroundShape("BackgroundShape","BackgroundShape",MassDT,DT_Background);
//RooGenericPdf BackgroundShape("BackgroundShape","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
//							RooArgSet(MassDT,BG_DT_ParB,BG_DT_ParC,BG_DT_ParD));
//BackgroundShape.fitTo(DT_Background,Range(5,15));

TCanvas* cRbg = new TCanvas("cRbg","cRbg",750,550) ;
cRbg->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_BG_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_Background.plotOn(plot_BG_DT); 
BackgroundShape.plotOn(plot_BG_DT,LineStyle(2),LineWidth(3),LineColor(1));
DT_Background.plotOn(plot_BG_DT,MarkerSize(1.1),MarkerColor(kGreen+2),LineColor(kGreen+2));
plot_BG_DT->GetYaxis()->SetTitleOffset(1.4);
plot_BG_DT->GetYaxis()->SetTitleSize(0.04);
plot_BG_DT->GetXaxis()->SetTitleSize(0.04);
plot_BG_DT->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_BG_DT->Draw();


RooRealVar BG_DT_ParBfix("BG_DT_ParBfix","BG_DT_ParBfix",BG_DT_ParB.getVal(),BG_DT_ParB.getVal(),BG_DT_ParB.getVal());
RooRealVar BG_DT_ParCfix("BG_DT_ParCfix","BG_DT_ParCfix",BG_DT_ParC.getVal(),BG_DT_ParC.getVal(),BG_DT_ParC.getVal());
RooRealVar BG_DT_ParDfix("BG_DT_ParDfix","BG_DT_ParDfix",BG_DT_ParD.getVal(),BG_DT_ParD.getVal(),BG_DT_ParD.getVal());

//RooGenericPdf BackgroundShape_DT("BackgroundShape","TMath::Power(MassDT,BG_DT_ParBfix)/TMath::Power(1+MassDT/BG_DT_ParCfix,BG_DT_ParDfix)",
//							RooArgSet(MassDT,BG_DT_ParBfix,BG_DT_ParCfix,BG_DT_ParDfix));

//RooGenericPdf BackgroundShape_DT("BackgroundShape_DT","TMath::Power(MassDT,BG_DT_ParB)/TMath::Power(1+MassDT/BG_DT_ParC,BG_DT_ParD)",
//							MassDT,RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal()));
//							RooArgSet(MassDT,RooConst(BG_DT_ParA.getVal()),RooConst(BG_DT_ParB.getVal()),RooConst(BG_DT_ParC.getVal()),RooConst(BG_DT_ParD.getVal())));


//RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S_DT,CrystalBall_3S_DT,BackgroundShape_DT,CBShape_DT),
//							   RooArgList(nSignal_1S,nSignal_2S,nSignal_3S,nBackground,nCB));
RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_1S_DT,CrystalBall_2S3S_DT,BackgroundShape,CBShapeLS),
							   RooArgList(nSignal_1S,nSignal_2S3S,nBackground,nCB));
DT_FitFunction.fitTo(DT_UnBin);

TCanvas* cDT = new TCanvas("cDT3","cDT3",750,550) ;
cDT->cd(); gPad->SetLeftMargin(0.15);
RooPlot* plot_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_UnBin.plotOn(plot_DT); 
DT_FitFunction.plotOn(plot_DT,Components(BackgroundShape),LineStyle(2),LineWidth(3),LineColor(1),RooFit::Name("bg3"));
DT_FitFunction.plotOn(plot_DT,Components(CBShapeLS),LineStyle(2),LineWidth(3),LineColor(kBlue),RooFit::Name("bg2"));
DT_FitFunction.plotOn(plot_DT,LineWidth(2),LineColor(1),RooFit::Name("fiit"));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_1S_DT),LineStyle(1),LineColor(myDarkRed),LineWidth(3));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_2S3S_DT),LineStyle(1),LineColor(myDarkGreen),LineWidth(3));
DT_UnBin.plotOn(plot_DT,MarkerSize(1.1),MarkerColor(kRed),LineColor(kRed));
DT_UnBinLS.plotOn(plot_DT,MarkerSize(1.1),MarkerColor(kBlue),LineColor(kBlue));
plot_DT->GetYaxis()->SetTitleOffset(1.4);
plot_DT->GetYaxis()->SetTitleSize(0.04);
plot_DT->GetXaxis()->SetTitleSize(0.04);
plot_DT->GetYaxis()->SetTitle(TString::Format("Counts/%.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));

plot_DT->Draw();
RooChi2Var chi2ndf("chi2ndf","chi2ndf",DT_FitFunction,DT_hist);

TLegend *myLegend_DT = new TLegend(0.064,0.588,0.547,0.87);
myLegendSetUp(myLegend_DT,0.03,1);
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{1S} = %.1f #pm %.1f",nSignal_1S.getVal(),nSignal_1S.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{2S+3S} = %.1f #pm %.1f",nSignal_2S3S.getVal(),nSignal_2S3S.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{cbbg} = %.1f #pm %.1f",nCB2.getVal(),nCB2.getError())," ");
myLegend_DT->AddEntry((TObject*)0,TString::Format("N_{corrbg} = %.1f #pm %.1f",nBackground.getVal(),nBackground.getError())," ");
//myLegend_DT->AddEntry((TObject*)0,TString::Format("#chi^{2}/ndf = %.1f / %i",chi2ndf.getVal(),nBins)," ");
myLegend_DT->Draw();

TLegend *myLegend_BG = new TLegend(0.57,0.60,0.94,0.87);
myLegendSetUp(myLegend_BG,0.03,1);
myLegend_BG->AddEntry(plot_DT->findObject("fiit"), "sig. + comb. bg + corr. bg", "l");
myLegend_BG->AddEntry(plot_DT->findObject("bg2"), "comb. bg", "l");
myLegend_BG->AddEntry(plot_DT->findObject("bg3"), "corr. bg", "l");
TString textstr = Form("#bf{%i-%i %}",centValues[centralityMax],centValues[centralityMin]);
cout << "string is " << textstr << endl;
myLegend_BG->AddEntry((TObject*)0,textstr," ");
/*TLegend *myLegend_BG = new TLegend(0.46,0.55,0.94,0.87);
myLegendSetUp(myLegend_BG,0.03,1);
myLegend_BG->AddEntry((TObject*)0,"BG = A*x^B/(1+x/C)^D"," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("A = %.1f #pm %.1f",BG_DT_ParA.getVal(),BG_DT_ParA.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("B = %.1f #pm %.1f",BG_DT_ParB.getVal(),BG_DT_ParB.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("C = %.1f #pm %.1f",BG_DT_ParC.getVal(),BG_DT_ParC.getError())," ");
myLegend_BG->AddEntry((TObject*)0,TString::Format("D = %.1f #pm %.1f",BG_DT_ParD.getVal(),BG_DT_ParD.getError())," ");
*/myLegend_BG->Draw();

cDT->SaveAs(TString::Format("ME_DataFit_Cent%d_%d.png",centralityMin,centralityMax));
#endif

}


void FitUpsilonMassAll(Int_t centralityMin = 7, Int_t centralityMax = 9) {

	FitUpsilonMass(centralityMin, centralityMax);
	//FitUpsilonMassUSLS(centralityMin, centralityMax);
	//FitUpsilonMassME(centralityMin, centralityMax);

   	// SAVE YIELDS IN FILE
   	ofstream outfile;
  	outfile.open("tmp.txt");
  	outfile << 0.5*(yields[0]+yields2[0]) << endl;
  	outfile << 0.5*(yields[1]+yields2[1]) << endl;
  	outfile << 0.5*(yields[2]+yields2[2]) << endl;
  	outfile << 0.5*(yields[3]+yields2[3]) << endl;
	outfile.close();

}
