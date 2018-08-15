// rootlogon.C
// Christopher Powell
using namespace std;

void style()
{// Add my own options here:
  TStyle* myStyle = new TStyle("myStyle","Chris's Root Styles");
  myStyle->SetPalette(1,0); // avoid horrible default color scheme
  //myStyle->SetOptStat(0);
  myStyle->SetOptTitle(1);
  myStyle->SetOptDate(0);
  myStyle->SetStatColor(10);
  myStyle->SetStatFontSize(0.05);
  myStyle->SetStatH(0.26);
  myStyle->SetStatW(0.26);
  myStyle->SetTitleFont(42,"xyz"); // font option 
  myStyle->SetLabelFont(42,"xyz");
  myStyle->SetLabelSize(0.051,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.053,"xz"); // size of axis title font
  myStyle->SetTitleSize(0.053,"y"); // size of axis title font
  myStyle->SetTitleOffset(1.1,"x");
  myStyle->SetTitleOffset(1.2,"y");
  myStyle->SetTitleOffset(1.0,"z");
  myStyle->SetNdivisions(10, "x");
  myStyle->SetNdivisions(10, "y");
  myStyle->SetPadBottomMargin(0.13); //margins...
  myStyle->SetPadTopMargin(0.12);
  myStyle->SetPadLeftMargin(0.16);
  myStyle->SetPadRightMargin(0.13);
  myStyle->SetTitleFillColor(10);
  myStyle->SetLineWidth(2);
  myStyle->SetHistLineWidth(2);
  // default canvas options
  myStyle->SetCanvasDefW(700);
  myStyle->SetCanvasDefH(600);
  //myStyle->SetCanvasColor(10);
  myStyle->SetCanvasColor(0);// canvas...
  myStyle->SetCanvasBorderMode(0);
  //myStyle->SetCanvasBorderMode(-1);
  myStyle->SetCanvasBorderSize(0);
  //myStyle->SetCanvasBorderSize(1);
  myStyle->SetPadColor(0);
  myStyle->SetPadBorderSize(1);
  myStyle->SetPadBorderMode(-1);
  myStyle->SetPadGridX(0); // grids, tickmarks
  myStyle->SetPadGridY(0);
  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameBorderMode(1);
  //myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameFillColor(0);
  //myStyle->SetFrameFillColor(10);
  myStyle->SetFrameLineWidth(1.2);
  myStyle->SetPaperSize(20,24); // US letter size
  gROOT->SetStyle("myStyle");
  cout << "Styles are Set!" << endl;
  return; 
}

