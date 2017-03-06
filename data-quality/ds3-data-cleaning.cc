#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEntryList.h"
#include "GATDataSet.hh"

using namespace std;

void ExamplePlot();

int main()
{
  gStyle->SetOptStat(0);
  // gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
  // gROOT->ForceStyle();
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;"); // suppresses "file creation" messages

  ExamplePlot();
}

char theCut[1000];
char basicCut[1000] = "trapENFDBSGCal < 200 && trapENFDBSGCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0";

"gat_variable < xxx && gat_variable > xxxx"

void ExamplePlot()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS3/gatrev_37365928/*.root");
  cout << "I found " << skimTree->GetEntries() << " entries.\n";

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogy(1);

  // basic cut spectrum
  TH1D *h1 = new TH1D("h1","h1",1000,2,200); // bins,x-lo,x-hi
  sprintf(theCut,"%s",basicCut);
  skimTree->Project("h1","trapENFDBSGCal",basicCut);
  h1->GetXaxis()->SetTitle("Energy (trapENFDBSGCal)");
  h1->SetLineColor(kBlue);
  h1->Draw();
  double tot = h1->Integral(h1->GetXaxis()->FindBin(0.8),h1->GetXaxis()->FindBin(10));
  cout << "There are " << tot << " counts in this range.\n";

  c1->Print("./plots/ds3-ExamplePlot.pdf");
}