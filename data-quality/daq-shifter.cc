// daq-shifter.cc
// Used to look at the M2 behavior after the re-bias on 11/10/2016.
// C. Wiseman

#include <iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMath.h"
#include "GATDataSet.hh"

using namespace std;

char theCut[1000];

void energySpec1()
{
  vector<int> runs = {60002399,19502,19515};
  for (auto i : runs) {
    GATDataSet ds(i);
    TChain *gatChain = ds.GetGatifiedChain();
    TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
    c1->SetLogy();
  	TH1D *h1 = new TH1D("h1","h1",1000,200,9000);
  	gatChain->Project("h1","trapE",TString::Format("channel%2==1"));
    h1->SetLineColor(kBlue);
    h1->Draw();
    TH1D *h2 = new TH1D("h2","h2",1000,200,9000);
    gatChain->Project("h2","trapE",TString::Format("channel%2==0"));
    h2->SetLineColor(kRed);
    h2->Draw("same");
  	h1->GetXaxis()->SetTitle("Energy (trapE)");
  	TLegend* l1 = new TLegend(0.7,0.7,0.87,0.92);
  	l1->AddEntry(h1,"HG","l");
    l1->AddEntry(h2,"LG","l");
  	l1->Draw("SAME");
  	c1->Update();
  	c1->Print(TString::Format("energySpec1-%i.pdf",i));
  }
}

void energySpec2(GATDataSet &ds)
{
  int lo = ds.GetRunNumber(0), hi=0, nruns=(int)ds.GetNRuns();
  char name[200];
  sprintf(name,"energySpec2-%i.pdf",lo);
  if (nruns>1) {
    hi = ds.GetRunNumber(nruns-1);
    sprintf(name,"energySpec2-%i-%i.pdf",lo,hi);
  }
  TChain *gatChain = ds.GetGatifiedChain();
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogy();
	TH1D *h1 = new TH1D("h1","h1",100,200,3000);
	gatChain->Project("h1","trapE",TString::Format("channel%2==1"));
  h1->SetLineColor(kBlue);
  h1->Draw();
	h1->GetXaxis()->SetTitle("Energy (trapE)");
	TLegend* l1 = new TLegend(0.7,0.7,0.87,0.92);
	l1->AddEntry(h1,"HG","l");
	l1->Draw("SAME");
	c1->Update();
	c1->Print(name);
}

void channelEnergyVsTime(GATDataSet &ds)
{
  TChain *gatChain = ds.GetGatifiedChain();
  int cts = gatChain->Draw("timestamp*1e-8","","GOFF");
  double firstTime = gatChain->GetV1()[0];
  double lastTime = gatChain->GetV1()[cts-1];
  double bins = (lastTime-firstTime)/10;
  printf("%i  %.2f  %.2f  %.2f\n",cts,firstTime,lastTime,lastTime-firstTime);
  int firstRun = ds.GetRunNumber(0);
  MJTChannelMap *map = ds.GetChannelMap();
  MJTChannelSettings *set = ds.GetChannelSettings();
  vector<uint32_t> enabledChannels = set->GetEnabledIDList();
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->cd();
  TVirtualPad *p1 = c1->GetPad(0);
  p1->SetLeftMargin(0.15);
  gatChain->SetMarkerStyle(kFullDotMedium);
  gatChain->SetMarkerColor(kRed);
  gatChain->SetLineColor(kBlue);
  for (auto i : enabledChannels)
  {
    string detPos = map->GetDetectorPos(i);
    if (detPos == "D") continue;
    double maxE = 100., minE = 100.;
    int ctsE = gatChain->Draw("trapE","GOFF");
    for(int i=0; i<ctsE; i++) {
      if(gatChain->GetV1()[i] > maxE) maxE = gatChain->GetV1()[i];
      if(gatChain->GetV1()[i] < minE) minE = gatChain->GetV1()[i];
    }
    double rangeB = maxE-minE;
    if(rangeB < 10) rangeB = 50;
    maxE = TMath::Ceil(maxE + rangeB*0.1);
    minE = TMath::Floor(minE - rangeB*0.1);
    TH2D *h1 = new TH2D("h1","h1", (int)bins,0,lastTime-firstTime, 10, minE, maxE);
    h1->GetXaxis()->SetTitle("Time After Start (sec)");
    h1->GetYaxis()->SetTitle(TString::Format("%s trapE",detPos.c_str()));
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->Draw();
    gatChain->Draw(TString::Format("trapE:(timestamp*1e-8 - %f)",firstTime),TString::Format("channel==%i",i),"LP same");
    c1->Print(TString::Format("eVsT-%i-%i.pdf",i,firstRun));
    delete h1;
  }
}

void energyVsTime(GATDataSet &ds)
{
  // right now this is edited to give all baselines vs. time

  TChain *gatChain = ds.GetGatifiedChain();
  int lo = ds.GetRunNumber(0), hi=0, nruns=(int)ds.GetNRuns();
  char name[200];
  // sprintf(name,"enVsT-%i.png",lo);
  sprintf(name,"blVsT-%i.png",lo);
  if (nruns>1) {
    hi = ds.GetRunNumber(nruns-1);
    // sprintf(name,"enVsT-%i-%i.png",lo,hi);
    sprintf(name,"blVsT-%i-%i.png",lo,hi);
  }
  MJTChannelMap *theMap = ds.GetChannelMap();
  vector<uint32_t> pulserChans = theMap->GetPulserChanList();
  string pulserChanCut = "(";
  for (auto i : pulserChans){
    pulserChanCut += "channel!="+to_string(i)+" &&";
  }
  pulserChanCut.pop_back(); pulserChanCut.pop_back(); // remove the last "&&"
  pulserChanCut += ")";

  vector<double> *ts=0;
  gatChain->SetBranchAddress("timestamp",&ts);
  gatChain->GetEntry(0);
  double firstTime=ts->at(0)*1e-8;
  gatChain->GetEntry(gatChain->GetEntries()-1);
  double lastTime=ts->at(0)*1e-8;
  double bins = (lastTime-firstTime)/10;
  printf("%.2f  %.2f  %.2f\n",firstTime,lastTime,lastTime-firstTime);

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",1400,800);
  c1->cd();
  TVirtualPad *p1 = c1->GetPad(0);
  p1->SetLeftMargin(0.15);
  // gatChain->SetMarkerStyle(kFullDotMedium);
  gatChain->SetMarkerStyle(kDot);
  gatChain->SetMarkerColor(kRed);
  gatChain->SetLineColor(kBlue);
  double maxE = 100., minE = 100.;
  // int ctsE = gatChain->Draw("trapE","channel%2==0","GOFF");
  int ctsE = gatChain->Draw("RawWFblOffset","channel%2==0 && (RawWFblOffset > -200 && RawWFblOffset < 300)","GOFF");

  for(int i=0; i<ctsE; i++) {
    if(gatChain->GetV1()[i] > maxE) maxE = gatChain->GetV1()[i];
    if(gatChain->GetV1()[i] < minE) minE = gatChain->GetV1()[i];
  }
  double rangeB = maxE-minE;
  if(rangeB < 10) rangeB = 50;
  maxE = TMath::Ceil(maxE + rangeB*0.1);
  minE = TMath::Floor(minE - rangeB*0.1);
  TH2D *h1 = new TH2D("h1","h1", 10,0,lastTime-firstTime, 10,minE,maxE);
  h1->GetXaxis()->SetTitle("Time After Start (sec)");
  h1->GetYaxis()->SetTitle("Baseline (ADC)");
  h1->GetYaxis()->SetTitleOffset(1.4);
  h1->Draw();
  // gatChain->Draw(TString::Format("trapE:(timestamp*1e-8 - %f)",firstTime),TString::Format("channel%2==0 && %s",pulserChanCut.c_str()),"LP same");
  gatChain->Draw(TString::Format("RawWFblOffset:(timestamp*1e-8 - %f)",firstTime),TString::Format("channel%2==0 && (RawWFblOffset > -200 && RawWFblOffset < 300) && %s",pulserChanCut.c_str()),"LP same");
  c1->Print(name);
}

void waveforms()
{
  // save ~50 random waveforms for each run into a ROOT file.

  // 19502 - 19515 - night of 11/10/2016
  // 19615 - 19624 - night of 11/11/2016

  vector<int> runs = {19502};
  for (auto i : runs) {
    GATDataSet ds(i);
    TChain *gatChain = ds.GetGatifiedChain();
    TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  }
}


int main()
{
  gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");

  // energySpec1();

  // GATDataSet ds1(19502,19515);
  // energySpec2(ds1);
  // GATDataSet ds2(19436,19450);
  // energySpec2(ds2);

  // full run list
  vector<int> runs = {19502,19503,19504,19505,19506,19507,19508,19509,19510,19511,19512,19513,19514,19515};
  GATDataSet ds;
  for (auto i : runs) ds.AddRunNumber(i);
  energyVsTime(ds);

  // waveforms();


}


