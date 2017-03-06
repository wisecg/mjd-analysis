// DS-1 data cleaning.
// Creates plots to justify TCuts.
// C. Wiseman, 12/11/2016

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

void LargeFileCuts();
void TCutSkimmer();
void FindNoisyRuns();
void IntegrateSpec(TH1D *h, double &tot, double &lower, double &upper);
void CompareCuts();
void ListActiveDetectors();
void ListActiveRuns();
void CalculateLivetime();
void StandardCutSpectrum();
void TuneTOverECut();
void TOverECutSpectrum();
void EnrVsNatSpectrum();
void DS0DS1EnrichedSpectra();
void CheckEnrichedSpectra();
void CheckWaveSkim();
void CheckNoiseTag();
void CheckCalibFFTCut();
void CheckMCMCCut();

int main()
{
  /*
  Source files:
  - DS1 background (PDSF):  ~/ds1-low/skimDS1_complete.root
  - Calibration file:  ./data/skimDS1_calib_13550.root
  - tcut-skimmed background basic+TETM cut:  ./data/skimDS1_basicCuts.root
  - tcut-skimmed calibration basic+TETM cut:  ./data/skimDS1_calib_basicCuts.root
  - ds-0 enriched & natural spectra: ./data/m1DiagnosticTree4Aug2016.root
  */
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
  gROOT->ForceStyle();
  // gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");

  // LargeFileCuts();         // 1. PDSF: Make basic + TETM cut spectra and save histos to a root file
  // TCutSkimmer();           // 2. Use a TEntryList to apply the basic+TETM cut to a new file
  // FindNoisyRuns();         // 3. Plot rate under 5 keV, make a list of runs w/ rates over 100Hz
  // CompareCuts();           // 4. Rank the reduction from different data cleaning cuts
  // ListActiveDetectors();   // 5. Combine all "standard" cuts and list active detectors
  // ListActiveRuns();        // 6. Combine all "standard" cuts and list active runs
  // CalculateLivetime();     // 7. PDSF: Calculate the raw and "standard cut" livetime
  // StandardCutSpectrum();   // 8. Look at the reduction between the basic cut and standard cut.
  // TuneTOverECut();         // 9. Tune the T/E cut for DS-1
  // TOverECutSpectrum();     // 10. Look at the T/E reduction from the standard cut.
  // EnrVsNatSpectrum();      // 11. Enriched versus natural spectra (with & without T/E)
  // DS0DS1EnrichedSpectra(); // 12. Compare enriched spectra using Kris's ROOT file.
  // CheckEnrichedSpectra();  // 13. Compare DS-1 enriched spectra and integrate counts.
  // CheckWaveSkim();         // 14. Used TCutSkimmer + standardCut, now check file
  // CheckNoiseTag();         // 15. Investigate using "d2wf0MHzTo50MHzPower" to eliminate noise
  // CheckCalibFFTCut();      // 16. PDSF: Check the FFT cut on BG and calibration data
  CheckMCMCCut();          // 17. Check MCMC cut
}

char theCut[1000];
char basicCut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0";
char burstCut[1000] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3)";
// char noisyRunsCut[1000] = "run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";
char noisyRunsCut[1000] = "run!=9648 && run!=10663 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004 && run!=13205 && run!=13306 && run!=13307 && run!=13330 && run!=13400";
// char noisyRunsOnly[1000] = "(run==9648 || run==10663 || run==10745 || run==11175 || run==12445 || run==12723 || run==12735 || run==12745 || run==12746 || run==12765 || run==12766 || run==12767 || run==13004)";
char noisyRunsOnly[1000] = "(run==9648 || run==10663 || run==11175 || run==12445 || run==12723 || run==12735 || run==12746 || run==12765 || run==12766 || run==12767 || run==13004 || run==13205 || run==13306 || run==13307 || run==13330 || run==13400)";
char standardCut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run!=9648 && run!=10663 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004 && run!=13205 && run!=13306 && run!=13307 && run!=13330 && run!=13400 && trapETailMin < 0 && mH==1 && !isLNFill1 && isGood && !muVeto && !wfDCBits";
char ds0_toeCut[1000] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";
char ds1_toeCut[1000] = "(((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";
char calibCut[1000] = "trapENFCal < 500 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && channel%2==0 && trapETailMin < 0 && mH==1 && channel!=616 && channel!=594";

void LargeFileCuts()
{
  // string inFile = "~/ds1-low/skimDS1_complete.root";
  // string outFile = "./data/skimDS1_basicCutSpectra.root";
  string inFile = "~/ds1-low/skimDS1_calib_13550-13556.root";
  string outFile = "./data/skimDS1_calib_basicCutSpectra.root";

  TFile *f = new TFile(outFile.c_str(),"RECREATE");

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add(inFile.c_str());
  cout << "Found " << skimTree->GetEntries() << " entries.\n";

  // 1 basic cut - raw HG energy spectrum
  TH1D *h1 = new TH1D("basicCut","h1",250,0,50);
  sprintf(theCut,"%s",basicCut);
  int cts = skimTree->Project("basicCut","trapENFCal",theCut);
  h1->Write();

  // 2 trapETailMin cut - significantly cuts down on file size.
  TH1D *h2 = new TH1D("basicTETMCut","h2",250,0,50);
  sprintf(theCut,"%s && trapETailMin < 0",basicCut);
  int cts2 = skimTree->Project("basicTETMCut","trapENFCal",theCut);
  h2->Write();

  cout << "basic cut: " << cts << "  basic+TETM cut: " << cts2 << endl;

  f->Close();
}

void TCutSkimmer()
{
  // string inFile = "~/ds1-low/skimDS1_complete.root";
  // string outFile = "./data/skimDS1_basicCuts.root";
  // string inFile = "~/ds1-low/skimDS1_calib_13550-13556.root";
  string inFile = "./data/skimDS1_calib_13550.root";
  string outFile = "./data/skimDS1_calibCut.root";
  // string inFile = "./data/skimDS1_basicCuts.root";
  // string outFile = "./data/skimDS1_standardCut.root";
  // string outFile = "./data/skimDS1_standardCut13078.root"; // for debugging wave-skim

  // sprintf(theCut,"%s && trapETailMin < 0",basicCut);
  sprintf(theCut,"%s",calibCut);
  // sprintf(theCut,"%s && run > 13078",standardCut);
  cout << "Skimming file " << inFile << " using this cut: " << theCut << "\n\n";

  TChain *skim = new TChain("skimTree");
  skim->Add(inFile.c_str());
  skim->Draw(">>entList",theCut,"entrylist GOFF");
  TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
  skim->SetEntryList(entList);
  TFile *f2 = new TFile(outFile.c_str(),"recreate");
  TTree *small = skim->CopyTree("");
  small->Write();
  cout << "Wrote " << small->GetEntries() << " entries.\n";
  TNamed thisCut("cutUsedHere",theCut);	// save the cut used into the file.
  thisCut.Write();
  f2->Close();
}

void FindNoisyRuns()
{
  double maxRate = 100;

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  c1->SetLogy();

  sprintf(theCut,"%s && %s && trapENFCal < 5",basicCut,burstCut);
  TH1D *h0 = new TH1D("h0","h0",4970,9420,14390);
  skimTree->Project("h0","run",theCut);
  h0->GetXaxis()->SetTitle("Run");
  h0->GetYaxis()->SetTitle("Rate < 5keV (cts/run)");
  h0->GetXaxis()->SetNdivisions(509);
  h0->Draw();
  vector<int> noisyRuns;
  for (int i = 0; i < 4970; i++)
  {
    double run = h0->GetXaxis()->GetBinCenter(i) - 0.5;
    double rate = h0->GetBinContent(i);
    if (rate > maxRate) {
      // cout << run << "  " << rate << endl;
      noisyRuns.push_back(run);
    }
  }
  c1->Print("./plots/FindNoisyRuns.pdf");

  char noisyCut[1000];
  char noisyOnly[1000];
  string cut1, cut2 = "(";
  for (int i = 0; i < (int)noisyRuns.size()-1; i++){
    cut1 += "run!="+to_string(noisyRuns[i])+" && ";
    cut2 += "run=="+to_string(noisyRuns[i])+ " || ";
  }
  cut1 += "run!="+to_string(noisyRuns[noisyRuns.size()-1]);
  cut2 += "run=="+to_string(noisyRuns[noisyRuns.size()-1])+")";
  strcpy(noisyCut, cut1.c_str());
  strcpy(noisyOnly, cut2.c_str());
  cout << cut1 << endl;
  cout << cut2 << endl;

  // basic cut spectrum
  TH1D *h1 = new TH1D("h1","h1",150,0,10);
  sprintf(theCut,"%s",basicCut);
  skimTree->Project("h1","trapENFCal",basicCut);
  h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h1->SetLineColor(kBlack);
  h1->Draw();
  double tot = h1->Integral(h1->GetXaxis()->FindBin(0.8),h1->GetXaxis()->FindBin(10));

  // remove noisy runs
  TH1D *h2 = new TH1D("h2","h2",150,0,10);
  sprintf(theCut,"%s && %s",basicCut,noisyCut);
  skimTree->Project("h2","trapENFCal",theCut);
  h2->SetLineColor(kBlue);
  h2->Draw("same");
  double fracLoNoise = 100*((h2->Integral(h2->GetXaxis()->FindBin(0.8),h2->GetXaxis()->FindBin(10)))/tot);

  // noisy runs only
  TH1D *h3 = new TH1D("h3","h3",150,0,10);
  sprintf(theCut,"%s && %s",basicCut,noisyOnly);
  skimTree->Project("h3","trapENFCal",theCut);
  h3->SetLineColor(kRed);
  h3->Draw("same");
  double fracNoisy = 100*((h3->Integral(h3->GetXaxis()->FindBin(0.8),h3->GetXaxis()->FindBin(10)))/tot);

  TLegend* leg1 = new TLegend(0.3,0.65,0.87,0.92);
  leg1->AddEntry(h1,"Standard DC Cut","l");
  leg1->AddEntry(h2,TString::Format("With Run Cut (%.1f%% persist)",100-fracLoNoise),"l");
  leg1->AddEntry(h3,TString::Format("Noisy Runs Only (%.1f%% removed)",100-fracNoisy),"l");
  leg1->Draw("SAME");
  c1->Update();
  c1->Print("./plots/noisyCut.pdf");

  // channels without noisy runs
  c1->SetLogy(0);
  c1->SetLogz(1);
  TH2D *h4 = new TH2D("h4","h4",150,0,10,125,575,700);
  sprintf(theCut,"%s && %s && trapENFCal < 10",basicCut,noisyCut);
  skimTree->Project("h4","channel:trapENFCal",theCut);
  h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h4->GetYaxis()->SetTitle("Channel");
  h4->Draw("COLZ");
  c1->Print("./plots/noisyRuns_EnergyVsChannel.pdf");

  // channels w/ noisy runs only.
  TH2D *h5 = new TH2D("h5","h5",150,0,10,125,575,700);
  sprintf(theCut,"%s && %s && trapENFCal < 10",basicCut,noisyOnly);
  skimTree->Project("h5","channel:trapENFCal",theCut);
  h5->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h5->GetYaxis()->SetTitle("Channel");
  h5->Draw("COLZ");
  c1->Print("./plots/noisyRunsOnly_EnergyVsChannel.pdf");
}

void IntegrateSpec(TH1D *h1, double &tot, double &lower, double &upper)
{
  TAxis *ax1 = h1->GetXaxis();
  tot = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(50));
  lower = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5));
  upper = h1->Integral(ax1->FindBin(5),ax1->FindBin(50));
}

void CompareCuts()
{
  TFile *f1 = new TFile("./data/skimDS1_basicCutSpectra.root");
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  double tot, lower, upper, totPrev, lowerPrev, upperPrev;

  // 1 basic cut
  TH1D *h1 = (TH1D*)f1->Get("basicCut");
  IntegrateSpec(h1,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f  [0.8 - 5]: %.0f  [5 - 50]: %.0f \n","Basic Cut",tot,lower,upper);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 2 basic cut + trapETailMin
  TH1D *h2 = (TH1D*)f1->Get("basicTETMCut");
  IntegrateSpec(h2,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "Basic+TETM", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracTETM = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 3 low-e noisy runs cut
  sprintf(theCut,"%s && trapETailMin < 0 && %s",basicCut,noisyRunsCut);
  TH1D *h3 = new TH1D("h3","h3",250,0,50);
  skimTree->Project("h3","trapENFCal",theCut);
  IntegrateSpec(h3,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "Noisy Runs", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracNoisyRuns = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 4 single-detector (mH==1)
  sprintf(theCut,"%s && trapETailMin < 0 && mH==1 && %s",basicCut,noisyRunsCut);
  TH1D *h4 = new TH1D("h4","h4",250,0,50);
  skimTree->Project("h4","trapENFCal",theCut);
  IntegrateSpec(h4,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "mH==1", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracMH = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 5 !isLNFill1
  sprintf(theCut,"%s && trapETailMin < 0 && !isLNFill1 && %s && mH==1",basicCut,noisyRunsCut);
  TH1D *h5 = new TH1D("h5","h5",250,0,50);
  skimTree->Project("h5","trapENFCal",theCut);
  IntegrateSpec(h5,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "!isLNFill1", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracLN = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 6 isGood
  sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill1",basicCut,noisyRunsCut);
  TH1D *h6 = new TH1D("h6","h6",250,0,50);
  skimTree->Project("h6","trapENFCal",theCut);
  IntegrateSpec(h6,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "isGood", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracGood = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 7 !muVeto
  sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill1 && !muVeto",basicCut,noisyRunsCut);
  TH1D *h7 = new TH1D("h7","h7",250,0,50);
  skimTree->Project("h7","trapENFCal",theCut);
  IntegrateSpec(h7,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "!muVeto", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracVeto = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 8 burst cut
  sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill1 && !muVeto && %s",basicCut,noisyRunsCut,burstCut);
  TH1D *h8 = new TH1D("h8","h8",250,0,50);
  skimTree->Project("h8","trapENFCal",theCut);
  IntegrateSpec(h8,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "!burst", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracBurst = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // 9 !wfDCBits
  sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill1 && !muVeto && %s && !wfDCBits",basicCut,noisyRunsCut,burstCut);
  TH1D *h9 = new TH1D("h9","h9",250,0,50);
  skimTree->Project("h9","trapENFCal",theCut);
  IntegrateSpec(h9,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "!wfDCBits", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracWFDC = 100*(tot/totPrev);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  // Now that we've established the proper order of the cuts,
  // draw them all together to see the reduction.

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  c1->SetLogy();
  h1->SetLineColor(1);
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  h3->SetLineColor(3);
  h3->Draw("same");
  h4->SetLineColor(4);
  h4->Draw("same");
  h5->SetLineColor(6);
  h5->Draw("same");
  h6->SetLineColor(7);
  h6->Draw("same");
  h7->SetLineColor(kOrange);
  h7->Draw("same");
  h8->SetLineColor(kMagenta+3);
  h8->Draw("same");
  h9->SetLineColor(kMagenta+3);
  h9->Draw("same");
  TLegend* leg1 = new TLegend(0.5,0.3,0.87,0.92);
  leg1->AddEntry(h1,"HG Only (% reduction)","l");
  leg1->AddEntry(h2,TString::Format("+TETM (%.1f%%)",100-fracTETM),"l");
  leg1->AddEntry(h3,TString::Format("+!NoisyRuns (%.1f%%)",100-fracNoisyRuns),"l");
  leg1->AddEntry(h4,TString::Format("+mH==1 (%.1f%%)",100-fracMH),"l");
  leg1->AddEntry(h5,TString::Format("+!isLNFill1 (%.1f%%)",100-fracLN),"l");
  leg1->AddEntry(h6,TString::Format("+isGood (%.1f%%)",100-fracGood),"l");
  leg1->AddEntry(h7,TString::Format("+!muVeto (%.1f%%)",100-fracVeto),"l");
  leg1->AddEntry(h8,TString::Format("+!burst (%.1f%%)",100-fracBurst),"l");
  leg1->AddEntry(h9,TString::Format("+!wfDCBits (%.1f%%)",100-fracWFDC),"l");
  leg1->Draw("SAME");
  c1->Update();
  c1->Print("./plots/CompareCuts.pdf");
}

void ListActiveDetectors()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");
  GATDataSet ds(13335);
  MJTChannelMap *chanMap = ds.GetChannelMap();

  vector<int> chanList;
  sprintf(theCut,"%s",standardCut);
  int cts = skimTree->Draw("channel",theCut,"GOFF");
  cout << "Found " << cts << " entries\n";
  for (int i = 0; i < cts; i++)
  {
    int chan = skimTree->GetV1()[i];
    bool newChan = true;
    for (auto j : chanList) {
      if (chan == j) { newChan = false; break; }
    }
    if (newChan) chanList.push_back(chan);
  }
  cout << "Final channel list: \n";
  for (auto i : chanList)
    cout << i << " " << chanMap->GetDetectorPos(i) << endl;
}

void ListActiveRuns()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  vector<int> runList;
  sprintf(theCut,"%s",standardCut);
  int cts = skimTree->Draw("run",theCut,"GOFF");
  cout << "Found " << cts << " entries\n";
  for (int i = 0; i < cts; i++)
  {
    int run = skimTree->GetV1()[i];
    if (!(find(runList.begin(), runList.end(), run) != runList.end()))
      runList.push_back(run);
  }
  ofstream runFile("./runs/ds1-standardCutRuns.txt");
  cout << "Final run list: \n";
  for (auto i : runList) {
    cout << i << endl;
    runFile << i << endl;
  }
  runFile.close();
}

void CalculateLivetime()
{
  int run;
  cout << "Finding complete livetime ... ";
  ifstream runFile1("./runs/ds1-complete.txt");
  GATDataSet *ds1 = new GATDataSet();
  while (runFile1 >> run) ds1->AddRunNumber(run);
  cout << ds1->GetRunTime()/1e9 << " seconds.\n";
  delete ds1;

  cout << "Finding livetime after standard cuts ...";
  ifstream runFile2("./runs/ds1-standardCutRuns.txt");
  GATDataSet *ds2 = new GATDataSet();
  while (runFile2 >> run) ds2->AddRunNumber(run);
  cout << ds2->GetRunTime()/1e9 << " seconds.\n";
  delete ds2;
}

void StandardCutSpectrum()
{
  TFile *f1 = new TFile("./data/skimDS1_basicCutSpectra.root");
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  double tot, lower, upper, totPrev, lowerPrev, upperPrev;

  TH1D *h1 = (TH1D*)f1->Get("basicCut");
  IntegrateSpec(h1,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f  [0.8 - 5]: %.0f  [5 - 50]: %.0f \n","Basic Cut",tot,lower,upper);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  sprintf(theCut,"%s",standardCut);
  TH1D *h2 = new TH1D("h2","h2",250,0,50);
  skimTree->Project("h2","trapENFCal",theCut);
  IntegrateSpec(h2,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "Standard Cut", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracStandard = 100*(tot/totPrev);

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  c1->SetLogy();
  h1->SetLineColor(1);
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  TLegend* l1 = new TLegend(0.4,0.7,0.87,0.92);
  l1->AddEntry(h1,"HG Only (% reduction)","l");
  l1->AddEntry(h2,TString::Format("+Standard Cut (%.1f%%)",100-fracStandard),"l");
  l1->Draw("SAME");
  c1->Update();

  c1->Print("./plots/StandardCutReduction.pdf");

}

void TuneTOverECut()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  // channel list is from ListActiveDetectors().  P2D1-692 and P7D1-600 are BEGes.
  int chans[16] = {582, 580, 578, 692, 648, 640, 610, 608, 664, 672, 632, 626, 690, 600, 598, 592};

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",2600,800);
  c1->Divide(3,1,0,0);
  char hname[100];
  TH1D *htoe[16];
  TH1D *htoe2[16];
  TH2D *htoe2D[16];
  for (int i=0; i<16; i++)
  {
    c1->cd(1);
    TVirtualPad *p1 = c1->GetPad(1);
    p1->SetTopMargin(0.05);
    p1->SetLogx(0);
    sprintf(hname,"h_%i",chans[i]);
    htoe[i] = new TH1D(hname,hname,200,0,10);
    sprintf(theCut,"%s && channel==%i",standardCut,chans[i]);
    skimTree->Project(hname,"kvorrT/trapENFCal",theCut);
    htoe[i]->GetXaxis()->SetTitle(TString::Format("kvorrT/trapENFCal (%i)",chans[i]));
    htoe[i]->GetXaxis()->SetTitleOffset(1.2);
    htoe[i]->Draw();

    c1->cd(2);
    TVirtualPad *p2 = c1->GetPad(2);
    p2->SetTopMargin(0.05);
    p2->SetLeftMargin(0.1);
    p2->SetRightMargin(0.1);
    p2->SetLogx(1);
    p2->SetLogz(1);
    sprintf(hname,"h_2d_%i",chans[i]);
    htoe2D[i] = new TH2D(hname,hname,1000,0,200,200,0,10);
    skimTree->Project(hname,"kvorrT/trapENFCal:trapENFCal",theCut);
    htoe2D[i]->GetXaxis()->SetTitle(TString::Format("trapENFCal (%i)",chans[i]));
    htoe2D[i]->GetXaxis()->SetTitleOffset(1.2);
    htoe2D[i]->GetYaxis()->SetTitle("T/E");
    htoe2D[i]->SetMinimum(1);
    if (htoe2D[i]->GetMaximum() < 100) htoe2D[i]->SetMaximum(100);
    htoe2D[i]->Draw("COLZ");

    c1->cd(3);
    TVirtualPad *p3 = c1->GetPad(3);
    p3->SetTopMargin(0.05);
    p3->SetLeftMargin(0.05);
    p3->SetRightMargin(0.05);
    p3->SetLogx(0);
    sprintf(hname,"h_%i",chans[i]);
    htoe2[i] = new TH1D(hname,hname,500,0.,2.5);
    sprintf(theCut,"%s && channel==%i",standardCut,chans[i]);
    skimTree->Project(hname,"kvorrT/trapENFCal",theCut);
    htoe2[i]->GetXaxis()->SetTitle(TString::Format("kvorrT/trapENFCal (%i)",chans[i]));
    htoe2[i]->GetXaxis()->SetTitleOffset(1.2);
    double ymax = htoe2[i]->GetMaximum();
    double binmax = htoe2[i]->GetMaximumBin();
    double xmax = htoe2[i]->GetBinCenter(binmax);
    htoe2[i]->GetXaxis()->SetRange(binmax-80,binmax+80);
    htoe2[i]->Draw();

    cout << chans[i] << "  " << ymax << "  " << binmax << "  " << xmax << endl;

    sprintf(hname,"./plots/toe_%i.pdf",chans[i]);
    c1->Print(hname);
  }
}

void TOverECutSpectrum()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  double tot, lower, upper, totPrev, lowerPrev, upperPrev;

  sprintf(theCut,"%s && isEnr",standardCut);
  TH1D *h1 = new TH1D("h1","h1",250,0,50);
  skimTree->Project("h1","trapENFCal",theCut);
  IntegrateSpec(h1,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f  [0.8 - 5]: %.0f  [5 - 50]: %.0f \n","Basic Cut",tot,lower,upper);
  totPrev=tot; lowerPrev=lower; upperPrev=upper;

  sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
  TH1D *h2 = new TH1D("h2","h2",250,0,50);
  skimTree->Project("h2","trapENFCal",theCut);
  IntegrateSpec(h2,tot,lower,upper);
  printf("%-10s [0.8 - 50]: %.0f (%.2f%%)  [0.8 - 5]: %.0f (%.2f%%)  [5 - 50]: %.0f (%.2f%%)\n", "DS-1 T/E Cut", tot,100*(tot/totPrev), lower,100*(lower/lowerPrev), upper,100*(upper/upperPrev));
  double fracTETM = 100*(tot/totPrev);

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  // c1->SetLogy();
  h1->SetLineColor(1);
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h1->SetMinimum(0.001);
  h1->SetMaximum(40);
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  TLegend* l1 = new TLegend(0.3,0.7,0.87,0.92);
  l1->AddEntry(h1,"isEnr & Standard (% reduction)","l");
  l1->AddEntry(h2,TString::Format("+DS-1 T/E Cut (%.1f%%)",100-fracTETM),"l");
  l1->Draw("SAME");
  c1->Update();

  c1->Print("./plots/TOverECutReduction.pdf");
}

void EnrVsNatSpectrum()
{
  double ds0enrExp = 478;  // taken from RS&DC DS-0 document
  double ds0natExp = 195;
  double ds1enrExp = 662;  // calculated in the above code
  double ds1natExp = 69;
  double binsPerKeV = 5;
  double ene_low = 1;
  double ene_hi = 30;

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  c1->SetLogy();

  // DS-0
  TFile *f1 = new TFile("./data/m1DiagnosticTree4Aug2016.root");
  TTree *t1 = (TTree*)f1->Get("diagTree");

  TH1D *ds0nat = new TH1D("ds0nat","ds0nat",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  t1->Project("ds0nat","calENF","!enr");
  ds0nat->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  ds0nat->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  ds0nat->GetYaxis()->SetTitleOffset(1.1);
  ds0nat->SetLineColor(kBlue);
  ds0nat->Scale(binsPerKeV/ds0natExp);
  ds0nat->SetMinimum(0.01);
  ds0nat->SetMaximum(5);
  ds0nat->Draw();

  TH1D *ds0enr = new TH1D("ds0enr","ds0enr",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  t1->Project("ds0enr","calENF","enr");
  ds0enr->SetLineColor(kRed);
  ds0enr->Scale(binsPerKeV/ds0enrExp);
  ds0enr->Draw("same");

  TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
  leg1->AddEntry(ds0nat,TString::Format("Natural: %.0f kg-d",ds0natExp),"l");
  leg1->AddEntry(ds0enr,TString::Format("Enriched: %.0f kg-d",ds0enrExp),"l");
  leg1->Draw("SAME");
  c1->Update();
  c1->Print("./plots/ds0EnrVsNat-log.pdf");

  // DS-1
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  TH1D *ds1nat = new TH1D("ds1nat","ds1nat",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  sprintf(theCut,"%s && %s && !isEnr",standardCut,ds1_toeCut);
  skimTree->Project("ds1nat","trapENFCal",theCut);
  ds1nat->GetXaxis()->SetTitle("Energy (keV)");
  ds1nat->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  ds1nat->GetYaxis()->SetTitleOffset(1.1);
  ds1nat->SetLineColor(kBlue);
  ds1nat->Scale(binsPerKeV/ds1natExp);
  ds1nat->SetMinimum(0.01);
  ds1nat->SetMaximum(5);
  ds1nat->Draw();

  TH1D *ds1enr = new TH1D("ds1enr","ds1enr",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
  skimTree->Project("ds1enr","trapENFCal",theCut);
  ds1enr->SetLineColor(kRed);
  ds1enr->Scale(binsPerKeV/ds1enrExp);
  ds1enr->Draw("same");

  TLegend* leg2 = new TLegend(0.45,0.75,0.87,0.92);
  leg2->AddEntry(ds1nat,TString::Format("Natural: %.0f kg-d",ds1natExp),"l");
  leg2->AddEntry(ds1enr,TString::Format("Enriched: %.0f kg-d",ds1enrExp),"l");
  leg2->Draw("SAME");
  c1->Update();
  c1->Print("./plots/ds1EnrVsNat-log.pdf");
}

void DS0DS1EnrichedSpectra()
{
  double ds0enrExp = 478;  // taken from RS&DC DS-0 document
  double ds1enrExp = 662;  // calculated in the above code
  double binsPerKeV = 5;
  double ene_low = 1;
  double ene_hi = 30;

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  // c1->SetLogy();

  TFile *f1 = new TFile("./data/m1DiagnosticTree4Aug2016.root");
  TTree *t1 = (TTree*)f1->Get("diagTree");
  TH1D *ds0enr = new TH1D("ds0enr","ds0enr",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  t1->Project("ds0enr","calENF","enr");
  ds0enr->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  ds0enr->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  ds0enr->GetYaxis()->SetTitleOffset(1.1);
  ds0enr->SetLineColor(kBlue);
  ds0enr->Scale(binsPerKeV/ds0enrExp);
  // ds0enr->SetMinimum(0.01);
  ds0enr->SetMaximum(0.3);
  ds0enr->Draw();

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");
  TH1D *ds1enr = new TH1D("ds1enr","ds1enr",(int)binsPerKeV*(ene_hi-ene_low),ene_low,ene_hi);
  sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
  skimTree->Project("ds1enr","trapENFCal",theCut);
  ds1enr->GetXaxis()->SetTitle("Energy (keV)");
  ds1enr->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  ds1enr->GetYaxis()->SetTitleOffset(1.1);
  ds1enr->SetLineColor(kRed);
  ds1enr->Scale(binsPerKeV/ds1enrExp);
  // ds1enr->SetMinimum(0.01);
  // ds1enr->SetMaximum(5);
  ds1enr->Draw("same");

  TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
  leg1->AddEntry(ds0enr,TString::Format("DS0 Enriched: %.0f kg-d",ds0enrExp),"l");
  leg1->AddEntry(ds1enr,TString::Format("DS1 Enriched: %.0f kg-d",ds1enrExp),"l");
  leg1->Draw("SAME");
  c1->Update();
  c1->Print("./plots/DS0EnrVsDS1Enr.pdf");
}

void CheckEnrichedSpectra()
{
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
  // c1->SetLogy();

  double ds1enrExp = 662;  // calculated in the above code
  double binsPerKeV = 5;
  double lo = 1, hi = 50;

  // Standard Cut
  TH1D *h1 = new TH1D("h1","h1",(int)binsPerKeV*(hi-lo),lo,hi);
  sprintf(theCut,"%s && isEnr",standardCut);
  skimTree->Project("h1","trapENFCal",theCut);
  h1->GetXaxis()->SetTitle("Energy (keV)");
  h1->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->SetLineColor(kBlue);
  h1->Scale(binsPerKeV/ds1enrExp);
  // h1->Draw();

  // With T/E Cut - "DNP Cut"
  TH1D *h2 = new TH1D("h2","h2",(int)binsPerKeV*(hi-lo),lo,hi);
  sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
  skimTree->Project("h2","trapENFCal",theCut);
  h2->SetLineColor(kBlue);
  h2->Scale(binsPerKeV/ds1enrExp);
  h2->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h2->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  h2->GetYaxis()->SetTitleOffset(1.1);
  // h2->Draw("same");
  h2->Draw();

  // TLegend* leg2 = new TLegend(0.3,0.7,0.87,0.92);
  // leg2->AddEntry(h1,TString::Format("DS1 Enriched (%.0f kg-d) (no T/E) ",ds1enrExp),"l");
  // leg2->AddEntry(h2,"DS1 Enriched (with T/E)","l");
  // leg2->Draw("SAME");
  c1->Update();
  c1->Print("./plots/DS1EnrichedSpectrum.pdf");

  // Clint, this computes the TOTAL counts between "lo" and "hi" keV, so the
	// result is in cts/kg-d. Then to get the average, divide by the number of bins.
	double tot1 = h1->Integral(h1->GetXaxis()->FindBin(lo),h1->GetXaxis()->FindBin(hi-1),"width");
	double tot2 = h2->Integral(h2->GetXaxis()->FindBin(lo),h2->GetXaxis()->FindBin(hi-1),"width");
	cout << "cts/kg-d                : standard: " << tot1 << "  with T/E: " << tot2 << "  reduction: " << tot1/tot2 << endl;
	cout << "cts/kg-d-keV            : standard: " << tot1/(hi-lo) << "  with T/E: " << tot2/(hi-lo) << "  reduction: " << tot1/tot2 << endl;

}

void CheckWaveSkim()
{
  TChain *waveTree = new TChain("waveTree");
  waveTree->Add("./data/waveSkimDS1_standardCut.root");

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_standardCut.root");

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogy();

  TH1D *h1 = new TH1D("h1","h1",250,0,50);
  int cts = skimTree->Project("h1","trapENFCal",standardCut);
  h1->Draw();
  cout << "skimTree: " << cts << " entries passed cuts.\n";

  TH1D *h2 = new TH1D("h2","h2",250,0,50);
  int cts2 = waveTree->Project("h2","trapENFCal","");
  h2->SetLineColor(kRed);
  h2->Draw("same");
  cout << "waveTree: " << cts2 << " entries in the file.\n";

  // note:
  // skimTree: 12420 entries passed cuts.
  // waveTree: 12416 entries in the file.  (ok, close enough.)

  TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
  leg1->AddEntry(h1,"Standard Cut","l");
  leg1->AddEntry(h2,"waveSkim","l");
  leg1->Draw("SAME");
  c1->Update();
  c1->Print("./plots/waveSkimCheck.pdf");
}

void CheckNoiseTag()
{
  TChain *waveTree = new TChain("waveTree");
  waveTree->Add("./data/waveSkimDS1_standardCut.root");

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas", 800, 600);
  c1->SetLogz();

  sprintf(theCut,"trapENFCal < 10 && trapENFCal > 0.8");
  TH2D *h1 = new TH2D("h1","h1",100,0,10,225,1000,10000);
	double cts = waveTree->Project("h1","d2wf0MHzTo50MHzPower:trapENFCal",theCut);
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("d2wf0MHzTo50MHzPower");
  h1->Draw("COLZ");
  h1->SetMinimum(1);
  c1->Print("./plots/FFTNoiseTagCut_10kev.C");

  // Noise blob: 4100 < d2wf < 6000, 0.8 < trapENF < 1.5
  TAxis *ax = h1->GetXaxis(); // trapENFCal
  TAxis *ay = h1->GetYaxis(); // d2wf0MHzTo50MHzPower
  double blob = h1->Integral(ax->FindBin(0.8),ax->FindBin(1.5),ay->FindBin(4100),ay->FindBin(6000));
  printf("blob %.0f  cts %.0f  remaining %.0f\n",blob,cts,cts-blob);

  // Now do it vs run and try to identify where the two splotches come from.
  c1->SetLogz(0);
  sprintf(theCut,"trapENFCal < 10 && trapENFCal > 0.8");
  TH2D *h2 = new TH2D("h2","h2",200,9420,14390,45,1000,10000); // 1 bin/run: 4970,9420,14390
  waveTree->Project("h2","d2wf0MHzTo50MHzPower:run",theCut);
  h2->GetXaxis()->SetTitle("Run Number");
	h2->GetYaxis()->SetTitle("d2wf0MHzTo50MHzPower");
  h2->Draw("COLZ");
  h2->SetMinimum(1);
  c1->Print("./plots/FFTNoiseTagVsRun.C");



  // Then try to look at the waveforms in the splotches (use wave-view)
}

void CheckCalibFFTCut()
{
  double binsPerKeV = 5;

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  // c1->SetLogx();
  // c1->SetLogy();
  c1->SetLogz();

  TChain *waveTree = new TChain("waveTree");
  waveTree->Add("./data/waveSkimDS1_standardCut.root");

  GATDataSet ds(13550,13556);
  TChain *gatChain = ds.GetGatifiedChain();
  sprintf(theCut,"%s",calibCut);
  TH2D *h1 = new TH2D("h1","h1",(int)binsPerKeV*25,0.8,50,500,1000,15000);
	int cts = gatChain->Project("h1","d2wf0MHzTo50MHzPower:trapENFCal",theCut);
  cout << "Found " << cts << " entries.\n";
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("d2wf0MHzTo50MHzPower");

  // TVirtualPad *p1 = c1->GetPad(0);
  // p1->SetLeftMargin(0.2);
  // p1->SetLogz(1);
  // h1->GetXaxis()->SetTitleOffset(1.5); // move the damn label over

  // int h2binx1 = h2->GetXaxis()->FindBin(10);
  // int h2binx2 = h2->GetXaxis()->FindBin(200);
  // int h2biny1 = h2->GetYaxis()->FindBin(0.);
  // int h2biny2 = h2->GetYaxis()->FindBin(150);
  // double h2result = h2->Integral(h2binx1,h2binx2,h2biny1,h2biny2);
  //
  // printf("h1result: %.1f  h2result: %.1f  frac: %.2f\n",h1result,h2result,h2result/h1result);

	h1->Draw("COLZ");
	c1->Print("./plots/FFTPwrVsEnergy.C");
}

void CheckMCMCCut()
{
  // TFile *f = new TFile("./data/waveFitDS1_extendedCut.root");
  TFile *f = new TFile("./data/wave-standardCutDS1-mcmc.root");  // this is actually from October data
  TTree *waveTree = (TTree*)f->Get("waveTreeMCMC");

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);

  TH2D *h1 = new TH2D("h1","h1",250,0,50,100,0,100);
  double tot = waveTree->Project("h1","slowness:trapENFCal","");
  h1->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h1->GetYaxis()->SetTitle("Slowness");
  h1->Draw("COLZ");

  c1->Print("./plots/MCMCSlownessVsEnergy.pdf");

  TAxis *ax = h1->GetXaxis(); // trapENFCal
  TAxis *ay = h1->GetYaxis(); // slowness
  double sptot = h1->Integral(ax->FindBin(0.8),ax->FindBin(50),ay->FindBin(10),ay->FindBin(100));
  printf("slow pulses %.0f  tot %.0f  remaining %.0f\n",sptot,tot,tot-sptot);

  // -----------------

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add("./data/skimDS1_basicCuts.root");

  double ds1enrExp = 662;  // calculated in the above code
  double binsPerKeV = 2;
  double lo = 5, hi = 50;

  // DS-1 T/E Cut
  TH1D *h2 = new TH1D("h2","h2",(int)binsPerKeV*(hi-lo),lo,hi);
  sprintf(theCut,"%s && isEnr",standardCut);
  skimTree->Project("h2","trapENFCal",theCut);
  h2->GetXaxis()->SetTitle("Hit Energy (trapENFCal)");
  h2->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
  h2->GetYaxis()->SetTitleOffset(1.1);
  h2->SetLineColor(kBlue);
  h2->Scale(binsPerKeV/ds1enrExp);
  h2->Draw();

  // Slowness Cut
  TH1D *h3 = new TH1D("h3","h3",(int)binsPerKeV*(hi-lo),lo,hi);
  sprintf(theCut,"slowness < 10");
  waveTree->Project("h3","trapENFCal",theCut);
  h3->SetLineColor(kRed);
  h3->Scale(binsPerKeV/ds1enrExp);
  h3->Draw("same");

  TLegend* leg2 = new TLegend(0.3,0.7,0.87,0.92);
  leg2->AddEntry(h2,TString::Format("DS1 Enriched (%.0f kg-d) (with T/E) ",ds1enrExp),"l");
  leg2->AddEntry(h3,"DS1 Enriched (Slowness < 10)","l");
  leg2->Draw("SAME");
  c1->Update();
  c1->Print("./plots/SlownessVsTECut.pdf");

  // Clint, this computes the TOTAL counts between "lo" and "hi" keV, so the
	// result is in cts/kg-d. Then to get the average, divide by the number of bins.
	double tot1 = h2->Integral(h2->GetXaxis()->FindBin(lo),h2->GetXaxis()->FindBin(hi-1),"width");
	double tot2 = h3->Integral(h3->GetXaxis()->FindBin(lo),h3->GetXaxis()->FindBin(hi-1),"width");
	cout << "cts/kg-d                : T/E: " << tot1 << " slowness: " << tot2 << "  reduction: " << tot1/tot2 << endl;
	cout << "cts/kg-d-keV            : T/E: " << tot1/(hi-lo) << " slowness: " << tot2/(hi-lo) << "  reduction: " << tot1/tot2 << endl;
}