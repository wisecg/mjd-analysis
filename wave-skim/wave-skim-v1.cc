// This identifies events in skim or gatified data
// passing a given TCut, and saves the corresponding
// waveforms from built data into a new file.
//
// Typical memory usage on a big skim: virt 1238m, res 448m
//
// Clint Wiseman, USC
// 12/16/2016

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
#include "TTreeReaderArray.h"
#include "GATDataSet.hh"
#include "GATWaveformBrowser.hh"
#include "MGTWaveform.hh"
#include "MGTEventList.hh"

using namespace std;

void SkimWaveforms(string theCut, string skimLoc, string outFile);
void GATWaveforms(string runListFile, string theCut, string outFile);

int main(int argc, char** argv)
{
  string outputDir = "./";
  bool skimFiles=true, gatFiles=false;
  vector<string> opt(argc);
  for (int i=0; i<argc; i++) opt[i]=argv[i];
  if (find(opt.begin(), opt.end(), "-g") != opt.end()) {
    gatFiles=true;
    skimFiles=false;
  }
  if (find(opt.begin(), opt.end(), "-o") != opt.end()) {
    int pos = find(opt.begin(), opt.end(), "-o") - opt.begin();
    outputDir = opt[pos+1]+"/";
  }
  if (find(opt.begin(), opt.end(), "-h") != opt.end()) {
    cout << "Usage: ./wave-skim [-h (print this help)]\n"
         << "                   [-o [directory] (optional - specify output location)]\n"
         << "                   [-g (optional - directly scan gatified files)]\n";
    return 1;
  }
  string skimLoc = "./data/skimDS1_basicCuts.root";
  // string skimLoc = "./data/skimDS1_standardCut.root";
  // string skimLoc = "./data/skimDS1_calibCut.root";
  // string outFile = "./data/waveSkimDS1_standardCut_kvorrT.root";
  string outFile = "./data/waveSkimDS1_test.root";
  // string outFile = "./data/wave-m1Cal.root";
  string runList = "./runs/test.txt";

  // string theCut = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run!=9648 && run!=10663 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004 && run!=13205 && run!=13306 && run!=13307 && run!=13330 && run!=13400 && trapETailMin < 0 && mH==1 && !isLNFill1 && isGood && !muVeto && !wfDCBits && run!=13099";
  string theCut = "channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH>1 && isGood && !isLNFill1 && !wfDCBits && Entry$ < 5000";
  // string theCut = "trapENFCal > 1591 && trapENFCal < 1593 && !wfDCBits && channel%2==0";
  // string theCut = "trapENFCal > 0.8 && trapENFCal < 100 && mH > 1 && channel%2==0 && Entry$ < 2000";

  // ----------------------------------------------------------------------------
  // 1. Skim file waveform skimmer (default)
  if (skimFiles) SkimWaveforms(theCut, skimLoc, outFile);

  // 2. GAT-based waveform skimmer (-g option)
  // Take a list of runs from a text file, and apply a cut based on gatified parameters.
  if (gatFiles) GATWaveforms(runList, theCut, outFile);
  // ----------------------------------------------------------------------------
}

void SkimWaveforms(string theCut, string skimLoc, string outFile)
{
  // Make an event list from the skim tree
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add(skimLoc.c_str());
  size_t n = skimTree->Draw("run:iEvent:iHit:channel",theCut.c_str(),"GOFF");
  if (n == 0) {
    cout << "No events found passing cuts.  Exiting ...\n";
    return;
  }
  double* runList = skimTree->GetV1();
  double* iEventList = skimTree->GetV2();
  double* iHitList = skimTree->GetV3();
  double* channelList = skimTree->GetV4();
  vector<double> runVec, iEventVec, iHitVec;
  for (size_t i = 0; i < n; i++)
  {
    cout << Form("skim -> iWF %lu  run %.0f  iEvent %.0f  iHit %.0f  chan %.0f\n", i,runList[i],iEventList[i],iHitList[i],channelList[i]);
    runVec.push_back(runList[i]);
    iEventVec.push_back(iEventList[i]);
    iHitVec.push_back(iHitList[i]);
  }
  int nWF = runVec.size(), iWF = 0;

  // Create a smaller "cut" tree and write it to a new output file
  skimTree->Draw(">>elist",theCut.c_str(), "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  skimTree->SetEntryList(elist);
  TFile *output = new TFile(outFile.c_str(),"RECREATE");
  TTree *cutTree = skimTree->CopyTree("");
  cutTree->Write("",TObject::kOverwrite);
  TNamed thisCut("theCut",theCut);	// save the cut used into the file.
  thisCut.Write();
  cout << "Wrote " << cutTree->GetEntries() << " entries to the cut tree.\n";

  // Add waveforms to the cut tree, keeping only one run in memory at a time
  int run=0, iEvent=0;
  cutTree->SetBranchAddress("run",&run);
  cutTree->SetBranchAddress("iEvent",&iEvent);
  vector<MGTWaveform*> *waveVector=0;
  vector<int> *waveIdx=0;
  TBranch *waveBranch = cutTree->Branch("MGTWaveforms","vector<MGTWaveform*>",&waveVector,32000,0);
  TBranch *waveIndexBranch = cutTree->Branch("waveformIdx","vector<int>",&waveIdx);
  GATDataSet *ds = new GATDataSet();
  string runPath = "";
  TChain *built = new TChain("MGTree");
  TChain *gat = new TChain("mjdTree");
  TTreeReader bReader(built);
  TTreeReader gReader(gat);
  TTreeReaderValue<TClonesArray> wfBranch(bReader,"fWaveforms");
  TTreeReaderArray<double> wfChan(gReader,"channel");
  int prevRun=0;
  int exceptionCount=0;
  cout << "Adding waveforms to the cut tree ...\n";
  for (int i = 0; i < cutTree->GetEntries(); i++)
  {
    cutTree->GetEntry(i);
    if (run != prevRun) {
      built->Reset();
      gat->Reset();
      runPath = ds->GetPathToRun(run,GATDataSet::kBuilt);
      built->Add(runPath.c_str());
      bReader.SetTree(built);
      runPath = ds->GetPathToRun(run,GATDataSet::kGatified);
      gat->Add(runPath.c_str());
      gReader.SetTree(gat);
    }
    waveVector->resize(0);
    iWF = find(runVec.begin(), runVec.end(), run) - runVec.begin();
    while (iEventVec[iWF] != iEvent) iWF++;
    if ((long)iWF >= nWF) {
      cout << "Hit the end of the WF list, breaking.\n";
      break;
    }
    while (iEventVec[iWF]==iEvent) {
      size_t iEvt = (size_t)iEventVec[iWF];
      size_t jEvt = (size_t)iHitVec[iWF];
      bReader.SetEntry(iEvt);
      gReader.SetEntry(iEvt);
      try {
        MGTWaveform *wf = dynamic_cast<MGTWaveform*>((*wfBranch).At(jEvt)->Clone());
        waveVector->push_back(wf);
        waveIdx->push_back(jEvt);
        cout << Form("wave -> iWF %i  entry %i/%i  run %i  iEvent(i) %lu  iHit(j) %lu  chan %.0f\n", iWF,i,(int)cutTree->GetEntriesFast(),run,iEvt,jEvt,wfChan[jEvt]);
        iWF++;
      }
      catch (...) { exceptionCount++; }
    }
    waveBranch->Fill();
    waveIndexBranch->Fill();
    if (i%500==0 && i!=0) {
      cout << i << " done.\n";
      cutTree->Write("", TObject::kOverwrite);
    }
    prevRun = run;
  }
  cutTree->Write("",TObject::kOverwrite);
  output->Close();
  delete ds;
  cout << "Wrote file: " << outFile << endl;

  if (exceptionCount > 0)
    cout << Form("WARNING: Missed %i of %i waveforms.\n",exceptionCount,nWF);
}


void GATWaveforms(string runListFile, string theCut, string outFile)
{
  // Get input chains
  GATDataSet *ds = new GATDataSet();
  ifstream runFile(runListFile);
  int rundummy;
  while (runFile >> rundummy) ds->AddRunNumber(rundummy);
  runFile.close();
  TChain *gatChain = ds->GetGatifiedChain(false);

  // Apply the cut to the chain and save tree to a new file
  TFile *output = new TFile(outFile.c_str(),"RECREATE");
  gatChain->Draw(">>elist",theCut.c_str(), "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  gatChain->SetEntryList(elist);
  TTree *cutTree = gatChain->CopyTree("");
  cutTree->Write("",TObject::kOverwrite);
  TNamed thisCut("theCut",theCut);	// save the cut used into the file.
  thisCut.Write();
  cout << "Wrote " << cutTree->GetEntries() << " entries to the cut tree.\n";

  // Add a waveform branch to the cut tree
  cout << "Adding waveforms to cut tree ...\n";
  vector<MGTWaveform*> waveVector;
  TBranch *waveBranch = cutTree->Branch("MGTWaveforms","vector<MGTWaveform*>",&waveVector,32000,0);
  TChain *builtChain = ds->GetBuiltChain(false);
  TTreeReader fReader(builtChain);
  TTreeReaderValue<TClonesArray> wfBranch(fReader,"fWaveforms");
  TTreeReaderValue<UInt_t> wfID(fReader,"fID");
  size_t nWF = gatChain->Draw("Entry$:Iteration$",theCut.c_str(),"GOFF");
  if (nWF == 0) {
    cout << "No events found passing cuts.  Exiting ...\n";
    output->Close();
    return;
  }
  double* iEventList = gatChain->GetV1();
  double* iHitList = gatChain->GetV2();
  double prevEvent = iEventList[0];
  int nEvts=0;
  for (size_t i = 0; i < nWF; i++)
  {
    fReader.SetEntry(iEventList[i]);
    MGTWaveform *wf = dynamic_cast<MGTWaveform*>((*wfBranch).At(iHitList[i])->Clone());
    waveVector.push_back(wf);
    cout << Form("i %lu  wvsize %lu\n",i,waveVector.size());
    if (iEventList[i] != prevEvent){
      waveBranch->Fill();
      cout << Form("pushed back event %i (gat entry %.0f).  wvsize: %lu\n",nEvts,iEventList[i],waveVector.size());
      waveVector.resize(0);
      nEvts++;
    }
    if (i%1000==0 && i!=0){
      cout << i << " done.\n";
      cutTree->Write("", TObject::kOverwrite);
    }
    prevEvent = iEventList[i];
  }
  waveBranch->Fill();
  cout << Form("pushed back event %i (gat entry %.0f).  wvsize: %lu\n",nEvts,iEventList[nWF-1],waveVector.size());
  cutTree->Write("",TObject::kOverwrite);
  cout << "Added " << waveBranch->GetEntries() << " entries to the waveform branch.\n";
  output->Close();
  cout << "Wrote file: " << outFile << endl;
  delete ds;
}