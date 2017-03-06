// Identify events in skim (or gat) data passing
// a given TCut, and append the corresponding waveforms
// from built data into a new file.
// C. Wiseman, 1/18/2017
//         v.2 3/06/2017

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "GATDataSet.hh"
#include "MGTWaveform.hh"

using namespace std;

void TCutSkimmer(string theCut, string inFile, string outFile);
void SkimWaveforms(string theCut, string inFile, string outFile);
void GATWaveforms(string theCut, string runFile, string outFile);

int main(int argc, char** argv)
{
  // Set cut and file I/O
  string theCut = "trapENFCal > 0 && trapENFCal < 200";
  theCut += " && gain==0 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0";
  theCut += " && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336";
  string runList = "./runs/test.txt";
  string inFile = "./data/skimDS1_basicCuts.root";
  string outFile = "waveSkimDS1_test.root";

  // Get user inputs
  bool s=1, g=0, t=0;
  vector<string> opt(argv, argv+argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-g") { s=0; g=1; }
    if (opt[i] == "-t") { s=0; t=1; }
    if (opt[i] == "-o") { outFile = opt[i+1] + "/" + outFile; }
    if (opt[i] == "-h") {
      cout << "Usage: ./wave-skim [-h (print this help)]\n"
         << "                     [-o [directory] (specify output location)]\n"
         << "                     [-g (directly scan gatified files)]\n"
         << "                     [-s (skim a skim file w/ a TCut)]\n";
    }
    return 1;
  }

  // Go running
  if (t) TCutSkimmer(theCut, inFile, outFile);
  if (s) SkimWaveforms(theCut, inFile, outFile);
  if (g) GATWaveforms(theCut, runList, outFile);
}


void TCutSkimmer(string theCut, string inFile, string outFile)
{
  sprintf(theCut,"%s",theCut.c_str());
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


void SkimWaveforms(string theCut, string inFile, string outFile)
{
  // From the skim tree, create a smaller "trimmed" tree
  // and write it to a new output file.
  // NOTE: The copied vectors are NOT resized to contain only entries passing cuts.

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add(inFile.c_str());
  size_t n = skimTree->Draw("run:iEvent",theCut.c_str(),"GOFF");
  if (n == 0) {
    cout << "No events found passing cuts.  Exiting ...\n";
    return;
  }
  skimTree->Draw(">>elist",theCut.c_str(), "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  skimTree->SetEntryList(elist);
  TFile *output = new TFile(outFile.c_str(),"RECREATE");
  TTree *cutTree = skimTree->CopyTree("");
  cutTree->Write("",TObject::kOverwrite);
  TNamed thisCut("theCut",theCut);	// save the cut used into the file.
  thisCut.Write();
  cout << Form("Using this cut:  \n%s  \nWrote %lli entries to the cut tree.\n",theCut.c_str(),cutTree->GetEntries());

  // Add waveforms to the cut tree, keeping only one run in memory at a time
  int run=0, iEvent=0;
  vector<double> *channel=0;
  cutTree->SetBranchAddress("run",&run);
  cutTree->SetBranchAddress("iEvent",&iEvent);
  cutTree->SetBranchAddress("channel",&channel);
  vector<MGTWaveform*> *waveVector=0;
  TBranch *waveBranch = cutTree->Branch("MGTWaveforms","vector<MGTWaveform*>",&waveVector,32000,0);
  GATDataSet *ds = new GATDataSet();
  string runPath = "";
  TChain *built = new TChain("MGTree");
  TChain *gat = new TChain("mjdTree");
  TTreeReader bReader(built);
  TTreeReader gReader(gat);
  TTreeReaderValue<TClonesArray> wfBranch(bReader,"fWaveforms");
  TTreeReaderArray<double> wfChan(gReader,"channel");
  int prevRun=0;
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
    bReader.SetEntry(iEvent);
    gReader.SetEntry(iEvent);
    waveVector->resize(0);

    int nWF = (*wfBranch).GetEntriesFast();
    int skimChans = channel->size();
    cout << Form("Entry %i -- run %i  iEvent %-10i(nWF skimChans) = (%i %i)\n" ,i,run,iEvent,nWF,skimChans);

    // Figure out which hits made it into the skim file (some are cut by data cleaning)
    int numPass = cutTree->Draw("channel","","GOFF",1,i);
    double *channelList = cutTree->GetV1();
    vector<double> chanVec;
    for (int j = 0; j < numPass; j++) chanVec.push_back(channelList[j]);

    // Fill the waveform branch
    for (int iWF = 0; iWF < nWF; iWF++) {
      if ( find(chanVec.begin(), chanVec.end(), wfChan[iWF]) != chanVec.end() ) {
        MGTWaveform *wave = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF)->Clone());
        waveVector->push_back(wave);
        cout << Form("   wave -> iWF %i  channel %.0f\n", iWF,wfChan[iWF]);
      }
    }
    waveBranch->Fill();


    // update progress and save
    if (i%500==0 && i!=0) {
      cout << i << " done.\n";
      cutTree->Write("", TObject::kOverwrite);
    }
    // save for next entry
    prevRun = run;
  }

  // Save and quit
  cutTree->Write("",TObject::kOverwrite);
  output->Close();
  delete ds;
  cout << "Wrote file: " << outFile << endl;
}

void GATWaveforms(string runFile, string theCut, string outFile)
{
  GATDataSet *ds = new GATDataSet();
  ifstream runFile(runFile);
  int rundummy;
  while (runFile >> rundummy) ds->AddRunNumber(rundummy);
  runFile.close();
  TChain *gatChain = ds->GetGatifiedChain(false);
  TChain *builtChain = ds->GetBuiltChain(false);

  // Apply the cut to the chain and save tree to a new file
  TFile *output = new TFile(outFile.c_str(),"RECREATE");
  gatChain->Draw(">>elist",theCut.c_str(), "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  gatChain->SetEntryList(elist);
  TTree *cutTree = gatChain->CopyTree("");
  cutTree->Write("",TObject::kOverwrite);
  TNamed thisCut("theCut",theCut);	// save the cut used into the file.
  thisCut.Write();
  cout << Form("Using this cut:  \n%s  \nWrote %lli entries to the cut tree.\n",theCut.c_str(),cutTree->GetEntries());

  // Now add a waveform branch
  cout << "Adding waveforms to the cut tree ...\n";
  vector<MGTWaveform*> *waveVector=0;
  TBranch *waveBranch = cutTree->Branch("MGTWaveforms","vector<MGTWaveform*>",&waveVector,32000,0);
  TTreeReader bReader(builtChain);
  TTreeReaderValue<TClonesArray> wfBranch(bReader,"fWaveforms");
  TTreeReader gReader(gatChain);
  TTreeReaderArray<double> wfChan(gReader,"channel");
  TTreeReaderValue<double> gRun(gReader,"run");
  for (int i = 0; i < cutTree->GetEntries(); i++)
  {
    cutTree->GetEntry(i);
    bReader.SetEntry(i);
    gReader.SetEntry(i);

    // Fill the waveform branch
    waveVector->resize(0);
    int nWF = (*wfBranch).GetEntriesFast();
    int nCh = wfChan.GetSize();
    if (nWF != nCh) {
      // Hopefully this will never happen because it will screw up the matching of the vector indexes
      cout << Form("Warning! At entry %i,  num WF's (%i) != num channels (%i).  Skipping this event to avoid a segfault.\n",i,nWF,nCh);
      // waveBranch->Fill();
      continue;
    }
    for (int iWF = 0; iWF < nWF; iWF++) {
      MGTWaveform *wave = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF)->Clone());
      waveVector->push_back(wave);
      cout << Form("   wave -> iWF %i  channel %.0f\n", iWF,wfChan[iWF]);
    }
    waveBranch->Fill();
    cout << Form("Entry %i -- run %-10.0f (nWF nCh vec) = (%i %i %lu)\n",i,*gRun,nWF,nCh,waveVector->size());

    // update progress and save
    if (i%500==0 && i!=0) {
      cout << i << " done.\n";
      cutTree->Write("", TObject::kOverwrite);
    }
  }
  // Save and quit
  cutTree->Write("",TObject::kOverwrite);
  output->Close();
  delete ds;
  cout << "Wrote file: " << outFile << endl;
}