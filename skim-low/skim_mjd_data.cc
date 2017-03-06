#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TCut.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TEntryList.h"
#include "TROOT.h"

#include "GATDataSet.hh"
#include "MJTChannelMap.hh"
#include "TClonesArray.h"
#include "MGTEvent.hh"
#include "GATUtils.hh"

using namespace std;
using namespace CLHEP;

void LoadDataSet(GATDataSet& ds, int dsNumber, size_t iRunSeq);
void LoadRun(GATDataSet& ds, size_t iRunSeq);
void LoadActiveMasses(map<int,double>& activeMassForDetID_g, int dsNumber);
double GetAOverENorm(int channel, double a50, double a100, double a200, double trapENF, double trapENFCal, int dsNumber);
double GetAOverENorm85(int channel, double a50, double a100, double a200, double trapENF, double trapENFCal, int dsNumber);
double GetAvsE(int channel, double TSCurrent50nsMax, double TSCurrent100nsMax, double TSCurrent200nsMax, double trapENF, double trapENFCal, int dsNumber);
void LoadLNFillTimes(vector<double>& lnFillTimes, int dsNumber);
double GetDCRraw(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);
double GetDCR85(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);
double GetDCR90(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);
double GetDCR95(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);
double GetDCR98(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);
double GetDCR99(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber);

int main(int argc, const char** argv)
{
  if(argc < 3 || argc > 7) {
    cout << "For raw DCR add flag -r " << endl;
    cout << "For minimal skim file add flag -m " << endl;
    cout << "For extensive skim file (multiple DCR and aenorm) add flag -e " << endl;
    cout << "Usage for single file: " << argv[0] << " -f [runNum] (output path)" << endl;
    cout << "Usage for data sets: " << argv[0] << " [dataset number] [runseq] (output path)" << endl;
    return 1;
  }

  GATDataSet ds;
  string outputPath = "";
  int noDS = -999999;
  int dsNumber = noDS;
  int runSeq;
  bool singleFile = false;
  bool writeRawDCR = false;
  bool smallOutput = false;
  bool extendedOutput = false;
  TString flags ="";
  int i = 1;
  while(!TString(argv[i]).IsDigit()){
    flags += TString(argv[i]);
    i++;
  }
  if(flags.Contains("f")){
    singleFile = true;
    runSeq = atoi(argv[i]);
    i++;
    if(runSeq >= 2335 && runSeq <= 8183) dsNumber = 0;
      else if(runSeq >= 8722 && runSeq < 600000) dsNumber = 1;
      else {
        cout << "Error: I don't know what dataset run " << runSeq << " is from." << endl;
        return 1;
      }
    cout << "loading run " << runSeq << endl;
    LoadRun(ds, runSeq);
  }
  else{
    dsNumber = atoi(argv[i]);
    i++;
    runSeq = atoi(argv[i]);
    i++;
    cout << "loading dataset " << dsNumber << " run sequence " << runSeq << endl;
    LoadDataSet(ds, dsNumber, runSeq);
  }

  if(flags.Contains("r")){
    writeRawDCR = true;
    cout<<"Raw DCR option selected."<<endl;
  }
  if(flags.Contains("m")){
    smallOutput = true;
    cout<<"Minimal skim file option selected."<<endl;
  }
  if(flags.Contains("e")){
    extendedOutput = true;
    cout<<"Extended skim file option selected."<<endl;
  }

  if(argc > i) outputPath += string(argv[i]);

  // set up dataset
  TChain* gatChain = ds.GetGatifiedChain(false);
  TTreeReader gatReader(gatChain);

  // set up input chain value readers

  // run level variables and indices
  TTreeReaderValue<unsigned int> gatrevIn(gatReader, "gatrev");
  TTreeReaderValue<double> runIn(gatReader, "run");
  double runSave = -1;

  // ID variables
  TTreeReaderValue< vector<double> > channelIn(gatReader, "channel");
  TTreeReaderValue< vector<int> > detIDIn(gatReader, "detID");
  TTreeReaderValue< vector<int> > posIn(gatReader, "P");
  TTreeReaderValue< vector<int> > detIn(gatReader, "D");
  TTreeReaderValue< vector<int> > mageIDIn(gatReader, "mageID");
  TTreeReaderValue< vector<string> > detNameIn(gatReader, "detName");
  TTreeReaderValue< vector<bool> > isEnrIn(gatReader, "isEnr");
  TTreeReaderValue< vector<bool> > isNatIn(gatReader, "isNat");

  // time variables
  TTreeReaderValue<double> startTimeIn(gatReader, "startTime");
  TTreeReaderValue<double> stopTimeIn(gatReader, "stopTime");
  TTreeReaderValue< vector<double> > timestampIn(gatReader, "timestamp");
  TTreeReaderValue< vector<double> > timeMTIn(gatReader, "timeMT");
  TTreeReaderValue< vector<int> > dateMTIn(gatReader, "dateMT");

  // energy variables
  TTreeReaderValue< vector<double> > trapENFIn(gatReader, "trapENF");
  TTreeReaderValue< vector<double> > trapENFCalIn(gatReader, "trapENFCal");
  TTreeReaderValue< vector<double> > trapECalIn(gatReader, "trapECal");
  TTreeReaderValue< vector<double> > trap4usMaxIn(gatReader, "trap4usMax");
  TTreeReaderValue< vector<double> > energyIn(gatReader, "energy");

  // analysis cut variables
  TTreeReaderValue< vector<double> > blrwfFMR1In(gatReader, "blrwfFMR1");
  TTreeReaderValue< vector<double> > blrwfFMR50In(gatReader, "blrwfFMR50");
  TTreeReaderValue< vector<double> > tsCurrent50nsMaxIn(gatReader, "TSCurrent50nsMax");
  TTreeReaderValue< vector<double> > tsCurrent100nsMaxIn(gatReader, "TSCurrent100nsMax");
  TTreeReaderValue< vector<double> > tsCurrent200nsMaxIn(gatReader, "TSCurrent200nsMax");
  TTreeReaderValue< vector<double> > trirt100nsft10nsMaxIn(gatReader, "trirt100nsft10nsMax");
  TTreeReaderValue< vector<double> > toeIn(gatReader, "toe");
  TTreeReaderValue< vector<double> > dcrSlopeIn(gatReader, "nlcblrwfSlope");

  // data cleaning variables
  TTreeReaderValue<unsigned int> eventDC1BitsIn(gatReader, "EventDC1Bits");
  // pulser tag channel seems to miss a lot of pulsers! So let's just use
  // Pinghan's (only).
  //const int kDoublePulserMask = (0x1 << 1) + 0x1; // pinghan pulsers + pulser tag channels
  const int kPinghanPulserMask = 0x1 << 1; // pinghan pulsers
  TTreeReaderValue< vector<unsigned int> > wfDCBitsIn(gatReader, "wfDCBits");
  TTreeReaderValue< vector<double> > trapETailMinIn(gatReader, "trapETailMin");

  // removed by clint, was causing segfaults (probably wasn't present in all DS1 data at the time)
  // TTreeReaderValue< vector<double> > nRisingXIn(gatReader, "fastTrapNLCWFsnRisingX");

  // muon variables
  vector<int> muRuns;
  vector<double> muRunTStarts_s;
  vector<double> muTimes_s;
  vector<int> muTypes;
  vector<bool> badScalers;
  cout << "Loading muon data..." << endl;
  string muFileName = TString::Format("MuonList_DS%d_aug.txt", dsNumber).Data();
  ifstream muFile(muFileName.c_str());
  if(!muFile.good()) {
    cout << "Couldn't open " << muFileName << endl;
    return 0;
  }
  int rundummy, typedummy;
  double rstartdummy, timedummy;
  bool baddummy;
  while(muFile >> rundummy >> rstartdummy >> timedummy >> typedummy >> baddummy) {
    muRuns.push_back(rundummy);
    muRunTStarts_s.push_back(rstartdummy);
    muTimes_s.push_back(timedummy);
    muTypes.push_back(typedummy);
    badScalers.push_back(baddummy);
  }
  muFile.close();
  size_t iMu = 0;
  size_t nMu = muTimes_s.size();
  if(nMu == 0) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }

  // set up output file and tree
  string filename = TString::Format("skimDS%d_", dsNumber).Data();
  char runStr[10];
  sprintf(runStr, "%d", runSeq);
  if(singleFile) filename += "run";
  filename += runStr;
  if(writeRawDCR) filename += "_rawDCR";
  if(smallOutput) filename += "_small";
  if(extendedOutput) filename += "_ext";
  filename += ".root";
  if(outputPath != "") filename = outputPath + "/" + filename;
  TFile *fOut = TFile::Open(filename.c_str(), "recreate");
  TTree* skimTree = new TTree("skimTree", "skimTree");

  // set up output chain branches

  // run level variables and indices
  unsigned int skimgatrev = strtol(GATUtils::GetGATRevision(), NULL, 16);
  skimTree->Branch("skimgatrev", &skimgatrev, "skimgatrev/i");
  unsigned int gatrev = 0;
  skimTree->Branch("gatrev", &gatrev, "gatrev/i");
  int run = 0;
  skimTree->Branch("run", &run, "run/I");
  int iEvent = 0;
  skimTree->Branch("iEvent", &iEvent, "iEvent/I");
  vector<int> iHit;
  skimTree->Branch("iHit", &iHit);

  // ID variables
  vector<int> channel;
  skimTree->Branch("channel", &channel);
  vector<int> pos;
  skimTree->Branch("P", &pos);
  vector<int> det;
  skimTree->Branch("D", &det);
  vector<int> gain;
  skimTree->Branch("gain", &gain);
  vector<int> mageID;
  skimTree->Branch("mageID", &mageID);
  vector<int> detID;
  skimTree->Branch("detID", &detID);
  vector<string> detName;
  skimTree->Branch("detName", &detName);
  vector<bool> isEnr;
  skimTree->Branch("isEnr", &isEnr);
  vector<bool> isNat;
  skimTree->Branch("isNat", &isNat);
  map<int, double> actM4Det_g;
  LoadActiveMasses(actM4Det_g, dsNumber);
  vector<double> mAct_g;
  skimTree->Branch("mAct_g", &mAct_g);
  vector<bool> isGood;
  skimTree->Branch("isGood", &isGood);

  // total mass variables
  double mAct_Total_kg = 0;
  double mAct_enr_kg = 0;
  double mAct_nat_kg = 0;
  if(dsNumber == 0) {
    mAct_Total_kg = 14.60;
    mAct_enr_kg = 10.69;
    mAct_nat_kg = 3.91;
  }
  else if(dsNumber == 1) {
    mAct_Total_kg = 12.43;
    mAct_enr_kg = 11.31;
    mAct_nat_kg = 1.12;
  }
  skimTree->Branch("mAct_Total_kg", &mAct_Total_kg, "mAct_Total_kg/D");
  skimTree->Branch("mAct_enr_kg", &mAct_enr_kg, "mAct_enr_kg/D");
  skimTree->Branch("mAct_nat_kg", &mAct_nat_kg, "mAct_nat_kg/D");

  // time variables
  double runTime_s = 0;
  double startTime = 0;
  double startTime0 = 0;
  double stopTime = 0;
  if(dsNumber == 0) {
    runTime_s = 4121110;
    startTime0 = 1435687000; // start time of run 2580
  }
  else if(dsNumber == 1) {
    runTime_s = 5278680; //was 4728790; before blinding
    startTime0 = 1452643100; // start time of run 9422
  }
  skimTree->Branch("startTime", &startTime, "startTime/D");
  skimTree->Branch("startTime0", &startTime0, "startTime0/D");
  skimTree->Branch("runTime_s", &runTime_s, "runTime_s/D");
  skimTree->Branch("stopTime", &stopTime, "stopTime/D");
  vector<double> tloc_s;
  skimTree->Branch("tloc_s", &tloc_s);
  vector<double> time_s;
  skimTree->Branch("time_s", &time_s);
  vector<double> timeMT;
  skimTree->Branch("timeMT", &timeMT);
  vector<int> dateMT;
  skimTree->Branch("dateMT", &dateMT);

  // energy variables
  vector<double> trapECal;
  vector<double> trap4usMax;
  vector<double> onBoardE;
  if(!smallOutput){
    skimTree->Branch("trapECal", &trapECal);
    skimTree->Branch("trap4usMax", &trap4usMax);
    skimTree->Branch("onBoardE", &onBoardE);
  }
  vector<double> trapENFCal;
  skimTree->Branch("trapENFCal", &trapENFCal);
  double sumEH = 0;
  skimTree->Branch("sumEH", &sumEH, "sumEH/D");
  double sumEL = 0;
  skimTree->Branch("sumEL", &sumEL, "sumEL/D");
  double sumEHClean = 0;
  skimTree->Branch("sumEHClean", &sumEHClean, "sumEHClean/D");
  double sumELClean = 0;
  skimTree->Branch("sumELClean", &sumELClean, "sumELClean/D");

  // analysis cut variables
  int mH = 0;
  skimTree->Branch("mH", &mH, "mH/I");
  int mL = 0;
  skimTree->Branch("mL", &mL, "mL/I");
  int mHClean = 0;
  skimTree->Branch("mHClean", &mHClean, "mHClean/I");
  int mLClean = 0;
  skimTree->Branch("mLClean", &mLClean, "mLClean/I");
  vector<double> aenorm;
  skimTree->Branch("aenorm", &aenorm);
  vector<double> avse;
  skimTree->Branch("avse", &avse);
  vector<double> t150;
  vector<double> kvorrT;
  vector<double> toe;
  vector<double> aenorm85;
  if(!smallOutput){
    skimTree->Branch("t150", &t150);
    skimTree->Branch("kvorrT", &kvorrT);
    skimTree->Branch("toe", &toe);
    if(extendedOutput) skimTree->Branch("aenorm85", &aenorm85);
  }
  vector<double> rawDCR;
  vector<double> dcrSlope85;
  vector<double> dcrSlope90;
  vector<double> dcrSlope95;
  vector<double> dcrSlope98;
  vector<double> dcrSlope99;
  if(writeRawDCR) skimTree->Branch("rawDCR", &rawDCR);
  else{
    if(!smallOutput && extendedOutput) {
      skimTree->Branch("dcrSlope85", &dcrSlope85);
      skimTree->Branch("dcrSlope95", &dcrSlope95);
      skimTree->Branch("dcrSlope98", &dcrSlope98);
      skimTree->Branch("dcrSlope99", &dcrSlope99);
    }
    skimTree->Branch("dcrSlope90", &dcrSlope90);
  }
 // data cleaning variables
  unsigned int eventDC1Bits = 0;
  skimTree->Branch("EventDC1Bits", &eventDC1Bits, "eventDC1Bits/i");
  vector<unsigned int> wfDCBits;
  skimTree->Branch("wfDCBits", &wfDCBits);
  vector<double> lnFillTimes;
  LoadLNFillTimes(lnFillTimes, dsNumber);
  vector<bool> isLNFill;
  skimTree->Branch("isLNFill", &isLNFill);
  vector<int> nX;
  skimTree->Branch("nX", &nX);
  vector<double> trapETailMin;
  if(!smallOutput) skimTree->Branch("trapETailMin", &trapETailMin);

  // veto variables
  vector<double> dtmu_s;
  skimTree->Branch("dtmu_s", &dtmu_s);
  vector<int> muType;
  vector<bool> badScaler;
  if(!smallOutput) {
    skimTree->Branch("muType", &muType);
    skimTree->Branch("badScaler", &badScaler);
  }
  vector<bool> muVeto;
  skimTree->Branch("muVeto", &muVeto);
  map<int,bool> detIDIsBad;
  if(dsNumber == 0) {
    detIDIsBad[28474] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28480] = true;
    detIDIsBad[1426980] = true;
    detIDIsBad[1426620] = true;
    detIDIsBad[1425370] = true;
  }
  if(dsNumber == 1) {
    detIDIsBad[1426981] = true;
    detIDIsBad[1426622] = true;
    detIDIsBad[28455] = true;
    detIDIsBad[28470] = true;
    detIDIsBad[28463] = true;
    detIDIsBad[28465] = true;
    detIDIsBad[28469] = true;
    detIDIsBad[28477] = true;
    detIDIsBad[1425751] = true;
    detIDIsBad[1425731] = true;
    detIDIsBad[1426611] = true;
  }

  map<int,bool> detIDIsVetoOnly;
  if(dsNumber == 0) {
    detIDIsVetoOnly[1425381] = true;
    detIDIsVetoOnly[1425742] = true;
  }
  if(dsNumber == 1) {
    detIDIsVetoOnly[28480] = true;
    detIDIsVetoOnly[1426621] = true;
  }

  // start loop over all events
  while(gatReader.Next()) {

    // stuff to do on run boundaries
    if(runSave != *runIn) {
      runSave = *runIn;
      cout << "Processing run " << *runIn << ", "
           << skimTree->GetEntries() << " entries saved so far"
           << endl;
      skimTree->Write("", TObject::kOverwrite);
    }

    // Skip this event if it is a pulser event as identified by Pinghan
    if(*eventDC1BitsIn & kPinghanPulserMask) continue;
    iEvent = gatChain->GetTree()->GetReadEntry();

    // copy the event-level info to the output fields
    gatrev = *gatrevIn;
    run = int(*runIn);
    startTime = *startTimeIn;
    stopTime = *stopTimeIn;
    eventDC1Bits = *eventDC1BitsIn;
    sumEH = 0;
    sumEL = 0;
    mH = 0;
    mL = 0;
    sumEHClean = 0;
    sumELClean = 0;
    mHClean = 0;
    mLClean = 0;

    // clear all hit-level info fields
    iHit.resize(0);
    if(!smallOutput){
      trapECal.resize(0);
      trap4usMax.resize(0);
      onBoardE.resize(0);
      t150.resize(0);
      kvorrT.resize(0);
      aenorm85.resize(0);
      toe.resize(0);
      trapETailMin.resize(0);
      muType.resize(0);
      badScaler.resize(0);
    }
    dtmu_s.resize(0);
    trapENFCal.resize(0);
    channel.resize(0);
    tloc_s.resize(0);
    time_s.resize(0);
    timeMT.resize(0);
    dateMT.resize(0);
    pos.resize(0);
    det.resize(0);
    gain.resize(0);
    mageID.resize(0);
    detID.resize(0);
    detName.resize(0);
    isEnr.resize(0);
    isNat.resize(0);
    mAct_g.resize(0);
    isGood.resize(0);
    wfDCBits.resize(0);
    aenorm.resize(0);
    avse.resize(0);
    if(writeRawDCR) rawDCR.resize(0);
    else {
      if(!smallOutput){
        dcrSlope85.resize(0);
        dcrSlope95.resize(0);
        dcrSlope98.resize(0);
      }
      dcrSlope90.resize(0);
      dcrSlope99.resize(0);
    }
    isLNFill.resize(0);
    nX.resize(0);
    muVeto.resize(0);
    // loop over hits
    bool skipMe = false;
    size_t nHits = trapENFCalIn->size();
    for(size_t i=0; i<nHits; i++) {

      // skip all hits with E_H < 2 keV, E_L < 10 keV in -both- trapE and trapENF
      // For small skim files, skip all hits with E_H and E_L < 200 keV in trapE and trapENF
      double hitENFCal = (*trapENFCalIn)[i];
      double hitENF = (*trapENFIn)[i];
      double hitEMax = (*trapECalIn)[i];
      double hit4usMax = (*trap4usMaxIn)[i];
      int hitCh = (*channelIn)[i];

	  // mod by clint - play limbo with the energy floor
	  //if(!smallOutput && hitCh%2 == 0 && hitENFCal <  2. && hitEMax < 2.) continue;
	  if(!smallOutput && hitCh%2 == 0 && hitENFCal <  0.8 && hitEMax < 0.8) continue;
	//   if(!smallOutput && hitCh%2 == 0 && hitENFCal <  0.1 && hitEMax < 0.1) continue;

	  if(!smallOutput && hitCh%2 == 1 && hitENFCal < 10. && hitEMax < 10.) continue;
      if(smallOutput && hitCh%2 == 0 && hitENFCal <  200. && hitEMax < 200.) continue;
      if(smallOutput && hitCh%2 == 1 && hitENFCal < 200. && hitEMax < 200.) continue;
      // skip hits from totally "bad" detectors (not biased, etc), or from
      // use-for-veto-only detectors if E < 10 keV
      int hitDetID = (*detIDIn)[i];
      if(detIDIsBad[hitDetID] || (detIDIsVetoOnly[hitDetID] && hitEMax < 10.)) continue;
      // copy over hit info
      iHit.push_back(i);
      trapENFCal.push_back(hitENFCal);
      channel.push_back(hitCh);
      double hitT_s = (*timestampIn)[i]*1.e-8;
      tloc_s.push_back(hitT_s);
      time_s.push_back( (startTime - startTime0) + hitT_s );
      timeMT.push_back((*timeMTIn)[i]);
      dateMT.push_back((*dateMTIn)[i]);
      pos.push_back((*posIn)[i]);
      det.push_back((*detIn)[i]);
      gain.push_back(hitCh % 2);
      double a50 = (*tsCurrent50nsMaxIn)[i];
      double a100 = (*tsCurrent100nsMaxIn)[i];
      double a200 = (*tsCurrent200nsMaxIn)[i];
      aenorm.push_back(GetAOverENorm(hitCh, a50, a100, a200, hitENF, hitENFCal, dsNumber));
      avse.push_back(GetAvsE(hitCh, a50, a100, a200, hitENF, hitENFCal, dsNumber));
      if(writeRawDCR) rawDCR.push_back(GetDCRraw(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
      else {
        dcrSlope90.push_back(GetDCR90(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
        dcrSlope99.push_back(GetDCR99(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
        if(!smallOutput){
          dcrSlope85.push_back(GetDCR85(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
          dcrSlope95.push_back(GetDCR95(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
          dcrSlope98.push_back(GetDCR98(hitCh, (*dcrSlopeIn)[i], hit4usMax, dsNumber));
        }
      }

      mageID.push_back((*mageIDIn)[i]);
      detID.push_back(hitDetID);
      detName.push_back((*detNameIn)[i]);
      isEnr.push_back((*detNameIn)[i][0] == 'P');
      isNat.push_back((*detNameIn)[i][0] == 'B');
      mAct_g.push_back(actM4Det_g[hitDetID]);
      isGood.push_back(!detIDIsVetoOnly[hitDetID]);
      wfDCBits.push_back((*wfDCBitsIn)[i]);
      // nX.push_back((*nRisingXIn)[i]); // removed by clint - see above note
      if(!smallOutput){
        trapECal.push_back(hitEMax);
        trap4usMax.push_back(hit4usMax);
        onBoardE.push_back((*energyIn)[i]);
        aenorm85.push_back(GetAOverENorm85(hitCh, a50, a100, a200, hitENF, hitENFCal, dsNumber));
        double t1 = (*blrwfFMR1In)[i];
        double t50 = (*blrwfFMR50In)[i];
        t150.push_back(t50-t1);
        kvorrT.push_back((*trirt100nsft10nsMaxIn)[i]);
        toe.push_back((*toeIn)[i] / hitEMax);
        trapETailMin.push_back((*trapETailMinIn)[i]);
      }

      // sum energies and multiplicities
      if(hitCh%2 == 0) {
        mH++;
        if(!detIDIsVetoOnly[hitDetID]) sumEH += hitENFCal;
      }
      else {
        mL++;
        if(!detIDIsVetoOnly[hitDetID]) sumEL += hitENFCal;
      }
      if(hitCh%2 == 0 && ~((~0x010) & wfDCBits[wfDCBits.size()-1])) {
        mHClean++;
        if(!detIDIsVetoOnly[hitDetID]) sumEHClean += hitENFCal;
      }
      if (hitCh%2 == 1 && ~((~0x010) & wfDCBits[wfDCBits.size()-1])) {
        mLClean++;
        if(!detIDIsVetoOnly[hitDetID]) sumELClean += hitENFCal;
      }

      // Calculate mu vetos
      double dtGoodMu_s = 0.0002;
      double dtBadMu_s = 8;
      while(1) {
        if(iMu >= nMu-1) break;
        double dt = muRunTStarts_s[iMu+1]-startTime;
        dt += muTimes_s[iMu+1] - hitT_s;
        if(badScalers[iMu+1]) dt -= dtBadMu_s;
        else dt -= dtGoodMu_s;
        if(dt > 0) break;
        iMu++;
      }
      dtmu_s.push_back((startTime -muRunTStarts_s[iMu]) + hitT_s - muTimes_s[iMu]);
      if(!smallOutput){
        muType.push_back(muTypes[iMu]);
        badScaler.push_back(badScalers[iMu]);
      }
      if(badScalers[iMu] && fabs(dtmu_s[dtmu_s.size()-1]) < dtBadMu_s) muVeto.push_back(true);
      else muVeto.push_back(dtmu_s[dtmu_s.size()-1] > -0.2e-3 && dtmu_s[dtmu_s.size()-1] < 1);

      // tag LN fills
      double utctime = startTime + hitT_s;
      isLNFill.push_back(false);
      for(size_t i=0; i<lnFillTimes.size(); i++) {
        if(lnFillTimes[i]+300 < utctime) continue;
        if(lnFillTimes[i]-900 > utctime) break;
        isLNFill[isLNFill.size()-1] = true;
        break;
      }
    }

    // If no good hits in the event or skipped for some other reason, don't
    // write this event to the output tree.
    if(trapENFCal.size() == 0 || skipMe) continue;

    // finally, fill the tree for this event
    skimTree->Fill();
  }

  // write output tree to output file
  cout << "Closing out skim file..." << endl;
  skimTree->Write("", TObject::kOverwrite);
  fOut->Close();

  return 0;
}

void LoadRun(GATDataSet& ds, size_t i)
{
  ds.AddRunNumber(i);
}

void LoadDataSet(GATDataSet& ds, int dsNumber, size_t i)
{
  // DS0, aka P3JDY
  if(dsNumber == 0) {
    if(i==0) {
      ds.AddRunRange(2580,2580); //1
      ds.AddRunRange(2582,2612); //31
    }
    else if (i==1) {
      ds.AddRunRange(2614,2629); //16
      ds.AddRunRange(2644,2649);//6
      ds.AddRunRange(2658,2673);//16
    }
    else if (i==2) ds.AddRunRange(2689,2715);//28
    else if (i==3) ds.AddRunRange(2717,2750);//33
    else if (i==4) ds.AddRunRange(2751,2784);//33
    else if (i==5) ds.AddRunRange(2785,2820);//35
    else if (i==6) ds.AddRunRange(2821,2855);//34
    else if (i==7) ds.AddRunRange(2856,2890);//34
    else if (i==8) {
      ds.AddRunRange(2891,2907);//17
      ds.AddRunRange(2909,2920);//12
    }
    else if (i==9) ds.AddRunRange(3137,3166);//30
    else if (i==10) ds.AddRunRange(3167,3196);//30
    else if (i==11) ds.AddRunRange(3197,3226);//30
    else if (i==12) ds.AddRunRange(3227,3256);//30
    else if (i==13) {
      ds.AddRunRange(3257,3271);//14
      ds.AddRunRange(3293,3310);//18
    }
    else if (i==14) ds.AddRunRange(3311,3340);//30
    else if (i==15) ds.AddRunRange(3341,3370);//30
    else if (i==16) ds.AddRunRange(3371,3400);//30
    else if (i==17) ds.AddRunRange(3401,3432);//32
    else if (i==18) {
      ds.AddRunRange(3461,3462);//2
      ds.AddRunRange(3464,3500);//36
      //ds.AddRunRange(3464,3580);//117
    }
    else if (i==19) ds.AddRunRange(3501,3530);//30
    else if (i==20) ds.AddRunRange(3531,3560);//30
    else if (i==21) {
      ds.AddRunRange(3561,3580);//20
      ds.AddRunRange(3596,3610);//15
    }
    else if (i==22) ds.AddRunRange(3611,3645);//50
    else if (i==23) {
      ds.AddRunRange(4034,4035);//2
      ds.AddRunRange(4038,4040);//3
      ds.AddRunRange(4045,4074);//30
    }
    else if (i==24) ds.AddRunRange(4075,4104);//30
    else if (i==25) ds.AddRunRange(4105,4134);//30
    else if (i==26) {
      ds.AddRunRange(4239,4245);//7
      ds.AddRunRange(4248,4254);//7
      ds.AddRunRange(4256,4268);//13
    }
    else if (i==27) {
      ds.AddRunRange(4270,4271);//2
      ds.AddRunRange(4273,4283);//11
    }
    else if (i==28) ds.AddRunRange(4285,4311);//27
    else if (i==29) {
      ds.AddRunRange(4313,4318);//6
      ds.AddRunRange(4320,4320);//1
      ds.AddRunRange(4322,4326);//5
      ds.AddRunRange(4328,4336);//9
    }
    else if (i==30) ds.AddRunRange(4338,4361);//24
    else if (i==31) ds.AddRunRange(4363,4382);//20
    else if (i==32) ds.AddRunRange(4384,4401);//18
    else if (i==33) ds.AddRunRange(4403,4428);//26
    else if (i==34) ds.AddRunRange(4436,4454);//19
    else if (i==35) ds.AddRunRange(4457,4489);//33
    else if (i==36) {
      ds.AddRunRange(4491,4493);//3
      ds.AddRunRange(4497,4503);//7
      ds.AddRunRange(4505,4518);//14
    }
    else if (i==37) ds.AddRunRange(4549,4590);//32
    else if (i==38) ds.AddRunRange(4591,4624);//35
    else if (i==39) ds.AddRunRange(4625,4654);//30
    else if (i==40) ds.AddRunRange(4655,4684);//30
    else if (i==41) ds.AddRunRange(4685,4714);//30
    else if (i==42) ds.AddRunRange(4715,4744);//30
    else if (i==43) ds.AddRunRange(4745,4777);//33
    else if (i==44) {
      ds.AddRunRange(4789,4797);//10
      ds.AddRunRange(4800,4831);//32
    }
    else if (i==45) ds.AddRunRange(4854,4872);//19
    else if (i==46) {
      ds.AddRunRange(4874,4883);//10
      ds.AddRunRange(4885,4907);//23
    }
    else if (i==47) ds.AddRunRange(4938,4960);//23
    else if (i==48) {
      ds.AddRunRange(4962,4968);//7
      ds.AddRunRange(4970,4980);//11
    }
    else if (i==49) ds.AddRunRange(5007,5038);//32
    else if (i==50) ds.AddRunRange(5040,5061);//22
    else if (i==51) ds.AddRunRange(5090,5118);//29
    else if (i==52) ds.AddRunRange(5125,5154);//30
    else if (i==53) ds.AddRunRange(5155,5184);//30
    else if (i==54) ds.AddRunRange(5185,5224);//40
    else if (i==55) ds.AddRunRange(5225,5252);//28
    else if (i==56) ds.AddRunRange(5277,5300);//34
    else if (i==57) ds.AddRunRange(5301,5330);//30
    else if (i==58) {
      ds.AddRunRange(5372,5393);//22
      ds.AddRunRange(5405,5414);//10
    }
    else if (i==59) ds.AddRunRange(5449,5479);//30
    else if (i==60) {
      ds.AddRunRange(5480,5501);//23
      ds.AddRunRange(5525,5527);//3
      ds.AddRunRange(5531,5534);//4
    }
    else if (i==61) ds.AddRunRange(5555,5589);//35
    else if (i==62) ds.AddRunRange(5591,5608);//18
    else if (i==63) ds.AddRunRange(5610,5639);//30
    else if (i==64) ds.AddRunRange(5640,5669);//30
    else if (i==65) ds.AddRunRange(5670,5699);//30
    else if (i==66) ds.AddRunRange(5700,5729);//30
    else if (i==67) {
      ds.AddRunRange(5730,5751);//22
      ds.AddRunRange(5753,5764);//12
    }
    else if (i==68) ds.AddRunRange(5766,5795);//30
    else if (i==69) ds.AddRunRange(5796,5822);//27
    else if (i==70) ds.AddRunRange(5826,5850);//24
    else if (i==71) {
      ds.AddRunRange(5889,5890);//2
      ds.AddRunRange(5894,5902);//9
    }
    else if (i==72) {
      ds.AddRunRange(6553,6577);//25
      ds.AddRunRange(6775,6775);//1
    }
    else if (i==73) ds.AddRunRange(6776,6809);//34
    else if (i==74) ds.AddRunRange(6811,6830);//20
    else if (i==75) ds.AddRunRange(6834,6853);//20
    else if (i==76) {
      ds.AddRunRange(6887,6903); //17
      ds.AddRunRange(6957,6963);//7
    }
    else cout << "unknown run sequence (" << i << ") for DS" << dsNumber << endl;
  }

  // DS1, aka P3KJR
  else if(dsNumber == 1) {
    // 52 sets total
    if(i == 0)  ds.AddRunRange(9422, 9440); // 18
    else if(i == 1) {
      ds.AddRunRange(9471, 9487); // 16
      ds.AddRunRange(9492, 9492); // 0
    }
    else if(i == 2) ds.AddRunRange(9536, 9565); // 29
    else if(i == 3)  {
      ds.AddRunRange(9638, 9648); // 10
      ds.AddRunRange(9650, 9668); // 18
    }
    else if(i == 4) {
      ds.AddRunRange(9674, 9676); // 2
      ds.AddRunRange(9678, 9678); // 0
      ds.AddRunRange(9711, 9727); // 16
    }
    else if(i == 5) ds.AddRunRange(9763, 9780); // 17
    else if(i == 6) {
      ds.AddRunRange(9815, 9821); // 6
      ds.AddRunRange(9823, 9832); // 9
      ds.AddRunRange(9848, 9849); // 1
      ds.AddRunRange(9851, 9854); // 3
    }
    else if(i == 7) ds.AddRunRange(9856, 9912); // 56
    else if(i == 8) {
      //ds.AddRunRange(9856, 9912); // 26
      ds.AddRunRange(9928, 9928); // 0
    }
    else if(i == 9) {
      ds.AddRunRange(9952, 9966); // 14
      ds.AddRunRange(10019, 10035); // 16
    }
    else if(i == 10) {
      ds.AddRunRange(10074, 10090); // 16
      ds.AddRunRange(10114, 10125); // 11
    }
    else if(i == 11) ds.AddRunRange(10129, 10149); // 20
    else if(i == 12) ds.AddRunRange(10150, 10171); // 21
    else if(i == 13) ds.AddRunRange(10173, 10203); // 30
    else if(i == 14) ds.AddRunRange(10204, 10231); // 27
    else if(i == 15) {
      ds.AddRunRange(10262, 10278); // 16
      ds.AddRunRange(10298, 10299); // 1
      ds.AddRunRange(10301, 10301); // 0
      ds.AddRunRange(10304, 10308); // 4
    }
    else if(i == 16) ds.AddRunRange(10312, 10342); // 30
    else if(i == 17) {
      ds.AddRunRange(10344, 10350); // 6
      ds.AddRunRange(10378, 10394); // 16
      ds.AddRunRange(10552, 10558); // 6
    }
    else if(i == 18) ds.AddRunRange(10608, 10648); // 40
    else if(i == 19) ds.AddRunRange(10651, 10677); // 26
    else if(i == 20) ds.AddRunRange(10679, 10717); // 38
    else if(i == 21) {
      ds.AddRunRange(10745, 10761); // 16
      ds.AddRunRange(10788, 10803); // 15
    }
    else if(i == 22) {
      ds.AddRunRange(10830, 10845); // 15
      ds.AddRunRange(10963, 10976); // 13
    }
    else if(i == 23) {
      ds.AddRunRange(11002, 11008); // 6
      ds.AddRunRange(11010, 11019); // 9
      ds.AddRunRange(11046, 11066); // 20
    }
    else if(i == 24) ds.AddRunRange(11083, 11113); // 30
    else if(i == 25) ds.AddRunRange(11114, 11144); // 30
    else if(i == 26) ds.AddRunRange(11145, 11175); // 30
    else if(i == 27) {
      ds.AddRunRange(11176, 11200); // 24
      ds.AddRunRange(11350, 11350); // 0
      ds.AddRunRange(11403, 11410); // 7
    }
    else if(i == 28) {
      ds.AddRunRange(11414, 11417); // 3
      ds.AddRunRange(11419, 11426); // 7
      ds.AddRunRange(11428, 11432); // 4
      ds.AddRunRange(11434, 11444); // 10
      ds.AddRunRange(11446, 11451); // 5
    }
    else if(i == 29) {
      ds.AddRunRange(11453, 11453); // 0
      ds.AddRunRange(11455, 11458); // 3
      ds.AddRunRange(11466, 11476); // 10
      ds.AddRunRange(11477, 11483); // 6
      ds.AddRunRange(12445, 12445); // 0
      ds.AddRunRange(12466, 12467); // 1
      ds.AddRunRange(12477, 12483); // 6
      ds.AddRunRange(12486, 12493); // 7
    }
    else if(i == 30) ds.AddRunRange(12521, 12550); // 29
    else if(i == 31) ds.AddRunRange(12551, 12580); // 60
    else if(i == 32) {
      ds.AddRunRange(12607, 12625); // 18
      ds.AddRunRange(12636, 12647); // 11
      ds.AddRunRange(12652, 12653); // 1
    }
    else if(i == 33) ds.AddRunRange(12664, 12675); // 11
    else if(i == 34) ds.AddRunRange(12677, 12724); // 47
    else if(i == 35) ds.AddRunRange(12735, 12765); // 33
    else if(i == 36) ds.AddRunRange(12766, 12798); // 32
    else if(i == 37) {
      ds.AddRunRange(12816, 12816); // 0
      ds.AddRunRange(12819, 12819); // 1
      ds.AddRunRange(12821, 12821); // 0
      ds.AddRunRange(12823, 12824); // 1
      ds.AddRunRange(12827, 12831); // 5
      ds.AddRunRange(12833, 12838); // 6
      ds.AddRunRange(12842, 12842); // 0
      ds.AddRunRange(12843, 12861); // 18
      ds.AddRunRange(12875, 12875); // 0
    }
    else if(i == 38) ds.AddRunRange(13000, 13028); // 28
    else if(i == 39) {
      ds.AddRunRange(13029, 13053); // 24
      ds.AddRunRange(13055, 13056); // 1
    }
    else if(i == 40) {
      ds.AddRunRange(13066, 13072); // 6
      ds.AddRunRange(13074, 13074); // 1
      ds.AddRunRange(13076, 13092); // 16
      ds.AddRunRange(13094, 13096); // 2
    }
    else if(i == 41) {
      ds.AddRunRange(13100, 13115); // 15
      ds.AddRunRange(13117, 13119); // 2
      ds.AddRunRange(13123, 13137); // 14
    }
    else if(i == 42) {
      ds.AddRunRange(13148, 13150); // 2
      ds.AddRunRange(13154, 13156); // 2
      ds.AddRunRange(13186, 13189); // 3
      ds.AddRunRange(13191, 13211); // 20
    }
    else if(i == 43) ds.AddRunRange(13212, 13242); // 30
    else if(i == 44) ds.AddRunRange(13243, 13275); // 32
    else if(i == 45) {
      ds.AddRunRange(13276, 13287); // 11
      ds.AddRunRange(13304, 13304); // 0
      ds.AddRunRange(13306, 13311); // 5
      ds.AddRunRange(13313, 13325); // 12
    }
    else if(i == 46) {
      ds.AddRunRange(13326, 13350); // 24
      ds.AddRunRange(13362, 13368); // 6
    }
    else if(i == 47) {
      ds.AddRunRange(13369,13383); // 15
      ds.AddRunRange(13395,13411); // 17
    }
    else if(i == 48) ds.AddRunRange(13519,13548); // 30
    else if(i == 49) {
      ds.AddRunRange(13572,13573); // 2
      ds.AddRunRange(13667,13688); // 22
      ds.AddRunRange(13699,13704); // 6
      ds.AddRunRange(13715,13719); // 5
    }
    else if(i == 50) {
      ds.AddRunRange(14010,14040); // 31
      ds.AddRunRange(14041,14041); // 1
    }
    else if(i == 51) {
      ds.AddRunRange(14342,14372); // 31
      ds.AddRunRange(14386,14387); // 2
    }
    else cout << "unknown run sequence (" << i << ") for DS" << dsNumber << endl;
  }

  else cout << "LoadDataSet(): unknown dataset number DS" << dsNumber << endl;
}

void LoadActiveMasses(map<int,double>& activeMassForDetID_g, int dsNumber)
{
  if(dsNumber == 0 || dsNumber == 1) {
    activeMassForDetID_g[1426981] = 509.9;
    activeMassForDetID_g[1425750] = 978.8;
    activeMassForDetID_g[1426612] = 811.3;
    activeMassForDetID_g[1425380] = 967.9;
    activeMassForDetID_g[28474] = 587.0;
    activeMassForDetID_g[1426640] = 722.9;
    activeMassForDetID_g[1426650] = 659.1;
    activeMassForDetID_g[1426622] = 688.6;
    activeMassForDetID_g[28480] = 577.6;
    activeMassForDetID_g[1426980] = 886.3;
    activeMassForDetID_g[1425381] = 949.0;
    activeMassForDetID_g[1425730] = 1023.8;
    activeMassForDetID_g[28455] = 584.2;
    activeMassForDetID_g[28470] = 590.8;
    activeMassForDetID_g[28463] = 594.5;
    activeMassForDetID_g[28465] = 571.1;
    activeMassForDetID_g[28469] = 593.3;
    activeMassForDetID_g[28477] = 579.5;
    activeMassForDetID_g[1425751] = 730.6;
    activeMassForDetID_g[1426610] = 632.0;
    activeMassForDetID_g[1425731] = 982.0;
    activeMassForDetID_g[1425742] = 731.6;
    activeMassForDetID_g[1426611] = 675.5;
    activeMassForDetID_g[1425740] = 701.4;
    activeMassForDetID_g[1426620] = 572.3;
    activeMassForDetID_g[28482] = 588.0;
    activeMassForDetID_g[1425741] = 709.9;
    activeMassForDetID_g[1426621] = 590.9;
    activeMassForDetID_g[1425370] = 964.3;
  }

  else cout << "LoadActiveMasses(): unknown dataset number DS" << dsNumber << endl;
}

double GetAOverENorm(int channel, double a50, double a100, double a200, double trapENF, double trapENFCal, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel==692) return (fabs(a50/trapENF)-(0.00000014892843*(trapENFCal)))/0.00690959;
    if(channel==693) return (fabs(a50/trapENF)-(0.00000014892843*(trapENFCal*0.299765)))/0.00690959;
    if(channel==690) return (fabs(a200/trapENF)-(-0.00000000619575*(trapENFCal)))/0.00461375;
    if(channel==691) return (fabs(a200/trapENF)-(-0.00000000619575*(trapENFCal*0.296368)))/0.00461375;
    if(channel==688) return (fabs(a50/trapENF)-(-0.00000002716361*(trapENFCal)))/0.00684412;
    if(channel==689) return (fabs(a50/trapENF)-(-0.00000002716361*(trapENFCal*0.296717)))/0.00684412;
    if(channel==640) return (fabs(a100/trapENF)-(-0.00000002167902*(trapENFCal)))/0.00620691;
    if(channel==641) return (fabs(a100/trapENF)-(-0.00000002167902*(trapENFCal*0.29473)))/0.00620691;
    if(channel==674) return (fabs(a200/trapENF)-(-0.00000001286328*(trapENFCal)))/0.00476108;
    if(channel==675) return (fabs(a200/trapENF)-(-0.00000001286328*(trapENFCal*0.298476)))/0.00476108;
    if(channel==576) return (fabs(a100/trapENF)-(-0.00000001202658*(trapENFCal)))/0.00660908;
    if(channel==577) return (fabs(a100/trapENF)-(-0.00000001202658*(trapENFCal*0.295277)))/0.00660908;
    if(channel==614) return (fabs(a100/trapENF)-(0.00000000600596*(trapENFCal)))/0.00658321;
    if(channel==615) return (fabs(a100/trapENF)-(0.00000000600596*(trapENFCal*0.294667)))/0.00658321;
    if(channel==610) return (fabs(a100/trapENF)-(-0.00000005247095*(trapENFCal)))/0.00421871;
    if(channel==611) return (fabs(a100/trapENF)-(-0.00000005247095*(trapENFCal*0.298349)))/0.00421871;
    if(channel==608) return (fabs(a200/trapENF)-(-0.00000001402004*(trapENFCal)))/0.00529826;
    if(channel==609) return (fabs(a200/trapENF)-(-0.00000001402004*(trapENFCal*0.293618)))/0.00529826;
    if(channel==598) return (fabs(a50/trapENF)-(0.00000003987425*(trapENFCal)))/0.00609923;
    if(channel==599) return (fabs(a50/trapENF)-(0.00000003987425*(trapENFCal*0.298569)))/0.00609923;
    if(channel==600) return (fabs(a200/trapENF)-(-0.00000001350690*(trapENFCal)))/0.00486026;
    if(channel==601) return (fabs(a200/trapENF)-(-0.00000001350690*(trapENFCal*0.300375)))/0.00486026;
    if(channel==594) return (fabs(a200/trapENF)-(-0.00000000603499*(trapENFCal)))/0.00496379;
    if(channel==595) return (fabs(a200/trapENF)-(-0.00000000603499*(trapENFCal*0.299192)))/0.00496379;
    if(channel==592) return (fabs(a200/trapENF)-(-0.00000000562414*(trapENFCal)))/0.00485191;
    if(channel==593) return (fabs(a200/trapENF)-(-0.00000000562414*(trapENFCal*0.300634)))/0.00485191;
    if(channel==664) return (fabs(a200/trapENF)-(-0.00000001597636*(trapENFCal)))/0.00477151;
    if(channel==665) return (fabs(a200/trapENF)-(-0.00000001597636*(trapENFCal*0.300183)))/0.00477151;
    if(channel==662) return (fabs(a200/trapENF)-(-0.00000001024009*(trapENFCal)))/0.00480596;
    if(channel==663) return (fabs(a200/trapENF)-(-0.00000001024009*(trapENFCal*0.149383)))/0.00480596;
    if(channel==656) return (fabs(a100/trapENF)-(-0.00000001140632*(trapENFCal)))/0.0061491;
    if(channel==657) return (fabs(a100/trapENF)-(-0.00000001140632*(trapENFCal*0.294046)))/0.0061491;
    if(channel==696) return (fabs(a200/trapENF)-(-0.00000001341208*(trapENFCal)))/0.00466965;
    if(channel==697) return (fabs(a200/trapENF)-(-0.00000001341208*(trapENFCal*0.297445)))/0.00466965;
    if(channel==628) return (fabs(a100/trapENF)-(-0.00000005461271*(trapENFCal)))/0.00618007;
    if(channel==629) return (fabs(a100/trapENF)-(-0.00000005461271*(trapENFCal*0.302283)))/0.00618007;
    if(channel==626) return (fabs(a100/trapENF)-(-0.00000002099260*(trapENFCal)))/0.00632098;
    if(channel==627) return (fabs(a100/trapENF)-(-0.00000002099260*(trapENFCal*0.296453)))/0.00632098;
    if(channel==624) return (fabs(a100/trapENF)-(-0.00000001254000*(trapENFCal)))/0.00589336;
    if(channel==625) return (fabs(a100/trapENF)-(-0.00000001254000*(trapENFCal*0.292326)))/0.00589336;
    if(channel==646) return (fabs(a200/trapENF)-(-0.00000001108327*(trapENFCal)))/0.00459846;
    if(channel==647) return (fabs(a200/trapENF)-(-0.00000001108327*(trapENFCal*0.29658)))/0.00459846;
    if(channel==644) return (fabs(a100/trapENF)-(0.00000000533074*(trapENFCal)))/0.00625633;
    if(channel==645) return (fabs(a100/trapENF)-(0.00000000533074*(trapENFCal*0.294694)))/0.00625633;
    if(channel==642) return (fabs(a200/trapENF)-(-0.00000000600096*(trapENFCal)))/0.0048162;
    if(channel==643) return (fabs(a200/trapENF)-(-0.00000000600096*(trapENFCal*0.298424)))/0.0048162    ;
  }

  else if(dsNumber == 1) {
    if(channel == 582) return (fabs(a200/trapENF)-(-0.00000000819708*(trapENFCal)))/0.00464055;
    if(channel == 583) return (fabs(a200/trapENF)-(-0.00000000819708*(trapENFCal*.301385)))/0.00464055;
    if(channel == 580) return (fabs(a50/trapENF)-(-0.00000002519509*(trapENFCal)))/0.00458294;
    if(channel == 581) return (fabs(a50/trapENF)-(-0.00000002519509*(trapENFCal*.298659)))/0.00458294;
    if(channel == 578) return (fabs(a100/trapENF)-(-0.00000002481901*(trapENFCal)))/0.00636773;
    if(channel == 579) return (fabs(a100/trapENF)-(-0.00000002481901*(trapENFCal*.297716)))/0.00636773;
    if(channel == 692) return (fabs(a200/trapENF)-(-0.00000000924144*(trapENFCal)))/0.00487448;
    if(channel == 693) return (fabs(a200/trapENF)-(-0.00000000924144*(trapENFCal*.2996647)))/0.00487448;
    if(channel == 648) return (fabs(a200/trapENF)-(-0.00000001531963*(trapENFCal)))/0.00489646;
    if(channel == 649) return (fabs(a200/trapENF)-(-0.00000001531963*(trapENFCal*.29774)))/0.00489646;
    if(channel == 640) return (fabs(a50 /trapENF)-(-0.00000003320082*(trapENFCal)))/0.00734893;
    if(channel == 641) return (fabs(a50 /trapENF)-(-0.00000003320082*(trapENFCal*.294796)))/0.00734893;
    if(channel == 616) return (fabs(a200/trapENF)-(-0.00000007284703*(trapENFCal)))/0.00447279;
    if(channel == 617) return (fabs(a200/trapENF)-(-0.00000007284703*(trapENFCal*.29391)))/0.00447279;
    if(channel == 610) return (fabs(a100 /trapENF)-(-0.00000001561029*(trapENFCal)))/0.00636434;
    if(channel == 611) return (fabs(a100 /trapENF)-(-0.00000001561029*(trapENFCal*.298432)))/0.00636434;
    if(channel == 608) return (fabs(a50 /trapENF)-(-0.00000002710169*(trapENFCal)))/0.00732955;
    if(channel == 609) return (fabs(a50 /trapENF)-(-0.00000002710169*(trapENFCal*.293385)))/0.00732955;
    if(channel == 664) return (fabs(a200/trapENF)-(-0.00000002234744*(trapENFCal)))/0.00438536;
    if(channel == 665) return (fabs(a200/trapENF)-(-0.00000002234744*(trapENFCal*.299588)))/0.00438536;
    if(channel == 672) return (fabs(a100/trapENF)-(-0.00000002853387*(trapENFCal)))/0.00644415;
    if(channel == 673) return (fabs(a100/trapENF)-(-0.00000002853387*(trapENFCal*.294691)))/0.00644415;
    if(channel == 632) return (fabs(a50/trapENF)-(-0.00000001627841*(trapENFCal)))/0.00690461;
    if(channel == 633) return (fabs(a50/trapENF)-(-0.00000001627841*(trapENFCal*.299157)))/0.00690461;
    if(channel == 626) return (fabs(a100/trapENF)-(-0.00000000133761*(trapENFCal)))/0.00615713;
    if(channel == 627) return (fabs(a100/trapENF)-(-0.00000000133761*(trapENFCal*.29593)))/0.00615713;
    if(channel == 690) return (fabs(a100/trapENF)-(-0.00000001916973*(trapENFCal)))/0.00618773;
    if(channel == 691) return (fabs(a100/trapENF)-(-0.00000001916973*(trapENFCal*.29633)))/0.00618773;
    if(channel == 600) return (fabs(a100/trapENF)-(-0.00000002458613*(trapENFCal)))/0.0059823;
    if(channel == 601) return (fabs(a100/trapENF)-(-0.00000002458613*(trapENFCal*.295671)))/0.0059823;
    if(channel == 598) return (fabs(a100/trapENF)-(-0.00000001621977*(trapENFCal)))/0.00620608;
    if(channel == 599) return (fabs(a100/trapENF)-(-0.00000001621977*(trapENFCal*.296562)))/0.00620608;
    if(channel == 594) return (fabs(a100 /trapENF)-(-0.00000000679516*(trapENFCal)))/0.00612278;
    if(channel == 595) return (fabs(a100 /trapENF)-(-0.00000000679516*(trapENFCal*.294345)))/0.00612278;
    if(channel == 592) return (fabs(a200 /trapENF)-(-0.00000004463564*(trapENFCal)))/0.00496273;
    if(channel == 593) return (fabs(a200 /trapENF)-(-0.00000004463564*(trapENFCal*.303767)))/0.00496273;
  }

  else cout << "GetAOverENorm(): unknown dataset number DS" << dsNumber << endl;
  return 0.0;
}

double GetAOverENorm85(int channel, double a50, double a100, double a200, double trapENF, double trapENFCal, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel==692) return (fabs(a100/trapENF)-(-0.00000004546863*(trapENFCal)))/0.00654001;
    if(channel==693) return (fabs(a100/trapENF)-(-0.00000004546863*(trapENFCal*0.299765)))/0.00654001;
    if(channel==690) return (fabs(a50/trapENF)-(0.00000002275037*(trapENFCal)))/0.00639616;
    if(channel==691) return (fabs(a50/trapENF)-(0.00000002275037*(trapENFCal*0.296368)))/0.00639616;
    if(channel==688) return (fabs(a100/trapENF)-(-0.00000001410593*(trapENFCal)))/0.00639616;
    if(channel==689) return (fabs(a100/trapENF)-(-0.00000001410593*(trapENFCal*0.296717)))/0.00639616;
    if(channel==640) return (fabs(a100/trapENF)-(-0.00000002167902*(trapENFCal)))/0.00640769;
    if(channel==641) return (fabs(a100/trapENF)-(-0.00000002167902*(trapENFCal*0.29473)))/0.00640769;
    if(channel==674) return (fabs(a200/trapENF)-(-0.00000001286328*(trapENFCal)))/0.00490712;
    if(channel==675) return (fabs(a200/trapENF)-(-0.00000001286328*(trapENFCal*0.298476)))/0.00490712;
    if(channel==576) return (fabs(a100/trapENF)-(-0.00000001202658*(trapENFCal)))/0.00665333;
    if(channel==577) return (fabs(a100/trapENF)-(-0.00000001202658*(trapENFCal*0.295277)))/0.00665333;
    if(channel==614) return (fabs(a100/trapENF)-(0.00000000600596*(trapENFCal)))/0.00659255;
    if(channel==615) return (fabs(a100/trapENF)-(0.00000000600596*(trapENFCal*0.294667)))/0.00659255;
    if(channel==610) return (fabs(a200/trapENF)-(-0.00000001907038*(trapENFCal)))/0.00436024;
    if(channel==611) return (fabs(a200/trapENF)-(-0.00000001907038*(trapENFCal*0.298349)))/0.00436024;
    if(channel==608) return (fabs(a100/trapENF)-(-0.00000002840734*(trapENFCal)))/0.00692363;
    if(channel==609) return (fabs(a100/trapENF)-(-0.00000002840734*(trapENFCal*0.293618)))/0.00692363;
    if(channel==598) return (fabs(a50/trapENF)-(0.00000003987425*(trapENFCal)))/0.00619092;
    if(channel==599) return (fabs(a50/trapENF)-(0.00000003987425*(trapENFCal*0.298569)))/0.00619092;
    if(channel==600) return (fabs(a200/trapENF)-(-0.00000001350690*(trapENFCal)))/0.00486836;
    if(channel==601) return (fabs(a200/trapENF)-(-0.00000001350690*(trapENFCal*0.300375)))/0.00486836;
    if(channel==594) return (fabs(a50/trapENF)-(0.00000000589302*(trapENFCal)))/0.00695222;
    if(channel==595) return (fabs(a50/trapENF)-(0.00000000589302*(trapENFCal*0.299192)))/0.00695222;
    if(channel==592) return (fabs(a100/trapENF)-(0.00000001978263*(trapENFCal)))/0.00642204;
    if(channel==593) return (fabs(a100/trapENF)-(0.00000001978263*(trapENFCal*0.300634)))/0.00642204;
    if(channel==664) return (fabs(a200/trapENF)-(-0.00000001597636*(trapENFCal)))/0.00477869;
    if(channel==665) return (fabs(a200/trapENF)-(-0.00000001597636*(trapENFCal*0.300183)))/0.00477869;
    if(channel==662) return (fabs(a100/trapENF)-(-0.00000001803842*(trapENFCal)))/0.00613487;
    if(channel==663) return (fabs(a100/trapENF)-(-0.00000001803842*(trapENFCal*0.149383)))/0.00613487;
    if(channel==656) return (fabs(a100/trapENF)-(-0.00000001140632*(trapENFCal)))/0.00618347;
    if(channel==657) return (fabs(a100/trapENF)-(-0.00000001140632*(trapENFCal*0.294046)))/0.00618347;
    if(channel==696) return (fabs(a200/trapENF)-(-0.00000001341208*(trapENFCal)))/0.0046871;
    if(channel==697) return (fabs(a200/trapENF)-(-0.00000001341208*(trapENFCal*0.297445)))/0.0046871;
    if(channel==628) return (fabs(a200/trapENF)-(-0.00000001274462*(trapENFCal)))/0.00480479;
    if(channel==629) return (fabs(a200/trapENF)-(-0.00000001274462*(trapENFCal*0.302283)))/0.00480479;
    if(channel==626) return (fabs(a100/trapENF)-(-0.00000002099260*(trapENFCal)))/0.00634179;
    if(channel==627) return (fabs(a100/trapENF)-(-0.00000002099260*(trapENFCal*0.296453)))/0.00634179;
    if(channel==624) return (fabs(a50/trapENF)-(-0.00000000973262*(trapENFCal)))/0.00649791;
    if(channel==625) return (fabs(a50/trapENF)-(-0.00000000973262*(trapENFCal*0.292326)))/0.00649791;
    if(channel==646) return (fabs(a200/trapENF)-(-0.00000001108327*(trapENFCal)))/0.00460577;
    if(channel==647) return (fabs(a200/trapENF)-(-0.00000001108327*(trapENFCal*0.29658)))/0.00460577;
    if(channel==644) return (fabs(a100/trapENF)-(0.00000000533074*(trapENFCal)))/0.00627598;
    if(channel==645) return (fabs(a100/trapENF)-(0.00000000533074*(trapENFCal*0.294694)))/0.00627598;
    if(channel==642) return (fabs(a100/trapENF)-(-0.00000001156293*(trapENFCal)))/0.00617827;
    if(channel==643) return (fabs(a100/trapENF)-(-0.00000001156293*(trapENFCal*0.298424)))/0.00617827       ;
  }

  else if(dsNumber == 1) {
    if(channel == 582) return (fabs(a200/trapENF)-(-8.19708E-09*(trapENFCal)))/0.00464859;
    if(channel == 583) return (fabs(a200/trapENF)-(-8.19708E-09*(trapENFCal*.301385)))/0.00464859;
    if(channel == 580) return (fabs(a50/trapENF)-(-2.47886E-08*(trapENFCal)))/0.00459782;
    if(channel == 581) return (fabs(a50/trapENF)-(-2.47886E-08*(trapENFCal*.298659)))/0.00459782;
    if(channel == 578) return (fabs(a50/trapENF)-(-2.2312E-08*(trapENFCal)))/0.00705478;
    if(channel == 579) return (fabs(a50/trapENF)-(-2.2312E-08*(trapENFCal*.297716)))/0.00705478;
    if(channel == 692) return (fabs(a100/trapENF)-(-1.6445E-08*(trapENFCal)))/0.00623492;
    if(channel == 693) return (fabs(a100/trapENF)-(-1.6445E-08*(trapENFCal*.2996647)))/0.00623492;
    if(channel == 648) return (fabs(a200/trapENF)-(-1.32103E-08*(trapENFCal)))/0.00492125;
    if(channel == 649) return (fabs(a200/trapENF)-(-1.32103E-08*(trapENFCal*.29774)))/0.00492125;
    if(channel == 640) return (fabs(a100 /trapENF)-(-1.04414E-08*(trapENFCal)))/0.00661752;
    if(channel == 641) return (fabs(a100 /trapENF)-(-1.04414E-08*(trapENFCal*.294796)))/0.00661752;
    if(channel == 616) return (fabs(a200/trapENF)-(-1.00872E-07*(trapENFCal)))/0.00451685;
    if(channel == 617) return (fabs(a200/trapENF)-(-1.00872E-07*(trapENFCal*.29391)))/0.00451685;
    if(channel == 610) return (fabs(a50 /trapENF)-(-2.47924E-08*(trapENFCal)))/0.00710022;
    if(channel == 611) return (fabs(a50 /trapENF)-(-2.47924E-08*(trapENFCal*.298432)))/0.00710022;
    if(channel == 608) return (fabs(a50 /trapENF)-(-1.80863E-08*(trapENFCal)))/0.00731791;
    if(channel == 609) return (fabs(a50 /trapENF)-(-1.80863E-08*(trapENFCal*.293385)))/0.00731791;
    if(channel == 664) return (fabs(a100/trapENF)-(-6.10461E-08*(trapENFCal)))/0.00554796;
    if(channel == 665) return (fabs(a100/trapENF)-(-6.10461E-08*(trapENFCal*.299588)))/0.00554796;
    if(channel == 672) return (fabs(a200/trapENF)-(-9.28838E-09*(trapENFCal)))/0.00496773;
    if(channel == 673) return (fabs(a200/trapENF)-(-9.28838E-09*(trapENFCal*.294691)))/0.00496773;
    if(channel == 632) return (fabs(a100/trapENF)-(-1.71507E-08*(trapENFCal)))/0.00622425;
    if(channel == 633) return (fabs(a100/trapENF)-(-1.71507E-08*(trapENFCal*.299157)))/0.00622425;
    if(channel == 626) return (fabs(a50/trapENF)-(-1.89126E-08*(trapENFCal)))/0.00679347;
    if(channel == 627) return (fabs(a50/trapENF)-(-1.89126E-08*(trapENFCal*.29593)))/0.00679347;
    if(channel == 690) return (fabs(a50/trapENF)-(-2.87825E-08*(trapENFCal)))/0.00689086;
    if(channel == 691) return (fabs(a50/trapENF)-(-2.87825E-08*(trapENFCal*.29633)))/0.00689086;
    if(channel == 600) return (fabs(a100/trapENF)-(-2.45861E-08*(trapENFCal)))/0.00601021;
    if(channel == 601) return (fabs(a100/trapENF)-(-2.45861E-08*(trapENFCal*.295671)))/0.00601021;
    if(channel == 598) return (fabs(a100/trapENF)-(-1.62198E-08*(trapENFCal)))/0.00622601;
    if(channel == 599) return (fabs(a100/trapENF)-(-1.62198E-08*(trapENFCal*.296562)))/0.00622601;
    if(channel == 594) return (fabs(a50 /trapENF)-(-1.95205E-09*(trapENFCal)))/0.00679792;
    if(channel == 595) return (fabs(a50 /trapENF)-(-1.95205E-09*(trapENFCal*.294345)))/0.00679792;
    if(channel == 592) return (fabs(a200 /trapENF)-(-1.71877E-08*(trapENFCal)))/0.00496213         ;
    if(channel == 593) return (fabs(a200 /trapENF)-(-1.71877E-08*(trapENFCal*.303767)))/0.00496213         ;
  }

  else cout << "GetAOverENorm85(): unknown dataset number DS" << dsNumber << endl;
  return 0.0;
}

double GetAvsE(int channel, double TSCurrent50nsMax, double TSCurrent100nsMax, double TSCurrent200nsMax, double trapENF, double trapENFCal, int dsNumber)
{
  if(dsNumber == 1) {

    if    (channel==582) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.01854576347230)-(0.00468276622549*(trapENFCal))-(-0.00000001407301*(trapENFCal)*(trapENFCal)))/-0.0702263);
    if    (channel==583) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.16369855142045)-(0.00628861548784*(trapENFCal))-(0.00000009715465*(trapENFCal)*(trapENFCal)))/-0.208138);
    if    (channel==580) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.05674199580730)-(0.00425035340190*(trapENFCal))-(0.00000000171536*(trapENFCal)*(trapENFCal)))/-0.118208);
    if    (channel==581) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.32060502828949)-(0.00423312877746*(trapENFCal))-(0.00000014226922*(trapENFCal)*(trapENFCal)))/-0.222847);
    if    (channel==578) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.08927516275763)-(0.00703937089571*(trapENFCal))-(0.00000003604893*(trapENFCal)*(trapENFCal)))/-0.248755);
    if    (channel==579) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.29254857541595)-(0.00689197923341*(trapENFCal))-(0.00000015350572*(trapENFCal)*(trapENFCal)))/-0.216434);
    if    (channel==692) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.06259010350846)-(0.00483586165925*(trapENFCal))-(0.00000001627416*(trapENFCal)*(trapENFCal)))/-0.0820596);
    if    (channel==693) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.12134467443495)-(0.00691895832301*(trapENFCal))-(-0.00000000470471*(trapENFCal)*(trapENFCal)))/-0.224861);
    if    (channel==648) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.06437389210040)-(0.00719101008714*(trapENFCal))-(-0.00000005202873*(trapENFCal)*(trapENFCal)))/-0.504072);
    if    (channel==649) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.19367715373649)-(0.00704134129392*(trapENFCal))-(0.00000000395276*(trapENFCal)*(trapENFCal)))/-0.47942);
    if    (channel==640) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10070555974961)-(0.00742514512250*(trapENFCal))-(0.00000002369359*(trapENFCal)*(trapENFCal)))/-0.387218);
    if    (channel==641) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.07526220050901)-(0.00667997978625*(trapENFCal))-(0.00000002407041*(trapENFCal)*(trapENFCal)))/-0.270941);
    if    (channel==616) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02486114215672)-(0.00455738079802*(trapENFCal))-(-0.00000010811773*(trapENFCal)*(trapENFCal)))/-0.129789);
    if    (channel==617) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02787931415456)-(0.00458580492241*(trapENFCal))-(-0.00000012117977*(trapENFCal)*(trapENFCal)))/-0.124577);
    if    (channel==610) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.05649949694814)-(0.00716511550872*(trapENFCal))-(-0.00000003055894*(trapENFCal)*(trapENFCal)))/-0.199234);
    if    (channel==611) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.14399987715688)-(0.00714778413265*(trapENFCal))-(0.00000008102384*(trapENFCal)*(trapENFCal)))/-0.310912);
    if    (channel==608) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.05524811515552)-(0.00743711815084*(trapENFCal))-(-0.00000003756238*(trapENFCal)*(trapENFCal)))/-0.281594);
    if    (channel==609) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.15876502286021)-(0.00732122138428*(trapENFCal))-(0.00000004729350*(trapENFCal)*(trapENFCal)))/-0.430407);
    if    (channel==664) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.07779551311681)-(0.00607567539278*(trapENFCal))-(-0.00000004657088*(trapENFCal)*(trapENFCal)))/-0.250382);
    if    (channel==665) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.07304968321561)-(0.00564943039142*(trapENFCal))-(-0.00000008192775*(trapENFCal)*(trapENFCal)))/-0.182482);
    if    (channel==672) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.06848460163308)-(0.00725365163400*(trapENFCal))-(-0.00000002128673*(trapENFCal)*(trapENFCal)))/-0.290611);
    if    (channel==673) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.05831439575329)-(0.00500538585273*(trapENFCal))-(0.00000000830948*(trapENFCal)*(trapENFCal)))/-0.137654);
    if    (channel==632) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.07917762952496)-(0.00694792043183*(trapENFCal))-(-0.00000000758383*(trapENFCal)*(trapENFCal)))/-0.181945);
    if    (channel==633) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.06954722595344)-(0.00626199652534*(trapENFCal))-(-0.00000001155030*(trapENFCal)*(trapENFCal)))/-0.133663);
    if    (channel==626) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.01810678070787)-(0.00488771705012*(trapENFCal))-(-0.00000000804786*(trapENFCal)*(trapENFCal)))/-0.086283);
    if    (channel==627) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.02984410185282)-(0.00496024053732*(trapENFCal))-(-0.00000000937571*(trapENFCal)*(trapENFCal)))/-0.102878);
    if    (channel==690) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.02821723031971)-(0.00621139001811*(trapENFCal))-(0.00000000902572*(trapENFCal)*(trapENFCal)))/-0.154616);
    if    (channel==691) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.06022352707069)-(0.00629275630451*(trapENFCal))-(0.00000005072433*(trapENFCal)*(trapENFCal)))/-0.152819);
    if    (channel==600) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10055600376558)-(0.00656624540158*(trapENFCal))-(0.00000001755970*(trapENFCal)*(trapENFCal)))/-0.195797);
    if    (channel==601) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.07060566139558)-(0.00462456482179*(trapENFCal))-(0.00000003567736*(trapENFCal)*(trapENFCal)))/-0.059421);
    if    (channel==598) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.08363810522446)-(0.00697336661343*(trapENFCal))-(-0.00000004853185*(trapENFCal)*(trapENFCal)))/-0.174753);
    if    (channel==599) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.10050448277932)-(0.00620078423930*(trapENFCal))-(-0.00000000963453*(trapENFCal)*(trapENFCal)))/-0.132879);
    if    (channel==594) return (-1*((TSCurrent100nsMax*(trapENFCal)/(trapENF))-(0.03637253831151)-(0.00614898547205*(trapENFCal))-(0.00000002174726*(trapENFCal)*(trapENFCal)))/-0.16588);
    if    (channel==595) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.13382021226886)-(0.00684308040220*(trapENFCal))-(0.00000009806745*(trapENFCal)*(trapENFCal)))/-0.2612);
    if    (channel==592) return (-1*((TSCurrent200nsMax*(trapENFCal)/(trapENF))-(0.00308216986982)-(0.00513834129293*(trapENFCal))-(-0.00000006730542*(trapENFCal)*(trapENFCal)))/-0.231685);
    if    (channel==593) return (-1*((TSCurrent50nsMax*(trapENFCal)/(trapENF))-(0.10854259868427)-(0.00695707196683*(trapENFCal))-(0.00000002953986*(trapENFCal)*(trapENFCal)))/-0.961633);
}

else cout << "GetAvsE(): unknown dataset number DS" << dsNumber << endl;
  return 0.0;
}

void LoadLNFillTimes(vector<double>& lnFillTimes, int dsNumber)
{
  // we don't really need to make DS-specific lists, but look-up
  // time is shorter if the lists are shorter.

  if(dsNumber == 0) {
    lnFillTimes.push_back(1435870160);
    lnFillTimes.push_back(1436020533);
    lnFillTimes.push_back(1436168622);
    lnFillTimes.push_back(1436316362);
    lnFillTimes.push_back(1436463827);
    lnFillTimes.push_back(1436614832);
    lnFillTimes.push_back(1436759402);
    lnFillTimes.push_back(1436908966);
    lnFillTimes.push_back(1437058884);
    lnFillTimes.push_back(1437060687);
    lnFillTimes.push_back(1437204895);
    lnFillTimes.push_back(1437350621);
    lnFillTimes.push_back(1437638022);
    lnFillTimes.push_back(1437781652);
    lnFillTimes.push_back(1437926383);
    lnFillTimes.push_back(1438070819);
    lnFillTimes.push_back(1438121812);
    lnFillTimes.push_back(1438212611);
    lnFillTimes.push_back(1438352837);
    lnFillTimes.push_back(1438496004);
    lnFillTimes.push_back(1438639835);
    lnFillTimes.push_back(1438782282);
    lnFillTimes.push_back(1438932169);
    lnFillTimes.push_back(1439080850);
    lnFillTimes.push_back(1439226624);
    lnFillTimes.push_back(1439371591);
    lnFillTimes.push_back(1439514291);
    lnFillTimes.push_back(1439657247);
    lnFillTimes.push_back(1439801649);
    lnFillTimes.push_back(1439944622);
    lnFillTimes.push_back(1440080091);
    lnFillTimes.push_back(1440192295);
    lnFillTimes.push_back(1440335863);
    lnFillTimes.push_back(1440477294);
    lnFillTimes.push_back(1440618411);
    lnFillTimes.push_back(1440759782);
    lnFillTimes.push_back(1440902658);
    lnFillTimes.push_back(1441046730);
    lnFillTimes.push_back(1441187019);
    lnFillTimes.push_back(1441325878);
    lnFillTimes.push_back(1441462000);
    lnFillTimes.push_back(1441600116);
    lnFillTimes.push_back(1441741779);
    lnFillTimes.push_back(1441883940);
    lnFillTimes.push_back(1442027368);
    lnFillTimes.push_back(1442169713);
    lnFillTimes.push_back(1442312599);
    lnFillTimes.push_back(1442453920);
    lnFillTimes.push_back(1442595578);
    lnFillTimes.push_back(1442737259);
    lnFillTimes.push_back(1442879000);
    lnFillTimes.push_back(1443021647);
  }

  if(dsNumber == 1) {
    lnFillTimes.push_back(1452519860);
    lnFillTimes.push_back(1452655867);
    lnFillTimes.push_back(1452791415);
    lnFillTimes.push_back(1453541032);
    lnFillTimes.push_back(1453670314);
    lnFillTimes.push_back(1453800137);
    lnFillTimes.push_back(1453929779);
    lnFillTimes.push_back(1453937178);
    lnFillTimes.push_back(1453998125);
    lnFillTimes.push_back(1454000591);
    lnFillTimes.push_back(1454002456);
    lnFillTimes.push_back(1454014971);
    lnFillTimes.push_back(1454107981);
    lnFillTimes.push_back(1454219392);
    lnFillTimes.push_back(1454332307);
    lnFillTimes.push_back(1454447212);
    lnFillTimes.push_back(1454559953);
    lnFillTimes.push_back(1454679496);
    lnFillTimes.push_back(1454769079);
    lnFillTimes.push_back(1454882301);
    lnFillTimes.push_back(1454946492);
    lnFillTimes.push_back(1454951907);
    lnFillTimes.push_back(1454954012);
    lnFillTimes.push_back(1454955958);
    lnFillTimes.push_back(1455225161);
    lnFillTimes.push_back(1455228289);
    lnFillTimes.push_back(1455440690);
    lnFillTimes.push_back(1455568585);
    lnFillTimes.push_back(1455696789);
    lnFillTimes.push_back(1455822149);
    lnFillTimes.push_back(1455952573);
    lnFillTimes.push_back(1456082166);
    lnFillTimes.push_back(1456206057);
    lnFillTimes.push_back(1456333236);
    lnFillTimes.push_back(1456460237);
    lnFillTimes.push_back(1456588495);
    lnFillTimes.push_back(1456717776);
    lnFillTimes.push_back(1456846882);
    lnFillTimes.push_back(1456943983);
    lnFillTimes.push_back(1456981934);
    lnFillTimes.push_back(1457110918);
    lnFillTimes.push_back(1457238098);
    lnFillTimes.push_back(1457365179);
    lnFillTimes.push_back(1457491997);
    lnFillTimes.push_back(1457619662);
    lnFillTimes.push_back(1457747884);
    lnFillTimes.push_back(1457874684);
    lnFillTimes.push_back(1458001143);
    lnFillTimes.push_back(1458130495);
    lnFillTimes.push_back(1458259402);
    lnFillTimes.push_back(1458387930);
    lnFillTimes.push_back(1458515657);
    lnFillTimes.push_back(1458639685);
    lnFillTimes.push_back(1458767518);
    lnFillTimes.push_back(1458896260);
    lnFillTimes.push_back(1459023862);
    lnFillTimes.push_back(1459150863);
    lnFillTimes.push_back(1459275989);
    lnFillTimes.push_back(1459402548);
    lnFillTimes.push_back(1459529663);
    lnFillTimes.push_back(1459654930);
    lnFillTimes.push_back(1459785212);
    lnFillTimes.push_back(1459912507);
    lnFillTimes.push_back(1460042708);
    lnFillTimes.push_back(1460169429);
    lnFillTimes.push_back(1460297779);
    lnFillTimes.push_back(1460426527);
    lnFillTimes.push_back(1460552742);
    lnFillTimes.push_back(1460678641);
    lnFillTimes.push_back(1460808916);
    lnFillTimes.push_back(1460939150);
    lnFillTimes.push_back(1461064619);
  }
}

double GetDCRraw(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05);
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05);
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05);
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05);
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05);
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05);
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05);
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05);
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05);
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05);
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05);
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05);
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05);
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05);
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05);
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05);
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05);
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05);
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05);
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05);
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05);
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05);
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05);
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05);
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05);
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05);
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05);
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05);
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05);
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05);
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05);
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05);
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05);
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05);
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05);
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05);
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05);
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05);
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05);
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05);
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05);
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05);
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05);
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05);
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05);
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05);
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  if(dsNumber == 1) {//updated 17 Aug 2016 using runs 11507-11592
    if(channel == 578) return nlcblrwfSlope-(4.620083E-05+trap4usMax*-1.306E-05) ;
    if(channel == 579) return nlcblrwfSlope-(1.163801E-05+trap4usMax*-1.304E-05) ;
    if(channel == 580) return nlcblrwfSlope-(2.604215E-05+trap4usMax*-1.291E-05) ;
    if(channel == 581) return nlcblrwfSlope-(2.994274E-05+trap4usMax*-1.290E-05) ;
    if(channel == 582) return nlcblrwfSlope-(5.175215E-05+trap4usMax*-1.304E-05) ;
    if(channel == 583) return nlcblrwfSlope-(1.902729E-05+trap4usMax*-1.302E-05) ;
    if(channel == 592) return nlcblrwfSlope-(5.145000E-05+trap4usMax*-1.297E-05) ;
    if(channel == 593) return nlcblrwfSlope-(4.298575E-05+trap4usMax*-1.297E-05) ;
    if(channel == 594) return nlcblrwfSlope-(7.298870E-05+trap4usMax*-1.353E-05) ;
    if(channel == 595) return nlcblrwfSlope-(6.178016E-06+trap4usMax*-1.351E-05) ;
    if(channel == 598) return nlcblrwfSlope-(3.497582E-05+trap4usMax*-1.298E-05) ;
    if(channel == 599) return nlcblrwfSlope-(9.553672E-06+trap4usMax*-1.298E-05) ;
    if(channel == 600) return nlcblrwfSlope-(5.649999E-05+trap4usMax*-1.287E-05) ;
    if(channel == 601) return nlcblrwfSlope-(5.482757E-05+trap4usMax*-1.288E-05) ;
    if(channel == 608) return nlcblrwfSlope-(4.586047E-05+trap4usMax*-1.284E-05) ;
    if(channel == 609) return nlcblrwfSlope-(1.016108E-05+trap4usMax*-1.282E-05) ;
    if(channel == 610) return nlcblrwfSlope-(2.209244E-05+trap4usMax*-1.298E-05) ;
    if(channel == 611) return nlcblrwfSlope-(3.159581E-05+trap4usMax*-1.299E-05) ;
    if(channel == 616) return nlcblrwfSlope-(7.722590E-05+trap4usMax*-1.316E-05) ;
    if(channel == 617) return nlcblrwfSlope-(6.568340E-05+trap4usMax*-1.309E-05) ;
    if(channel == 626) return nlcblrwfSlope-(1.786671E-05+trap4usMax*-1.266E-05) ;
    if(channel == 627) return nlcblrwfSlope-(8.580088E-06+trap4usMax*-1.264E-05) ;
    if(channel == 632) return nlcblrwfSlope-(9.481784E-05+trap4usMax*-1.307E-05) ;
    if(channel == 633) return nlcblrwfSlope-(6.909768E-06+trap4usMax*-1.304E-05) ;
    if(channel == 640) return nlcblrwfSlope-(7.977354E-05+trap4usMax*-1.287E-05) ;
    if(channel == 641) return nlcblrwfSlope-(3.680967E-05+trap4usMax*-1.287E-05) ;
    if(channel == 648) return nlcblrwfSlope-(9.035880E-05+trap4usMax*-1.312E-05) ;
    if(channel == 649) return nlcblrwfSlope-(6.199340E-06+trap4usMax*-1.308E-05) ;
    if(channel == 664) return nlcblrwfSlope-(5.873798E-05+trap4usMax*-1.321E-05) ;
    if(channel == 665) return nlcblrwfSlope-(2.416351E-05+trap4usMax*-1.320E-05) ;
    if(channel == 672) return nlcblrwfSlope-(2.815562E-05+trap4usMax*-1.328E-05) ;
    if(channel == 673) return nlcblrwfSlope-(4.809700E-05+trap4usMax*-1.334E-05) ;
    if(channel == 690) return nlcblrwfSlope-(4.148262E-05+trap4usMax*-1.308E-05) ;
    if(channel == 691) return nlcblrwfSlope-(1.779754E-05+trap4usMax*-1.310E-05) ;
    if(channel == 692) return nlcblrwfSlope-(4.737009E-05+trap4usMax*-1.309E-05) ;
    if(channel == 693) return nlcblrwfSlope-(2.395867E-05+trap4usMax*-1.312E-05) ;

    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCRraw(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}

double GetDCR85(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05)-1.404733E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05)-5.853623E-05;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05)-1.969379E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05)-6.711741E-05;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05)-1.273623E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05)-5.400800E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05)-7.336924E-05;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05)-4.684072E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05)-1.193772E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05)-5.622187E-05;
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05)-7.634619E-05;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05)-3.847147E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05)-9.726243E-05;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05)-4.177841E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05)-8.493529E-05;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05)-3.811032E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05)-4.101792E-05;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05)-2.216033E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05)-1.110490E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05)-5.854669E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05)-8.393719E-05;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05)-3.169238E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05)-1.104404E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05)-5.060284E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05)-1.857343E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05)-7.794238E-05;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05)-9.827360E-05;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05)-4.216143E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05)-8.512630E-05;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05)-4.005612E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05)-1.161883E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05)-5.772356E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05)-1.079951E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05)-5.054795E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05)-8.726088E-05;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05)-4.125616E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05)-1.252195E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05)-4.269257E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05)-7.771228E-05;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05)-2.347411E-05;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05)-7.504308E-05;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05)-3.671951E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05)-2.332299E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05)-4.664937E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05)-7.504187E-05;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05)-3.670170E-05;

    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }


  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trap4usMax*-1.328E-05)-8.966599E-05;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trap4usMax*-1.334E-05)-4.815561E-05;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trap4usMax*-1.308E-05)-6.733372E-05;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trap4usMax*-1.310E-05)-3.688012E-05;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trap4usMax*-1.309E-05)-8.342106E-05;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trap4usMax*-1.312E-05)-4.434964E-05;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trap4usMax*-1.306E-05)-1.178722E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trap4usMax*-1.304E-05)-4.880046E-05;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trap4usMax*-1.291E-05)-2.066446E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trap4usMax*-1.290E-05)-6.389357E-05;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trap4usMax*-1.304E-05)-8.078366E-05;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trap4usMax*-1.302E-05)-4.854558E-05; //average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trap4usMax*-1.297E-05)-1.025187E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trap4usMax*-1.296E-05)-5.053970E-05;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trap4usMax*-1.352E-05)-8.236275E-05;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trap4usMax*-1.351E-05)-4.854558E-05; //average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trap4usMax*-1.298E-05)-1.014487E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trap4usMax*-1.298E-05)-4.279017E-05;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trap4usMax*-1.287E-05)-6.607836E-05;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trap4usMax*-1.288E-05)-3.617333E-05;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trap4usMax*-1.284E-05)-1.825068E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trap4usMax*-1.282E-05)-6.448547E-05;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trap4usMax*-1.298E-05)-9.937699E-05;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trap4usMax*-1.299E-05)-4.632804E-05;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trap4usMax*-1.316E-05)-8.402205E-05;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trap4usMax*-1.310E-05)-1.073512E-04;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trap4usMax*-1.265E-05)-8.080291E-05;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trap4usMax*-1.264E-05)-4.500772E-05;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trap4usMax*-1.307E-05)-1.300685E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trap4usMax*-1.304E-05)-6.070530E-05;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trap4usMax*-1.287E-05)-7.355431E-05;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trap4usMax*-1.287E-05)-4.635975E-05;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trap4usMax*-1.311E-05)-6.962783E-05;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trap4usMax*-1.308E-05)-5.785846E-05;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trap4usMax*-1.321E-05)-1.088147E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trap4usMax*-1.320E-05)-4.178785E-05;


    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCR85(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}

double GetDCR90(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05)-1.019405E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05)-4.913911E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05)-1.279785E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05)-5.424818E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05)-1.092631E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05)-4.862555E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05)-7.087017E-05;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05)-3.625546E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05)-1.345640E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05)-7.199881E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05)-1.075847E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05)-4.179749E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05)-1.458138E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05)-6.535173E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05)-2.365925E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05)-9.628056E-05;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05)-1.227120E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05)-5.426306E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05)-1.125461E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05)-5.166456E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05)-1.513889E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05)-7.207791E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05)-1.396431E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05)-6.263415E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05)-1.099235E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05)-5.128737E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05)-1.573303E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05)-5.677784E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05)-1.042533E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05)-5.084706E-05;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05)-9.601873E-05;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05)-4.681626E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05)-2.948478E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05)-6.062255E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05)-9.675150E-05;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05)-4.822259E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05)-2.173713E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05)-8.310702E-05;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05)-2.461283E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05)-8.321697E-05;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05)-1.720172E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05)-6.872953E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05)-1.009845E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05)-5.797197E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05)-1.537598E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05)-7.108055E-05;
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }


  if(dsNumber == 1) {
    //params from 11507-11592 Updated 29 Aug 2016
	if(channel == 578) return nlcblrwfSlope-(4.620083E-05+trap4usMax*-1.306E-05)-1.414597E-04 ;
	if(channel == 579) return nlcblrwfSlope-(1.163801E-05+trap4usMax*-1.304E-05)-6.024276E-05 ;
	if(channel == 580) return nlcblrwfSlope-(2.604215E-05+trap4usMax*-1.291E-05)-2.600445E-04 ;
	if(channel == 581) return nlcblrwfSlope-(2.994274E-05+trap4usMax*-1.290E-05)-7.895203E-05 ;
	if(channel == 582) return nlcblrwfSlope-(5.175215E-05+trap4usMax*-1.304E-05)-8.452177E-05 ;
	if(channel == 583) return nlcblrwfSlope-(1.902729E-05+trap4usMax*-1.302E-05)-5.807816E-05 ;
	if(channel == 592) return nlcblrwfSlope-(5.145000E-05+trap4usMax*-1.297E-05)-1.158769E-04 ;
	if(channel == 593) return nlcblrwfSlope-(4.298575E-05+trap4usMax*-1.297E-05)-6.062252E-05 ;
	if(channel == 594) return nlcblrwfSlope-(7.298870E-05+trap4usMax*-1.353E-05)-1.216112E-04 ;
	if(channel == 595) return nlcblrwfSlope-(6.178016E-06+trap4usMax*-1.351E-05)-6.158457E-05 ;
	if(channel == 598) return nlcblrwfSlope-(3.497582E-05+trap4usMax*-1.298E-05)-1.198058E-04 ;
	if(channel == 599) return nlcblrwfSlope-(9.553672E-06+trap4usMax*-1.298E-05)-5.338186E-05 ;
	if(channel == 600) return nlcblrwfSlope-(5.649999E-05+trap4usMax*-1.287E-05)-9.150926E-05 ;
	if(channel == 601) return nlcblrwfSlope-(5.482757E-05+trap4usMax*-1.288E-05)-4.738173E-05 ;
	if(channel == 608) return nlcblrwfSlope-(4.586047E-05+trap4usMax*-1.284E-05)-2.255628E-04 ;
	if(channel == 609) return nlcblrwfSlope-(1.016108E-05+trap4usMax*-1.282E-05)-8.085080E-05 ;
	if(channel == 610) return nlcblrwfSlope-(2.209244E-05+trap4usMax*-1.298E-05)-1.144211E-04 ;
	if(channel == 611) return nlcblrwfSlope-(3.159581E-05+trap4usMax*-1.299E-05)-5.614616E-05 ;
	if(channel == 616) return nlcblrwfSlope-(7.722590E-05+trap4usMax*-1.316E-05)-1.175098E-04 ;
	if(channel == 617) return nlcblrwfSlope-(6.568340E-05+trap4usMax*-1.309E-05)-3.279914E-05 ;
	if(channel == 626) return nlcblrwfSlope-(1.786671E-05+trap4usMax*-1.266E-05)-1.285227E-04 ;
	if(channel == 627) return nlcblrwfSlope-(8.580088E-06+trap4usMax*-1.264E-05)-5.442935E-05 ;
	if(channel == 632) return nlcblrwfSlope-(9.481784E-05+trap4usMax*-1.307E-05)-1.709287E-04 ;
	if(channel == 633) return nlcblrwfSlope-(6.909768E-06+trap4usMax*-1.304E-05)-8.019256E-05 ;
	if(channel == 640) return nlcblrwfSlope-(7.977354E-05+trap4usMax*-1.287E-05)-9.866131E-05 ;
	if(channel == 641) return nlcblrwfSlope-(3.680967E-05+trap4usMax*-1.287E-05)-5.716358E-05 ;
	if(channel == 648) return nlcblrwfSlope-(9.035880E-05+trap4usMax*-1.312E-05)-1.874037E-04 ;
	if(channel == 649) return nlcblrwfSlope-(6.199340E-06+trap4usMax*-1.308E-05)-8.963580E-05 ;
	if(channel == 664) return nlcblrwfSlope-(5.873798E-05+trap4usMax*-1.321E-05)-1.320418E-04 ;
	if(channel == 665) return nlcblrwfSlope-(2.416351E-05+trap4usMax*-1.320E-05)-5.207803E-05 ;
	if(channel == 672) return nlcblrwfSlope-(2.815562E-05+trap4usMax*-1.328E-05)-1.103071E-04 ;
	if(channel == 673) return nlcblrwfSlope-(4.809700E-05+trap4usMax*-1.334E-05)-5.966680E-05 ;
	if(channel == 690) return nlcblrwfSlope-(4.148262E-05+trap4usMax*-1.308E-05)-8.145841E-05 ;
	if(channel == 691) return nlcblrwfSlope-(1.779754E-05+trap4usMax*-1.310E-05)-4.578676E-05 ;
	if(channel == 692) return nlcblrwfSlope-(4.737009E-05+trap4usMax*-1.309E-05)-1.106426E-04 ;
	if(channel == 693) return nlcblrwfSlope-(2.395867E-05+trap4usMax*-1.312E-05)-5.510341E-05 ;
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCR90(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}



double GetDCR95(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05)-1.438638E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05)-6.627588E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05)-1.777946E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05)-7.428328E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05)-1.487208E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05)-6.560241E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05)-1.310455E-04;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05)-6.411484E-05;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05)-1.709638E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05)-9.198958E-05;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05)-1.461563E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05)-5.850014E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05)-2.135500E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05)-9.026296E-05;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05)-3.157850E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05)-1.252092E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05)-1.609082E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05)-7.434378E-05;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05)-1.518572E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05)-6.920177E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05)-2.094821E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05)-9.346685E-05;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05)-1.868514E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05)-8.219029E-05;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05)-1.453889E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05)-6.690158E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05)-2.107852E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05)-8.136305E-05;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05)-1.473133E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05)-3.482776E-04;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05)-1.287283E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05)-6.244994E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05)-3.979238E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05)-8.256436E-05;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05)-1.323383E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05)-6.696894E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05)-3.899950E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05)-1.327352E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05)-3.258232E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05)-1.088782E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05)-2.461961E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05)-9.352369E-05;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05)-1.456008E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05)-7.560564E-05;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05)-2.075377E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05)-9.389924E-05;
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trap4usMax*-1.328E-05)-1.418052E-04;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trap4usMax*-1.334E-05)-8.184476E-05;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trap4usMax*-1.308E-05)-1.221969E-04;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trap4usMax*-1.310E-05)-6.775107E-05;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trap4usMax*-1.309E-05)-1.519825E-04;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trap4usMax*-1.312E-05)-7.302743E-05;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trap4usMax*-1.306E-05)-1.928345E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trap4usMax*-1.304E-05)-8.146071E-05;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trap4usMax*-1.291E-05)-3.260552E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trap4usMax*-1.290E-05)-1.044327E-04;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trap4usMax*-1.304E-05)-1.747793E-04;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trap4usMax*-1.302E-05)-8.834911E-05; //average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trap4usMax*-1.297E-05)-2.718900E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trap4usMax*-1.296E-05)-1.013371E-04;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trap4usMax*-1.352E-05)-1.443377E-04;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trap4usMax*-1.351E-05)-8.834911E-05; //average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trap4usMax*-1.298E-05)-1.715410E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trap4usMax*-1.298E-05)-7.209174E-05;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trap4usMax*-1.287E-05)-1.289802E-04;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trap4usMax*-1.288E-05)-6.895361E-05;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trap4usMax*-1.284E-05)-3.084515E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trap4usMax*-1.282E-05)-1.090966E-04;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trap4usMax*-1.298E-05)-1.758732E-04;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trap4usMax*-1.299E-05)-7.925371E-05;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trap4usMax*-1.316E-05)-1.587131E-04;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trap4usMax*-1.310E-05)-4.362600E-04;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trap4usMax*-1.265E-05)-1.412070E-04;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trap4usMax*-1.264E-05)-7.404711E-05;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trap4usMax*-1.307E-05)-2.461187E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trap4usMax*-1.304E-05)-1.122029E-04;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trap4usMax*-1.287E-05)-1.427142E-04;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trap4usMax*-1.287E-05)-7.623847E-05;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trap4usMax*-1.311E-05)-1.863428E-04;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trap4usMax*-1.308E-05)-1.220285E-04;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trap4usMax*-1.321E-05)-1.822770E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trap4usMax*-1.320E-05)-7.736052E-05;

    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCR95(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}


double GetDCR98(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05)-2.128491E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05)-9.009598E-05;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05)-2.478023E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05)-9.960543E-05;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05)-2.021753E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05)-8.815796E-05;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05)-5.932375E-04;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05)-2.631672E-04;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05)-2.128397E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05)-1.156453E-04;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05)-2.108129E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05)-8.500168E-05;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05)-3.363783E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05)-1.291137E-04;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05)-4.345738E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05)-1.682398E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05)-2.108711E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05)-1.006748E-04;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05)-1.994223E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05)-9.152186E-05;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05)-2.816064E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05)-1.197040E-04;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05)-2.513343E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05)-1.063627E-04;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05)-1.911555E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05)-8.660089E-05;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05)-2.776511E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05)-1.153726E-04;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05)-2.211304E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05)-8.877251E-04;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05)-1.702642E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05)-8.161982E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05)-5.350548E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05)-1.100618E-04;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05)-1.847819E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05)-9.712239E-05;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05)-6.408777E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05)-2.113617E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05)-4.178157E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05)-1.380808E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05)-3.561727E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05)-1.279477E-04;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05)-2.265576E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05)-1.004911E-04;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05)-2.745397E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05)-1.215824E-04;
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trap4usMax*-1.328E-05)-1.815976E-04 ;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trap4usMax*-1.334E-05)-1.102766E-04 ;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trap4usMax*-1.308E-05)-1.861127E-04 ;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trap4usMax*-1.310E-05)-1.024712E-04 ;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trap4usMax*-1.309E-05)-2.062761E-04 ;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trap4usMax*-1.312E-05)-9.459848E-05 ;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trap4usMax*-1.306E-05)-2.471341E-04 ;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trap4usMax*-1.304E-05)-1.087355E-04 ;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trap4usMax*-1.291E-05)-4.108245E-04 ;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trap4usMax*-1.290E-05)-1.306187E-04 ;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trap4usMax*-1.304E-05)-2.544365E-04 ;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trap4usMax*-1.302E-05)-1.316353E-04 ;//average for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trap4usMax*-1.297E-05)-4.974235E-04 ;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trap4usMax*-1.296E-05)-1.586067E-04 ;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trap4usMax*-1.352E-05)-1.914537E-04 ;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trap4usMax*-1.351E-05)-1.316353E-04 ;//average for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trap4usMax*-1.298E-05)-2.391778E-04 ;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trap4usMax*-1.298E-05)-9.659031E-05 ;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trap4usMax*-1.287E-05)-2.038157E-04 ;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trap4usMax*-1.288E-05)-1.089046E-04 ;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trap4usMax*-1.284E-05)-4.175150E-04 ;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trap4usMax*-1.282E-05)-1.439053E-04 ;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trap4usMax*-1.298E-05)-2.325430E-04 ;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trap4usMax*-1.299E-05)-1.038957E-04 ;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trap4usMax*-1.316E-05)-2.999593E-04 ;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trap4usMax*-1.310E-05)-6.829416E-04 ;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trap4usMax*-1.265E-05)-1.895566E-04 ;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trap4usMax*-1.264E-05)-9.793286E-05 ;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trap4usMax*-1.307E-05)-3.559619E-04 ;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trap4usMax*-1.304E-05)-1.610833E-04 ;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trap4usMax*-1.287E-05)-2.227665E-04 ;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trap4usMax*-1.287E-05)-1.031861E-04 ;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trap4usMax*-1.311E-05)-4.997305E-04 ;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trap4usMax*-1.308E-05)-1.948269E-04 ;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trap4usMax*-1.321E-05)-2.549099E-04 ;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trap4usMax*-1.320E-05)-1.133329E-04 ;


    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCR98(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}


double GetDCR99(int channel, double nlcblrwfSlope, double trap4usMax, int dsNumber)
{
  if(dsNumber == 0) {
    if(channel == 576) return nlcblrwfSlope-( 3.673385E-05+trap4usMax*-1.279E-05)-2.943105E-04;
    if(channel == 577) return nlcblrwfSlope-( 3.122161E-05+trap4usMax*-1.278E-05)-1.116883E-04;
    if(channel == 592) return nlcblrwfSlope-( 1.371969E-04+trap4usMax*-1.287E-05)-3.338696E-04;
    if(channel == 593) return nlcblrwfSlope-( 4.661555E-05+trap4usMax*-1.285E-05)-1.257624E-04;
    if(channel == 594) return nlcblrwfSlope-( 7.684589E-05+trap4usMax*-1.294E-05)-2.709597E-04;
    if(channel == 595) return nlcblrwfSlope-( 5.906154E-05+trap4usMax*-1.294E-05)-1.215498E-04;
    if(channel == 598) return nlcblrwfSlope-( 1.842168E-04+trap4usMax*-1.315E-05)-1.895059E-03;
    if(channel == 599) return nlcblrwfSlope-( 7.434566E-05+trap4usMax*-1.315E-05)-7.290344E-04;
    if(channel == 600) return nlcblrwfSlope-( 5.645349E-05+trap4usMax*-1.264E-05)-2.532124E-04;
    if(channel == 601) return nlcblrwfSlope-( 1.644013E-05+trap4usMax*-1.260E-05)-1.343034E-04;
    if(channel == 608) return nlcblrwfSlope-( 9.091576E-05+trap4usMax*-1.298E-05)-3.866514E-04;
    if(channel == 609) return nlcblrwfSlope-( 1.319922E-05+trap4usMax*-1.296E-05)-1.362231E-04;
    if(channel == 610) return nlcblrwfSlope-( 3.708218E-05+trap4usMax*-1.331E-05)-4.881074E-04;
    if(channel == 611) return nlcblrwfSlope-( 1.691204E-05+trap4usMax*-1.331E-05)-1.643671E-04;
    if(channel == 614) return nlcblrwfSlope-( 7.716150E-05+trap4usMax*-1.272E-05)-5.700044E-04;
    if(channel == 615) return nlcblrwfSlope-( 1.479029E-05+trap4usMax*-1.271E-05)-2.195198E-04;
    if(channel == 624) return nlcblrwfSlope-( 3.891078E-05+trap4usMax*-1.262E-05)-2.539550E-04;
    if(channel == 625) return nlcblrwfSlope-( 4.344515E-06+trap4usMax*-1.262E-05)-1.206225E-04;
    if(channel == 626) return nlcblrwfSlope-( 2.578407E-05+trap4usMax*-1.298E-05)-2.380339E-04;
    if(channel == 627) return nlcblrwfSlope-( 1.913080E-05+trap4usMax*-1.296E-05)-1.086605E-04;
    if(channel == 628) return nlcblrwfSlope-( 3.266857E-05+trap4usMax*-1.297E-05)-3.433920E-04;
    if(channel == 629) return nlcblrwfSlope-( 2.171356E-05+trap4usMax*-1.295E-05)-1.375756E-04;
    if(channel == 640) return nlcblrwfSlope-( 7.818519E-05+trap4usMax*-1.296E-05)-3.114656E-04;
    if(channel == 641) return nlcblrwfSlope-( 3.668705E-05+trap4usMax*-1.296E-05)-1.278493E-04;
    if(channel == 642) return nlcblrwfSlope-( 5.146852E-05+trap4usMax*-1.373E-05)-2.291227E-04;
    if(channel == 643) return nlcblrwfSlope-( 2.370897E-05+trap4usMax*-1.375E-05)-1.034860E-04;
    if(channel == 644) return nlcblrwfSlope-( 6.066896E-05+trap4usMax*-1.308E-05)-3.411170E-04;
    if(channel == 645) return nlcblrwfSlope-(-1.188733E-05+trap4usMax*-1.305E-05)-1.401627E-04;
    if(channel == 646) return nlcblrwfSlope-( 7.972405E-05+trap4usMax*-1.291E-05)-7.431032E-04;
    if(channel == 647) return nlcblrwfSlope-( 2.540927E-05+trap4usMax*-1.286E-05)-1.122732E-03;
    if(channel == 656) return nlcblrwfSlope-( 1.262340E-05+trap4usMax*-1.337E-05)-2.086850E-04;
    if(channel == 657) return nlcblrwfSlope-( 1.086301E-05+trap4usMax*-1.342E-05)-9.781184E-05;
    if(channel == 662) return nlcblrwfSlope-( 1.485632E-04+trap4usMax*-1.460E-05)-6.434840E-04;
    if(channel == 663) return nlcblrwfSlope-( 2.742351E-06+trap4usMax*-1.457E-05)-1.311618E-04;
    if(channel == 664) return nlcblrwfSlope-( 5.110392E-05+trap4usMax*-1.321E-05)-2.597014E-04;
    if(channel == 665) return nlcblrwfSlope-( 6.499560E-06+trap4usMax*-1.319E-05)-1.339893E-04;
    if(channel == 674) return nlcblrwfSlope-( 9.983398E-05+trap4usMax*-1.310E-05)-8.513833E-04;
    if(channel == 675) return nlcblrwfSlope-( 3.757473E-05+trap4usMax*-1.309E-05)-2.744186E-04;
    if(channel == 688) return nlcblrwfSlope-( 2.461617E-05+trap4usMax*-1.287E-05)-4.855323E-04;
    if(channel == 689) return nlcblrwfSlope-( 3.018505E-05+trap4usMax*-1.287E-05)-1.607455E-04;
    if(channel == 690) return nlcblrwfSlope-( 9.944427E-05+trap4usMax*-1.300E-05)-4.550779E-04;
    if(channel == 691) return nlcblrwfSlope-(-4.737995E-06+trap4usMax*-1.300E-05)-1.582763E-04;
    if(channel == 692) return nlcblrwfSlope-( 8.146529E-05+trap4usMax*-1.354E-05)-3.431355E-04;
    if(channel == 693) return nlcblrwfSlope-( 1.839544E-05+trap4usMax*-1.356E-05)-1.244729E-04;
    if(channel == 696) return nlcblrwfSlope-( 2.093490E-05+trap4usMax*-1.308E-05)-3.268483E-04;
    if(channel == 697) return nlcblrwfSlope-( 8.938457E-06+trap4usMax*-1.310E-05)-1.413939E-04;
    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  if(dsNumber == 1) {
    //params from 11507-11592
    if(channel == 672) return nlcblrwfSlope-(2.833847E-05+trap4usMax*-1.328E-05)-2.196966E-04;
    if(channel == 673) return nlcblrwfSlope-(4.734834E-05+trap4usMax*-1.334E-05)-1.301283E-04;
    if(channel == 690) return nlcblrwfSlope-(3.690947E-05+trap4usMax*-1.308E-05)-3.892826E-04;
    if(channel == 691) return nlcblrwfSlope-(1.605341E-05+trap4usMax*-1.310E-05)-1.785248E-04;
    if(channel == 692) return nlcblrwfSlope-(5.681803E-05+trap4usMax*-1.309E-05)-2.678930E-04;
    if(channel == 693) return nlcblrwfSlope-(2.453714E-05+trap4usMax*-1.312E-05)-1.176745E-04;
    if(channel == 578) return nlcblrwfSlope-(4.417989E-05+trap4usMax*-1.306E-05)-2.913426E-04;
    if(channel == 579) return nlcblrwfSlope-(1.106159E-05+trap4usMax*-1.304E-05)-1.342151E-04;
    if(channel == 580) return nlcblrwfSlope-(2.577260E-05+trap4usMax*-1.291E-05)-4.716845E-04;
    if(channel == 581) return nlcblrwfSlope-(3.041374E-05+trap4usMax*-1.290E-05)-1.497926E-04;
    if(channel == 582) return nlcblrwfSlope-(5.396620E-05+trap4usMax*-1.304E-05)-3.266592E-04;
    if(channel == 583) return nlcblrwfSlope-(1.977025E-05+trap4usMax*-1.302E-05)-1.795207E-04;//average value for all LG channels
    if(channel == 592) return nlcblrwfSlope-(6.803550E-05+trap4usMax*-1.297E-05)-7.344548E-04;
    if(channel == 593) return nlcblrwfSlope-(4.309706E-05+trap4usMax*-1.296E-05)-2.255520E-04;
    if(channel == 594) return nlcblrwfSlope-(4.643130E-05+trap4usMax*-1.352E-05)-2.328330E-04;
    if(channel == 595) return nlcblrwfSlope-(1.550955E-06+trap4usMax*-1.351E-05)-1.795207E-04;//average value for all LG channels
    if(channel == 598) return nlcblrwfSlope-(2.897476E-05+trap4usMax*-1.298E-05)-3.320390E-04;
    if(channel == 599) return nlcblrwfSlope-(9.612110E-06+trap4usMax*-1.298E-05)-1.228940E-04;
    if(channel == 600) return nlcblrwfSlope-(5.842046E-05+trap4usMax*-1.287E-05)-8.579336E-04;
    if(channel == 601) return nlcblrwfSlope-(5.372889E-05+trap4usMax*-1.288E-05)-2.916546E-04;
    if(channel == 608) return nlcblrwfSlope-(4.693324E-05+trap4usMax*-1.284E-05)-5.116411E-04;
    if(channel == 609) return nlcblrwfSlope-(1.055481E-05+trap4usMax*-1.282E-05)-1.752132E-04;
    if(channel == 610) return nlcblrwfSlope-(1.216405E-05+trap4usMax*-1.298E-05)-2.729985E-04;
    if(channel == 611) return nlcblrwfSlope-(2.884038E-05+trap4usMax*-1.299E-05)-1.225513E-04;
    if(channel == 616) return nlcblrwfSlope-(8.295462E-05+trap4usMax*-1.316E-05)-6.543178E-04;
    if(channel == 617) return nlcblrwfSlope-(7.294466E-05+trap4usMax*-1.310E-05)-1.052964E-03;
    if(channel == 626) return nlcblrwfSlope-(1.860959E-06+trap4usMax*-1.265E-05)-2.306604E-04;
    if(channel == 627) return nlcblrwfSlope-(7.000642E-06+trap4usMax*-1.264E-05)-1.191646E-04;
    if(channel == 632) return nlcblrwfSlope-(9.473908E-05+trap4usMax*-1.307E-05)-4.958992E-04;
    if(channel == 633) return nlcblrwfSlope-(6.874599E-06+trap4usMax*-1.304E-05)-1.981503E-04;
    if(channel == 640) return nlcblrwfSlope-(8.088126E-05+trap4usMax*-1.287E-05)-3.634485E-04;
    if(channel == 641) return nlcblrwfSlope-(3.667507E-05+trap4usMax*-1.287E-05)-1.330648E-04;
    if(channel == 648) return nlcblrwfSlope-(8.317390E-05+trap4usMax*-1.311E-05)-9.261722E-04;
    if(channel == 649) return nlcblrwfSlope-(4.360668E-06+trap4usMax*-1.308E-05)-2.763007E-04;
    if(channel == 664) return nlcblrwfSlope-(5.916704E-05+trap4usMax*-1.321E-05)-4.001758E-04;
    if(channel == 665) return nlcblrwfSlope-(2.373893E-05+trap4usMax*-1.320E-05)-1.605581E-04;

    if(trap4usMax == 0) return 0;
    return nlcblrwfSlope/trap4usMax;
  }

  else cout << "GetDCR99(): unknown dataset number DS" << dsNumber << endl;
  return 0;
}