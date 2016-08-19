// waveSkim.cc
// Grabs entry numbers of cleaned data
// and pulls the corresponding built/gatified data.
// Clint Wiseman, USC.

///////////////////////////////////////////////////////////
//                                                       //
// =================== DS1 REFERENCE =================== //
//                                                       //
///////////////////////////////////////////////////////////

// Channels range from 575 to 700
//
// int enabled[36] = {640, 641, 648, 649, 664, 665, 672, 673, 690, 691, 692, 693,
// 				   578, 579, 580, 581, 582, 583, 592, 593, 594, 595, 598, 599,
// 				   600, 601, 608, 609, 610, 611, 616, 617, 626, 627, 632, 633};

// int highGain[18] = {640, 648, 664, 672, 690, 692,
// 				   	578, 580, 582, 592, 594, 598,
// 				   	600, 608, 610, 616, 626, 632};

// corresponds to the highGain array
// string hgPos[18] = {"P2D3","P2D2","P3D4","P5D3","P6D4","P2D1",
// 					"P1D4","P1D3","P1D2","P7D4","P7D3","P7D2",
// 					"P7D1","P3D3","P3D2","P3D1","P6D3","P6D1"};

// string detEnr[18] = {"enr","enr","enr","enr","enr","nat",
//					 "enr","enr","enr","enr","enr","enr",
// 		  			 "nat","enr","enr","nat","enr","enr"};

// int lowGain[18] = {641, 649, 665, 673, 691, 693,
// 		    	   579, 581, 583, 593, 594, 599,
// 				   601, 609, 611, 617, 627, 633};
//

// Drawing ---
// energy : "trapENFCal"
// T/E vs E : "kvorrT/trapENFCal:trapENFCal"
// T/E : "kvorrT/trapENFCal"
// "channel:trapENFCal"
// "channel:run"
//
// Cutting ---
// energy : "trapENFCal>0 && trapENFCal<100"
// channel : "channel==648"
// HG : "gain==0"
// SD : "mH==1"
// IsEnr: "isEnr"
// DCBit : "!wfDCBits"
// LNFill : "!isLNFill"
// muon : "(dtmu_s < -0.2e-3 || dtmu_s > 1)"
// TE : "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1"
// TailMin : "trapETailMin < 0"

#include "iostream"

#include "TFile.h"
#include "TChain.h"
#include "TEntryList.h"

#include "MGTRun.hh"
#include "MGTEvent.hh"
#include "MJTChannelData.hh"

#include "GATDataSet.hh"
#include "GATMultiplicityProcessor.hh"

using namespace std;

void skimmer(char cutstr[500], char skimloc[500]);
void grabber(char wsout[500]);

int main()
{
	// =========== THE BIG CUT ===========
	char cutstr[1000];
	char HGCut[500] = "channel%2==0";	// only look at high gain
	// char SDCut[500] = "mH==1";
	char SDCut[500] = "1";
	// char energyCut[500] = "trapENFCal > 0 && trapENFCal < 30";
	// char energyCut[500] = "trapENFCal > 1590 && trapENFCal < 1595";
	char energyCut[500] = "trapENFCal > 1550 && trapENFCal < 1650";
	// char TrapETailMinCut[500] = "trapETailMin < 0";
	char TrapETailMinCut[500] = "1";
	char DCCut[500] = "!wfDCBits";
	char vetoCut[500] = "!muVeto";
	char burstCut[500] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";
	char runCut[500] = "run!=13312 && run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863";
	char segfaultCut[500] = "run!=9422. && run!=9423 && run!=9425 && run!=9426 && run!=9427";
	// char segfaultCut[500] = "1";
	sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && %s"
		,HGCut,SDCut,energyCut,TrapETailMinCut,DCCut,vetoCut,burstCut,runCut,segfaultCut);

	// =========== Skim file location ===========
	// char skimloc[500] = "/Users/wisecg/dev/datasets/ds1/*.root";
	char skimloc[500] = "$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/*.root";

	// =========== waveSkim output file ===========
	// char wsout[500] = "./output/waveSkim_DS1.root";
	char wsout[500] = "./output/waveSkim-1550-1650.root";


	// =========== Ready? Go! ===========
	skimmer(cutstr,skimloc);
	grabber(wsout);

	cout << "Cletus codes good.\n";
}

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// =========================== ENTRY LIST GENERATOR =========================== //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

void skimmer(char cutstr[500], char skimloc[500])
{
	TChain *skim = new TChain("skimTree");
	skim->Add(skimloc);
	long entries = skim->GetEntries();
	printf("Scanning skim files ... \nFound %li entries.\n",entries);

	// skim file branches.
	// May 2016 version.
	// get the list from a raw file with:
	// root[0] skimTree->GetListOfBranches()->Print()
	unsigned int gatrev=0,EventDC1Bits=0;
	int run=0,iEvent=0,mH=0,mL=0;
	double startTime=0,stopTime=0,sumEH=0,sumEL=0;
	vector<int> *iHit=0,*channel=0,*P=0,*D=0,*gain=0,*mageID=0,*detID=0,*dateMT=0,*muType=0;
	vector<double> *mAct_g=0,*tloc_s=0,*time_s=0,*timeMT=0,*trapECal=0,*trapENFCal=0,*aenorm=0,*kvorrT=0,*toe=0,*dcrSlope95=0,*dcrSlope99=0,*trapETailMin=0,*dtmu_s=0,*trap4usMax=0,*t150=0;
	vector<bool> *isEnr=0,*isNat=0,*isGood=0,*isLNFill=0,*badScaler=0,*muVeto=0;
	// vector<string> *detName;
	vector<unsigned int> *wfDCBits=0;
	skim->SetBranchAddress("gatrev",&gatrev);
	skim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
	skim->SetBranchAddress("run",&run);
	skim->SetBranchAddress("iEvent",&iEvent);
	skim->SetBranchAddress("mH",&mH);
	skim->SetBranchAddress("mL",&mL);
	skim->SetBranchAddress("t150",&t150);
	skim->SetBranchAddress("startTime",&startTime);
	skim->SetBranchAddress("stopTime",&stopTime);
	skim->SetBranchAddress("sumEH",&sumEH);
	skim->SetBranchAddress("sumEL",&sumEL);
	skim->SetBranchAddress("iHit",&iHit);
	skim->SetBranchAddress("channel",&channel);
	skim->SetBranchAddress("P",&P);
	skim->SetBranchAddress("D",&D);
	skim->SetBranchAddress("gain",&gain);
	skim->SetBranchAddress("mageID",&mageID);
	skim->SetBranchAddress("detID",&detID);
	skim->SetBranchAddress("dateMT",&dateMT);
	skim->SetBranchAddress("muType",&muType);
	skim->SetBranchAddress("mAct_g",&mAct_g);
	skim->SetBranchAddress("tloc_s",&tloc_s);
	skim->SetBranchAddress("time_s",&time_s);
	skim->SetBranchAddress("timeMT",&timeMT);
	skim->SetBranchAddress("trapECal",&trapECal);
	skim->SetBranchAddress("trapENFCal",&trapENFCal);
	skim->SetBranchAddress("trap4usMax",&trap4usMax);
	skim->SetBranchAddress("aenorm",&aenorm);
	skim->SetBranchAddress("kvorrT",&kvorrT);  // in GAT: trirt100nsft10nsMax / trapECal
	skim->SetBranchAddress("toe",&toe);
	skim->SetBranchAddress("dcrSlope95",&dcrSlope95);
	skim->SetBranchAddress("dcrSlope99",&dcrSlope99);
	skim->SetBranchAddress("trapETailMin",&trapETailMin);	// events with positive trapETailMin are removed
	skim->SetBranchAddress("dtmu_s",&dtmu_s);
	skim->SetBranchAddress("isEnr",&isEnr);
	skim->SetBranchAddress("isNat",&isNat);
	skim->SetBranchAddress("isGood",&isGood);
	skim->SetBranchAddress("isLNFill",&isLNFill);
	skim->SetBranchAddress("badScaler",&badScaler);
	skim->SetBranchAddress("muVeto",&muVeto);
	// skim->SetBranchAddress("detName",&detName);
	skim->SetBranchAddress("wfDCBits",&wfDCBits);

	// create TEntryList from skim file using draw command
 	skim->Draw(">>myElist",cutstr,"entrylist");
 	TEntryList *myElist = (TEntryList*)gDirectory->Get("myElist");
	skim->SetEntryList(myElist);
	Long64_t listEntries=myElist->GetN();
	printf("Created TEntryList from skim with %lli entries.\n",listEntries);

	// set output file
	Char_t OutputFile[200] = "./output/skimTEL.root";
	TFile *RootFile = new TFile(OutputFile, "RECREATE");
	TTree *list = new TTree("list","Run entry list");
	list->Branch("run",&run);
	list->Branch("iEvent",&iEvent);
	list->Branch("mH",&mH);	// pass these skim-file params into the waveSkim output file, because they're handy
	list->Branch("mL",&mL);
	list->Branch("dtmu_s",&dtmu_s);	// remember, this is a vector here.

	// loop over skim file with entry list
	int treenum=0;
	for (Long64_t el = 0; el < listEntries; el++)
	{
		Long64_t treeEntry = myElist->GetEntryAndTree(el,treenum);
		Long64_t chainEntry = treeEntry+skim->GetTreeOffset()[treenum];
		// printf("el=%lld, treeEntry=%lld, chainEntry=%lld, treenum=%d\n",el,treeEntry,chainEntry,treenum);
		skim->GetEntry(chainEntry);
		// if (el < 10) printf("el %lli  skim entry %lli  run %i  gat entry %i  multHG %i  time %.2f  sumEH %.2f\n",el,chainEntry,run,iEvent,mH,tloc_s->at(0),sumEH);
		list->Fill();
	}

	list->Write();
	RootFile->Close();
	// delete skim;
	// delete RootFile;
	// delete list;
}

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// ========================== GAT/BLT DATA COLLECTOR ========================== //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

void grabber(char wsout[500])
{
	printf("\nGrabbing GAT/BLT data ...\n");

	// loop over gatified/built data using "skimmer" output file,
	// store results in vectors, and make a list of unique run numbers
	// to loop over with the gat/blt data.
	vector<int> runList, entryList;
	vector<int> mHList, mLList;
	vector<double> dtmu_sList;
	TFile *skimList = new TFile("./output/skimTEL.root");
	TTree *list = (TTree*)skimList->Get("list");
	int sRun=0,iEvent=0;
	int mHtemp=0, mLtemp=0;
	vector <double> *dtmu_sTemp=0;
	double dtmu_sec=0;
	list->SetBranchAddress("run",&sRun);
	list->SetBranchAddress("iEvent",&iEvent);
	list->SetBranchAddress("mH",&mHtemp);
	list->SetBranchAddress("mL",&mLtemp);
	list->SetBranchAddress("dtmu_s",&dtmu_sTemp);	// remember, this is a vector

	long listEntries=list->GetEntries();
	for (long i = 0; i < listEntries; i++)
	{
		list->GetEntry(i);
		runList.push_back(sRun);
		entryList.push_back(iEvent);
		mHList.push_back(mHtemp);
		mLList.push_back(mLtemp);
		dtmu_sList.push_back(dtmu_sTemp->at(0));	// take the 0 entry to get a single value
	}
	printf("Skim list has %lli entries.\n",list->GetEntries());
	vector<int> uniqueRuns;
	for(int i = 0; i < (int)runList.size(); i++)
	{
		if (uniqueRuns.empty()) uniqueRuns.push_back(runList[i]);
		else {
			bool newRun = true;
			for (int j = 0; j < (int)uniqueRuns.size(); j++) {
				if (runList[i] == uniqueRuns[j]) {
					newRun = false;
					break;
				}
			}
			if (newRun) uniqueRuns.push_back(runList[i]);
		}
	}
	printf("Found %lu unique runs: \n",uniqueRuns.size());
	// for (int i = 0; i < (int)uniqueRuns.size(); i++) cout << uniqueRuns[i] << "  ";
	// cout << endl;
	skimList->Close();
	delete skimList;

	// create waveSkim output file
	TFile *waveSkim = new TFile(wsout,"RECREATE");
	TTree *builtSkim = new TTree("builtSkim","skimmed built data");
	TTree *gatSkim = new TTree("gatSkim","skimmed gatified data");

	// built branches (april 2016 format)
	MGTRun *bRun = new MGTRun();
	MGTEvent *event = new MGTEvent();
	// TClonesArray *channelData = new TClonesArray("MJTChannelData");
	builtSkim->Branch("run",&bRun);
	builtSkim->Branch("event","MGTEvent",&event,32000,0);
	// builtSkim->Branch("channelData",&channelData);	// causes too many "TProcessID" objects to be written to output

	// gatified branches (may 2016 format) + a few skim/custom branches
	// NOTE: difficulties with GATMIJ cause me to not include it, and pass in mH and mL from the skim file.
	int mL64=0,mH64=0;
	// GATMIJ* mij = new GATMIJ();
	unsigned int gatrev=0,EventDC1Bits=0;
	double run=0,startTime=0,stopTime=0;
	// vector<string> *detName=0;
	vector<int> *dateMT=0,*detID=0,*C=0,*P=0,*D=0,*mageID=0;
	vector<bool> *isEnr=0,*isNat=0;
	vector<unsigned int> *wfDCBits=0;
	vector<double> *rawWFMax=0,*rawWFMin=0,*timeMT=0,
		*blrwfFMR0p1=0, *sgswfFMR0p1=0, *blrwfFMR1=0, *sgswfFMR1=0, *blrwfFMR3=0, *sgswfFMR3=0, *blrwfFMR10=0,
		*sgswfFMR10=0,*blrwfFMR20=0, *sgswfFMR20=0, *blrwfFMR50=0, *sgswfFMR50=0, *blrwfFMR80=0, *sgswfFMR80=0,
		*blrwfFMR90=0, *sgswfFMR90=0,*blrwfFMR97=0, *sgswfFMR97=0, *blrwfFMR99=0, *sgswfFMR99=0, *sgswft0=0,
		*TSCurrent50nsMax=0, *TSCurrent100nsMax=0, *TSCurrent200nsMax=0, *RawWFblSlope=0, *RawWFblOffset=0,
		*RawWFblSlopeUnc=0, *RawWFblOffsetUnc=0, *RawWFblChi2=0, *RawWFftSlope=0, *RawWFftSlopeUnc=0,
		*RawWFftOffset=0, *RawWFftOffsetUnc=0, *RawWFftChi2=0, *RawWFwholeWFLinFitSlope=0,
		*RawWFwholeWFLinFitSlopeUnc=0, *RawWFwholeWFLinFitOffset=0, *RawWFwholeWFLinFitOffsetUnc=0,
		*RawWFwholeWFLinFitChi2=0, *trapENF55us=0, *trapENF70us=0, *trapENF85us=0, *trapENF100us=0,
		*trapENF115us=0, *trapENF130us=0, *trapBL=0, *trap500nsMax=0, *trapENF500nsrt=0,
		*trap1usMax=0, *trapENF1usrt=0, *trap2usMax=0, *trapENF2usrt=0, *trap4usMax=0, *trapENF4usrt=0,
		*trap6usMax=0, *trapENF6usrt=0, *trap8usMax=0, *trapENF8usrt=0, *trapE    =0, *trapEMin =0,
		*longGapTrapMax=0, *trap1ust0=0, *trapENF=0, *trapBL1us=0, *trapETailMin=0, *trirt50nsft10nsMax=0,
		*trirt100nsft10nsMax=0, *trirt200nsft10nsMax=0, *triFilMin=0, *trirt100nsft10nsIntegralW=0, *blrwfIntegralW=0,
		*smoothTrirt100nsft10nsMax=0, *energyCal=0, *trapECal =0, *trapEMinCal=0, *trapENFCal=0, *nlcblrwfSlope=0,
		*RawWFdcblSlope=0, *RawWFdcblOffset=0, *RawWFdcblSlopeUnc=0, *RawWFdcblOffsetUnc=0, *RawWFdcblChi2=0,
		*RawWFdcftSlope=0, *RawWFdcftOffset=0, *RawWFdcftSlopeUnc=0, *RawWFdcftOffsetUnc=0, *blrwf2ndDiffPeakValue=0,
		*blrwfNorm2ndDiffPeakValue=0, *blrwfSSF =0, *BLMoverMerr=0, *FTMoverMerr=0, *delOffset=0, *toe=0,
		*wfDC_Bit_0_SSSpikeBL=0, *wfDC_Bit_1_SSSpikePhys=0, *wfDC_Bit_2_EarlyTrigger=0, *wfDC_Bit_3_LateTrigger=0,
		*wfDC_Bit_4_PosSaturatedWFs=0, *wfDC_Bit_5_NegSaturatedWFs=0,
		*energy=0, *channel=0, *timestamp=0,
		*d2wfnoiseTagNorm=0,*d2wfDCPower=0,*d2wf0MHzTo2MHzPower=0,*d2wf2MHzTo10MHzPower=0,*d2wf10MHzTo20MHzPower=0,
		*d2wf30MHzTo35MHzPower=0,*d2wf48MHzTo50MHzPower=0,*d2wf0MHzTo50MHzPower=0;
	if(1){
		gatSkim->Branch("iEvent",&iEvent);	// added by clint
		// gatSkim->Branch("mij","GATMIJ",&mij,32000,0);
		// gatSkim->Branch("m",&m64);
		gatSkim->Branch("mL",&mL64);	// passed in from skim file
		gatSkim->Branch("mH",&mH64);
		gatSkim->Branch("dtmu_sec",&dtmu_sec);	// skim file, downgraded to a single value
		// gatSkim->Branch("i",&i);
		// gatSkim->Branch("j",&j);
		// gatSkim->Branch("iH",&iH);
		// gatSkim->Branch("jH",&jH);
		// gatSkim->Branch("iL",&iL);
		// gatSkim->Branch("jL",&jL);
		gatSkim->Branch("timeMT",&timeMT);
		gatSkim->Branch("dateMT",&dateMT);
		// gatSkim->Branch("detName",&detName);
		gatSkim->Branch("detID",&detID);
		gatSkim->Branch("C",&C);
		gatSkim->Branch("P",&P);
		gatSkim->Branch("D",&D);
		gatSkim->Branch("mageID",&mageID);
		gatSkim->Branch("isEnr",&isEnr);
		gatSkim->Branch("isNat",&isNat);
		gatSkim->Branch("rawWFMax",&rawWFMax);
		gatSkim->Branch("rawWFMin",&rawWFMin);
		gatSkim->Branch("blrwfFMR0p1",&blrwfFMR0p1);
		gatSkim->Branch("sgswfFMR0p1",&sgswfFMR0p1);
		gatSkim->Branch("blrwfFMR1",&blrwfFMR1);
		gatSkim->Branch("sgswfFMR1",&sgswfFMR1);
		gatSkim->Branch("blrwfFMR3",&blrwfFMR3);
		gatSkim->Branch("sgswfFMR3",&sgswfFMR3);
		gatSkim->Branch("blrwfFMR10",&blrwfFMR10);
		gatSkim->Branch("sgswfFMR10",&sgswfFMR10);
		gatSkim->Branch("blrwfFMR20",&blrwfFMR20);
		gatSkim->Branch("sgswfFMR20",&sgswfFMR20);
		gatSkim->Branch("blrwfFMR50",&blrwfFMR50);
		gatSkim->Branch("sgswfFMR50",&sgswfFMR50);
		gatSkim->Branch("blrwfFMR80",&blrwfFMR80);
		gatSkim->Branch("sgswfFMR80",&sgswfFMR80);
		gatSkim->Branch("blrwfFMR90",&blrwfFMR90);
		gatSkim->Branch("sgswfFMR90",&sgswfFMR90);
		gatSkim->Branch("blrwfFMR97",&blrwfFMR97);
		gatSkim->Branch("sgswfFMR97",&sgswfFMR97);
		gatSkim->Branch("blrwfFMR99",&blrwfFMR99);
		gatSkim->Branch("sgswfFMR99",&sgswfFMR99);
		gatSkim->Branch("sgswft0",&sgswft0);
		gatSkim->Branch("TSCurrent50nsMax",&TSCurrent50nsMax);
		gatSkim->Branch("TSCurrent100nsMax",&TSCurrent100nsMax);
		gatSkim->Branch("TSCurrent200nsMax",&TSCurrent200nsMax);
		gatSkim->Branch("RawWFblSlope",&RawWFblSlope);
		gatSkim->Branch("RawWFblOffset",&RawWFblOffset);
		gatSkim->Branch("RawWFblSlopeUnc",&RawWFblSlopeUnc);
		gatSkim->Branch("RawWFblOffsetUnc",&RawWFblOffsetUnc);
		gatSkim->Branch("RawWFblChi2",&RawWFblChi2);
		gatSkim->Branch("RawWFftSlope",&RawWFftSlope);
		gatSkim->Branch("RawWFftSlopeUnc",&RawWFftSlopeUnc);
		gatSkim->Branch("RawWFftOffset",&RawWFftOffset);
		gatSkim->Branch("RawWFftOffsetUnc",&RawWFftOffsetUnc);
		gatSkim->Branch("RawWFftChi2",&RawWFftChi2);
		gatSkim->Branch("RawWFwholeWFLinFitSlope",&RawWFwholeWFLinFitSlope);
		gatSkim->Branch("RawWFwholeWFLinFitSlopeUnc",&RawWFwholeWFLinFitSlopeUnc);
		gatSkim->Branch("RawWFwholeWFLinFitOffset",&RawWFwholeWFLinFitOffset);
		gatSkim->Branch("RawWFwholeWFLinFitOffsetUnc",&RawWFwholeWFLinFitOffsetUnc);
		gatSkim->Branch("RawWFwholeWFLinFitChi2",&RawWFwholeWFLinFitChi2);
		// gatSkim->Branch("trapENF55us",&trapENF55us);
		gatSkim->Branch("trapENF70us",&trapENF70us);
		// gatSkim->Branch("trapENF85us",&trapENF85us);
		gatSkim->Branch("trapENF100us",&trapENF100us);
		// gatSkim->Branch("trapENF115us",&trapENF115us);
		gatSkim->Branch("trapENF130us",&trapENF130us);
		// gatSkim->Branch("trapBL",&trapBL);
		gatSkim->Branch("trap500nsMax",&trap500nsMax);
		gatSkim->Branch("trapENF500nsrt",&trapENF500nsrt);
		gatSkim->Branch("trap1usMax",&trap1usMax);
		gatSkim->Branch("trapENF1usrt",&trapENF1usrt);
		gatSkim->Branch("trap2usMax",&trap2usMax);
		gatSkim->Branch("trapENF2usrt",&trapENF2usrt);
		gatSkim->Branch("trap4usMax",&trap4usMax);
		gatSkim->Branch("trapENF4usrt",&trapENF4usrt);
		gatSkim->Branch("trap6usMax",&trap6usMax);
		gatSkim->Branch("trapENF6usrt",&trapENF6usrt);
		gatSkim->Branch("trap8usMax",&trap8usMax);
		gatSkim->Branch("trapENF8usrt",&trapENF8usrt);
		gatSkim->Branch("trapE",&trapE);
		gatSkim->Branch("trapEMin",&trapEMin);
		gatSkim->Branch("longGapTrapMax",&longGapTrapMax);
		gatSkim->Branch("trap1ust0",&trap1ust0);
		gatSkim->Branch("trapENF",&trapENF);
		gatSkim->Branch("trapBL1us",&trapBL1us);
		gatSkim->Branch("trapETailMin",&trapETailMin);
		gatSkim->Branch("trirt50nsft10nsMax",&trirt50nsft10nsMax);
		gatSkim->Branch("trirt100nsft10nsMax",&trirt100nsft10nsMax);
		gatSkim->Branch("trirt200nsft10nsMax",&trirt200nsft10nsMax);
		gatSkim->Branch("triFilMin",&triFilMin);
		gatSkim->Branch("trirt100nsft10nsIntegralW",&trirt100nsft10nsIntegralW);
		gatSkim->Branch("blrwfIntegralW",&blrwfIntegralW);
		gatSkim->Branch("smoothTrirt100nsft10nsMax",&smoothTrirt100nsft10nsMax);
		gatSkim->Branch("energyCal",&energyCal);
		gatSkim->Branch("trapECal",&trapECal);
		gatSkim->Branch("trapEMinCal",&trapEMinCal);
		gatSkim->Branch("trapENFCal",&trapENFCal);
		gatSkim->Branch("nlcblrwfSlope",&nlcblrwfSlope);
		gatSkim->Branch("RawWFdcblSlope",&RawWFdcblSlope);
		gatSkim->Branch("RawWFdcblOffset",&RawWFdcblOffset);
		gatSkim->Branch("RawWFdcblSlopeUnc",&RawWFdcblSlopeUnc);
		gatSkim->Branch("RawWFdcblOffsetUnc",&RawWFdcblOffsetUnc);
		gatSkim->Branch("RawWFdcblChi2",&RawWFdcblChi2);
		gatSkim->Branch("RawWFdcftSlope",&RawWFdcftSlope);
		gatSkim->Branch("RawWFdcftOffset",&RawWFdcftOffset);
		gatSkim->Branch("RawWFdcftSlopeUnc",&RawWFdcftSlopeUnc);
		gatSkim->Branch("RawWFdcftOffsetUnc",&RawWFdcftOffsetUnc);
		gatSkim->Branch("blrwf2ndDiffPeakValue",&blrwf2ndDiffPeakValue);
		gatSkim->Branch("blrwfNorm2ndDiffPeakValue",&blrwfNorm2ndDiffPeakValue);
		gatSkim->Branch("blrwfSSF",&blrwfSSF);
		gatSkim->Branch("BLMoverMerr",&BLMoverMerr);
		gatSkim->Branch("FTMoverMerr",&FTMoverMerr);
		gatSkim->Branch("delOffset",&delOffset);
		gatSkim->Branch("toe",&toe);
		gatSkim->Branch("wfDCBits",&wfDCBits);
		gatSkim->Branch("wfDC_Bit_0_SSSpikeBL",&wfDC_Bit_0_SSSpikeBL);
		gatSkim->Branch("wfDC_Bit_1_SSSpikePhys",&wfDC_Bit_1_SSSpikePhys);
		gatSkim->Branch("wfDC_Bit_2_EarlyTrigger",&wfDC_Bit_2_EarlyTrigger);
		gatSkim->Branch("wfDC_Bit_3_LateTrigger",&wfDC_Bit_3_LateTrigger);
		gatSkim->Branch("wfDC_Bit_4_PosSaturatedWFs",&wfDC_Bit_4_PosSaturatedWFs);
		gatSkim->Branch("wfDC_Bit_5_NegSaturatedWFs",&wfDC_Bit_5_NegSaturatedWFs);
		gatSkim->Branch("run",&run);
		gatSkim->Branch("startTime",&startTime);
		gatSkim->Branch("stopTime",&stopTime);
		gatSkim->Branch("energy",&energy);
		gatSkim->Branch("channel",&channel);
		gatSkim->Branch("timestamp",&timestamp);
		gatSkim->Branch("gatrev",&gatrev);
		gatSkim->Branch("EventDC1Bits",&EventDC1Bits);
		gatSkim->Branch("d2wfnoiseTagNorm",&d2wfnoiseTagNorm);
		gatSkim->Branch("d2wfDCPower",&d2wfDCPower);
		gatSkim->Branch("d2wf0MHzTo2MHzPower",&d2wf0MHzTo2MHzPower);
		gatSkim->Branch("d2wf2MHzTo10MHzPower",&d2wf2MHzTo10MHzPower);
		gatSkim->Branch("d2wf10MHzTo20MHzPower",&d2wf10MHzTo20MHzPower);
		gatSkim->Branch("d2wf30MHzTo35MHzPower",&d2wf30MHzTo35MHzPower);
		gatSkim->Branch("d2wf48MHzTo50MHzPower",&d2wf48MHzTo50MHzPower);
		gatSkim->Branch("d2wf0MHzTo50MHzPower",&d2wf0MHzTo50MHzPower);
	}

	// Loop over unique runs
	int entryCount = 0;
	int runCount = (int)uniqueRuns.size();
	// cout << "0\n";
	for (int q = 0; q < (int)uniqueRuns.size(); q++)
	{
		// Load data
		int tempRun = uniqueRuns[q];
		// cout << "1  tempRun: " << tempRun << "\n";
		GATDataSet ds(tempRun);
		// GATDataSet *ds = new GATDataSet(tempRun);
		TChain *gatChain = ds.GetGatifiedChain();
		TChain *bltChain = ds.GetBuiltChain();
		// cout << "2\n";
		long entries = gatChain->GetEntries();
		printf("[%.2f%%] Scanning run %i, %li entries. ... ",((double)q/runCount)*100,uniqueRuns[q],entries);

		// set branches
		bltChain->SetBranchAddress("run",&bRun);
		bltChain->SetBranchAddress("event",&event);
		// bltChain->SetBranchAddress("channelData",&channelData);
		// cout << "3\n";
		if(1){
			// gatChain->SetBranchAddress("mij",&mij);	// this breaks the loop for some weird reason.
			gatChain->SetBranchAddress("timeMT",&timeMT);
			gatChain->SetBranchAddress("dateMT",&dateMT);
			// gatChain->SetBranchAddress("detName",&detName);
			gatChain->SetBranchAddress("detID",&detID);
			gatChain->SetBranchAddress("C",&C);
			gatChain->SetBranchAddress("P",&P);
			gatChain->SetBranchAddress("D",&D);
			gatChain->SetBranchAddress("mageID",&mageID);
			gatChain->SetBranchAddress("isEnr",&isEnr);
			gatChain->SetBranchAddress("isNat",&isNat);
			gatChain->SetBranchAddress("rawWFMax",&rawWFMax);
			gatChain->SetBranchAddress("rawWFMin",&rawWFMin);
			gatChain->SetBranchAddress("blrwfFMR0p1",&blrwfFMR0p1);
			gatChain->SetBranchAddress("sgswfFMR0p1",&sgswfFMR0p1);
			gatChain->SetBranchAddress("blrwfFMR1",&blrwfFMR1);
			gatChain->SetBranchAddress("sgswfFMR1",&sgswfFMR1);
			gatChain->SetBranchAddress("blrwfFMR3",&blrwfFMR3);
			gatChain->SetBranchAddress("sgswfFMR3",&sgswfFMR3);
			gatChain->SetBranchAddress("blrwfFMR10",&blrwfFMR10);
			gatChain->SetBranchAddress("sgswfFMR10",&sgswfFMR10);
			gatChain->SetBranchAddress("blrwfFMR20",&blrwfFMR20);
			gatChain->SetBranchAddress("sgswfFMR20",&sgswfFMR20);
			gatChain->SetBranchAddress("blrwfFMR50",&blrwfFMR50);
			gatChain->SetBranchAddress("sgswfFMR50",&sgswfFMR50);
			gatChain->SetBranchAddress("blrwfFMR80",&blrwfFMR80);
			gatChain->SetBranchAddress("sgswfFMR80",&sgswfFMR80);
			gatChain->SetBranchAddress("blrwfFMR90",&blrwfFMR90);
			gatChain->SetBranchAddress("sgswfFMR90",&sgswfFMR90);
			gatChain->SetBranchAddress("blrwfFMR97",&blrwfFMR97);
			gatChain->SetBranchAddress("sgswfFMR97",&sgswfFMR97);
			gatChain->SetBranchAddress("blrwfFMR99",&blrwfFMR99);
			gatChain->SetBranchAddress("sgswfFMR99",&sgswfFMR99);
			gatChain->SetBranchAddress("sgswft0",&sgswft0);
			gatChain->SetBranchAddress("TSCurrent50nsMax",&TSCurrent50nsMax);
			gatChain->SetBranchAddress("TSCurrent100nsMax",&TSCurrent100nsMax);
			gatChain->SetBranchAddress("TSCurrent200nsMax",&TSCurrent200nsMax);
			gatChain->SetBranchAddress("RawWFblSlope",&RawWFblSlope);
			gatChain->SetBranchAddress("RawWFblOffset",&RawWFblOffset);
			gatChain->SetBranchAddress("RawWFblSlopeUnc",&RawWFblSlopeUnc);
			gatChain->SetBranchAddress("RawWFblOffsetUnc",&RawWFblOffsetUnc);
			gatChain->SetBranchAddress("RawWFblChi2",&RawWFblChi2);
			gatChain->SetBranchAddress("RawWFftSlope",&RawWFftSlope);
			gatChain->SetBranchAddress("RawWFftSlopeUnc",&RawWFftSlopeUnc);
			gatChain->SetBranchAddress("RawWFftOffset",&RawWFftOffset);
			gatChain->SetBranchAddress("RawWFftOffsetUnc",&RawWFftOffsetUnc);
			gatChain->SetBranchAddress("RawWFftChi2",&RawWFftChi2);
			gatChain->SetBranchAddress("RawWFwholeWFLinFitSlope",&RawWFwholeWFLinFitSlope);
			gatChain->SetBranchAddress("RawWFwholeWFLinFitSlopeUnc",&RawWFwholeWFLinFitSlopeUnc);
			gatChain->SetBranchAddress("RawWFwholeWFLinFitOffset",&RawWFwholeWFLinFitOffset);
			gatChain->SetBranchAddress("RawWFwholeWFLinFitOffsetUnc",&RawWFwholeWFLinFitOffsetUnc);
			gatChain->SetBranchAddress("RawWFwholeWFLinFitChi2",&RawWFwholeWFLinFitChi2);
			// gatChain->SetBranchAddress("trapENF55us",&trapENF55us);
			gatChain->SetBranchAddress("trapENF70us",&trapENF70us);
			// gatChain->SetBranchAddress("trapENF85us",&trapENF85us);
			gatChain->SetBranchAddress("trapENF100us",&trapENF100us);
			// gatChain->SetBranchAddress("trapENF115us",&trapENF115us);
			gatChain->SetBranchAddress("trapENF130us",&trapENF130us);
			// gatChain->SetBranchAddress("trapBL",&trapBL);
			gatChain->SetBranchAddress("trap500nsMax",&trap500nsMax);
			gatChain->SetBranchAddress("trapENF500nsrt",&trapENF500nsrt);
			gatChain->SetBranchAddress("trap1usMax",&trap1usMax);
			gatChain->SetBranchAddress("trapENF1usrt",&trapENF1usrt);
			gatChain->SetBranchAddress("trap2usMax",&trap2usMax);
			gatChain->SetBranchAddress("trapENF2usrt",&trapENF2usrt);
			gatChain->SetBranchAddress("trap4usMax",&trap4usMax);
			gatChain->SetBranchAddress("trapENF4usrt",&trapENF4usrt);
			gatChain->SetBranchAddress("trap6usMax",&trap6usMax);
			gatChain->SetBranchAddress("trapENF6usrt",&trapENF6usrt);
			gatChain->SetBranchAddress("trap8usMax",&trap8usMax);
			gatChain->SetBranchAddress("trapENF8usrt",&trapENF8usrt);
			gatChain->SetBranchAddress("trapE",&trapE);
			gatChain->SetBranchAddress("trapEMin",&trapEMin);
			gatChain->SetBranchAddress("longGapTrapMax",&longGapTrapMax);
			gatChain->SetBranchAddress("trap1ust0",&trap1ust0);
			gatChain->SetBranchAddress("trapENF",&trapENF);
			gatChain->SetBranchAddress("trapBL1us",&trapBL1us);
			gatChain->SetBranchAddress("trapETailMin",&trapETailMin);
			gatChain->SetBranchAddress("trirt50nsft10nsMax",&trirt50nsft10nsMax);
			gatChain->SetBranchAddress("trirt100nsft10nsMax",&trirt100nsft10nsMax);
			gatChain->SetBranchAddress("trirt200nsft10nsMax",&trirt200nsft10nsMax);
			gatChain->SetBranchAddress("triFilMin",&triFilMin);
			gatChain->SetBranchAddress("trirt100nsft10nsIntegralW",&trirt100nsft10nsIntegralW);
			gatChain->SetBranchAddress("blrwfIntegralW",&blrwfIntegralW);
			gatChain->SetBranchAddress("smoothTrirt100nsft10nsMax",&smoothTrirt100nsft10nsMax);
			gatChain->SetBranchAddress("energyCal",&energyCal);
			gatChain->SetBranchAddress("trapECal",&trapECal);
			gatChain->SetBranchAddress("trapEMinCal",&trapEMinCal);
			gatChain->SetBranchAddress("trapENFCal",&trapENFCal);
			gatChain->SetBranchAddress("nlcblrwfSlope",&nlcblrwfSlope);
			gatChain->SetBranchAddress("RawWFdcblSlope",&RawWFdcblSlope);
			gatChain->SetBranchAddress("RawWFdcblOffset",&RawWFdcblOffset);
			gatChain->SetBranchAddress("RawWFdcblSlopeUnc",&RawWFdcblSlopeUnc);
			gatChain->SetBranchAddress("RawWFdcblOffsetUnc",&RawWFdcblOffsetUnc);
			gatChain->SetBranchAddress("RawWFdcblChi2",&RawWFdcblChi2);
			gatChain->SetBranchAddress("RawWFdcftSlope",&RawWFdcftSlope);
			gatChain->SetBranchAddress("RawWFdcftOffset",&RawWFdcftOffset);
			gatChain->SetBranchAddress("RawWFdcftSlopeUnc",&RawWFdcftSlopeUnc);
			gatChain->SetBranchAddress("RawWFdcftOffsetUnc",&RawWFdcftOffsetUnc);
			gatChain->SetBranchAddress("blrwf2ndDiffPeakValue",&blrwf2ndDiffPeakValue);
			gatChain->SetBranchAddress("blrwfNorm2ndDiffPeakValue",&blrwfNorm2ndDiffPeakValue);
			gatChain->SetBranchAddress("blrwfSSF",&blrwfSSF);
			gatChain->SetBranchAddress("BLMoverMerr",&BLMoverMerr);
			gatChain->SetBranchAddress("FTMoverMerr",&FTMoverMerr);
			gatChain->SetBranchAddress("delOffset",&delOffset);
			gatChain->SetBranchAddress("toe",&toe);
			gatChain->SetBranchAddress("wfDCBits",&wfDCBits);
			gatChain->SetBranchAddress("wfDC_Bit_0_SSSpikeBL",&wfDC_Bit_0_SSSpikeBL);
			gatChain->SetBranchAddress("wfDC_Bit_1_SSSpikePhys",&wfDC_Bit_1_SSSpikePhys);
			gatChain->SetBranchAddress("wfDC_Bit_2_EarlyTrigger",&wfDC_Bit_2_EarlyTrigger);
			gatChain->SetBranchAddress("wfDC_Bit_3_LateTrigger",&wfDC_Bit_3_LateTrigger);
			gatChain->SetBranchAddress("wfDC_Bit_4_PosSaturatedWFs",&wfDC_Bit_4_PosSaturatedWFs);
			gatChain->SetBranchAddress("wfDC_Bit_5_NegSaturatedWFs",&wfDC_Bit_5_NegSaturatedWFs);
			gatChain->SetBranchAddress("run",&run);
			gatChain->SetBranchAddress("startTime",&startTime);
			gatChain->SetBranchAddress("stopTime",&stopTime);
			gatChain->SetBranchAddress("energy",&energy);
			gatChain->SetBranchAddress("channel",&channel);
			gatChain->SetBranchAddress("timestamp",&timestamp);
			gatChain->SetBranchAddress("gatrev",&gatrev);
			gatChain->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
			gatChain->SetBranchAddress("d2wf0MHzTo2MHzPower",&d2wf0MHzTo2MHzPower);
			gatChain->SetBranchAddress("d2wf2MHzTo10MHzPower",&d2wf2MHzTo10MHzPower);
			gatChain->SetBranchAddress("d2wf10MHzTo20MHzPower",&d2wf10MHzTo20MHzPower);
			gatChain->SetBranchAddress("d2wf30MHzTo35MHzPower",&d2wf30MHzTo35MHzPower);
			gatChain->SetBranchAddress("d2wf48MHzTo50MHzPower",&d2wf48MHzTo50MHzPower);
			gatChain->SetBranchAddress("d2wf0MHzTo50MHzPower",&d2wf0MHzTo50MHzPower);
		}
		// cout << "4\n";

		// loop over entries for this run.
		int count = 0;
		for (int i = 0; i < (int)entryList.size(); i++)
		{
			if (runList[i] == uniqueRuns[q])
			{
				gatChain->GetEntry(entryList[i]);
				bltChain->GetEntry(entryList[i]);

				// Get variables from "skimmer"
				mH64 = mHList[i];
				mL64 = mLList[i];
				dtmu_sec = dtmu_sList[i];
				iEvent = entryList[i];

				// printf("run %.0f  skim entry %i  gat entry %lld  mH %i  dtmu_sec %.2f\n"
					// ,run,entryList[i],gatChain->GetEntryNumber(entryList[i]),mH64,dtmu_sec);

				count++;
				entryCount++;
				gatSkim->Fill();
				builtSkim->Fill();
			}
		}
		printf("Filled waveSkim with %i entries.\n",count);
		// delete gatChain;
		// delete bltChain;
		// delete ds;
	}
	printf("\nwaveSkim output has %i entries total.\n",entryCount);

	gatSkim->Write();
	builtSkim->Write();
	waveSkim->Close();
}