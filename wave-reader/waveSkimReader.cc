// waveSkimReader.cc
// Clint Wiseman, USC

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

#include "iostream"
#include <vector>

#include "TH1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTreeFormula.h"
#include "TLegend.h"

#include "MGTRun.hh"
#include "MGTEvent.hh"
#include "MJTChannelData.hh"
#include "MGTWaveform.hh"

#include "GATWaveformBrowser.hh"

using namespace std;

void waveSkimStats(char input[500]);
void plotter(char input[500]);
void wavePlots(char input[500]);
void entryPlots(char input[500]);

int main(int argc, char** argv)
{
	// char input[500] = "waveSkim_DS1LowE.root";
	// char input[500] = "waveSkim_DEP_DS1.root";
	char input[500] = "waveSkim_DS1.root";

	// waveSkimStats(input);
	// plotter(input);
	wavePlots(input);
	// entryPlots(input);
}

void waveSkimStats(char input[500])
{
	// initialize
	TFile *f = new TFile(input);
	// TTree *b = (TTree*)f->Get("builtSkim");
	TTree *g = (TTree*)f->Get("gatSkim");
	long entries = g->GetEntries();

	int mH64=0;
	vector<int> *P=0,*D=0;
	vector<double> *channel=0, *trapENFCal=0;
	g->SetBranchAddress("P",&P);
	g->SetBranchAddress("D",&D);
	g->SetBranchAddress("channel",&channel);
	g->SetBranchAddress("trapENFCal",&trapENFCal);
	g->SetBranchAddress("mH",&mH64);

	int cts[7][5] = {{0}};
	int highM = 0;
	int totCts = 0;
	for (long i = 0; i < entries; i++)
	{
		g->GetEntry(i);

		int chans = (int)channel->size();
		if (mH64 > 1) highM++;

		for (int j = 0; j < chans; j++)
		{
			int pos=P->at(j);
			int det=D->at(j);
			double en = trapENFCal->at(j);
			int chan = channel->at(j);

			if (en > 2 && chan%2==0) {
				cts[pos-1][det-1]++;
				totCts++;
			}

			if (chans > 1) printf("%li  vec size %i  j %i  chan %i  en %.2f  mH %i\n",i,chans,j,chan,en,mH64);
		}
	}

	printf("============ waveSkim Stats ============\n");
	printf("Total entries: %li.  Entries with mult > 1: %i\n",entries,highM);
	printf("Total detector counts > 2 keV: %i\n",totCts);
	printf("P1D1 %-6i  P1D2 %-6i  P1D3 %-6i  P1D4 %-6i\n",cts[0][0],cts[0][1],cts[0][2],cts[0][3]);
	printf("P2D1 %-6i  P2D2 %-6i  P2D3 %-6i  P2D4 %-6i\n",cts[1][0],cts[1][1],cts[1][2],cts[1][3]);
	printf("P3D1 %-6i  P3D2 %-6i  P3D3 %-6i  P3D4 %-6i\n",cts[2][0],cts[2][1],cts[2][2],cts[2][3]);
	printf("P4D1 %-6i  P4D2 %-6i  P4D3 %-6i  P4D4 %-6i  P4D5 %-6i\n",cts[3][0],cts[3][1],cts[3][2],cts[3][3],cts[3][4]);
	printf("P5D1 %-6i  P5D2 %-6i  P5D3 %-6i  P5D4 %-6i\n",cts[4][0],cts[4][1],cts[4][2],cts[4][3]);
	printf("P6D1 %-6i  P6D2 %-6i  P6D3 %-6i  P6D4 %-6i\n",cts[5][0],cts[5][1],cts[5][2],cts[5][3]);
	printf("P7D1 %-6i  P7D2 %-6i  P7D3 %-6i  P7D4 %-6i\n",cts[6][0],cts[6][1],cts[6][2],cts[6][3]);
}

void plotter(char input[500])
{
	// initialize
	TFile *f = new TFile(input);
	TTree *gatSkim = (TTree*)f->Get("gatSkim");
	TTree *builtSkim = (TTree*)f->Get("builtSkim");
	builtSkim->AddFriend(gatSkim);

	// built branches (april 2016 format)
	MGTRun *bRun = new MGTRun();
	MGTEvent *event = new MGTEvent();
	// TClonesArray *channelData = new TClonesArray("MJTChannelData");
	builtSkim->SetBranchAddress("run",&bRun);
	builtSkim->SetBranchAddress("event",&event);
	// builtSkim->SetBranchAddress("channelData",&channelData);	// causes too many "TProcessID" objects to be written to output

	// gatified branches (may 2016 format) + a few additions by Clint from waveSkim
	int iEvent=0;	// gat file entry number
	double dtmu_sec=0;	// time since last muon
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
		gatSkim->SetBranchAddress("iEvent",&iEvent);	// added by clint
		gatSkim->SetBranchAddress("dtmu_sec",&dtmu_sec);	// skim file, downgraded to a single value
		// gatSkim->SetBranchAddress("mij","GATMIJ",&mij,32000,0);
		// gatSkim->SetBranchAddress("m",&m64);
		gatSkim->SetBranchAddress("mL",&mL64);	// passed in from skim file
		gatSkim->SetBranchAddress("mH",&mH64);
		// gatSkim->SetBranchAddress("i",&i);
		// gatSkim->SetBranchAddress("j",&j);
		// gatSkim->SetBranchAddress("iH",&iH);
		// gatSkim->SetBranchAddress("jH",&jH);
		// gatSkim->SetBranchAddress("iL",&iL);
		// gatSkim->SetBranchAddress("jL",&jL);
		gatSkim->SetBranchAddress("timeMT",&timeMT);
		gatSkim->SetBranchAddress("dateMT",&dateMT);
		// gatSkim->SetBranchAddress("detName",&detName);
		gatSkim->SetBranchAddress("detID",&detID);
		gatSkim->SetBranchAddress("C",&C);
		gatSkim->SetBranchAddress("P",&P);
		gatSkim->SetBranchAddress("D",&D);
		gatSkim->SetBranchAddress("mageID",&mageID);
		gatSkim->SetBranchAddress("isEnr",&isEnr);
		gatSkim->SetBranchAddress("isNat",&isNat);
		gatSkim->SetBranchAddress("rawWFMax",&rawWFMax);
		gatSkim->SetBranchAddress("rawWFMin",&rawWFMin);
		gatSkim->SetBranchAddress("blrwfFMR0p1",&blrwfFMR0p1);
		gatSkim->SetBranchAddress("sgswfFMR0p1",&sgswfFMR0p1);
		gatSkim->SetBranchAddress("blrwfFMR1",&blrwfFMR1);
		gatSkim->SetBranchAddress("sgswfFMR1",&sgswfFMR1);
		gatSkim->SetBranchAddress("blrwfFMR3",&blrwfFMR3);
		gatSkim->SetBranchAddress("sgswfFMR3",&sgswfFMR3);
		gatSkim->SetBranchAddress("blrwfFMR10",&blrwfFMR10);
		gatSkim->SetBranchAddress("sgswfFMR10",&sgswfFMR10);
		gatSkim->SetBranchAddress("blrwfFMR20",&blrwfFMR20);
		gatSkim->SetBranchAddress("sgswfFMR20",&sgswfFMR20);
		gatSkim->SetBranchAddress("blrwfFMR50",&blrwfFMR50);
		gatSkim->SetBranchAddress("sgswfFMR50",&sgswfFMR50);
		gatSkim->SetBranchAddress("blrwfFMR80",&blrwfFMR80);
		gatSkim->SetBranchAddress("sgswfFMR80",&sgswfFMR80);
		gatSkim->SetBranchAddress("blrwfFMR90",&blrwfFMR90);
		gatSkim->SetBranchAddress("sgswfFMR90",&sgswfFMR90);
		gatSkim->SetBranchAddress("blrwfFMR97",&blrwfFMR97);
		gatSkim->SetBranchAddress("sgswfFMR97",&sgswfFMR97);
		gatSkim->SetBranchAddress("blrwfFMR99",&blrwfFMR99);
		gatSkim->SetBranchAddress("sgswfFMR99",&sgswfFMR99);
		gatSkim->SetBranchAddress("sgswft0",&sgswft0);
		gatSkim->SetBranchAddress("TSCurrent50nsMax",&TSCurrent50nsMax);
		gatSkim->SetBranchAddress("TSCurrent100nsMax",&TSCurrent100nsMax);
		gatSkim->SetBranchAddress("TSCurrent200nsMax",&TSCurrent200nsMax);
		gatSkim->SetBranchAddress("RawWFblSlope",&RawWFblSlope);
		gatSkim->SetBranchAddress("RawWFblOffset",&RawWFblOffset);
		gatSkim->SetBranchAddress("RawWFblSlopeUnc",&RawWFblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFblOffsetUnc",&RawWFblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFblChi2",&RawWFblChi2);
		gatSkim->SetBranchAddress("RawWFftSlope",&RawWFftSlope);
		gatSkim->SetBranchAddress("RawWFftSlopeUnc",&RawWFftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFftOffset",&RawWFftOffset);
		gatSkim->SetBranchAddress("RawWFftOffsetUnc",&RawWFftOffsetUnc);
		gatSkim->SetBranchAddress("RawWFftChi2",&RawWFftChi2);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlope",&RawWFwholeWFLinFitSlope);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlopeUnc",&RawWFwholeWFLinFitSlopeUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffset",&RawWFwholeWFLinFitOffset);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffsetUnc",&RawWFwholeWFLinFitOffsetUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitChi2",&RawWFwholeWFLinFitChi2);
		gatSkim->SetBranchAddress("trapENF55us",&trapENF55us);
		gatSkim->SetBranchAddress("trapENF70us",&trapENF70us);
		gatSkim->SetBranchAddress("trapENF85us",&trapENF85us);
		gatSkim->SetBranchAddress("trapENF100us",&trapENF100us);
		gatSkim->SetBranchAddress("trapENF115us",&trapENF115us);
		gatSkim->SetBranchAddress("trapENF130us",&trapENF130us);
		gatSkim->SetBranchAddress("trapBL",&trapBL);
		gatSkim->SetBranchAddress("trap500nsMax",&trap500nsMax);
		gatSkim->SetBranchAddress("trapENF500nsrt",&trapENF500nsrt);
		gatSkim->SetBranchAddress("trap1usMax",&trap1usMax);
		gatSkim->SetBranchAddress("trapENF1usrt",&trapENF1usrt);
		gatSkim->SetBranchAddress("trap2usMax",&trap2usMax);
		gatSkim->SetBranchAddress("trapENF2usrt",&trapENF2usrt);
		gatSkim->SetBranchAddress("trap4usMax",&trap4usMax);
		gatSkim->SetBranchAddress("trapENF4usrt",&trapENF4usrt);
		gatSkim->SetBranchAddress("trap6usMax",&trap6usMax);
		gatSkim->SetBranchAddress("trapENF6usrt",&trapENF6usrt);
		gatSkim->SetBranchAddress("trap8usMax",&trap8usMax);
		gatSkim->SetBranchAddress("trapENF8usrt",&trapENF8usrt);
		gatSkim->SetBranchAddress("trapE",&trapE);
		gatSkim->SetBranchAddress("trapEMin",&trapEMin);
		gatSkim->SetBranchAddress("longGapTrapMax",&longGapTrapMax);
		gatSkim->SetBranchAddress("trap1ust0",&trap1ust0);
		gatSkim->SetBranchAddress("trapENF",&trapENF);
		gatSkim->SetBranchAddress("trapBL1us",&trapBL1us);
		gatSkim->SetBranchAddress("trapETailMin",&trapETailMin);
		gatSkim->SetBranchAddress("trirt50nsft10nsMax",&trirt50nsft10nsMax);
		gatSkim->SetBranchAddress("trirt100nsft10nsMax",&trirt100nsft10nsMax);
		gatSkim->SetBranchAddress("trirt200nsft10nsMax",&trirt200nsft10nsMax);
		gatSkim->SetBranchAddress("triFilMin",&triFilMin);
		gatSkim->SetBranchAddress("trirt100nsft10nsIntegralW",&trirt100nsft10nsIntegralW);
		gatSkim->SetBranchAddress("blrwfIntegralW",&blrwfIntegralW);
		gatSkim->SetBranchAddress("smoothTrirt100nsft10nsMax",&smoothTrirt100nsft10nsMax);
		gatSkim->SetBranchAddress("energyCal",&energyCal);
		gatSkim->SetBranchAddress("trapECal",&trapECal);
		gatSkim->SetBranchAddress("trapEMinCal",&trapEMinCal);
		gatSkim->SetBranchAddress("trapENFCal",&trapENFCal);
		gatSkim->SetBranchAddress("nlcblrwfSlope",&nlcblrwfSlope);
		gatSkim->SetBranchAddress("RawWFdcblSlope",&RawWFdcblSlope);
		gatSkim->SetBranchAddress("RawWFdcblOffset",&RawWFdcblOffset);
		gatSkim->SetBranchAddress("RawWFdcblSlopeUnc",&RawWFdcblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcblOffsetUnc",&RawWFdcblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFdcblChi2",&RawWFdcblChi2);
		gatSkim->SetBranchAddress("RawWFdcftSlope",&RawWFdcftSlope);
		gatSkim->SetBranchAddress("RawWFdcftOffset",&RawWFdcftOffset);
		gatSkim->SetBranchAddress("RawWFdcftSlopeUnc",&RawWFdcftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcftOffsetUnc",&RawWFdcftOffsetUnc);
		gatSkim->SetBranchAddress("blrwf2ndDiffPeakValue",&blrwf2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfNorm2ndDiffPeakValue",&blrwfNorm2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfSSF",&blrwfSSF);
		gatSkim->SetBranchAddress("BLMoverMerr",&BLMoverMerr);
		gatSkim->SetBranchAddress("FTMoverMerr",&FTMoverMerr);
		gatSkim->SetBranchAddress("delOffset",&delOffset);
		gatSkim->SetBranchAddress("toe",&toe);
		gatSkim->SetBranchAddress("wfDCBits",&wfDCBits);
		gatSkim->SetBranchAddress("wfDC_Bit_0_SSSpikeBL",&wfDC_Bit_0_SSSpikeBL);
		gatSkim->SetBranchAddress("wfDC_Bit_1_SSSpikePhys",&wfDC_Bit_1_SSSpikePhys);
		gatSkim->SetBranchAddress("wfDC_Bit_2_EarlyTrigger",&wfDC_Bit_2_EarlyTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_3_LateTrigger",&wfDC_Bit_3_LateTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_4_PosSaturatedWFs",&wfDC_Bit_4_PosSaturatedWFs);
		gatSkim->SetBranchAddress("wfDC_Bit_5_NegSaturatedWFs",&wfDC_Bit_5_NegSaturatedWFs);
		gatSkim->SetBranchAddress("run",&run);
		gatSkim->SetBranchAddress("startTime",&startTime);
		gatSkim->SetBranchAddress("stopTime",&stopTime);
		gatSkim->SetBranchAddress("energy",&energy);
		gatSkim->SetBranchAddress("channel",&channel);
		gatSkim->SetBranchAddress("timestamp",&timestamp);
		gatSkim->SetBranchAddress("gatrev",&gatrev);
		gatSkim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
		gatSkim->SetBranchAddress("d2wfnoiseTagNorm",&d2wfnoiseTagNorm);
		gatSkim->SetBranchAddress("d2wfDCPower",&d2wfDCPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo2MHzPower",&d2wf0MHzTo2MHzPower);
		gatSkim->SetBranchAddress("d2wf2MHzTo10MHzPower",&d2wf2MHzTo10MHzPower);
		gatSkim->SetBranchAddress("d2wf10MHzTo20MHzPower",&d2wf10MHzTo20MHzPower);
		gatSkim->SetBranchAddress("d2wf30MHzTo35MHzPower",&d2wf30MHzTo35MHzPower);
		gatSkim->SetBranchAddress("d2wf48MHzTo50MHzPower",&d2wf48MHzTo50MHzPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo50MHzPower",&d2wf0MHzTo50MHzPower);
	}

	gROOT->ProcessLine(".x /Users/wisecg/dev/MJDClintPlotStyle.C");

	// drawing
	char drawstr[500];
	// char eDraw[50] = "trapENFCal";
	// char t_eDraw[50] = "kvorrT/trapENFCal:trapENFCal";
	// char t_eAccDraw[50] = "kvorrT/trapENFCal";
	// char e_chanDraw[50] = "channel:trapENFCal";
	// char chanVRunDraw[50] = "channel:run";

 	// cutting
	char cutstr[500];
	char HGCut[500] = "channel%2==0";
	char EFloor[500] = "trapENFCal>=2";
	// char energyCut[50] = "trapENFCal>0 && trapENFCal<100";
	// char channelCut[50] = "channel==648";
	// char HGCut[50] = "gain==0";
	// char SDCut[50] = "mH==1";
	// char IsEnr[50] = "isEnr";
	// char DCBitCut[50] = "!wfDCBits";
	// char LNFillCut[50] = "!isLNFill";
	// char muonCut[50] = "(dtmu_s < -0.2e-3 || dtmu_s > 1)";
	// char TECut[100] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";
	// char TMCut[50] = "trapETailMin < 0";

	// binning
	// kris's binning: 15 bins/kev.  Ex: 450,0,30
	// clint's binning: 2 bins/kev.  Ex: 60,0,30
	// alternate clint: 8 bins/kev.  Ex: 240,0,30

	///////////////////////////////////////////////////////////
	//                                                       //
	// ===================== 1-D PLOTS ===================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	// Verify the energy spectrum of the waveSkim file.
	if(0)
	{
		TCanvas *c1 = new TCanvas("c1", "Bob Ross's Canvas",800,600);
		// c1->SetLogy();

		TH1D *h1 = new TH1D("h1","",450,0,30);
		sprintf(drawstr,"trapENFCal");
		sprintf(cutstr,"%s && %s",HGCut,EFloor);
		gatSkim->Project("h1", drawstr, cutstr);

		h1->SetLineColor(kBlack);
		h1->GetXaxis()->SetTitle(drawstr);
		h1->GetYaxis()->SetTitle("Counts (15 bins / keV)");
		h1->Draw();

		c1->Update();
		c1->Print("output/LowE_Spectrum.pdf");
	}

	// Compare enriched and natural spectra.  Should be cleaner ...
	// Don't plot the waveSkim output below 2 keV.
	if(0)
	{
		TCanvas *c2 = new TCanvas("c2", "Bob Ross's Canvas",800,600);
		// c1->SetLogy();

		TH1D *h2 = new TH1D("h2","",240,0,30);
		sprintf(drawstr,"trapENFCal");
		sprintf(cutstr,"%s && %s && isEnr",HGCut,EFloor);
		gatSkim->Project("h2", drawstr, cutstr);

		TH1D *h3 = new TH1D("h3","",240,0,30);
		sprintf(drawstr,"trapENFCal");
		sprintf(cutstr,"%s && %s && isNat",HGCut,EFloor);
		gatSkim->Project("h3", drawstr, cutstr);

		h2->SetLineColor(kBlue);
		h2->GetXaxis()->SetTitle(drawstr);
		h2->GetYaxis()->SetTitle("Counts (8 bins / keV)");
		h2->Draw();

		h3->SetLineColor(kRed);
		h3->Draw("same");

		TLegend* leg1 = new TLegend(0.6,0.7,0.87,0.94);
		leg1->AddEntry(h2,"Enriched","l");
		leg1->AddEntry(h3,"Natural","l");
		leg1->Draw();

		c2->Update();
		c2->Print("output/LowE_EnrVsNat.pdf");
	}

	// You crazy bastard.
	// Plot ALL the gatified parameters for the waveforms matching the "HG + EFloor (+ SuperCut)?"
	// This is to have a reference when you try to remove shitty waveforms from the good list.
	// Save them in a ROOT file and do something smart with it.
	// NOTE: Think about the cut. Enriched vs natural?  detector by detector?
	// 	     You could store different parameters in different directories inside the ROOT file.
	if(0)
	{
		TFile *out = new TFile("WSR_GATParams.root","RECREATE");
		// TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
		// out->mkdir("GATParamHists");
		// out->cd("GATParamHists");

		// The Cut
		sprintf(cutstr,"%s && %s",HGCut,EFloor);

		gatSkim->Draw("rawWFMax>>h0",cutstr);
		TH1F *h0 = (TH1F*)gDirectory->Get("h0");
		h0->Write("rawWFMax",TObject::kOverwrite);

		gatSkim->Draw("rawWFMin>>h1",cutstr);
		TH1F *h1 = (TH1F*)gDirectory->Get("h1");
		h1->Write("rawWFMin",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR0p1>>h2",cutstr);
		TH1F *h2 = (TH1F*)gDirectory->Get("h2");
		h2->Write("blrwfFMR0p1",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR0p1>>h3",cutstr);
		TH1F *h3 = (TH1F*)gDirectory->Get("h3");
		h3->Write("sgswfFMR0p1",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR1>>h4",cutstr);
		TH1F *h4 = (TH1F*)gDirectory->Get("h4");
		h4->Write("blrwfFMR1",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR1>>h5",cutstr);
		TH1F *h5 = (TH1F*)gDirectory->Get("h5");
		h5->Write("sgswfFMR1",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR3>>h6",cutstr);
		TH1F *h6 = (TH1F*)gDirectory->Get("h6");
		h6->Write("blrwfFMR3",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR3>>h7",cutstr);
		TH1F *h7 = (TH1F*)gDirectory->Get("h7");
		h7->Write("sgswfFMR3",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR10>>h8",cutstr);
		TH1F *h8 = (TH1F*)gDirectory->Get("h8");
		h8->Write("blrwfFMR10",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR10>>h9",cutstr);
		TH1F *h9 = (TH1F*)gDirectory->Get("h9");
		h9->Write("sgswfFMR10",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR20>>h10",cutstr);
		TH1F *h10 = (TH1F*)gDirectory->Get("h10");
		h10->Write("blrwfFMR20",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR20>>h11",cutstr);
		TH1F *h11 = (TH1F*)gDirectory->Get("h11");
		h11->Write("sgswfFMR20",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR50>>h12",cutstr);
		TH1F *h12 = (TH1F*)gDirectory->Get("h12");
		h12->Write("blrwfFMR50",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR50>>h13",cutstr);
		TH1F *h13 = (TH1F*)gDirectory->Get("h13");
		h13->Write("sgswfFMR50",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR80>>h14",cutstr);
		TH1F *h14 = (TH1F*)gDirectory->Get("h14");
		h14->Write("blrwfFMR80",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR80>>h15",cutstr);
		TH1F *h15 = (TH1F*)gDirectory->Get("h15");
		h15->Write("sgswfFMR80",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR90>>h16",cutstr);
		TH1F *h16 = (TH1F*)gDirectory->Get("h16");
		h16->Write("blrwfFMR90",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR90>>h17",cutstr);
		TH1F *h17 = (TH1F*)gDirectory->Get("h17");
		h17->Write("sgswfFMR90",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR97>>h18",cutstr);
		TH1F *h18 = (TH1F*)gDirectory->Get("h18");
		h18->Write("blrwfFMR97",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR97>>h19",cutstr);
		TH1F *h19 = (TH1F*)gDirectory->Get("h19");
		h19->Write("sgswfFMR97",TObject::kOverwrite);

		gatSkim->Draw("blrwfFMR99>>h20",cutstr);
		TH1F *h20 = (TH1F*)gDirectory->Get("h20");
		h20->Write("blrwfFMR99",TObject::kOverwrite);

		gatSkim->Draw("sgswfFMR99>>h21",cutstr);
		TH1F *h21 = (TH1F*)gDirectory->Get("h21");
		h21->Write("sgswfFMR99",TObject::kOverwrite);

		gatSkim->Draw("sgswft0>>h22",cutstr);
		TH1F *h22 = (TH1F*)gDirectory->Get("h22");
		h22->Write("sgswft0",TObject::kOverwrite);

		gatSkim->Draw("TSCurrent50nsMax>>h23",cutstr);
		TH1F *h23 = (TH1F*)gDirectory->Get("h23");
		h23->Write("TSCurrent50nsMax",TObject::kOverwrite);

		gatSkim->Draw("TSCurrent100nsMax>>h24",cutstr);
		TH1F *h24 = (TH1F*)gDirectory->Get("h24");
		h24->Write("TSCurrent100nsMax",TObject::kOverwrite);

		gatSkim->Draw("TSCurrent200nsMax>>h25",cutstr);
		TH1F *h25 = (TH1F*)gDirectory->Get("h25");
		h25->Write("TSCurrent200nsMax",TObject::kOverwrite);

		gatSkim->Draw("RawWFblSlope>>h26",cutstr);
		TH1F *h26 = (TH1F*)gDirectory->Get("h26");
		h26->Write("RawWFblSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFblOffset>>h27",cutstr);
		TH1F *h27 = (TH1F*)gDirectory->Get("h27");
		h27->Write("RawWFblOffset",TObject::kOverwrite);

		gatSkim->Draw("RawWFblSlopeUnc>>h28",cutstr);
		TH1F *h28 = (TH1F*)gDirectory->Get("h28");
		h28->Write("RawWFblSlopeUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFblOffsetUnc>>h29",cutstr);
		TH1F *h29 = (TH1F*)gDirectory->Get("h29");
		h29->Write("RawWFblOffsetUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFblChi2>>h30",cutstr);
		TH1F *h30 = (TH1F*)gDirectory->Get("h30");
		h30->Write("RawWFblChi2",TObject::kOverwrite);

		gatSkim->Draw("RawWFftSlope>>h31",cutstr);
		TH1F *h31 = (TH1F*)gDirectory->Get("h31");
		h31->Write("RawWFftSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFftSlopeUnc>>h32",cutstr);
		TH1F *h32 = (TH1F*)gDirectory->Get("h32");
		h32->Write("RawWFftSlopeUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFftOffset>>h33",cutstr);
		TH1F *h33 = (TH1F*)gDirectory->Get("h33");
		h33->Write("RawWFftOffset",TObject::kOverwrite);

		gatSkim->Draw("RawWFftOffsetUnc>>h34",cutstr);
		TH1F *h34 = (TH1F*)gDirectory->Get("h34");
		h34->Write("RawWFftOffsetUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFftChi2>>h35",cutstr);
		TH1F *h35 = (TH1F*)gDirectory->Get("h35");
		h35->Write("RawWFftChi2",TObject::kOverwrite);

		gatSkim->Draw("RawWFwholeWFLinFitSlope>>h36",cutstr);
		TH1F *h36 = (TH1F*)gDirectory->Get("h36");
		h36->Write("RawWFwholeWFLinFitSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFwholeWFLinFitSlopeUnc>>h37",cutstr);
		TH1F *h37 = (TH1F*)gDirectory->Get("h37");
		h37->Write("RawWFwholeWFLinFitSlopeUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFwholeWFLinFitOffset>>h38",cutstr);
		TH1F *h38 = (TH1F*)gDirectory->Get("h38");
		h38->Write("RawWFwholeWFLinFitOffset",TObject::kOverwrite);

		gatSkim->Draw("RawWFwholeWFLinFitOffsetUnc>>h39",cutstr);
		TH1F *h39 = (TH1F*)gDirectory->Get("h39");
		h39->Write("RawWFwholeWFLinFitOffsetUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFwholeWFLinFitChi2>>h40",cutstr);
		TH1F *h40 = (TH1F*)gDirectory->Get("h40");
		h40->Write("RawWFwholeWFLinFitChi2",TObject::kOverwrite);

		gatSkim->Draw("trapENF55us>>h41",cutstr);
		TH1F *h41 = (TH1F*)gDirectory->Get("h41");
		h41->Write("trapENF55us",TObject::kOverwrite);

		gatSkim->Draw("trapENF70us>>h42",cutstr);
		TH1F *h42 = (TH1F*)gDirectory->Get("h42");
		h42->Write("trapENF70us",TObject::kOverwrite);

		gatSkim->Draw("trapENF85us>>h43",cutstr);
		TH1F *h43 = (TH1F*)gDirectory->Get("h43");
		h43->Write("trapENF85us",TObject::kOverwrite);

		gatSkim->Draw("trapENF100us>>h44",cutstr);
		TH1F *h44 = (TH1F*)gDirectory->Get("h44");
		h44->Write("trapENF100us",TObject::kOverwrite);

		gatSkim->Draw("trapENF115us>>h45",cutstr);
		TH1F *h45 = (TH1F*)gDirectory->Get("h45");
		h45->Write("trapENF115us",TObject::kOverwrite);

		gatSkim->Draw("trapENF130us>>h46",cutstr);
		TH1F *h46 = (TH1F*)gDirectory->Get("h46");
		h46->Write("trapENF130us",TObject::kOverwrite);

		gatSkim->Draw("trapBL>>h47",cutstr);
		TH1F *h47 = (TH1F*)gDirectory->Get("h47");
		h47->Write("trapBL",TObject::kOverwrite);

		gatSkim->Draw("trap500nsMax>>h48",cutstr);
		TH1F *h48 = (TH1F*)gDirectory->Get("h48");
		h48->Write("trap500nsMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF500nsrt>>h49",cutstr);
		TH1F *h49 = (TH1F*)gDirectory->Get("h49");
		h49->Write("trapENF500nsrt",TObject::kOverwrite);

		gatSkim->Draw("trap1usMax>>h50",cutstr);
		TH1F *h50 = (TH1F*)gDirectory->Get("h50");
		h50->Write("trap1usMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF1usrt>>h51",cutstr);
		TH1F *h51 = (TH1F*)gDirectory->Get("h51");
		h51->Write("trapENF1usrt",TObject::kOverwrite);

		gatSkim->Draw("trap2usMax>>h52",cutstr);
		TH1F *h52 = (TH1F*)gDirectory->Get("h52");
		h52->Write("trap2usMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF2usrt>>h53",cutstr);
		TH1F *h53 = (TH1F*)gDirectory->Get("h53");
		h53->Write("trapENF2usrt",TObject::kOverwrite);

		gatSkim->Draw("trap4usMax>>h54",cutstr);
		TH1F *h54 = (TH1F*)gDirectory->Get("h54");
		h54->Write("trap4usMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF4usrt>>h55",cutstr);
		TH1F *h55 = (TH1F*)gDirectory->Get("h55");
		h55->Write("trapENF4usrt",TObject::kOverwrite);

		gatSkim->Draw("trap6usMax>>h56",cutstr);
		TH1F *h56 = (TH1F*)gDirectory->Get("h56");
		h56->Write("trap6usMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF6usrt>>h57",cutstr);
		TH1F *h57 = (TH1F*)gDirectory->Get("h57");
		h57->Write("trapENF6usrt",TObject::kOverwrite);

		gatSkim->Draw("trap8usMax>>h58",cutstr);
		TH1F *h58 = (TH1F*)gDirectory->Get("h58");
		h58->Write("trap8usMax",TObject::kOverwrite);

		gatSkim->Draw("trapENF8usrt>>h59",cutstr);
		TH1F *h59 = (TH1F*)gDirectory->Get("h59");
		h59->Write("trapENF8usrt",TObject::kOverwrite);

		gatSkim->Draw("trapE>>h60",cutstr);
		TH1F *h60 = (TH1F*)gDirectory->Get("h60");
		h60->Write("trapE",TObject::kOverwrite);

		gatSkim->Draw("trapEMin>>h61",cutstr);
		TH1F *h61 = (TH1F*)gDirectory->Get("h61");
		h61->Write("trapEMin",TObject::kOverwrite);

		gatSkim->Draw("longGapTrapMax>>h62",cutstr);
		TH1F *h62 = (TH1F*)gDirectory->Get("h62");
		h62->Write("longGapTrapMax",TObject::kOverwrite);

		gatSkim->Draw("trap1ust0>>h63",cutstr);
		TH1F *h63 = (TH1F*)gDirectory->Get("h63");
		h63->Write("trap1ust0",TObject::kOverwrite);

		gatSkim->Draw("trapENF>>h64",cutstr);
		TH1F *h64 = (TH1F*)gDirectory->Get("h64");
		h64->Write("trapENF",TObject::kOverwrite);

		gatSkim->Draw("trapBL1us>>h65",cutstr);
		TH1F *h65 = (TH1F*)gDirectory->Get("h65");
		h65->Write("trapBL1us",TObject::kOverwrite);

		gatSkim->Draw("trapETailMin>>h66",cutstr);
		TH1F *h66 = (TH1F*)gDirectory->Get("h66");
		h66->Write("trapETailMin",TObject::kOverwrite);

		gatSkim->Draw("trirt50nsft10nsMax>>h67",cutstr);
		TH1F *h67 = (TH1F*)gDirectory->Get("h67");
		h67->Write("trirt50nsft10nsMax",TObject::kOverwrite);

		gatSkim->Draw("trirt100nsft10nsMax>>h68",cutstr);
		TH1F *h68 = (TH1F*)gDirectory->Get("h68");
		h68->Write("trirt100nsft10nsMax",TObject::kOverwrite);

		gatSkim->Draw("trirt200nsft10nsMax>>h69",cutstr);
		TH1F *h69 = (TH1F*)gDirectory->Get("h69");
		h69->Write("trirt200nsft10nsMax",TObject::kOverwrite);

		gatSkim->Draw("triFilMin>>h70",cutstr);
		TH1F *h70 = (TH1F*)gDirectory->Get("h70");
		h70->Write("triFilMin",TObject::kOverwrite);

		gatSkim->Draw("trirt100nsft10nsIntegralW>>h71",cutstr);
		TH1F *h71 = (TH1F*)gDirectory->Get("h71");
		h71->Write("trirt100nsft10nsIntegralW",TObject::kOverwrite);

		gatSkim->Draw("blrwfIntegralW>>h72",cutstr);
		TH1F *h72 = (TH1F*)gDirectory->Get("h72");
		h72->Write("blrwfIntegralW",TObject::kOverwrite);

		gatSkim->Draw("smoothTrirt100nsft10nsMax>>h73",cutstr);
		TH1F *h73 = (TH1F*)gDirectory->Get("h73");
		h73->Write("smoothTrirt100nsft10nsMax",TObject::kOverwrite);

		gatSkim->Draw("energyCal>>h74",cutstr);
		TH1F *h74 = (TH1F*)gDirectory->Get("h74");
		h74->Write("energyCal",TObject::kOverwrite);

		gatSkim->Draw("trapECal>>h75",cutstr);
		TH1F *h75 = (TH1F*)gDirectory->Get("h75");
		h75->Write("trapECal",TObject::kOverwrite);

		gatSkim->Draw("trapEMinCal>>h76",cutstr);
		TH1F *h76 = (TH1F*)gDirectory->Get("h76");
		h76->Write("trapEMinCal",TObject::kOverwrite);

		gatSkim->Draw("trapENFCal>>h77",cutstr);
		TH1F *h77 = (TH1F*)gDirectory->Get("h77");
		h77->Write("trapENFCal",TObject::kOverwrite);

		gatSkim->Draw("nlcblrwfSlope>>h78",cutstr);
		TH1F *h78 = (TH1F*)gDirectory->Get("h78");
		h78->Write("nlcblrwfSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcblSlope>>h79",cutstr);
		TH1F *h79 = (TH1F*)gDirectory->Get("h79");
		h79->Write("RawWFdcblSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcblOffset>>h80",cutstr);
		TH1F *h80 = (TH1F*)gDirectory->Get("h80");
		h80->Write("RawWFdcblOffset",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcblSlopeUnc>>h81",cutstr);
		TH1F *h81 = (TH1F*)gDirectory->Get("h81");
		h81->Write("RawWFdcblSlopeUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcblOffsetUnc>>h82",cutstr);
		TH1F *h82 = (TH1F*)gDirectory->Get("h82");
		h82->Write("RawWFdcblOffsetUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcblChi2>>h83",cutstr);
		TH1F *h83 = (TH1F*)gDirectory->Get("h83");
		h83->Write("RawWFdcblChi2",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcftSlope>>h84",cutstr);
		TH1F *h84 = (TH1F*)gDirectory->Get("h84");
		h84->Write("RawWFdcftSlope",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcftOffset>>h85",cutstr);
		TH1F *h85 = (TH1F*)gDirectory->Get("h85");
		h85->Write("RawWFdcftOffset",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcftSlopeUnc>>h86",cutstr);
		TH1F *h86 = (TH1F*)gDirectory->Get("h86");
		h86->Write("RawWFdcftSlopeUnc",TObject::kOverwrite);

		gatSkim->Draw("RawWFdcftOffsetUnc>>h87",cutstr);
		TH1F *h87 = (TH1F*)gDirectory->Get("h87");
		h87->Write("RawWFdcftOffsetUnc",TObject::kOverwrite);

		gatSkim->Draw("blrwf2ndDiffPeakValue>>h88",cutstr);
		TH1F *h88 = (TH1F*)gDirectory->Get("h88");
		h88->Write("blrwf2ndDiffPeakValue",TObject::kOverwrite);

		gatSkim->Draw("blrwfNorm2ndDiffPeakValue>>h89",cutstr);
		TH1F *h89 = (TH1F*)gDirectory->Get("h89");
		h89->Write("blrwfNorm2ndDiffPeakValue",TObject::kOverwrite);

		gatSkim->Draw("blrwfSSF>>h90",cutstr);
		TH1F *h90 = (TH1F*)gDirectory->Get("h90");
		h90->Write("blrwfSSF",TObject::kOverwrite);

		gatSkim->Draw("BLMoverMerr>>h91",cutstr);
		TH1F *h91 = (TH1F*)gDirectory->Get("h91");
		h91->Write("BLMoverMerr",TObject::kOverwrite);

		gatSkim->Draw("FTMoverMerr>>h92",cutstr);
		TH1F *h92 = (TH1F*)gDirectory->Get("h92");
		h92->Write("FTMoverMerr",TObject::kOverwrite);

		gatSkim->Draw("delOffset>>h93",cutstr);
		TH1F *h93 = (TH1F*)gDirectory->Get("h93");
		h93->Write("delOffset",TObject::kOverwrite);

		gatSkim->Draw("toe>>h94",cutstr);
		TH1F *h94 = (TH1F*)gDirectory->Get("h94");
		h94->Write("toe",TObject::kOverwrite);

		gatSkim->Draw("wfDCBits>>h95",cutstr);
		TH1F *h95 = (TH1F*)gDirectory->Get("h95");
		h95->Write("wfDCBits",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_0_SSSpikeBL>>h96",cutstr);
		TH1F *h96 = (TH1F*)gDirectory->Get("h96");
		h96->Write("wfDC_Bit_0_SSSpikeBL",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_1_SSSpikePhys>>h97",cutstr);
		TH1F *h97 = (TH1F*)gDirectory->Get("h97");
		h97->Write("wfDC_Bit_1_SSSpikePhys",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_2_EarlyTrigger>>h98",cutstr);
		TH1F *h98 = (TH1F*)gDirectory->Get("h98");
		h98->Write("wfDC_Bit_2_EarlyTrigger",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_3_LateTrigger>>h99",cutstr);
		TH1F *h99 = (TH1F*)gDirectory->Get("h99");
		h99->Write("wfDC_Bit_3_LateTrigger",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_4_PosSaturatedWFs>>h100",cutstr);
		TH1F *h100 = (TH1F*)gDirectory->Get("h100");
		h100->Write("wfDC_Bit_4_PosSaturatedWFs",TObject::kOverwrite);

		gatSkim->Draw("wfDC_Bit_5_NegSaturatedWFs>>h101",cutstr);
		TH1F *h101 = (TH1F*)gDirectory->Get("h101");
		h101->Write("wfDC_Bit_5_NegSaturatedWFs",TObject::kOverwrite);

		gatSkim->Draw("run>>h102",cutstr);
		TH1F *h102 = (TH1F*)gDirectory->Get("h102");
		h102->Write("run",TObject::kOverwrite);

		gatSkim->Draw("startTime>>h103",cutstr);
		TH1F *h103 = (TH1F*)gDirectory->Get("h103");
		h103->Write("startTime",TObject::kOverwrite);

		gatSkim->Draw("stopTime>>h104",cutstr);
		TH1F *h104 = (TH1F*)gDirectory->Get("h104");
		h104->Write("stopTime",TObject::kOverwrite);

		gatSkim->Draw("energy>>h105",cutstr);
		TH1F *h105 = (TH1F*)gDirectory->Get("h105");
		h105->Write("energy",TObject::kOverwrite);

		gatSkim->Draw("channel>>h106",cutstr);
		TH1F *h106 = (TH1F*)gDirectory->Get("h106");
		h106->Write("channel",TObject::kOverwrite);

		gatSkim->Draw("timestamp>>h107",cutstr);
		TH1F *h107 = (TH1F*)gDirectory->Get("h107");
		h107->Write("timestamp",TObject::kOverwrite);

		gatSkim->Draw("gatrev>>h108",cutstr);
		TH1F *h108 = (TH1F*)gDirectory->Get("h108");
		h108->Write("gatrev",TObject::kOverwrite);

		gatSkim->Draw("EventDC1Bits>>h109",cutstr);
		TH1F *h109 = (TH1F*)gDirectory->Get("h109");
		h109->Write("EventDC1Bits",TObject::kOverwrite);

		gatSkim->Draw("d2wfnoiseTagNorm>>h110",cutstr);
		TH1F *h110 = (TH1F*)gDirectory->Get("h110");
		h110->Write("d2wfnoiseTagNorm",TObject::kOverwrite);

		gatSkim->Draw("d2wfDCPower>>h111",cutstr);
		TH1F *h111 = (TH1F*)gDirectory->Get("h111");
		h111->Write("d2wfDCPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf0MHzTo2MHzPower>>h112",cutstr);
		TH1F *h112 = (TH1F*)gDirectory->Get("h112");
		h112->Write("d2wf0MHzTo2MHzPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf2MHzTo10MHzPower>>h113",cutstr);
		TH1F *h113 = (TH1F*)gDirectory->Get("h113");
		h113->Write("d2wf2MHzTo10MHzPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf10MHzTo20MHzPower>>h114",cutstr);
		TH1F *h114 = (TH1F*)gDirectory->Get("h114");
		h114->Write("d2wf10MHzTo20MHzPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf30MHzTo35MHzPower>>h115",cutstr);
		TH1F *h115 = (TH1F*)gDirectory->Get("h115");
		h115->Write("d2wf30MHzTo35MHzPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf48MHzTo50MHzPower>>h116",cutstr);
		TH1F *h116 = (TH1F*)gDirectory->Get("h116");
		h116->Write("d2wf48MHzTo50MHzPower",TObject::kOverwrite);

		gatSkim->Draw("d2wf0MHzTo50MHzPower>>h117",cutstr);
		TH1F *h117 = (TH1F*)gDirectory->Get("h117");
		h117->Write("d2wf0MHzTo50MHzPower",TObject::kOverwrite);

		out->Close();
	}
}

void wavePlots(char input[500])
{
	// initialize
	TFile *f = new TFile(input);
	TTree *gatSkim = (TTree*)f->Get("gatSkim");
	TTree *builtSkim = (TTree*)f->Get("builtSkim");
	builtSkim->AddFriend(gatSkim);

	// built branches (april 2016 format)
	MGTRun *bRun = new MGTRun();
	MGTEvent *event = new MGTEvent();
	// TClonesArray *channelData = new TClonesArray("MJTChannelData");
	builtSkim->SetBranchAddress("run",&bRun);
	builtSkim->SetBranchAddress("event",&event);
	// builtSkim->SetBranchAddress("channelData",&channelData);	// causes too many "TProcessID" objects to be written to output

	// gatified branches (may 2016 format) + a few additions by Clint from waveSkim
	int iEvent=0;	// gat file entry number
	double dtmu_sec=0;	// time since last muon
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
		gatSkim->SetBranchAddress("iEvent",&iEvent);	// added by clint
		gatSkim->SetBranchAddress("dtmu_sec",&dtmu_sec);	// skim file, downgraded to a single value
		// gatSkim->SetBranchAddress("mij","GATMIJ",&mij,32000,0);
		// gatSkim->SetBranchAddress("m",&m64);
		gatSkim->SetBranchAddress("mL",&mL64);	// passed in from skim file
		gatSkim->SetBranchAddress("mH",&mH64);
		// gatSkim->SetBranchAddress("i",&i);
		// gatSkim->SetBranchAddress("j",&j);
		// gatSkim->SetBranchAddress("iH",&iH);
		// gatSkim->SetBranchAddress("jH",&jH);
		// gatSkim->SetBranchAddress("iL",&iL);
		// gatSkim->SetBranchAddress("jL",&jL);
		gatSkim->SetBranchAddress("timeMT",&timeMT);
		gatSkim->SetBranchAddress("dateMT",&dateMT);
		// gatSkim->SetBranchAddress("detName",&detName);
		gatSkim->SetBranchAddress("detID",&detID);
		gatSkim->SetBranchAddress("C",&C);
		gatSkim->SetBranchAddress("P",&P);
		gatSkim->SetBranchAddress("D",&D);
		gatSkim->SetBranchAddress("mageID",&mageID);
		gatSkim->SetBranchAddress("isEnr",&isEnr);
		gatSkim->SetBranchAddress("isNat",&isNat);
		gatSkim->SetBranchAddress("rawWFMax",&rawWFMax);
		gatSkim->SetBranchAddress("rawWFMin",&rawWFMin);
		gatSkim->SetBranchAddress("blrwfFMR0p1",&blrwfFMR0p1);
		gatSkim->SetBranchAddress("sgswfFMR0p1",&sgswfFMR0p1);
		gatSkim->SetBranchAddress("blrwfFMR1",&blrwfFMR1);
		gatSkim->SetBranchAddress("sgswfFMR1",&sgswfFMR1);
		gatSkim->SetBranchAddress("blrwfFMR3",&blrwfFMR3);
		gatSkim->SetBranchAddress("sgswfFMR3",&sgswfFMR3);
		gatSkim->SetBranchAddress("blrwfFMR10",&blrwfFMR10);
		gatSkim->SetBranchAddress("sgswfFMR10",&sgswfFMR10);
		gatSkim->SetBranchAddress("blrwfFMR20",&blrwfFMR20);
		gatSkim->SetBranchAddress("sgswfFMR20",&sgswfFMR20);
		gatSkim->SetBranchAddress("blrwfFMR50",&blrwfFMR50);
		gatSkim->SetBranchAddress("sgswfFMR50",&sgswfFMR50);
		gatSkim->SetBranchAddress("blrwfFMR80",&blrwfFMR80);
		gatSkim->SetBranchAddress("sgswfFMR80",&sgswfFMR80);
		gatSkim->SetBranchAddress("blrwfFMR90",&blrwfFMR90);
		gatSkim->SetBranchAddress("sgswfFMR90",&sgswfFMR90);
		gatSkim->SetBranchAddress("blrwfFMR97",&blrwfFMR97);
		gatSkim->SetBranchAddress("sgswfFMR97",&sgswfFMR97);
		gatSkim->SetBranchAddress("blrwfFMR99",&blrwfFMR99);
		gatSkim->SetBranchAddress("sgswfFMR99",&sgswfFMR99);
		gatSkim->SetBranchAddress("sgswft0",&sgswft0);
		gatSkim->SetBranchAddress("TSCurrent50nsMax",&TSCurrent50nsMax);
		gatSkim->SetBranchAddress("TSCurrent100nsMax",&TSCurrent100nsMax);
		gatSkim->SetBranchAddress("TSCurrent200nsMax",&TSCurrent200nsMax);
		gatSkim->SetBranchAddress("RawWFblSlope",&RawWFblSlope);
		gatSkim->SetBranchAddress("RawWFblOffset",&RawWFblOffset);
		gatSkim->SetBranchAddress("RawWFblSlopeUnc",&RawWFblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFblOffsetUnc",&RawWFblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFblChi2",&RawWFblChi2);
		gatSkim->SetBranchAddress("RawWFftSlope",&RawWFftSlope);
		gatSkim->SetBranchAddress("RawWFftSlopeUnc",&RawWFftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFftOffset",&RawWFftOffset);
		gatSkim->SetBranchAddress("RawWFftOffsetUnc",&RawWFftOffsetUnc);
		gatSkim->SetBranchAddress("RawWFftChi2",&RawWFftChi2);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlope",&RawWFwholeWFLinFitSlope);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlopeUnc",&RawWFwholeWFLinFitSlopeUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffset",&RawWFwholeWFLinFitOffset);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffsetUnc",&RawWFwholeWFLinFitOffsetUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitChi2",&RawWFwholeWFLinFitChi2);
		gatSkim->SetBranchAddress("trapENF55us",&trapENF55us);
		gatSkim->SetBranchAddress("trapENF70us",&trapENF70us);
		gatSkim->SetBranchAddress("trapENF85us",&trapENF85us);
		gatSkim->SetBranchAddress("trapENF100us",&trapENF100us);
		gatSkim->SetBranchAddress("trapENF115us",&trapENF115us);
		gatSkim->SetBranchAddress("trapENF130us",&trapENF130us);
		gatSkim->SetBranchAddress("trapBL",&trapBL);
		gatSkim->SetBranchAddress("trap500nsMax",&trap500nsMax);
		gatSkim->SetBranchAddress("trapENF500nsrt",&trapENF500nsrt);
		gatSkim->SetBranchAddress("trap1usMax",&trap1usMax);
		gatSkim->SetBranchAddress("trapENF1usrt",&trapENF1usrt);
		gatSkim->SetBranchAddress("trap2usMax",&trap2usMax);
		gatSkim->SetBranchAddress("trapENF2usrt",&trapENF2usrt);
		gatSkim->SetBranchAddress("trap4usMax",&trap4usMax);
		gatSkim->SetBranchAddress("trapENF4usrt",&trapENF4usrt);
		gatSkim->SetBranchAddress("trap6usMax",&trap6usMax);
		gatSkim->SetBranchAddress("trapENF6usrt",&trapENF6usrt);
		gatSkim->SetBranchAddress("trap8usMax",&trap8usMax);
		gatSkim->SetBranchAddress("trapENF8usrt",&trapENF8usrt);
		gatSkim->SetBranchAddress("trapE",&trapE);
		gatSkim->SetBranchAddress("trapEMin",&trapEMin);
		gatSkim->SetBranchAddress("longGapTrapMax",&longGapTrapMax);
		gatSkim->SetBranchAddress("trap1ust0",&trap1ust0);
		gatSkim->SetBranchAddress("trapENF",&trapENF);
		gatSkim->SetBranchAddress("trapBL1us",&trapBL1us);
		gatSkim->SetBranchAddress("trapETailMin",&trapETailMin);
		gatSkim->SetBranchAddress("trirt50nsft10nsMax",&trirt50nsft10nsMax);
		gatSkim->SetBranchAddress("trirt100nsft10nsMax",&trirt100nsft10nsMax);
		gatSkim->SetBranchAddress("trirt200nsft10nsMax",&trirt200nsft10nsMax);
		gatSkim->SetBranchAddress("triFilMin",&triFilMin);
		gatSkim->SetBranchAddress("trirt100nsft10nsIntegralW",&trirt100nsft10nsIntegralW);
		gatSkim->SetBranchAddress("blrwfIntegralW",&blrwfIntegralW);
		gatSkim->SetBranchAddress("smoothTrirt100nsft10nsMax",&smoothTrirt100nsft10nsMax);
		gatSkim->SetBranchAddress("energyCal",&energyCal);
		gatSkim->SetBranchAddress("trapECal",&trapECal);
		gatSkim->SetBranchAddress("trapEMinCal",&trapEMinCal);
		gatSkim->SetBranchAddress("trapENFCal",&trapENFCal);
		gatSkim->SetBranchAddress("nlcblrwfSlope",&nlcblrwfSlope);
		gatSkim->SetBranchAddress("RawWFdcblSlope",&RawWFdcblSlope);
		gatSkim->SetBranchAddress("RawWFdcblOffset",&RawWFdcblOffset);
		gatSkim->SetBranchAddress("RawWFdcblSlopeUnc",&RawWFdcblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcblOffsetUnc",&RawWFdcblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFdcblChi2",&RawWFdcblChi2);
		gatSkim->SetBranchAddress("RawWFdcftSlope",&RawWFdcftSlope);
		gatSkim->SetBranchAddress("RawWFdcftOffset",&RawWFdcftOffset);
		gatSkim->SetBranchAddress("RawWFdcftSlopeUnc",&RawWFdcftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcftOffsetUnc",&RawWFdcftOffsetUnc);
		gatSkim->SetBranchAddress("blrwf2ndDiffPeakValue",&blrwf2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfNorm2ndDiffPeakValue",&blrwfNorm2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfSSF",&blrwfSSF);
		gatSkim->SetBranchAddress("BLMoverMerr",&BLMoverMerr);
		gatSkim->SetBranchAddress("FTMoverMerr",&FTMoverMerr);
		gatSkim->SetBranchAddress("delOffset",&delOffset);
		gatSkim->SetBranchAddress("toe",&toe);
		gatSkim->SetBranchAddress("wfDCBits",&wfDCBits);
		gatSkim->SetBranchAddress("wfDC_Bit_0_SSSpikeBL",&wfDC_Bit_0_SSSpikeBL);
		gatSkim->SetBranchAddress("wfDC_Bit_1_SSSpikePhys",&wfDC_Bit_1_SSSpikePhys);
		gatSkim->SetBranchAddress("wfDC_Bit_2_EarlyTrigger",&wfDC_Bit_2_EarlyTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_3_LateTrigger",&wfDC_Bit_3_LateTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_4_PosSaturatedWFs",&wfDC_Bit_4_PosSaturatedWFs);
		gatSkim->SetBranchAddress("wfDC_Bit_5_NegSaturatedWFs",&wfDC_Bit_5_NegSaturatedWFs);
		gatSkim->SetBranchAddress("run",&run);
		gatSkim->SetBranchAddress("startTime",&startTime);
		gatSkim->SetBranchAddress("stopTime",&stopTime);
		gatSkim->SetBranchAddress("energy",&energy);
		gatSkim->SetBranchAddress("channel",&channel);
		gatSkim->SetBranchAddress("timestamp",&timestamp);
		gatSkim->SetBranchAddress("gatrev",&gatrev);
		gatSkim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
		gatSkim->SetBranchAddress("d2wfnoiseTagNorm",&d2wfnoiseTagNorm);
		gatSkim->SetBranchAddress("d2wfDCPower",&d2wfDCPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo2MHzPower",&d2wf0MHzTo2MHzPower);
		gatSkim->SetBranchAddress("d2wf2MHzTo10MHzPower",&d2wf2MHzTo10MHzPower);
		gatSkim->SetBranchAddress("d2wf10MHzTo20MHzPower",&d2wf10MHzTo20MHzPower);
		gatSkim->SetBranchAddress("d2wf30MHzTo35MHzPower",&d2wf30MHzTo35MHzPower);
		gatSkim->SetBranchAddress("d2wf48MHzTo50MHzPower",&d2wf48MHzTo50MHzPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo50MHzPower",&d2wf0MHzTo50MHzPower);
	}

	// Waveform Output
	TFile *out = new TFile("WSR_WFOutput.root","RECREATE");
	TTree *waves = new TTree("waves","skimmed waveforms");
	MGTWaveform *wave = new MGTWaveform();
	string *cutstring=0;
	int chan=0;
	double e_ENFCal=0, e_ECal=0, t_blrwfFMR50=0;
	waves->Branch("cutstr",&cutstring);
	waves->Branch("wave",&wave);
	waves->Branch("run",&run);
	waves->Branch("iEvent",&iEvent);
	waves->Branch("channel",&chan);
	waves->Branch("trapENFCal",&e_ENFCal);
	waves->Branch("trapECal",&e_ECal);
	waves->Branch("blrwfFMR50",&t_blrwfFMR50);

	///////////////////////////////////////////////////////////
	//                                                       //
	// =================== PLOT WAVEFORMS ================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	gROOT->ProcessLine(".x /Users/wisecg/dev/MJDWavePlotStyle.C");

	// ========================== THE BIG CUT ==========================
	// (Should match the cut used to generate the waveSkim input file.)
	// =================================================================
	//
	char cutstr[1000];

	// low-energy cut
	char HGCut[500] = "channel%2==0";	// only look at high gain
	char SDCut[500] = "mH==1";
	char energyCut[500] = "trapENFCal > 0 && trapENFCal < 30";
	char TrapETailMinCut[500] = "trapETailMin < 0";
	char DCCut[500] = "!wfDCBits";
	char vetoCut[500] = "(dtmu_sec < -0.2e-3 || dtmu_sec > 1)";
	char burstCut[500] = "1";//!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";
	char runCut[500] = "run!=13312 && run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863";
	char segfaultCut[500] = "run!=9422. && run!=9423 && run!=9425 && run!=9426 && run!=9427";
	sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && %s"
		,HGCut,SDCut,energyCut,TrapETailMinCut,DCCut,vetoCut,burstCut,runCut,segfaultCut);

	// DEP cut
	// char HGCut[500] = "channel%2==0";	// only look at high gain
	// // char SDCut[500] = "mH==1";
	// char SDCut[500] = "1";
	// // char energyCut[500] = "trapENFCal > 0 && trapENFCal < 30";
	// char energyCut[500] = "trapENFCal > 1590 && trapENFCal < 1595";
	// // char TrapETailMinCut[500] = "trapETailMin < 0";
	// char TrapETailMinCut[500] = "1";
	// char DCCut[500] = "!wfDCBits";
	// char vetoCut[500] = "(dtmu_sec < -0.2e-3 || dtmu_sec > 1)";
	// char burstCut[500] = "1";//!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";
	// char runCut[500] = "run!=13312 && run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863";
	// // char segfaultCut[500] = "run!=9422. && run!=9423 && run!=9425 && run!=9426 && run!=9427";
	// char segfaultCut[500] = "1";
	// sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && %s"
	// 	,HGCut,SDCut,energyCut,TrapETailMinCut,DCCut,vetoCut,burstCut,runCut,segfaultCut);


	// Loop over waveSkim events
	//
	long entries = builtSkim->GetEntries();
	for (int i = 0; i < entries; i++)
	// for (int i = 0; i < 1; i++)
	{
		gatSkim->GetEntry(i);

		// progress
		if(i%100 == 0) printf("[%.1f%%]\n",((double)i/entries)*100);


		// Make a GWB with only one event so we can still apply the BIG CUT from waveSkim.
		// Otherwise, the other hits in the event will show up in the output.
		//
		char thisEntry[500];
		sprintf(thisEntry,"Entry$==%i && %s",i,cutstr);
		GATWaveformBrowser wb(builtSkim,thisEntry);
		int numWFs = (int)wb.GetNWaveforms();

		// Loop over waveforms in the event that pass cuts
		//
		for (int j = 0; j < (int)numWFs; j++)
		{

			MGTWaveform *twf = wb.GetWaveform(i,j);
			if(twf==NULL) { cout << "couldn't get waveform!\n"; return; }
			MGTWaveform wf = *(twf);
			TH1D *hist = wf.GimmeHist();

			// Pull some variables directly from the waveform object
			double ch=0, pos=0, det=0, enf=0;
			const char* wfInfo = "P:D:channel:trapENFCal";
			builtSkim->Draw(wfInfo, TString::Format("Iteration$ == %i", j), "GOFF", 1, i);
			if(builtSkim->GetV1()==NULL||builtSkim->GetV2()==NULL||builtSkim->GetV3()==NULL||builtSkim->GetV4()==NULL){
				cout << "Warning!  Variable doesn't exist!  Can't make waveform title!\n";
				break;
			}
			pos = builtSkim->GetV1()[0];
			det = builtSkim->GetV2()[0];
			ch = builtSkim->GetV3()[0];
			enf = builtSkim->GetV4()[0];
			TString title = TString::Format("P%gD%g(%g),%s=%g",pos,det,ch,builtSkim->GetVar4()->GetTitle(),enf);
			hist->SetTitle(title.Data());

			// printf("waveSkim entry %i  size %lu  run %.0f  mH %i  numWFs %i  iEvent %i\n"
			// 	,i,channel->size(),run,mH64,numWFs,iEvent);

			// PDF Output: run, channel, waveSkim event, gat event.
			//
			// TCanvas *c = new TCanvas("c","Bob Ross's Canvas",800,600);
			// char wfName[500];
			// sprintf(wfName,"./output/wave_%i_%.0f_%.0f_%i.pdf",i,run,ch,iEvent);
			// // hist->SetMaximum(300);
			// // hist->SetMinimum(0);
			// hist->Draw("L");
			// c->Print(wfName);
			// delete c;

			// Output file for Ben
			// waves variables: wave (filled here), run (already loaded), iEvent (already loaded),
			// chan (vec), e_ENFCal (vec), e_ECal (vec), t_blrwfFMR50 (vec) -> need to grab the individual values
			string str(cutstr);
			cutstring = &str;
			wave = &wf;
			int fullEventSize = (int)channel->size();
			for (int k = 0; k < fullEventSize; k++)
			{
				if (channel->at(k) == ch)
				{
					e_ENFCal = trapENFCal->at(k);
					e_ECal = trapECal->at(k);
					chan = channel->at(k);
					t_blrwfFMR50 = blrwfFMR50->at(k);
				}
			}
			waves->Fill();
		}
	}

	// Write output and close
	waves->Write();
	out->Close();
}

/*
void entryPlots(char input[500])
{
	// A list of entries from "wavePlots" we want to look at the corresponding GAT parameters of.
	// vector<int> entryList = {90};	// spike waveforms.  rawWFMin > 0 might cut some of them.
	// vector<int> entryList = {13,14};	// slow rise noisy events.  no obvious cut yet.

	// run 10185 has some weird baseline drift things without physics (events 243 - 256)

	// events 357 - 360 have some ringing events (different runs)

	// initialize
	TFile *f = new TFile(input);
	TTree *gatSkim = (TTree*)f->Get("gatSkim");
	TTree *builtSkim = (TTree*)f->Get("builtSkim");
	builtSkim->AddFriend(gatSkim);

	// built branches (april 2016 format)
	MGTRun *bRun = new MGTRun();
	MGTEvent *event = new MGTEvent();
	// TClonesArray *channelData = new TClonesArray("MJTChannelData");
	builtSkim->SetBranchAddress("run",&bRun);
	builtSkim->SetBranchAddress("event",&event);
	// builtSkim->SetBranchAddress("channelData",&channelData);	// causes too many "TProcessID" objects to be written to output

	// gatified branches (may 2016 format) + a few additions by Clint from waveSkim
	int iEvent=0;	// gat file entry number
	double dtmu_sec=0;	// time since last muon
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
		gatSkim->SetBranchAddress("iEvent",&iEvent);	// added by clint
		gatSkim->SetBranchAddress("dtmu_sec",&dtmu_sec);	// skim file, downgraded to a single value
		// gatSkim->SetBranchAddress("mij","GATMIJ",&mij,32000,0);
		// gatSkim->SetBranchAddress("m",&m64);
		gatSkim->SetBranchAddress("mL",&mL64);	// passed in from skim file
		gatSkim->SetBranchAddress("mH",&mH64);
		// gatSkim->SetBranchAddress("i",&i);
		// gatSkim->SetBranchAddress("j",&j);
		// gatSkim->SetBranchAddress("iH",&iH);
		// gatSkim->SetBranchAddress("jH",&jH);
		// gatSkim->SetBranchAddress("iL",&iL);
		// gatSkim->SetBranchAddress("jL",&jL);
		gatSkim->SetBranchAddress("timeMT",&timeMT);
		gatSkim->SetBranchAddress("dateMT",&dateMT);
		// gatSkim->SetBranchAddress("detName",&detName);
		gatSkim->SetBranchAddress("detID",&detID);
		gatSkim->SetBranchAddress("C",&C);
		gatSkim->SetBranchAddress("P",&P);
		gatSkim->SetBranchAddress("D",&D);
		gatSkim->SetBranchAddress("mageID",&mageID);
		gatSkim->SetBranchAddress("isEnr",&isEnr);
		gatSkim->SetBranchAddress("isNat",&isNat);
		gatSkim->SetBranchAddress("rawWFMax",&rawWFMax);
		gatSkim->SetBranchAddress("rawWFMin",&rawWFMin);
		gatSkim->SetBranchAddress("blrwfFMR0p1",&blrwfFMR0p1);
		gatSkim->SetBranchAddress("sgswfFMR0p1",&sgswfFMR0p1);
		gatSkim->SetBranchAddress("blrwfFMR1",&blrwfFMR1);
		gatSkim->SetBranchAddress("sgswfFMR1",&sgswfFMR1);
		gatSkim->SetBranchAddress("blrwfFMR3",&blrwfFMR3);
		gatSkim->SetBranchAddress("sgswfFMR3",&sgswfFMR3);
		gatSkim->SetBranchAddress("blrwfFMR10",&blrwfFMR10);
		gatSkim->SetBranchAddress("sgswfFMR10",&sgswfFMR10);
		gatSkim->SetBranchAddress("blrwfFMR20",&blrwfFMR20);
		gatSkim->SetBranchAddress("sgswfFMR20",&sgswfFMR20);
		gatSkim->SetBranchAddress("blrwfFMR50",&blrwfFMR50);
		gatSkim->SetBranchAddress("sgswfFMR50",&sgswfFMR50);
		gatSkim->SetBranchAddress("blrwfFMR80",&blrwfFMR80);
		gatSkim->SetBranchAddress("sgswfFMR80",&sgswfFMR80);
		gatSkim->SetBranchAddress("blrwfFMR90",&blrwfFMR90);
		gatSkim->SetBranchAddress("sgswfFMR90",&sgswfFMR90);
		gatSkim->SetBranchAddress("blrwfFMR97",&blrwfFMR97);
		gatSkim->SetBranchAddress("sgswfFMR97",&sgswfFMR97);
		gatSkim->SetBranchAddress("blrwfFMR99",&blrwfFMR99);
		gatSkim->SetBranchAddress("sgswfFMR99",&sgswfFMR99);
		gatSkim->SetBranchAddress("sgswft0",&sgswft0);
		gatSkim->SetBranchAddress("TSCurrent50nsMax",&TSCurrent50nsMax);
		gatSkim->SetBranchAddress("TSCurrent100nsMax",&TSCurrent100nsMax);
		gatSkim->SetBranchAddress("TSCurrent200nsMax",&TSCurrent200nsMax);
		gatSkim->SetBranchAddress("RawWFblSlope",&RawWFblSlope);
		gatSkim->SetBranchAddress("RawWFblOffset",&RawWFblOffset);
		gatSkim->SetBranchAddress("RawWFblSlopeUnc",&RawWFblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFblOffsetUnc",&RawWFblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFblChi2",&RawWFblChi2);
		gatSkim->SetBranchAddress("RawWFftSlope",&RawWFftSlope);
		gatSkim->SetBranchAddress("RawWFftSlopeUnc",&RawWFftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFftOffset",&RawWFftOffset);
		gatSkim->SetBranchAddress("RawWFftOffsetUnc",&RawWFftOffsetUnc);
		gatSkim->SetBranchAddress("RawWFftChi2",&RawWFftChi2);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlope",&RawWFwholeWFLinFitSlope);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitSlopeUnc",&RawWFwholeWFLinFitSlopeUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffset",&RawWFwholeWFLinFitOffset);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitOffsetUnc",&RawWFwholeWFLinFitOffsetUnc);
		gatSkim->SetBranchAddress("RawWFwholeWFLinFitChi2",&RawWFwholeWFLinFitChi2);
		gatSkim->SetBranchAddress("trapENF55us",&trapENF55us);
		gatSkim->SetBranchAddress("trapENF70us",&trapENF70us);
		gatSkim->SetBranchAddress("trapENF85us",&trapENF85us);
		gatSkim->SetBranchAddress("trapENF100us",&trapENF100us);
		gatSkim->SetBranchAddress("trapENF115us",&trapENF115us);
		gatSkim->SetBranchAddress("trapENF130us",&trapENF130us);
		gatSkim->SetBranchAddress("trapBL",&trapBL);
		gatSkim->SetBranchAddress("trap500nsMax",&trap500nsMax);
		gatSkim->SetBranchAddress("trapENF500nsrt",&trapENF500nsrt);
		gatSkim->SetBranchAddress("trap1usMax",&trap1usMax);
		gatSkim->SetBranchAddress("trapENF1usrt",&trapENF1usrt);
		gatSkim->SetBranchAddress("trap2usMax",&trap2usMax);
		gatSkim->SetBranchAddress("trapENF2usrt",&trapENF2usrt);
		gatSkim->SetBranchAddress("trap4usMax",&trap4usMax);
		gatSkim->SetBranchAddress("trapENF4usrt",&trapENF4usrt);
		gatSkim->SetBranchAddress("trap6usMax",&trap6usMax);
		gatSkim->SetBranchAddress("trapENF6usrt",&trapENF6usrt);
		gatSkim->SetBranchAddress("trap8usMax",&trap8usMax);
		gatSkim->SetBranchAddress("trapENF8usrt",&trapENF8usrt);
		gatSkim->SetBranchAddress("trapE",&trapE);
		gatSkim->SetBranchAddress("trapEMin",&trapEMin);
		gatSkim->SetBranchAddress("longGapTrapMax",&longGapTrapMax);
		gatSkim->SetBranchAddress("trap1ust0",&trap1ust0);
		gatSkim->SetBranchAddress("trapENF",&trapENF);
		gatSkim->SetBranchAddress("trapBL1us",&trapBL1us);
		gatSkim->SetBranchAddress("trapETailMin",&trapETailMin);
		gatSkim->SetBranchAddress("trirt50nsft10nsMax",&trirt50nsft10nsMax);
		gatSkim->SetBranchAddress("trirt100nsft10nsMax",&trirt100nsft10nsMax);
		gatSkim->SetBranchAddress("trirt200nsft10nsMax",&trirt200nsft10nsMax);
		gatSkim->SetBranchAddress("triFilMin",&triFilMin);
		gatSkim->SetBranchAddress("trirt100nsft10nsIntegralW",&trirt100nsft10nsIntegralW);
		gatSkim->SetBranchAddress("blrwfIntegralW",&blrwfIntegralW);
		gatSkim->SetBranchAddress("smoothTrirt100nsft10nsMax",&smoothTrirt100nsft10nsMax);
		gatSkim->SetBranchAddress("energyCal",&energyCal);
		gatSkim->SetBranchAddress("trapECal",&trapECal);
		gatSkim->SetBranchAddress("trapEMinCal",&trapEMinCal);
		gatSkim->SetBranchAddress("trapENFCal",&trapENFCal);
		gatSkim->SetBranchAddress("nlcblrwfSlope",&nlcblrwfSlope);
		gatSkim->SetBranchAddress("RawWFdcblSlope",&RawWFdcblSlope);
		gatSkim->SetBranchAddress("RawWFdcblOffset",&RawWFdcblOffset);
		gatSkim->SetBranchAddress("RawWFdcblSlopeUnc",&RawWFdcblSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcblOffsetUnc",&RawWFdcblOffsetUnc);
		gatSkim->SetBranchAddress("RawWFdcblChi2",&RawWFdcblChi2);
		gatSkim->SetBranchAddress("RawWFdcftSlope",&RawWFdcftSlope);
		gatSkim->SetBranchAddress("RawWFdcftOffset",&RawWFdcftOffset);
		gatSkim->SetBranchAddress("RawWFdcftSlopeUnc",&RawWFdcftSlopeUnc);
		gatSkim->SetBranchAddress("RawWFdcftOffsetUnc",&RawWFdcftOffsetUnc);
		gatSkim->SetBranchAddress("blrwf2ndDiffPeakValue",&blrwf2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfNorm2ndDiffPeakValue",&blrwfNorm2ndDiffPeakValue);
		gatSkim->SetBranchAddress("blrwfSSF",&blrwfSSF);
		gatSkim->SetBranchAddress("BLMoverMerr",&BLMoverMerr);
		gatSkim->SetBranchAddress("FTMoverMerr",&FTMoverMerr);
		gatSkim->SetBranchAddress("delOffset",&delOffset);
		gatSkim->SetBranchAddress("toe",&toe);
		gatSkim->SetBranchAddress("wfDCBits",&wfDCBits);
		gatSkim->SetBranchAddress("wfDC_Bit_0_SSSpikeBL",&wfDC_Bit_0_SSSpikeBL);
		gatSkim->SetBranchAddress("wfDC_Bit_1_SSSpikePhys",&wfDC_Bit_1_SSSpikePhys);
		gatSkim->SetBranchAddress("wfDC_Bit_2_EarlyTrigger",&wfDC_Bit_2_EarlyTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_3_LateTrigger",&wfDC_Bit_3_LateTrigger);
		gatSkim->SetBranchAddress("wfDC_Bit_4_PosSaturatedWFs",&wfDC_Bit_4_PosSaturatedWFs);
		gatSkim->SetBranchAddress("wfDC_Bit_5_NegSaturatedWFs",&wfDC_Bit_5_NegSaturatedWFs);
		gatSkim->SetBranchAddress("run",&run);
		gatSkim->SetBranchAddress("startTime",&startTime);
		gatSkim->SetBranchAddress("stopTime",&stopTime);
		gatSkim->SetBranchAddress("energy",&energy);
		gatSkim->SetBranchAddress("channel",&channel);
		gatSkim->SetBranchAddress("timestamp",&timestamp);
		gatSkim->SetBranchAddress("gatrev",&gatrev);
		gatSkim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
		gatSkim->SetBranchAddress("d2wfnoiseTagNorm",&d2wfnoiseTagNorm);
		gatSkim->SetBranchAddress("d2wfDCPower",&d2wfDCPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo2MHzPower",&d2wf0MHzTo2MHzPower);
		gatSkim->SetBranchAddress("d2wf2MHzTo10MHzPower",&d2wf2MHzTo10MHzPower);
		gatSkim->SetBranchAddress("d2wf10MHzTo20MHzPower",&d2wf10MHzTo20MHzPower);
		gatSkim->SetBranchAddress("d2wf30MHzTo35MHzPower",&d2wf30MHzTo35MHzPower);
		gatSkim->SetBranchAddress("d2wf48MHzTo50MHzPower",&d2wf48MHzTo50MHzPower);
		gatSkim->SetBranchAddress("d2wf0MHzTo50MHzPower",&d2wf0MHzTo50MHzPower);
	}



	// ========================== THE BIG CUT ==========================
	// (Should match the cut used to generate the waveSkim input file.)
	// =================================================================
	//
	char cutstr[1000];
	char HGCut[500] = "channel%2==0";	// only look at high gain
	char SDCut[500] = "mH==1";
	char energyCut[500] = "trapENFCal > 0 && trapENFCal < 30";
	char TrapETailMinCut[500] = "trapETailMin < 0";
	char DCCut[500] = "!wfDCBits";
	char vetoCut[500] = "(dtmu_sec < -0.2e-3 || dtmu_sec > 1)";
	char burstCut[500] = "1";//!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";
	char runCut[500] = "run!=13312 && run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863";
	char segfaultCut[500] = "run!=9422. && run!=9423 && run!=9425 && run!=9426 && run!=9427";
	sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && %s"
		,HGCut,SDCut,energyCut,TrapETailMinCut,DCCut,vetoCut,burstCut,runCut,segfaultCut);


	// gat branch list
	vector<string> gatName = {"rawWFMax", // 0
	"rawWFMin",  // 1
	"blrwfFMR0p1",  // 2
	"sgswfFMR0p1",  // 3
	"blrwfFMR1",  // 4
	"sgswfFMR1",  // 5
	"blrwfFMR3",  // 6
	"sgswfFMR3",  // 7
	"blrwfFMR10",  // 8
	"sgswfFMR10",  // 9
	"blrwfFMR20",  // 10
	"sgswfFMR20",  // 11
	"blrwfFMR50",  // 12
	"sgswfFMR50",  // 13
	"blrwfFMR80",  // 14
	"sgswfFMR80",  // 15
	"blrwfFMR90",  // 16
	"sgswfFMR90",  // 17
	"blrwfFMR97",  // 18
	"sgswfFMR97",  // 19
	"blrwfFMR99",  // 20
	"sgswfFMR99",  // 21
	"sgswft0",  // 22
	"TSCurrent50nsMax",  // 23
	"TSCurrent100nsMax",  // 24
	"TSCurrent200nsMax",  // 25
	"RawWFblSlope",  // 26
	"RawWFblOffset",  // 27
	"RawWFblSlopeUnc",  // 28
	"RawWFblOffsetUnc",  // 29
	"RawWFblChi2",  // 30
	"RawWFftSlope",  // 31
	"RawWFftSlopeUnc",  // 32
	"RawWFftOffset",  // 33
	"RawWFftOffsetUnc",  // 34
	"RawWFftChi2",  // 35
	"RawWFwholeWFLinFitSlope",  // 36
	"RawWFwholeWFLinFitSlopeUnc",  // 37
	"RawWFwholeWFLinFitOffset",  // 38
	"RawWFwholeWFLinFitOffsetUnc",  // 39
	"RawWFwholeWFLinFitChi2",  // 40
	"trapENF55us",  // 41
	"trapENF70us",  // 42
	"trapENF85us",  // 43
	"trapENF100us",  // 44
	"trapENF115us",  // 45
	"trapENF130us",  // 46
	"trapBL",  // 47
	"trap500nsMax",  // 48
	"trapENF500nsrt",  // 49
	"trap1usMax",  // 50
	"trapENF1usrt",  // 51
	"trap2usMax",  // 52
	"trapENF2usrt",  // 53
	"trap4usMax",  // 54
	"trapENF4usrt",  // 55
	"trap6usMax",  // 56
	"trapENF6usrt",  // 57
	"trap8usMax",  // 58
	"trapENF8usrt",  // 59
	"trapE",  // 60
	"trapEMin",  //61
	"longGapTrapMax", // 62
	"trap1ust0",  // 63
	"trapENF",  // 64
	"trapBL1us",  // 65
	"trapETailMin",  // 66
	"trirt50nsft10nsMax",  // 67
	"trirt100nsft10nsMax",  // 68
	"trirt200nsft10nsMax",  // 69
	"triFilMin",  // 70
	"trirt100nsft10nsIntegralW",  // 71
	"blrwfIntegralW",  // 72
	"smoothTrirt100nsft10nsMax",  // 73
	"energyCal",  // 74
	"trapECal",  // 75
	"trapEMinCal",  // 76
	"trapENFCal",  // 77
	"nlcblrwfSlope",  // 78
	"RawWFdcblSlope",  // 79
	"RawWFdcblOffset",  // 80
	"RawWFdcblSlopeUnc",  // 81
	"RawWFdcblOffsetUnc",  // 82
	"RawWFdcblChi2",  // 83
	"RawWFdcftSlope",  // 84
	"RawWFdcftOffset",  // 85
	"RawWFdcftSlopeUnc",  // 86
	"RawWFdcftOffsetUnc",  // 87
	"blrwf2ndDiffPeakValue",  // 88
	"blrwfNorm2ndDiffPeakValue",  // 89
	"blrwfSSF",  // 90
	"BLMoverMerr",  // 91
	"FTMoverMerr",  // 92
	"delOffset",  // 93
	"toe",  // 94
	"wfDCBits",  // 95
	"wfDC_Bit_0_SSSpikeBL",  // 96
	"wfDC_Bit_1_SSSpikePhys",  // 97
	"wfDC_Bit_2_EarlyTrigger",  // 98
	"wfDC_Bit_3_LateTrigger",  // 99
	"wfDC_Bit_4_PosSaturatedWFs",  // 100
	"wfDC_Bit_5_NegSaturatedWFs",  // 101
	"run",  // 102
	"startTime",  // 103
	"stopTime",  // 104
	"energy",  // 105
	"channel",  // 106
	"timestamp",  // 107
	"gatrev",  // 108
	"EventDC1Bits",  // 109
	"d2wfnoiseTagNorm",  // 110
	"d2wfDCPower",  // 111
	"d2wf0MHzTo2MHzPower",  // 112
	"d2wf2MHzTo10MHzPower",  // 113
	"d2wf10MHzTo20MHzPower",  // 114
	"d2wf30MHzTo35MHzPower",  // 115
	"d2wf48MHzTo50MHzPower",  // 116
	"d2wf0MHzTo50MHzPower"};  // 117

	// Load histogram file of GAT parameters made by "plotter"
	TFile *g = new TFile("WSR_GATParams.root");

	// two vectors of histograms that correspond to gatName.
	vector<TH1D*> gatHist;
	gatHist.resize(118);
	vector<TH1D*> gatSingle;
	gatSingle.resize(118);
	cout << "gat hist size: " << (int)gatHist.size() << endl;
	for(int i = 0; i < (int)gatHist.size(); i++)
	{
		gatHist[i] = (TH1D*)g->Get(gatName[i].c_str());
		double x_low = gatHist[i]->GetXaxis()->GetXmin();
		double x_hi = gatHist[i]->GetXaxis()->GetXmax();
		gatSingle[i] = new TH1D(gatName[i].c_str(),gatName[i].c_str(),gatHist[i]->GetNbinsX(),x_low,x_hi);
		gatSingle[i]->SetLineColor(kRed);
		gatSingle[i]->SetFillColor(kRed);
	}
	// g->Close();

	// Make an output file
	TFile *h = new TFile("WSR_GATSingleWaveformStats.root","RECREATE");
	gROOT->ProcessLine(".x /Users/wisecg/dev/MJDWavePlotStyle.C");

	// Pull out waveforms that are in the entry list.
	for (int i = 0; i < (int)entryList.size(); i++)
	{
		cout << "Getting entry " << entryList[i] << endl;
		gatSkim->GetEntry(entryList[i]);

		// Make a GWB with only one event so we can still apply the BIG CUT from waveSkim.
		// Otherwise, the other hits in the event will show up in the output.
		char thisEntry[500];
		sprintf(thisEntry,"Entry$==%i && %s",entryList[i],cutstr);
		GATWaveformBrowser wb(builtSkim,thisEntry);
		int numWFs = (int)wb.GetNWaveforms();

		// Loop over waveforms in the event that pass cuts
		for (int j = 0; j < (int)numWFs; j++)
		{
			MGTWaveform *twf = wb.GetWaveform(entryList[i],j);
			if(twf==NULL) { cout << "couldn't get waveform!\n"; return; }
			MGTWaveform wf = *(twf);
			TH1D *hist = wf.GimmeHist();

			// Pull some variables directly from the waveform object
			double ch=0, pos=0, det=0, enf=0;
			const char* wfInfo = "P:D:channel:trapENFCal";
			builtSkim->Draw(wfInfo, TString::Format("Iteration$ == %i", j), "GOFF", 1, entryList[i]);
			if(builtSkim->GetV1()==NULL||builtSkim->GetV2()==NULL||builtSkim->GetV3()==NULL||builtSkim->GetV4()==NULL){
				cout << "Warning!  Variable doesn't exist!  Can't make waveform title!\n";
				break;
			}
			pos = builtSkim->GetV1()[0];
			det = builtSkim->GetV2()[0];
			ch = builtSkim->GetV3()[0];
			enf = builtSkim->GetV4()[0];
			TString title = TString::Format("P%gD%g(%g),%s=%g",pos,det,ch,builtSkim->GetVar4()->GetTitle(),enf);

			// PDF Output: run, channel, waveSkim event, gat event.
			TCanvas *theWF = new TCanvas("waveform","Bob Ross's Canvas",800,600);
			char wfName[500];
			sprintf(wfName,"wave_%i_%.0f_%.0f_%i",entryList[i],run,ch,iEvent);
			// hist->SetMaximum(300);
			// hist->SetMinimum(0);
			hist->Draw("L");
			hist->SetTitle(title.Data());

			theWF->Write(wfName,TObject::kOverwrite);
			delete theWF;

			cout << "filling\n" << endl;

			int fullEventSize = (int)channel->size();
			for (int k = 0; k < fullEventSize; k++)
			{
				if (channel->at(k) == ch)
				{
					cout << "channel match\n" << endl;
					gatSingle[0]->Fill(rawWFMax->at(k));
					gatSingle[1]->Fill(rawWFMin->at(k));
					gatSingle[2]->Fill(blrwfFMR0p1->at(k));
					gatSingle[3]->Fill(sgswfFMR0p1->at(k));
					gatSingle[4]->Fill(blrwfFMR1->at(k));
					gatSingle[5]->Fill(sgswfFMR1->at(k));
					gatSingle[6]->Fill(blrwfFMR3->at(k));
					gatSingle[7]->Fill(sgswfFMR3->at(k));
					gatSingle[8]->Fill(blrwfFMR10->at(k));
					gatSingle[9]->Fill(sgswfFMR10->at(k));
					gatSingle[10]->Fill(blrwfFMR20->at(k));
					gatSingle[11]->Fill(sgswfFMR20->at(k));
					gatSingle[12]->Fill(blrwfFMR50->at(k));
					gatSingle[13]->Fill(sgswfFMR50->at(k));
					gatSingle[14]->Fill(blrwfFMR80->at(k));
					gatSingle[15]->Fill(sgswfFMR80->at(k));
					gatSingle[16]->Fill(blrwfFMR90->at(k));
					gatSingle[17]->Fill(sgswfFMR90->at(k));
					gatSingle[18]->Fill(blrwfFMR97->at(k));
					gatSingle[19]->Fill(sgswfFMR97->at(k));
					gatSingle[20]->Fill(blrwfFMR99->at(k));
					gatSingle[21]->Fill(sgswfFMR99->at(k));
					gatSingle[22]->Fill(sgswft0->at(k));
					gatSingle[23]->Fill(TSCurrent50nsMax->at(k));
					gatSingle[24]->Fill(TSCurrent100nsMax->at(k));
					gatSingle[25]->Fill(TSCurrent200nsMax->at(k));
					gatSingle[26]->Fill(RawWFblSlope->at(k));
					gatSingle[27]->Fill(RawWFblOffset->at(k));
					gatSingle[28]->Fill(RawWFblSlopeUnc->at(k));
					gatSingle[29]->Fill(RawWFblOffsetUnc->at(k));
					gatSingle[30]->Fill(RawWFblChi2->at(k));
					gatSingle[31]->Fill(RawWFftSlope->at(k));
					gatSingle[32]->Fill(RawWFftSlopeUnc->at(k));
					gatSingle[33]->Fill(RawWFftOffset->at(k));
					gatSingle[34]->Fill(RawWFftOffsetUnc->at(k));
					gatSingle[35]->Fill(RawWFftChi2->at(k));
					gatSingle[36]->Fill(RawWFwholeWFLinFitSlope->at(k));
					gatSingle[37]->Fill(RawWFwholeWFLinFitSlopeUnc->at(k));
					gatSingle[38]->Fill(RawWFwholeWFLinFitOffset->at(k));
					gatSingle[39]->Fill(RawWFwholeWFLinFitOffsetUnc->at(k));
					gatSingle[40]->Fill(RawWFwholeWFLinFitChi2->at(k));
					gatSingle[41]->Fill(trapENF55us->at(k));
					gatSingle[42]->Fill(trapENF70us->at(k));
					gatSingle[43]->Fill(trapENF85us->at(k));
					gatSingle[44]->Fill(trapENF100us->at(k));
					gatSingle[45]->Fill(trapENF115us->at(k));
					gatSingle[46]->Fill(trapENF130us->at(k));
					gatSingle[47]->Fill(trapBL->at(k));
					gatSingle[48]->Fill(trap500nsMax->at(k));
					gatSingle[49]->Fill(trapENF500nsrt->at(k));
					gatSingle[50]->Fill(trap1usMax->at(k));
					gatSingle[51]->Fill(trapENF1usrt->at(k));
					gatSingle[52]->Fill(trap2usMax->at(k));
					gatSingle[53]->Fill(trapENF2usrt->at(k));
					gatSingle[54]->Fill(trap4usMax->at(k));
					gatSingle[55]->Fill(trapENF4usrt->at(k));
					gatSingle[56]->Fill(trap6usMax->at(k));
					gatSingle[57]->Fill(trapENF6usrt->at(k));
					gatSingle[58]->Fill(trap8usMax->at(k));
					gatSingle[59]->Fill(trapENF8usrt->at(k));
					gatSingle[60]->Fill(trapE->at(k));
					gatSingle[61]->Fill(trapEMin->at(k));
					gatSingle[62]->Fill(longGapTrapMax->at(k));
					gatSingle[63]->Fill(trap1ust0->at(k));
					gatSingle[64]->Fill(trapENF->at(k));
					gatSingle[65]->Fill(trapBL1us->at(k));
					gatSingle[66]->Fill(trapETailMin->at(k));
					gatSingle[67]->Fill(trirt50nsft10nsMax->at(k));
					gatSingle[68]->Fill(trirt100nsft10nsMax->at(k));
					gatSingle[69]->Fill(trirt200nsft10nsMax->at(k));
					gatSingle[70]->Fill(triFilMin->at(k));
					gatSingle[71]->Fill(trirt100nsft10nsIntegralW->at(k));
					gatSingle[72]->Fill(blrwfIntegralW->at(k));
					gatSingle[73]->Fill(smoothTrirt100nsft10nsMax->at(k));
					gatSingle[74]->Fill(energyCal->at(k));
					gatSingle[75]->Fill(trapECal->at(k));
					gatSingle[76]->Fill(trapEMinCal->at(k));
					gatSingle[77]->Fill(trapENFCal->at(k));
					gatSingle[78]->Fill(nlcblrwfSlope->at(k));
					gatSingle[79]->Fill(RawWFdcblSlope->at(k));
					gatSingle[80]->Fill(RawWFdcblOffset->at(k));
					gatSingle[81]->Fill(RawWFdcblSlopeUnc->at(k));
					gatSingle[82]->Fill(RawWFdcblOffsetUnc->at(k));
					gatSingle[83]->Fill(RawWFdcblChi2->at(k));
					gatSingle[84]->Fill(RawWFdcftSlope->at(k));
					gatSingle[85]->Fill(RawWFdcftOffset->at(k));
					gatSingle[86]->Fill(RawWFdcftSlopeUnc->at(k));
					gatSingle[87]->Fill(RawWFdcftOffsetUnc->at(k));
					gatSingle[88]->Fill(blrwf2ndDiffPeakValue->at(k));
					gatSingle[89]->Fill(blrwfNorm2ndDiffPeakValue->at(k));
					gatSingle[90]->Fill(blrwfSSF->at(k));
					gatSingle[91]->Fill(BLMoverMerr->at(k));
					gatSingle[92]->Fill(FTMoverMerr->at(k));
					gatSingle[93]->Fill(delOffset->at(k));
					gatSingle[94]->Fill(toe->at(k));
					gatSingle[95]->Fill(wfDCBits->at(k));
					gatSingle[95]->Fill(wfDC_Bit_0_SSSpikeBL->at(k));
					gatSingle[97]->Fill(wfDC_Bit_1_SSSpikePhys->at(k));
					gatSingle[98]->Fill(wfDC_Bit_2_EarlyTrigger->at(k));
					gatSingle[99]->Fill(wfDC_Bit_3_LateTrigger->at(k));
					gatSingle[100]->Fill(wfDC_Bit_4_PosSaturatedWFs->at(k));
					gatSingle[101]->Fill(wfDC_Bit_5_NegSaturatedWFs->at(k));
					// gatSingle[102]->Fill((double)run->at(k));
					// gatSingle[103]->Fill((double)startTime->at(k));
					// gatSingle[104]->Fill((double)stopTime->at(k));
					gatSingle[105]->Fill(energy->at(k));
					gatSingle[106]->Fill(channel->at(k));
					gatSingle[107]->Fill(timestamp->at(k));
					// gatSingle[108]->Fill((double)gatrev->at(k));
					// gatSingle[109]->Fill((double)EventDC1Bits->at(k));
					// gatSingle[110]->Fill(d2wfnoiseTagNorm->at(k));	// empty in input file
					// gatSingle[111]->Fill(d2wfDCPower->at(k));	// empty in input file
					gatSingle[112]->Fill(d2wf0MHzTo2MHzPower->at(k));
					gatSingle[113]->Fill(d2wf2MHzTo10MHzPower->at(k));
					gatSingle[114]->Fill(d2wf10MHzTo20MHzPower->at(k));
					gatSingle[115]->Fill(d2wf30MHzTo35MHzPower->at(k));
					gatSingle[116]->Fill(d2wf48MHzTo50MHzPower->at(k));
					gatSingle[117]->Fill(d2wf0MHzTo50MHzPower->at(k));
				}
			}
		}
	}


	// Make some double canvases and put them in a ROOT file
	for (int i = 0; i < (int)gatSingle.size(); i++)
	{
		TCanvas *c = new TCanvas(gatName[i].c_str(),"Bob Ross's Canvas",800,600);
		c->SetLogy();
		gatHist[i]->Draw();
		gatSingle[i]->Draw("same");
		c->Write(gatName[i].c_str(),TObject::kOverwrite);
	}
	h->Close();

}
*/