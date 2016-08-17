// Generate a bunch of repeated code lines.
// Because I'm not a monkey.
// usage: clang++ -std=c++11 gatBranches.C -o gatBranches && ./gatBranches

#include <iostream>
#include <vector>

using namespace std;

void gatBranches();

int main()
{
	gatBranches();
}

void gatBranches()
{
	vector<string> gatList = {"rawWFMax",
	"rawWFMin",
	"blrwfFMR0p1",
	"sgswfFMR0p1",
	"blrwfFMR1",
	"sgswfFMR1",
	"blrwfFMR3",
	"sgswfFMR3",
	"blrwfFMR10",
	"sgswfFMR10",
	"blrwfFMR20",
	"sgswfFMR20",
	"blrwfFMR50",
	"sgswfFMR50",
	"blrwfFMR80",
	"sgswfFMR80",
	"blrwfFMR90",
	"sgswfFMR90",
	"blrwfFMR97",
	"sgswfFMR97",
	"blrwfFMR99",
	"sgswfFMR99",
	"sgswft0",
	"TSCurrent50nsMax",
	"TSCurrent100nsMax",
	"TSCurrent200nsMax",
	"RawWFblSlope",
	"RawWFblOffset",
	"RawWFblSlopeUnc",
	"RawWFblOffsetUnc",
	"RawWFblChi2",
	"RawWFftSlope",
	"RawWFftSlopeUnc",
	"RawWFftOffset",
	"RawWFftOffsetUnc",
	"RawWFftChi2",
	"RawWFwholeWFLinFitSlope",
	"RawWFwholeWFLinFitSlopeUnc",
	"RawWFwholeWFLinFitOffset",
	"RawWFwholeWFLinFitOffsetUnc",
	"RawWFwholeWFLinFitChi2",
	"trapENF55us",
	"trapENF70us",
	"trapENF85us",
	"trapENF100us",
	"trapENF115us",
	"trapENF130us",
	"trapBL",
	"trap500nsMax",
	"trapENF500nsrt",
	"trap1usMax",
	"trapENF1usrt",
	"trap2usMax",
	"trapENF2usrt",
	"trap4usMax",
	"trapENF4usrt",
	"trap6usMax",
	"trapENF6usrt",
	"trap8usMax",
	"trapENF8usrt",
	"trapE",
	"trapEMin",
	"longGapTrapMax",
	"trap1ust0",
	"trapENF",
	"trapBL1us",
	"trapETailMin",
	"trirt50nsft10nsMax",
	"trirt100nsft10nsMax",
	"trirt200nsft10nsMax",
	"triFilMin",
	"trirt100nsft10nsIntegralW",
	"blrwfIntegralW",
	"smoothTrirt100nsft10nsMax",
	"energyCal",
	"trapECal",
	"trapEMinCal",
	"trapENFCal",
	"nlcblrwfSlope",
	"RawWFdcblSlope",
	"RawWFdcblOffset",
	"RawWFdcblSlopeUnc",
	"RawWFdcblOffsetUnc",
	"RawWFdcblChi2",
	"RawWFdcftSlope",
	"RawWFdcftOffset",
	"RawWFdcftSlopeUnc",
	"RawWFdcftOffsetUnc",
	"blrwf2ndDiffPeakValue",
	"blrwfNorm2ndDiffPeakValue",
	"blrwfSSF",
	"BLMoverMerr",
	"FTMoverMerr",
	"delOffset",
	"toe",
	"wfDCBits",
	"wfDC_Bit_0_SSSpikeBL",
	"wfDC_Bit_1_SSSpikePhys",
	"wfDC_Bit_2_EarlyTrigger",
	"wfDC_Bit_3_LateTrigger",
	"wfDC_Bit_4_PosSaturatedWFs",
	"wfDC_Bit_5_NegSaturatedWFs",
	"run",
	"startTime",
	"stopTime",
	"energy",
	"channel",
	"timestamp",
	"gatrev",
	"EventDC1Bits",
	"d2wfnoiseTagNorm",
	"d2wfDCPower",
	"d2wf0MHzTo2MHzPower",
	"d2wf2MHzTo10MHzPower",
	"d2wf10MHzTo20MHzPower",
	"d2wf30MHzTo35MHzPower",
	"d2wf48MHzTo50MHzPower",
	"d2wf0MHzTo50MHzPower"};

	// cout << "size: " << gatList.size() << " entry 52: " << gatList[52] << endl;

	for (int i = 0; i < (int)gatList.size(); i++)
	{
		cout << "gatSkim->Draw(\"" << gatList[i] << ">>h" << i << "\",cutstr);\n";
		cout << "TH1F *h" << i << " = (TH1F*)gDirectory->Get(\"h" << i << "\");\n";
		cout << "h" << i << "->Write(\"" << gatList[i] << "\",TObject::kOverwrite);\n\n";
	}

}