//----------------------------------------------------------------------------
// mjd_channel_quality.cc
// Usage: ./mjd_channel_quality <run_number>
// Clint Wiseman, USC/Majorana
//------------------------------------------------------------------------------
// 
// CHANNEL QUALITY
// 0: true/pass, 1: false/fail
// Bit 0: Channel enabled in ORCA
// Bit 1: HV off
// Bit 2: HV different by more than 1% from voltage @ prev. calibration
// Bit 3: FWHM at 2.6 MeV (most recent calibration) different by more than 5% from prev. calibration.
// Bit 4: “Noise rate” greater than 100 Hz.
// Bit 5: Threshold greater than 100 keV
// Bit 6: Gain change from previous and next calibration greater than 1%
// Bit 7: t50 time
// Bit 8: Live detector (> 0 events)
// Bit 9: Unstable baseline (MJSCDoc)
// Bit 10: Pulser monitor channel
//------------------------------------------------------------------------------
// Pulser Monitor Channels (as of Nov. 2015)
// Monitor 	Pulses these MB's		Monitors these detectors
// 644		M1-1alpha, M1-1delta	P1D1 P1D2 P1D3 P1D4 P2D2 P2D3 P2D4
// 612		M1-1gamma, M1-1beta		P3D1 P3D2 P3D3 P3D4
// 596		M1-2beta, M1-2gamma		P2D1 P4D3 P4D4 P4D5 P5D1 P5D2 P5D3 P5D4 P6D4
// 676		M1-2alpha, M1-2delta	P4D1 P4D2 P6D1 P6D2 P6D3 P7D1 P7D2 P7D3 P7D4
//------------------------------------------------------------------------------

#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <getopt.h>
#include <bitset>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <TStyle.h>

#include "GATDataSet.hh"
#include "KTree.h"

using namespace std;

//------------------------------------------------------------------------------
// PROCESSING:
// Put all the stuff requiring a scan over the gatified data into the same scan
// to reduce computation load.  Returns multiple results in a vector.
//
vector< vector<double> > process_data(GATDataSet *ds, vector<uint32_t> list, int verbose = 0)
{
	TChain *g = ds->GetGatifiedChain();
	int duration = ds->GetRunTime()/CLHEP::second;
	long entries = g->GetEntries();
	vector<double>* trapECal = 0;
	vector<double>* chan = 0;
	vector<double>* blrwfFMR50 = 0;
	g->SetBranchAddress("trapECal",&trapECal);
	g->SetBranchAddress("channel",&chan);
	g->SetBranchAddress("blrwfFMR50",&blrwfFMR50);

	// csc 4
	double energyThresh = 100;	// keV.
	
	// Store count info for all channels in "list" in a 2d vector, "results"
	int numChans = (int)list.size();
	vector< vector<double> > results(numChans, vector<double>()); 
    for (int i = 0; i < numChans; i++)  {
    	results[i].push_back(list[i]);	// channel number: results[i][0]
    	results[i].push_back(0);		// low-energy count information: results[i][1]
    	results[i].push_back(0);		// noise rate: results[i][2]
    	results[i].push_back(0);		// average time point: results[i][3]
    	results[i].push_back(0);		// counts above 100 keV: results[i][4]
    }

	// loop over entries
	int ch = 0;
	double en = 0;
	double tp = 0;
	if (verbose) printf("Now scanning run %i, %li entries.\n",ds->GetRunNumber(0),entries);
	for (long i = 0; i < entries; i++)
	{
		g->GetEntry(i);

		// loop over channels in this entry
		for (int j = 0; j < (int)chan->size(); j++) { 
			
			ch = chan->at(j);
			en = trapECal->at(j);
			tp = blrwfFMR50->at(j);

			for (int k = 0; k < numChans; k++) {
				if (ch == results[k][0] && en < energyThresh) 
				{
					results[k][1]++;
				}
				if (ch == results[k][0] && en > energyThresh) 
				{
					results[k][4]++;
					results[k][3] += tp;
				}
			}
		}
	}

	// output
	for (int i = 0; i < numChans; i++){
		
		// low-energy count information
		results[i][2] = results[i][1]/duration;
		
		// average time point
		if (results[i][4] > 0) {
			results[i][3] = results[i][3]/results[i][4];
		}
		else results[i][3] = 0;

		if (verbose) 
		{
			for (int j = 0; j < 5; j++) printf("i:%i j:%i - %-4.0f  ",i,j,results[i][j]);
			cout << endl;
		}
	}
	return results;
}


//------------------------------------------------------------------------------
// Bit 0: Channel enabled in ORCA
//
int check_csc_0(MJTChannelMap *map, MJTChannelSettings *set, uint32_t channel, int verbose = 0)
{
	bool isEnabled = set->GetBool("Enabled",channel);
	int crate = map->GetInt(channel,"kVME");
	int card = map->GetInt(channel,"kCardSlot");

	int sigChan;
	if (map->IsHiGain(channel)) sigChan = map->GetInt(channel,"kChanHi");
	else sigChan = map->GetInt(channel,"kChanLo");

	if (verbose) printf("0. Crate: %i  Card: %i  Signal Chan: %i\n",crate,card,sigChan);
	if (isEnabled) {
		if (verbose) cout << "   Enabled in ORCA" << endl << "   csc0: " << 0 << endl;

		return 0;
	}
	else {
		if (verbose) cout << "   Not enabled in ORCA" << endl << "   csc1: " << 1 << endl;
		return 1;
	}
}

//------------------------------------------------------------------------------
// Bit 1: HV Off
//
int check_csc_1(MJTChannelMap *map, MJTChannelSettings *set, uint32_t channel, int verbose = 0)
{
	// get HV channel info
	int hvcrate = map->GetInt(channel,"kHVCrate");
	int hvcard = map->GetInt(channel,"kHVCard");
	int hvchannel = map->GetInt(channel,"kHVChan");

	// compare max voltage to actual voltage.
	int maxvoltage = map->GetInt(channel,"kMaxVoltage");
	int actualvoltage = set->GetInt("targets",hvcrate,hvcard,hvchannel,"OREHS8260pModel");
	
	if (verbose) {
		printf("1. Max HV: %i  Target HV: %i  Crate: %i  Card: %i  Channel: %i\n",
			maxvoltage,actualvoltage,hvcrate,hvcard,hvchannel);
	}

	double ratio = (double)maxvoltage/actualvoltage;
	if (ratio <= 1.02 && ratio >= 0.98 && actualvoltage != 0) {
		if (verbose) cout << "   HV within limits" << endl << "   csc1: " << 0 << endl;
		return 0;
	}
	else if (actualvoltage == 0){
		if (verbose) cout << "   Detector not biased" << endl << "   csc1: " << 1 << endl;
		return 1;
	}
	else if (ratio >= 1.02 || ratio <= 0.98) {
		if (verbose) cout << "   HV outside limits!" << endl << "   csc1: " << 1 << endl;
		return 1;
	}
	else return 1;
}

//------------------------------------------------------------------------------
// Bit 2: HV different by more than 1% from voltage @ prev. calibration
//
void check_csc_2() {}

//------------------------------------------------------------------------------
// Bit 3: FWHM at 2.6 MeV (most recent calibration) different by more than 5% 
// from prev. calibration.
//
void check_csc_3() {}

//------------------------------------------------------------------------------
// Bit 4: “Noise rate” greater than 100 Hz.
//
// Relies on the vector created by process_data.
int check_csc_4(vector< vector<double> > results, uint32_t channel, int verbose = 0)
{
	for (size_t i = 0; i < results.size(); i++)
	{
		if ((int)channel == (int)results[i][0]) {
			if (verbose) cout << "4. \"Noise\" rate (under 100 keV): " << results[i][2] << endl;
			if (results[i][2] > 100) {
				if (verbose) cout << "   csc4: " << 1 << endl;
				return 1;
			}
			else {
				if (verbose) cout << "   csc4: " << 0 << endl;
				return 0;
			}
		}
	}
	return -1; // channel not found
}

//------------------------------------------------------------------------------
// Bit 5: Threshold greater than 100 keV
// The TRAP threshold settings are in MJTChannelSettings, if only I knew the conversion ...
void check_csc_5() {}

//------------------------------------------------------------------------------
// Bit 6: Gain change from previous and next calibration greater than 1%
void check_csc_6() {}

//------------------------------------------------------------------------------
// Bit 7: t50 time
int check_csc_7(vector< vector<double> > results, uint32_t channel, int verbose = 0)
{
	for (size_t i = 0; i < results.size(); i++)
	{
		if ((int)channel == (int)results[i][0]) {
			if (verbose) printf("7. (T50)  Total counts: %.0f  Avg. T50: %.2f\n",results[i][4],results[i][3]);
		// 	if (results[i][2] > 100) {
		// 		cout << "   csc4: " << 1 << endl;
		// 		return 1;
		// 	}
		// 	else {
		// 		cout << "   csc4: " << 0 << endl;
		// 		return 0;
		// 	}
		}
	}
	return 1;
}

//------------------------------------------------------------------------------
// 8: Live detector (> 0 events)
//
// This is an UNOFFICIAL BIT.  We are not guaranteed to have a hit in every detector
// of every run.  This bit is "information-only", i.e. to be used for trend-spotting.
int check_csc_8(vector< vector<double> > results, uint32_t channel, int verbose = 0)
{
	for (size_t i = 0; i < results.size(); i++)
	{
		if ((int)channel == (int)results[i][0]) {
			int counts = (int)results[i][1];
			if (verbose) cout << "8. Events this run: " << counts << endl;
			if (counts == 0) {
				if (verbose) cout << "   csc8: " << 1 << endl;
				return 1;
			}
			else {
				if (verbose) cout << "   csc8: " << 0 << endl;
				return 0;
			}
		}
	}
	if (verbose) cout << "   Channel not found!" << endl << "   csc8: " << 1 << endl;
	return 1;
}

//------------------------------------------------------------------------------
// ? 9: Unstable baseline (MJSCDoc)
void check_csc_9() {}

//------------------------------------------------------------------------------
// ? 10: Pulser monitor channel
void check_csc_10(uint32_t channel) {}

//------------------------------------------------------------------------------
int main(int argc, char** argv) 
{
	// get run from input line
	if (argc < 2) {
		cout << "Usage: mjd_channel_quality  <run>  <verbose flag (optional, default=0)>" << endl;
		return 1;
	}
	
	int run = atoi(argv[1]);
	int verbose = 0;
	if (argc == 3) verbose = atoi(argv[2]);
	
	if (verbose) printf("Checking channel quality of run %i .... \n",run);
	
	// Grab the data and settings for this run.
	GATDataSet *ds = new GATDataSet(run);
	MJTChannelMap *map = ds->GetChannelMap();
	MJTChannelSettings *set = ds->GetChannelSettings();
	//set->DumpSettings();
	
	// Save a vector with all channels from the ChannelMap
	vector<uint32_t> list;
	int numDetectors = (int)map->NumberOfRows();
	for(int i = 0; i < numDetectors; i++) 
	{
    	list.push_back(map->GetRowIndex(i)["kIDHi"].AsLong());
    	list.push_back(map->GetRowIndex(i)["kIDLo"].AsLong());
	}
	
	// If any enabled channels (e.g. pulser monitors) are not on the 
	// channel map, add them to the list vector.
	vector<uint32_t> en = set->GetEnabledIDList(); 	
	bool foundExtraChannel = false;
	for (size_t m = 0; m < en.size(); m++)
	{
		foundExtraChannel = false;
		for (size_t n = 0; n < list.size(); n++)
		{
			if (en[m]==list[n]) {
				continue;
			}
			// else foundExtraChannel = true;
		}
		if (foundExtraChannel) list.push_back(en[m]);
	}

	// This 2D vector holds the quality bits for each channel.
	int numChannels = numDetectors*2;
	// add pulser monitors to this
	vector< vector<int> > channelQuality(numChannels, vector<int>());

	// Process data
	vector< vector<double> > processResults = process_data(ds,list);
	
	// Channel selection checks.
    int c = 0;	// channelQuality vector index.
	for(vector<uint32_t>::iterator it1 = list.begin(); it1!=list.end(); ++it1) {
		
		int channel = *it1;
		string name = map->GetDetectorName(channel);
		string pos = map->GetDetectorPos(channel);

		if (verbose) cout << "Detector: " << name << "  Channel: " << channel << "  Position: " << pos << endl;

		channelQuality[c].push_back(channel);	// save channel number as first entry
		
		channelQuality[c].push_back( check_csc_0(map,set,channel,verbose) );
		
		channelQuality[c].push_back( check_csc_1(map,set,channel,verbose) );

		channelQuality[c].push_back( check_csc_4(processResults,channel,verbose) );

		channelQuality[c].push_back( check_csc_7(processResults,channel,verbose) );

		channelQuality[c].push_back( check_csc_8(processResults,channel,verbose) );

		if (verbose) cout << endl;
		
		c++;
	}	

	// OUTPUT : Read back channelQuality, and condense the bits down to a QVal word.
	char OutputName[200];
	sprintf(OutputName,"ChannelQuality_run%i.txt",run);
	ofstream OutFile(OutputName);
	cout << "--- MJD Channel Quality (Run " << run << ") --- " << endl;
	for (int i = 0; i < numChannels; i++) 
	{
		cout << map->GetDetectorPos(channelQuality[i][0]) << " ";
		cout << map->GetDetectorName(channelQuality[i][0]) << " ";
		
		if (channelQuality[i][0] %2 == 0) cout << "HG ";
		else cout << "LG ";
		
		int QVal = 0;	
		for (size_t j = 0; j < channelQuality[i].size(); j++) 
		{
			cout << channelQuality[i][j] << " ";
			if (j==0) cout << "\t";
			if (j>0) QVal += channelQuality[i][j]*pow(2,j);
		}
		cout << "\tquality " << QVal << endl;

		OutFile << channelQuality[i][0] << " " << QVal << endl;
	}
	cout << "-------------" << endl;

	OutFile.close();
	delete ds;
	return 0;
}
