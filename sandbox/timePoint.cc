/*
	Un-sync'd Signal Card 11 detectors:
	P1D1 
	P1D2
	P1D3
	P5D4
*/

#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TChain.h>

#include "GATDataSet.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

void histVectTest();
void timePoint();

// int main(int argc, char** argv) 
int main()
{
	// histVectTest();
	timePoint();
}

void histVectTest()
{
	cout << "is this working?\n";
	
	// create the vector of histograms
	char hname[200];
	vector<TH1D*> histos;
	for (int i=0; i<3; i++) 
	{
		sprintf(hname,"hist_%i",i);
		histos.reserve(histos.size()+1);
		histos.push_back(new TH1D(hname,hname,10,0,10));
	}

	// fill the vector with data
	for (int i=0; i<3; i++) histos[i]->FillRandom("gaus",100);

	// draw, print, and delete canvas
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
	histos[0]->Draw();
	c1->Print("hist.C");
	cout << "this is working\n";
	delete c1;

	for (int i=0; i<3; i++) delete histos[i];
	histos.clear();
}

void timePoint() 
{
	// Input a run list
	// ifstream InputList("Debug.txt");	
	// ifstream InputList("LockFail.txt");	// 6220-6317
	// ifstream InputList("LockPass.txt");	// 5610-5710
	// if(!InputList.good()) {
    	// cout << "Couldn't open file!" << endl;
    	// return;
    // }

	string Name = "timePointRef_10984_10990";
	GATDataSet ds;
	cout << "Creating GatDS\n";
	ds.AddRunRange(10984,10990);

    // output a list for the DC framework
   	ofstream tps;
   	string outputFile = Name+".txt";
   	tps.open(outputFile.c_str());

	// int run = 0;
 	// while(!InputList.eof())
	// {
	// 	InputList >> run;
	// 	ds.AddRunNumber(run);
	// }
	TChain *g = ds.GetGatifiedChain();
	long entries = g->GetEntries();
	vector<double>* trapECal = 0;
	vector<double>* chan = 0;
	vector<double>* blrwfFMR50 = 0;
	g->SetBranchAddress("trapECal",&trapECal);
	g->SetBranchAddress("channel",&chan);
	g->SetBranchAddress("blrwfFMR50",&blrwfFMR50);

	// Which channel map does it use?
	MJTChannelMap *map = ds.GetChannelMap();

	cout << "Initialized. Found " << entries << " entries.\n";

	
    // Write histos into a ROOT file
    // Set up output files
	// Name.erase(Name.find_last_of("."),string::npos);
	// Name.erase(0,Name.find_last_of("\\/")+1);
	Char_t OutputFile[200];
	sprintf(OutputFile,"%s.root",Name.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE"); 	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."


	// Save a vector with all possible channels.
	vector<uint32_t> list;
	int numDetectors = (int)map->NumberOfRows();
	for(int i = 0; i < numDetectors; i++) 
	{
    	list.push_back(map->GetRowIndex(i)["kIDHi"].AsLong());
    	list.push_back(map->GetRowIndex(i)["kIDLo"].AsLong());
	}

	// Create a 2D vector of results.  1st entry is the channel number.
	int numChans = (int)list.size();
	vector< vector<double> > results(numChans, vector<double>()); 
    for (int i = 0; i < numChans; i++)  
    {
    	results[i].push_back(list[i]);	// channel number: results[i][0]
    	results[i].push_back(0);		// average time point: results[i][1]
    	results[i].push_back(0);		// counts above 100 keV: results[i][2]
    }

    // vector of histograms.
    char hname[200];
	vector<TH1D*> t50Hists;
	for (int i=0; i<numChans; i++) 
	{
		sprintf(hname,"hist_%i",i);
		t50Hists.reserve(t50Hists.size()+1);
		t50Hists.push_back(new TH1D(hname,hname,10000,0,20000));
	}

	// Loop over entries and calculate average time point for 
	// events over 100 keV
	int energyThresh = 100; // keV
	int ch = 0;
	double en = 0;
	double tp = 0;

	cout << "Scanning DS\n";
	for (long q = 0; q < entries; q++)
	{
		g->GetEntry(q);

		if (q % 50000 == 0) {
			printf("%.1f done\n",100*(double)q/entries);
		}

		// loop over channels in this entry
		for (int j = 0; j < (int)chan->size(); j++) 
		{ 
			ch = chan->at(j);
			en = trapECal->at(j);
			tp = blrwfFMR50->at(j);

			for (int k = 0; k < numChans; k++) 
			{
				if (ch == results[k][0] && en > energyThresh) 
				{
					results[k][1] += tp;
					results[k][2]++;
					t50Hists[k]->Fill(tp);
				}
			}
		}
	}

	// output
	for (int i = 0; i < numChans; i++)
	{
		// average time point
		if (results[i][2] > 0) results[i][1] = results[i][1]/results[i][2];
		else results[i][1] = -1;

		cout << map->GetDetectorPos(results[i][0]) << " ";
		cout << map->GetDetectorName(results[i][0]) << " ";
		cout << "\tCounts >100keV: " << results[i][2] << "\tAvg TP: " << results[i][1] << endl;
	}

	vector<TF1*> fits;
	printf("\n  Time Point Information:\n  Channel blrwfFMR50 3-sigma_hi 3-sigma_lo\n");
	char buffer[300];
	for (size_t j=0; j<t50Hists.size(); j++)
	{	
		fits.push_back(t50Hists[j]->GetFunction("gaus"));

		// try to fit if we have entries in this channel
		if (t50Hists[j]->GetEntries() > 0) {
			t50Hists[j]->Fit("gaus","q");
			fits[j] = t50Hists[j]->GetFunction("gaus");
			double hi = fits[j]->GetParameter(1)+3*fits[j]->GetParameter(2);
			double lo = fits[j]->GetParameter(1)-3*fits[j]->GetParameter(2);
			sprintf(buffer,"%.0f blrwfFMR50 %.2f %.2f\n",
				results[j][0],hi,lo);
		}
		else {
			sprintf(buffer,"%.0f blrwfFMR50 -1 -1\n",results[j][0]);
		}
		tps << buffer;
		
		sprintf(hname,"t50_%.0f",results[j][0]);		
		t50Hists[j]->Write(hname,TObject::kOverwrite);
	}	    	

	RootFile->Close();
	tps.close();
}