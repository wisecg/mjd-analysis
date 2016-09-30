// data-quality.cc
// T. Caldwell and C. Wiseman, 9/10/2016
//
// Generates near-term analysis plots from built data
// and puts them into a ROOT file.  A wrapper program,
// data-report.py, can be used to make a PDF of the output.

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"

#include "MJTRun.hh"
#include "MGTEvent.hh"
#include "GATDataSet.hh"

using namespace std;

void halve_array(vector<double>& v, bool ave)
{
	for(unsigned i=0; i<v.size(); i+=2)
	{
		if(ave) v[i/2] = (v[i] + v[i+1]) / 2;
		else v[i/2] = v[i];
	}
	v.resize(v.size() / 2);
}

TGraph* graph(string name, vector<double>& x, vector<double>& y, string sx, string sy, int color)
{
	if(x.size() == 0 || x.size() != y.size())
	return NULL;
	TGraph* g = new TGraph(x.size(), &x.front(), &y.front());
	g->SetName(name.c_str());
	g->SetTitle(("#font[132]{" + name + "}").c_str());
	g->GetXaxis()->SetTitle(("#font[132]{" + sx + "}").c_str());
	g->GetYaxis()->SetTitle(("#font[132]{" + sy + "}").c_str());
	g->GetXaxis()->SetLabelFont(132);
	g->GetYaxis()->SetTitleFont(132);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.2);
	g->SetLineColor(color);
	g->SetMarkerColor(color);
	return g;
}

int main(int argc, char* argv[])
{
	if (argc < 2){
		cout << "Usage:\n ./data-quality [run number] (Optional: [last run] [first event] [last event]) \n\n";
		return 1;
	}
	int fRun = stoi(argv[1]);
	int lRun=0, fEvent=0, lEvent=0;
	if (argc > 2) lRun = stoi(argv[2]);
	if (argc > 4) {
		fEvent = stoi(argv[3]);
		lEvent = stoi(argv[4]);
	}
	printf("\n\ninputs: %i  %i  %i  %i\n",fRun,lRun,fEvent,lEvent);

	GATDataSet ds;
	if (lRun == 0) ds.AddRunNumber(fRun);
	else if (lRun != 0) ds.AddRunRange(fRun,lRun);

	// NOTE: analysis machine version may only have access to built data,
	//       PDSF version will have both built & reconstructed.
	TChain *built = ds.GetBuiltChain(false); // false: don't friend the recon chain
	// TChain *recon = ds.GetGatifiedChain();
	printf("Found %lli entries.\n",built->GetEntries());
	if (lEvent==0) lEvent = (int)built->GetEntries();

	MGTEvent *event = new MGTEvent();
	MJTRun *mjRun = new MJTRun();
    built->SetBranchAddress("event", &event);
	built->SetBranchAddress("run",&mjRun);

	MJTChannelMap *map = ds.GetChannelMap();
	MJTChannelSettings *set = ds.GetChannelSettings();
	vector<uint32_t> enab = set->GetEnabledIDList(); // gets ID's for the first run only

	const int ndet = map->NumberOfRows();	// all possible detectors in channel map
	std::map<int, int> index_map;  // map the channel number to channel map column number (ndet)
	for(int i=1; i<2000; i++)
		if(map->GetInt(i, "kIDHi") == i)
			index_map[i] = map->GetInt(i, "kSegmentNumber");
	// for (std::map<int,int>::iterator it = index_map.begin(); it != index_map.end(); ++it)
    	// cout << it->first << " , " << it->second << "\n";

	const int nbase = 200;	// number of baseline samples per waveform
	const unsigned max_array_size = 2 << 16;

	vector< vector<double> > t(ndet);
    vector< vector<double> > v(ndet);
    vector< vector<double> > vbase(ndet);
	vector<int> wf_prescale(ndet, 1);

	vector< vector<double> > t1(ndet);
	vector< vector<double> > vmaxima(ndet);
	vector< vector<double> > tdecimate(ndet);
	vector< vector<double> > vdecimate(ndet);

	// Figure out time offsets s/t the first event requested is at t = 0.
	double localTime = 0; 	// 1. time since first event in this run ("local run time")
	double chainTime = 0;	// 2. time since first event in the chain ("global run time")
	double unixTime = 0;	// 3. absolute unix time of event
	double timeSinceFirst = 0;  // 4. time since first event specified ("t=0 time")
	built->GetEntry(0);
	int iRun = mjRun->GetRunNumber();
	double chainRunStart = mjRun->GetStartTime();
	double chainEventOffset = event->GetTime()/1e9;
	double thisRunStart = chainRunStart;
	double thisEventOffset = chainEventOffset;
	built->GetEntry(fEvent);
	double firstRunStart = mjRun->GetStartTime();
	double firstEventOffset = event->GetTime()/1e9;
	double firstChainTime = (firstRunStart-chainRunStart)+firstEventOffset-chainEventOffset;

	for(int iev = fEvent; iev < lEvent; iev++)
	{
		if(iev % 10000 == 0) printf("%i\t%i\n", iev, (int) built->GetEntries());
		if(iev >= built->GetEntries()) { cout << "user's event range is beyond the end of the run!\n"; break; }
		built->GetEntry(iev);

		if (mjRun->GetRunNumber() != iRun){
			thisRunStart = mjRun->GetStartTime();
			thisEventOffset = event->GetTime()/1e9;
			iRun = mjRun->GetRunNumber();
		}
		localTime = event->GetTime()/1e9 - thisEventOffset;
		chainTime = (thisRunStart - chainRunStart) + localTime;
		unixTime = thisRunStart + localTime;
		timeSinceFirst = chainTime - firstChainTime;
		// printf("loc %.2f  chain %.2f  unix %.0f  tsf %.2f\n", localTime ,chainTime, unixTime, timeSinceFirst);
		double time = timeSinceFirst;

		for(int ich=0; ich<(int)event->GetNWaveforms(); ich++)
		{
			MGTWaveform *wf = event->GetWaveform(ich);
			int ch_id = wf->GetID();
			if (ch_id%2==1) continue;	// skip LG channels
			int ch_index = index_map[ch_id];  // detector number

			// get string and detector positions from channel map and ignore pulser tag channels
			string strName = map->GetString(ch_id,"kStringName");
			string detName = map->GetString(ch_id,"kDetectorName");
			char pos = strName.back();
			int strNum = pos - '0';
			int det = map->GetInt(ch_id,"kDetectorPosition"); // if det = 0, pulser tag channel
			string detector = strName+"D"+to_string(det);
			if (det == 0) continue;
			// printf("ch_index %i  chan %i  detector %s  detName %s\n",ch_index,ch_id,detector.c_str(),detName.c_str());

			TH1D* h = wf->GimmeHist(detector.c_str());

			for(int i=0; i<4; i++) // replace last 4 samples with 4 next-to-last
				h->SetBinContent(h->GetNbinsX()-i, h->GetBinContent(h->GetNbinsX()-4));

			double offset = wf->GetTOffset()/1e9; // time offset (in sec) of waveforms that are grouped into the same event
			double base = 0.0;
			double rms = 0.0;
			for(int i=4; i<nbase; i++) base += h->GetBinContent(i);
			base /= nbase;
			vbase[ch_index].push_back(base);

			for(int i=4; i<nbase; i++) rms += pow(h->GetBinContent(i) - base, 2);
			rms = sqrt(rms) / base;

			// max and min values of waveform (useful to have both when we trigger on pos/neg pulses)
			double wmax = -1.0e6;
			double wmin =  1.0e6;
			bool rail = false;

			// looping over the waveform
			for(int i = 1; i <= h->GetNbinsX(); i += wf_prescale[ch_index])
			{
				double wfTime = time + offset + (h->GetBinCenter(i))/1e9;

				t[ch_index].push_back(wfTime);
				v[ch_index].push_back(h->GetBinContent(i));

				if(t[ch_index].size() == max_array_size)
				{
					halve_array(t[ch_index], false);
					halve_array(v[ch_index], false);
					wf_prescale[ch_index] *= 2;
					tdecimate[ch_index].push_back(wfTime);
					vdecimate[ch_index].push_back(1);
					printf("halving arrays, detector %s.  prescale now %i\n",detector.c_str(),wf_prescale[ch_index]);
				}
				wmax = std::max(h->GetBinContent(i), wmax);
				wmin = std::min(h->GetBinContent(i), wmin);
			}
			delete h;

			// this gets the max value found during sampling, not the true max/min value, but maybe close enough.
			t1[ch_index].push_back(time);
			if(TMath::Abs(wmax) > TMath::Abs(wmin))
				vmaxima[ch_index].push_back(wmax);
			else
				vmaxima[ch_index].push_back(wmin);
		}
	}


	// ==================== OUTPUT ====================
	string runs = to_string(fRun);
	if (lRun != 0) runs += "_" + to_string(lRun);
	string fname = "data-quality_" + runs;
	TFile* outfile = new TFile((fname + ".root").c_str(), "recreate");
    outfile->cd();

	for(int i = 0; i < ndet; i++)
	{
		// get detector position and name for plots
		int dpos = 0;
		string detName = "";
		for (std::map<int,int>::iterator it2 = index_map.begin(); it2 != index_map.end(); ++it2)
		{
			if (it2->second == i) {
				int ch_id = it2->first;
				string strName = map->GetString(ch_id,"kStringName");
				char pos = strName.back();
				int strNum = pos - '0';
				dpos = map->GetInt(ch_id,"kDetectorPosition");
				detName = strName+"D"+to_string(dpos);
				break;
			}
		}
		printf("detector %s  size %lu\n",detName.c_str(),t[i].size());
		if(t[i].size() == 0) continue;

		int color = (dpos + 1) * 2;
		if(dpos > 3) color = 1;

		TGraph *gwf = graph(detName + "_gwf", t[i], v[i],"Time (sec)", "ADC", color);
		gwf->Write();

		TGraph *gmax = graph(detName + "_gmax", t1[i], vmaxima[i], "Time (sec)", "ADC", kRed);
		gmax->Write();

		TGraph *gbase = graph(detName + "_gbase", t1[i], vbase[i], "Time (sec)", "Baseline ADC", kBlue);
		gbase->Write();

		TGraph *decTimes = graph(detName + "_dec", tdecimate[i], vdecimate[i], "Time (sec)", "ADC", color);
		decTimes->SetMarkerStyle(kFullDotLarge);
		decTimes->SetMarkerColor(kRed);
		decTimes->Write();

	}
	outfile->Close();
    return 1;
}