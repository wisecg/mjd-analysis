// LN fill tagging app.
// Clint Wiseman, USC/Majorana
//
// 3/1/2016 - 1st version
// 9/18/2016 - Updated to work with M1 and M2

#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "MJSlowControlsDoc.hh"
#include "GATDataSet.hh"
#include "MJTRun.hh"

using namespace std;
using namespace MJDB;

void GetHistory(int mod, string start, string stop, string zone, string levelHistory, string runList);
vector<long> TagFills(int mod, string levelHistory, string tagFile);
int *FindMatchingRuns(string runList, long fill, int *arr);
void GeRate(string geFile, long fill, int FirstRun, int LastRun);

int main(int argc, char** argv)
{
	string start = "2016/08/24 00:00:00";
	string stop = "2016/09/15 00:00:00";
	string zone = "MDT";

	if (argc < 2) {
		cout << "Usage: ./ln_fill_tag [module number] [options]\n"
			 << "     -H: Get new DB History file\n"
			 << "     -R: Do Ge Rate comparisions\n";
		return 1;
	}
	int mod = 0;
	bool getHistory = false, compareRate = false;
	string opt1 = "", opt2 = "";
	if (argc >= 2)
		mod = stoi(argv[1]);
	if (argc == 3)
		opt1 = argv[2];
	if (argc == 4) {
		opt1 = argv[2];
		opt2 = argv[3];
	}
	if (opt1 == "-H" || opt2 == "-H") getHistory = true;
	if (opt1 == "-R" || opt2 == "-R") compareRate = true;

	string levelHistory = "Cryo"+to_string(mod)+"LN.root"; 	// GetHistory and TagFills
	string runList = "RunList_M"+to_string(mod)+".txt";		// GetHistory and FindMatchingRuns
	string tagFile = "LNTags_M"+to_string(mod)+".root";		// TagFills
	string tagList = "LNFills_M"+to_string(mod)+".txt";		// ** LN fill unix times **
	string geFile = "GeRate_M"+to_string(mod)+".root";		// GeRate

	if (getHistory) GetHistory(mod, start, stop, zone, levelHistory, runList);
	vector<long> fills = TagFills(mod, levelHistory, tagFile);
	ofstream LNTags(tagList.c_str());
	for (int i=0; i < (int)fills.size(); i++)
	{
		int results[4];
		FindMatchingRuns(runList,fills[i],results);
		int FirstRun = results[0];
		int FirstStart = results[1];
		int LastRun = results[2];
		int LastStop = results[3];
	 	long loFill = fills[i]-30*60;
		long hiFill = fills[i]+10*60;
		LNTags << fills[i] << " " << loFill << " " << hiFill << " "
			   << FirstRun << " " << FirstStart << " " << LastRun << " " << LastStop << endl;
		cout << fills[i] << endl;	// only unix times are needed by skim_mjd_data
		if (compareRate) GeRate(geFile, fills[i], FirstRun, LastRun);
	}
}


void GetHistory(int mod, string start, string stop, string zone, string levelHistory, string runList)
{
	cout << "Downloading Module " << mod << " LN level sensor history ...\n";
	MJSlowControlsDoc doc;
	string variable;
	if (mod==1) variable = "Cryo1,Davis LN2";
	if (mod==2) variable = "Thermosyphon,Davis LN2";
	doc.SetDatabase(kHDB_SCM,kFeresaValues,variable,start,stop,zone,kIncDocs);
	doc.GrabHistory();
	doc.RootifyDocument(levelHistory,variable);

	cout << "Downloading Module " << mod << " list of corresponding runs ...\n";
	MJSlowControlsDoc doc2;
	string whichDB;
	if (mod==1) whichDB = kHDB_DAQ2;
	if (mod==2) whichDB = kHDB_DAQ1;
	doc2.SetDatabase(whichDB,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc2.GrabHistory();
	doc2.CreateRunList(runList,true);	// true: include time info
}


vector<long> TagFills(int mod, string levelHistory, string tagFile)
{
	cout << "Tagging LN fills in file : " << levelHistory << endl;
	TFile *f = new TFile(levelHistory.c_str());
	TTree *t;
	if (mod == 1) t = (TTree*)f->Get("Cryo1");
	if (mod == 2) t = (TTree*)f->Get("Thermosyphon");
	int entries = t->GetEntries();
	int unixtime;
	double value;
	bool gapInData;
	t->SetBranchAddress("unixtime",&unixtime);
	t->SetBranchAddress("value",&value);
	t->SetBranchAddress("gapInData",&gapInData);
	t->GetEntry(0);
	vector<double> timestamp;
	vector<double> LNLevel;
	vector<double> derivative;
	double deriv;
	int step = 120; // default.
	int prevtime = 0;
	double prevalue = 0;
	for (int i = 0; i < entries; i++) {
		t->GetEntry(i);
		timestamp.push_back((double)unixtime);
		LNLevel.push_back(value);
		step = unixtime - prevtime; // generally 60-80 seconds.
		deriv = (value - prevalue)/(double)step;
		derivative.push_back(deriv);
		prevalue = value;
		prevtime = unixtime;
	}
	f->Close();

	TFile *f1 = new TFile(tagFile.c_str(),"RECREATE");
	TGraph *g = new TGraph(entries,&(timestamp[0]),&(LNLevel[0]));
	TGraph *d = new TGraph(entries,&(timestamp[0]),&(derivative[0]));
	TCanvas *c1 = new TCanvas("LevelAndDeriv","Bob Ross's Canvas",800,600);
	TCanvas *c2 = new TCanvas("TaggedDeriv","Bob Ross's Canvas",800,600);
	c1->Divide(0,2);
	c1->cd(2);
	g->Draw();
	c1->cd(1);
	d->Draw();
	c2->cd();
	d->Draw();

	// Tag the LN fills and return a vector with the times.
	// If the value jumps by more than 10% over 10 minutes, assume it's an LN fill.
	long lastFillTime = 0;
	int fills = 0;
	vector<long> fillTimes;
	MJSlowControlsDoc doc;
	if (derivative.size()!=timestamp.size())
	{
		printf("Sizes not equal! deriv: %lu  timestamp: %lu\n",derivative.size(),timestamp.size());
		return vector<long>();
	}
	for (int i = 0; i < (int)derivative.size(); i++)
	{
		if (derivative[i] > 0.04 && timestamp[i] > lastFillTime + 30*60)
		{
			cout << (long)timestamp[i] << "  " << derivative[i]
				 << "\t" << doc.GetGMTString((long)timestamp[i]) << endl;
			lastFillTime = (long)timestamp[i];
			fills++;
			fillTimes.push_back(timestamp[i]);
			TLine *line = new TLine(timestamp[i],-0.02,timestamp[i],0.1);
			line->SetLineColor(kRed);
			line->SetLineWidth(1.0);
			c1->cd(1);
			line->Draw();
			c2->cd();
			line->Draw();
		}
	}
	cout << "Found " << fills << " fills." << endl;
	g->Write();
	d->Write();
	c1->Write();
	c2->Write();
	f1->Close();
	return fillTimes;
}


int* FindMatchingRuns(string runList, long fill, int *arr)
{
	int results[4];
	ifstream RunList(runList.c_str());
	if(!RunList.good()) {
    	cout << "Couldn't open " << runList << endl;
    	return NULL;
    }
	int run=0, FirstRun=0, LastRun=0;
	long startTime=0, stopTime=0, FirstStart=0, FirstStop=0, LastStart=0, LastStop=0;
	double duration=0;
	string startDate, startTimeString, stopDate, stopTimeString;
	while (!RunList.eof())
	{
		RunList >> run >> duration >> startDate >> startTimeString
				>> startTime >> stopDate >> stopTimeString >> stopTime;
		if (startTime < fill-60*100) // 100 mins before fill
		{
			FirstRun = run;
			FirstStart = startTime;
			FirstStop = stopTime;
		}
		if (stopTime < fill+60*100) // 30 minutes after fill
		{
			LastRun = run;
			LastStart = startTime;
			LastStop = stopTime;
		}
	}
	// cout << "Fill time: " << fill << endl;
	// printf("  First run: %i  start: %li  stop: %li\n",FirstRun,FirstStart,FirstStop);
	// cout << "  Time before fill: " << (fill-FirstStart)/(double)60 << " minutes" << endl;
	// printf("  Last run: %i  start: %li  stop: %li\n",LastRun,LastStart,LastStop);
	// cout << "  Time after fill: " << (LastStop - fill)/(double)60 << " minutes" << endl;
	results[0]=FirstRun;
	results[1]=FirstStart;
	results[2]=LastRun;
	results[3]=LastStop;
	memcpy(arr,results,sizeof(results));
	return arr;
}


void GeRate(string geFile, long fill, int FirstRun, int LastRun)
{
	GATDataSet ds(FirstRun,LastRun);

	TFile *f = new TFile(geFile.c_str(),"RECREATE");
	long loFill = fill-15*60;
	long hiFill = fill+5*60;
	long noFillTime = fill-25*60;
	long loNoFill = noFillTime-15*60;
	long hiNoFill = noFillTime+5*60;

	TChain *g = ds.GetGatifiedChain();
	long entries = g->GetEntries();
	vector<double>* trapENFCal = 0;
	vector<double>* chan = 0;
	vector<double>* timestamp = 0;
	double run = 0;
	double startTime = 0;
	g->SetBranchAddress("run",&run);
	g->SetBranchAddress("trapENFCal",&trapENFCal);
	g->SetBranchAddress("channel",&chan);
	g->SetBranchAddress("startTime",&startTime);
	g->SetBranchAddress("timestamp",&timestamp);
	cout << "Scanning GATDataSet with " << entries << " entries.\n";

	TChain *built = ds.GetBuiltChain();
	long bEntries = built->GetEntries();
	MJTRun *bRun = new MJTRun();
	built->SetBranchAddress("run",&bRun);
	built->GetEntry(0);
	time_t loUnix = bRun->GetStartTime();
	built->GetEntry(bEntries-1);
	time_t hiUnix = bRun->GetStopTime();

	int numSeconds = hiUnix - loUnix;
	TH1D *rate = new TH1D("rate","Ge event rate",numSeconds,loUnix,hiUnix);
	TH1D *rate_5keVCut = new TH1D("rate_5keVCut","Ge event rate (5 keV cut)",numSeconds,loUnix,hiUnix);
 	TH2D *rateByChan = new TH2D("rateByChan","rateByChan",numSeconds,loUnix,hiUnix,125,575,700);
 	TH1D *fullSpec = new TH1D("fullSpectrum","HG energy spectrum",1000,0,3000);
 	TH1D *fillSpec = new TH1D("fillSpectrum","HG energy spectrum",1000,0,3000);
 	TH1D *noFillSpec = new TH1D("noFillSpectrum","HG energy spectrum",1000,0,3000);

 	long geTime = 0;
 	int ch = 0;
 	double en = 0, tm = 0;
	for (int i = 0; i < entries; i++)
	{
		g->GetEntry(i);
		tm = timestamp->at(0);
		geTime = (long)(tm/1E8) + (long)startTime;	// unix timestamp (seconds)
		if (i % 50000 == 0)
			printf("%.1f done.\n",100*(double)i/entries);
		// printf("run: %.0f  timestamp: %.0f  timestamp/1E8: %.2f  unix: %li\n",run,tm,tm/1E8,geTime);

		// loop over HG channels
		for (int j = 0; j < (int)chan->size(); j++)
		{
			ch = chan->at(j);
			en = trapENFCal->at(j);
			if (ch%2 == 0) // no energy floor
			{
				rate->Fill(geTime);
				if (en > 5)
					rate_5keVCut->Fill(geTime);
				rateByChan->Fill(geTime,ch);	// x,y
				fullSpec->Fill(en);
				if (geTime > loFill && geTime < hiFill)
					fillSpec->Fill(en);
				if (geTime > loNoFill && geTime < hiNoFill)
					noFillSpec->Fill(en);
			}
		}
	}
	rate->Write();
	rate_5keVCut->Write();
	rateByChan->Write();
	fullSpec->Write();
	fillSpec->Write();
	noFillSpec->Write();

	// rate plot with fill tag window
	char name[200];
	sprintf(name,"rateDates");
	TCanvas *c1 = new TCanvas(name,"Bob Ross's Canvas",800,600);
	rate->GetXaxis()->SetTimeOffset(0,"gmt");
	rate->GetXaxis()->SetTimeFormat("%m/%d-%H:%M");
	rate->GetXaxis()->SetTimeDisplay(1);
	rate->GetXaxis()->SetNdivisions(-505);	//-503
	rate->GetXaxis()->SetLabelOffset(0.02);
	rate->GetXaxis()->SetTitleOffset(1.5);
	rate->Draw();
	double ymax = rate->GetMaximum();
	cout << "Max rate was " << ymax << endl;

	TLine *theFill = new TLine(fill,0,fill,ymax+10);
	theFill->SetLineColor(kRed);
	theFill->SetLineWidth(2.0);
	theFill->Draw();
	TLine *loFill_L = new TLine(loFill,0,loFill,ymax+10);
	loFill_L->SetLineColor(kGreen);
	loFill_L->SetLineWidth(2.0);
	loFill_L->Draw();
	TLine *hiFill_L = new TLine(hiFill,0,hiFill,ymax+10);
	hiFill_L->SetLineColor(kGreen);
	hiFill_L->SetLineWidth(2.0);
	hiFill_L->Draw();
 	TLine *loNoFill_L = new TLine(loNoFill,0,loNoFill,ymax+10);
	loNoFill_L->SetLineColor(kMagenta);
	loNoFill_L->SetLineWidth(2.0);
	// loNoFill_L->Draw();
	TLine *hiNoFill_L = new TLine(hiNoFill,0,hiNoFill,ymax+10);
	hiNoFill_L->SetLineColor(kMagenta);
	hiNoFill_L->SetLineWidth(2.0);
	// hiNoFill_L->Draw();
	c1->Write();

	// rate plot with fill and no fill windows drawn
	sprintf(name,"rateDates_bothWindows");
	TCanvas *c2 = new TCanvas(name,"Bob Ross's Canvas",800,600);
	rate->GetXaxis()->SetTimeOffset(0,"gmt");
	rate->GetXaxis()->SetTimeFormat("%m/%d-%H:%M");
	rate->GetXaxis()->SetTimeDisplay(1);
	rate->GetXaxis()->SetNdivisions(-505);	//-503
	rate->GetXaxis()->SetLabelOffset(0.02);
	rate->GetXaxis()->SetTitleOffset(1.5);
	rate->Draw();
	theFill->Draw();
	loFill_L->Draw();
	hiFill_L->Draw();
	loNoFill_L->Draw();
	hiNoFill_L->Draw();
	c2->Write();

	// rate > 5kev rate plot with fill and no fill windows drawn
	sprintf(name,"5kev_rateDates_bothWindows");
	TCanvas *c3 = new TCanvas(name,"Bob Ross's Canvas",800,600);
	rate_5keVCut->GetXaxis()->SetTimeOffset(0,"gmt");
	rate_5keVCut->GetXaxis()->SetTimeFormat("%m/%d-%H:%M");
	rate_5keVCut->GetXaxis()->SetTimeDisplay(1);
	rate_5keVCut->GetXaxis()->SetNdivisions(-505);	//-503
	rate_5keVCut->GetXaxis()->SetLabelOffset(0.02);
	rate_5keVCut->GetXaxis()->SetTitleOffset(1.5);
	rate_5keVCut->Draw();
	theFill->Draw();
	loFill_L->Draw();
	hiFill_L->Draw();
	loNoFill_L->Draw();
	hiNoFill_L->Draw();
	c3->Write();

	// rate by channel with fill window only
	sprintf(name,"rateByChan_fill");
	TCanvas *c4 = new TCanvas(name,"Bob Ross's Canvas",800,600);
	rateByChan->GetXaxis()->SetTimeOffset(0,"gmt");
	rateByChan->GetXaxis()->SetTimeFormat("%m/%d-%H:%M");
	rateByChan->GetXaxis()->SetTimeDisplay(1);
	rateByChan->GetXaxis()->SetNdivisions(-505);	//-503
	rateByChan->GetXaxis()->SetLabelOffset(0.02);
	rateByChan->GetXaxis()->SetTitleOffset(1.5);
	rateByChan->Draw("COLZ");
	TLine *theFill_chan = new TLine(fill,575,fill,700);
	theFill_chan->SetLineColor(kRed);
	theFill_chan->SetLineWidth(2.0);
	theFill_chan->Draw();
	TLine *loFill_chan = new TLine(loFill,575,loFill,700);
	loFill_chan->SetLineColor(kGreen);
	loFill_chan->SetLineWidth(2.0);
	loFill_chan->Draw();
	TLine *hiFill_chan = new TLine(hiFill,575,hiFill,700);
	hiFill_chan->SetLineColor(kGreen);
	hiFill_chan->SetLineWidth(2.0);
	hiFill_chan->Draw();
	c4->Write();
	f->Close();
}
