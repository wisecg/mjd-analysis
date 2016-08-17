// LN Fill Tagging app.
// Clint Wiseman, USC/Majorana
//
// To use, edit the list of strings at the top of "main".
// Run with ./lnFills (-H, -R)
// A list of LN fill times with windows is written to "tagFile".
//
// 3/1/2016

#include "MJSlowControlsDoc.hh"
#include "GATDataSet.hh"

#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;    
using namespace MJDB;

// fill tag
void GetHistory(string start, string stop, string zone, string levelHistory, string runList);
std::vector<long> TagFills(string levelHistory, string tagFile, string tagList);

// analysis of fill tag
int *FindMatchingRuns(string runList, long fill, int *arr);
void CompareGeRate(string geFile, int fillNumber, long fill, long loUnix, long hiUnix, GATDataSet *ds);

static const char Usage[] =
"Usage: ./lnFills [options]\n"
"     -H: Get new DB History file\n"
"     -R: Do Ge Rate comparisions\n";

//===============================================================
int main(int argc, char** argv)
{
	string start = "2016/02/01 08:00:00";
	string stop = "2016/02/28 15:00:00";  
	string zone = "MDT";
	string levelHistory = "Cryo1LN_Feb2016.root";
	string runList = "RunList_Feb2016.txt";
	string tagFile = "LNTags_Feb2016.root";
	string tagList = "LNFills_Feb2016.txt";
	string geFile = "GeRate_Feb2016.root";

	// Since getting history and comparing rate can take a long time,
	// these options allow you to choose whether you do them.
	// 
	bool getHistory = false, compareRate = false;
	int c;
	while ((c = getopt (argc, argv, "HRh")) != -1)
    switch (c)
	{
	case 'H':
		getHistory=1;
		break;
	case 'R':
		compareRate=1;
		break;
	case 'h':
		cout << Usage;
		return 0;
	default:
		abort ();
	}

	// Start analysis
	//
	if (getHistory) GetHistory(start,stop,zone,levelHistory,runList);
	vector<long> fills = TagFills(levelHistory,tagFile,tagList);
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
		
		// Write to the tag file.
		LNTags << fills[i] << " " << loFill << " " << hiFill << " " 
			   << FirstRun << " " << FirstStart << " " << LastRun << " " << LastStop << endl;

		if (compareRate) {
			GATDataSet ds(FirstRun,LastRun);
			cout << "Created GATDataSet.\n";
			CompareGeRate(geFile,i,fills[i],FirstStart,LastStop,&ds);	
		}
	}
}

//===============================================================
void GetHistory(string start, string stop, string zone, string levelHistory, string runList)
{
	// LN level sensor history
	MJSlowControlsDoc doc;
	string variable = "Cryo1,Davis LN2";
	doc.SetDatabase(kHDB_SCM,kFeresaValues,variable,start,stop,zone,kIncDocs);
	doc.GrabHistory();	
	doc.RootifyDocument(levelHistory,variable);

	// Grab a list of run numbers that happened during these dates
	MJSlowControlsDoc doc2;
	doc2.SetDatabase(kHDB_DAQ2,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc2.GrabHistory();
	doc2.CreateRunList(runList,true);	// true: include time info
}

// If the value jumps by more than 10% over 10 minutes, 
// assume it's an LN fill.
//
std::vector<long> TagFills(string levelHistory, string tagFile, string tagList)
{
	cout << "Getting file" << endl;
	TFile *f = new TFile(levelHistory.c_str());	// read-only
	TTree *t = (TTree*)f->Get("Cryo1");
	int entries = t->GetEntries();
	int unixtime;
	double value;
	bool gapInData;
	t->SetBranchAddress("unixtime",&unixtime);
	t->SetBranchAddress("value",&value);
	t->SetBranchAddress("gapInData",&gapInData);
	t->GetEntry(0);
	cout << "set branches" << endl;
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

		// step is generally 60-80 seconds.
		step = unixtime - prevtime;
		deriv = (value - prevalue)/(double)step;
		derivative.push_back(deriv);

		prevalue = value;
		prevtime = unixtime;
	}
	cout << "closing ..." << endl;
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
	long lastFillTime = 0;
	int fills = 0;
	vector<long> fillTimes;
	MJSlowControlsDoc doc;
	if (derivative.size()!=timestamp.size()) {
		printf("Sizes not equal! deriv: %lu  timestamp: %lu\n",derivative.size(),timestamp.size());
		return vector<long>();
	}
	for (int i = 0; i < (int)derivative.size(); i++ ) 
	{
		// Fill tag condition.
		if (derivative[i] > 0.04 && timestamp[i] > lastFillTime + 30*60)	
		{  
			cout << (long)timestamp[i] << "  " << derivative[i] << "\t" << doc.GetGMTString((long)timestamp[i]) << endl;
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

	// Input a list of run numbers 
	ifstream RunList(runList.c_str());
	if(!RunList.good()) {
    	cout << "Couldn't open " << runList << endl;
    	return NULL;
    }

	int run;
	double duration;
	string startDate;
	string startTimeString;
	long startTime;
	string stopDate;
	string stopTimeString;
	long stopTime;

	int FirstRun = 0;
	long FirstStart = 0;
	long FirstStop = 0;
	int LastRun = 0;
	long LastStart = 0;
	long LastStop = 0;

	while (!RunList.eof()) {
		RunList >> run >> duration >> startDate >> startTimeString >> startTime >> stopDate >> stopTimeString >> stopTime;
		
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

	cout << "Fill time: " << fill << endl;
	printf("  First run: %i  start: %li  stop: %li\n",FirstRun,FirstStart,FirstStop);
	cout << "  Time before fill: " << (fill-FirstStart)/(double)60 << " minutes" << endl;

	printf("  Last run: %i  start: %li  stop: %li\n",LastRun,LastStart,LastStop);
	cout << "  Time after fill: " << (LastStop - fill)/(double)60 << " minutes" << endl;

	results[0]=FirstRun;
	results[1]=FirstStart;
	results[2]=LastRun;
	results[3]=LastStop;
	memcpy(arr,results,sizeof(results));
	
	return arr;
}

void CompareGeRate(string geFile, int fillNumber, long fill, long loUnix, long hiUnix, GATDataSet *d)
{
	cout << "Comparing Ge Rate for the fill at " << fill << endl;

	long loFill = fill-30*60;
	long hiFill = fill+10*60;

 	long noFillTime = fill+40*60;
 	long loNoFill = noFillTime-10*60;
 	long hiNoFill = noFillTime+30*60;

 	if (loFill < loUnix || hiNoFill > hiUnix) 
 	{
 		printf("Out of range ... loFill-loUnix = %li  hiNoFill-hiUnix = %li\n",loFill-loUnix,hiNoFill-hiUnix);
 		return;
 	}

 	GATDataSet *ds = d;
	TChain *g = ds->GetGatifiedChain();
	// int duration = ds->GetRunTime()/CLHEP::second;
	long entries = g->GetEntries();
	vector<double>* trapECal = 0;
	vector<double>* chan = 0;
	vector<double>* timestamp = 0;
	double startTime = 0;
	g->SetBranchAddress("trapECal",&trapECal);
	g->SetBranchAddress("channel",&chan);
	g->SetBranchAddress("startTime",&startTime);
	g->SetBranchAddress("timestamp",&timestamp);
	cout << "Scanning GATDataSet with " << entries << " entries.\n";

	// rate spectra (channels go from 578 - 696)
	int numSeconds = hiUnix-loUnix;
	char hname[200];
	sprintf(hname,"rate_%i",fillNumber);
	TH1D *rate = new TH1D(hname,"Events over threshold",numSeconds,loUnix,hiUnix);
	sprintf(hname,"rateByChan_%i",fillNumber);
 	TH2D *rateByChan = new TH2D(hname,"rateByChan",numSeconds,loUnix,hiUnix,125,575,700); 

 	// Energy spectra (3 keV/bin)
 	sprintf(hname,"fullSpectrum_%i",fillNumber);
 	TH1D *fullSpec = new TH1D(hname,"HG energy spectrum",1000,0,3000);
 	sprintf(hname,"fillSpectrum_%i",fillNumber);
 	TH1D *fillSpec = new TH1D(hname,"HG energy spectrum",1000,0,3000);	
 	sprintf(hname,"noFillSpectrum_%i",fillNumber);
 	TH1D *noFillSpec = new TH1D(hname,"HG energy spectrum",1000,0,3000);

 	// Loop over entries, fill histograms
 	long geTime = 0;
 	int ch = 0;
 	double en = 0;
 	double tm = 0;

	for (int i = 0; i < entries; i++)
	{
		g->GetEntry(i);

		tm = timestamp->at(0);
		geTime = (long)(tm/1E8) + (long)startTime;	// unix timestamp (seconds)
		
		if (i % 50000 == 0) printf("%.1f done\n",100*(double)i/entries);

		// loop over HG channels
		for (int j = 0; j < (int)chan->size(); j++) 
		{ 
			ch = chan->at(j);
			en = trapECal->at(j);

			if (ch%2 == 0 && en > 5) // 5keV floor
			{
				rate->Fill(geTime);
				rateByChan->Fill(geTime,ch);	// x,y
				fullSpec->Fill(en);

				if (geTime > loFill && geTime < hiFill) 
					fillSpec->Fill(en);

				if (geTime > loNoFill && geTime < hiNoFill) 
					noFillSpec->Fill(en);
			} 
		}
	}

	// output
	TFile *f = new TFile(geFile.c_str(),"UPDATE");
	rate->Write();

	char name[200];
	sprintf(name,"rateDates_%i",fillNumber);
	TCanvas *c1 = new TCanvas(name,"Bob Ross's Canvas",800,600);
	rate->GetXaxis()->SetTimeOffset(0,"gmt");
	rate->GetXaxis()->SetTimeFormat("%m/%d-%H:%m");
  	rate->GetXaxis()->SetTimeDisplay(1);
  	rate->GetXaxis()->SetNdivisions(-503);	//-503
  	rate->GetXaxis()->SetLabelOffset(0.02);
  	rate->GetXaxis()->SetTitleOffset(1.5);
  	rate->Draw();

  	double ymax = rate->GetMaximum();
  	cout << "Max rate was " << ymax << endl;

  	// Draw time window boundaries
  	//
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
	loNoFill_L->Draw();		

	TLine *hiNoFill_L = new TLine(hiNoFill,0,hiNoFill,ymax+10);
	hiNoFill_L->SetLineColor(kMagenta);
	hiNoFill_L->SetLineWidth(2.0);
	hiNoFill_L->Draw();		

  	c1->Write();
	fullSpec->Write();
	fillSpec->Write();
	noFillSpec->Write();
	rateByChan->Write();
	
	f->Close();

	delete c1;
	delete f;
	delete theFill;
	delete loFill_L;
	delete hiFill_L;
	delete loNoFill_L;
	delete hiNoFill_L;
	delete fullSpec;
	delete fillSpec;
	delete noFillSpec;
	delete rate;
	delete rateByChan;

	return;
}