//----------------------------------------------------------------------------
// LN fill tagging app.
// To use, edit the list of strings at the top of "main".
// A list of LN fill times with windows is written to "tagList".
//
// Usage: ./ln_fill_tag (-H -R)
// Clint Wiseman, USC/Majorana
// 3/1/2016
//-----------------------------------------------------------------------------

#include <unistd.h>

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
std::vector<long> TagFills(string levelHistory, string tagFile);

// analysis of fill tag
int *FindMatchingRuns(string runList, long fill, int *arr);
void CompareGeRate(string geFile, int fillNumber, long fill, long loUnix, long hiUnix, GATDataSet *ds);
void GeRate(string geFile);

static const char Usage[] =
"Usage: ./ln_fill_tag [options]\n"
"     -H: Get new DB History file\n"
"     -R: Do Ge Rate comparisions\n"
"     -T: Make LN fill list from existing DB History file\n"
"\n";

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
	/*
		1st DS0 run: 2580 Jun 30th, 2015
		Las DS0 run: 6963 Sep 23rd, 2015

		1st DS1 run: 9455  2016-1-13-P3KJR_Run9455
		Last DS1 run: 13368  2016-4-14-P3KJR_Run13368

		Ralph's fill: April 8th, 20:37 (M(S/D)T)
	*/

	// DS-0
	// string start = "2015/06/30 08:00:00";
	// string stop = "2015/09/23 23:00:00";
	// string zone = "MDT";
	// string levelHistory = "Cryo1LN_DS0.root";
	// string runList = "RunList_DS0.txt";
	// string tagFile = "LNTags_DS0.root";
	// string tagList = "LNFills_DS0.txt";
	// string geFile = "GeRate_DS0.root";

	// DS-1
	string start = "2016/01/12 00:00:00";
	string stop = "2016/04/15 00:00:00";
	string zone = "MDT";
	string levelHistory = "Cryo1LN_DS1.root";
	string runList = "RunList_DS1.txt";
	string tagFile = "LNTags_DS1.root";
	string tagList = "LNFills_DS1.txt";
	string geFile = "GeRate_DS1.root";

	// Since getting history and comparing rate can take a long time,
	// these options allow you to choose whether you do them.
	//
	bool getHistory=0,compareRate=0,justTagFills=0,DoGeRate=0;
	int c;
	while ((c = getopt (argc, argv, "HRTGh")) != -1)
    switch (c)
	{
	case 'H':
		getHistory=1;
		break;
	case 'T':
		justTagFills=1;
		break;
	case 'R':
		compareRate=1;
		break;
	case 'G':
		DoGeRate=1;
		break;
	case 'h':
		cout << Usage;
		return 0;

	default:
		abort ();
	}

	// Start analysis
	//
	// if (getHistory) GetHistory(start,stop,zone,levelHistory,runList);
	// vector<long> fills = TagFills(levelHistory,tagFile);
	// ofstream LNTags(tagList.c_str());
	// for (int i=1; i < (int)fills.size(); i++)
	// // for (int i = 1; i < 2; i++)
	// {
	// 	int results[4];
	// 	FindMatchingRuns(runList,fills[i],results);
	//
	// 	int FirstRun = results[0];
	// 	int FirstStart = results[1];
	// 	int LastRun = results[2];
	// 	int LastStop = results[3];
	// 	long loFill = fills[i]-15*60;
	// 	long hiFill = fills[i]+5*60;
	//
	// 	printf("FMR results:\n fill %i  FirstRun %i  FirstStart %i  LastRun %i  LastStop %i  loFill %li  hiFill %li\n",
	// 		fills[i],FirstRun,FirstStart,LastRun,LastStop,loFill,hiFill);
	//
	// 	// Write to the tag file.
	// 	// LNTags << fills[i] << " " << loFill << " " << hiFill << " "
	// 		   // << FirstRun << " " << FirstStart << " " << LastRun << " " << LastStop << endl;
	//
	// 	// Format is unix timestamps: fill time, 15 minutes before, 5 minutes after.
	// 	LNTags << fills[i] << " " << loFill << " " << hiFill << endl;
	//
	// 	if (compareRate) {
	// 		GATDataSet ds(FirstRun,LastRun);
	// 		cout << "Created GATDataSet.\n";
	// 		CompareGeRate(geFile,i,fills[i],FirstStart,LastStop,&ds);
	// 	}
	// }

	// This is a shortcut plotting function.  See below.
	if (DoGeRate) GeRate(geFile);
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// If the value jumps by more than 10% over 10 minutes,
// assume it's an LN fill.
//
std::vector<long> TagFills(string levelHistory, string tagFile)
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

	// write the bare graphs
	g->Write();
	d->Write();

	TCanvas *c1 = new TCanvas("LevelAndDeriv","Bob Ross's Canvas",800,600);
	TCanvas *c2 = new TCanvas("TaggedDeriv","Bob Ross's Canvas",800,600);
	c1->Divide(0,2);
	c1->cd(1);

	// LN level
	g->GetXaxis()->SetTimeOffset(0,"gmt");
	g->GetXaxis()->SetTimeFormat("%m/%d/%y");
  	g->GetXaxis()->SetTimeDisplay(1);
  	g->GetXaxis()->SetNdivisions(-506);	//-503
  	g->GetXaxis()->SetLabelOffset(0.02);
  	g->GetXaxis()->SetTitleOffset(1.5);
  	g->SetTitle("LN Level");
  	g->GetXaxis()->SetLabelSize(0.06);
  	g->GetYaxis()->SetLabelSize(0.06);


	g->Draw();
	c1->cd(2);

	// Derivative
	d->GetXaxis()->SetTimeOffset(0,"gmt");
	d->GetXaxis()->SetTimeFormat("%m/%d/%y");
  	d->GetXaxis()->SetTimeDisplay(1);
  	d->GetXaxis()->SetNdivisions(-506);	//-503
  	d->GetXaxis()->SetLabelOffset(0.02);
  	d->GetXaxis()->SetTitleOffset(1.5);
  	d->SetTitle("LN Derivative");
  	d->GetXaxis()->SetLabelSize(0.06);
  	d->GetYaxis()->SetLabelSize(0.06);

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
			// cout << (long)timestamp[i] << "  " << derivative[i] << "\t" << doc.GetGMTString((long)timestamp[i]) << endl;
			lastFillTime = (long)timestamp[i];
			fills++;
			fillTimes.push_back(timestamp[i]);
			TLine *line = new TLine(timestamp[i],-0.02,timestamp[i],0.1);
			line->SetLineColor(kRed);
			line->SetLineWidth(1.0);

			c1->cd(2);
			line->Draw();
			c2->cd();
			line->Draw();

		}
	}
	cout << "Found " << fills << " fills." << endl;

	c1->Write();
	c2->Write();
	f1->Close();

	return fillTimes;
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
void CompareGeRate(string geFile, int fillNumber, long fill, long loUnix, long hiUnix, GATDataSet *d)
{
	cout << "Comparing Ge Rate for the fill at " << fill << endl;

	// FirstRun 13203  FirstStart 1460161157  LastRun 13205  LastStop 1460171959

 	// fill 1460169429  FirstRun 13203  FirstStart 1460161157  LastRun 13205  LastStop 1460171959  loFill 1460168529  hiFill 1460169729

 	// fill: 1460169429
 	// loUnix = FirstStart = 1460161157
 	// hiUnix = LastStop = 1460171959

	long loFill = fill-14*60;
	long hiFill = fill+4*60;

 	long noFillTime = fill+10*60;
 	long loNoFill = noFillTime-10*60;
 	long hiNoFill = noFillTime+20*60;

	printf("loFill %li  loUnix %li  hiNoFill %li  hiUnix %li\n",loFill,loUnix,hiNoFill,hiUnix);


 	if ((loFill < loUnix) || (hiNoFill > hiUnix))
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

void GeRate(string geFile)
{


	// This made the plot for Clara.
	// FirstRun 13203  FirstStart 1460161157  LastRun 13205  LastStop 1460171959
	// long fill = 1460169429;
	// GATDataSet *ds = new GATDataSet(13204,13205);
	// long loUnix=1460164758;
	// long hiUnix=1460171959;

	// This is looking at other ones for vince
	// A few random DS-1 fills.
	// fill 1454014971  FirstRun 10112  FirstStart 1454008578  LastRun 10112  LastStop 1454008651
	// fill 1454107981  FirstRun 10139  FirstStart 1454101428  LastRun 10141  LastStop 1454112234
	// fill 1454219392  FirstRun 10170  FirstStart 1454212245  LastRun 10173  LastStop 1454223215
	// done fill 1454332307  FirstRun 10202  FirstStart 1454324094  LastRun 10204  LastStop 1454334900
	// fill 1454447212  FirstRun 10238  FirstStart 1454440940  LastRun 10260  LastStop 1454452363
	// fill 1454559953  FirstRun 10278  FirstStart 1454511200  LastRun 10278  LastStop 1454514800
	// fill 1454679496  FirstRun 10292  FirstStart 1454672984  LastRun 10294  LastStop 1454683787
	long fill = 1454107981;
	GATDataSet *ds = new GATDataSet(10139,10141);
	long loUnix=1454101428;
	long hiUnix=1454112234;

	// output
	TFile *f = new TFile("DS1GeRatePlot_02.root","RECREATE");

	long loFill = fill-15*60;
	long hiFill = fill+5*60;
	long noFillTime = fill-25*60;
	long loNoFill = noFillTime-15*60;
	long hiNoFill = noFillTime+5*60;

	TChain *g = ds->GetGatifiedChain();
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

	// histograms
	int numSeconds = hiUnix-loUnix;
	char hname[200];

	sprintf(hname,"rate");
	TH1D *rate = new TH1D(hname,"Ge event rate",numSeconds,loUnix,hiUnix);

	sprintf(hname,"rate_5keVCut");
	TH1D *rate_5keVCut = new TH1D(hname,"Ge event rate (5 keV cut)",numSeconds,loUnix,hiUnix);

	sprintf(hname,"rateByChan");
 	TH2D *rateByChan = new TH2D(hname,"rateByChan",numSeconds,loUnix,hiUnix,125,575,700);

 	sprintf(hname,"fullSpectrum");
 	TH1D *fullSpec = new TH1D(hname,"HG energy spectrum",1000,0,3000);

 	sprintf(hname,"fillSpectrum");
 	TH1D *fillSpec = new TH1D(hname,"HG energy spectrum",1000,0,3000);

 	sprintf(hname,"noFillSpectrum");
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

		if (i % 50000 == 0) printf("%.1f done.\n",100*(double)i/entries);

		// verify timestamps (OK!)
		// printf("run: %.0f  timestamp: %.0f  timestamp/1E8: %.2f  unix: %li\n",run,tm,tm/1E8,geTime);

		// // loop over HG channels
		for (int j = 0; j < (int)chan->size(); j++)
		{
			ch = chan->at(j);
			en = trapENFCal->at(j);

			if (ch%2 == 0) // no energy floor
			{
				rate->Fill(geTime);

				if (en > 5) rate_5keVCut->Fill(geTime);

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

	// ==============================
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

	// ==============================
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

	// ==============================
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

	// ==============================
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
