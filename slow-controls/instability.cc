#include "MJSlowControlsDoc.hh"
using namespace std;
using namespace MJDB;

void CompareGeRate();

int main(int argc, char** argv)
{
	string start = "2016/02/09 15:00:00"; 
	string stop = "2016/02/09 17:00:00";  
	string zone = "MDT";
	string variable = "";

	// Find runs matching this date
	MJSlowControlsDoc doc;
	doc.SetDatabase(kHDB_DAQ2,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc.GrabHistory();
	doc.CreateRunList("MatchingRuns.txt",false);	// if true, includes time info
	
	GATDataSet *ds = new GATDataSet();
	ifstream InputList("MatchingRuns.txt");
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }
    while(!InputList.eof())
    {
		InputList >> run;
		ds->AddRunNumber(run);
	}

	CompareGeRate(ds);

	// //P4D1	B8455
	// MJSlowControlsDoc doc1;
	// variable = "Baseline Voltage,B8455";
	// doc1.SetDatabase(kHDB_DetHist,kFeresaValues,variable,start,stop,zone,kIncDocs);
	// doc1.GrabHistory();	
	// doc1.RootifyDocument("BaselineInfo.root",variable,"P4D1_B8455");
	
	// //P4D5	B8469
	// MJSlowControlsDoc doc2;
	// variable = "Baseline Voltage,B8469";
	// doc2.SetDatabase(kHDB_DetHist,kFeresaValues,variable,start,stop,zone,kIncDocs);
	// doc2.GrabHistory();	
	// doc2.RootifyDocument("BaselineInfo.root",variable,"P4D5_B8469",2);
}


void CompareGeRate(GATDataSet *d)
{
 	GATDataSet *ds = d;
	TChain *g = ds->GetGatifiedChain();
	int duration = ds->GetRunTime()/CLHEP::second;
	long entries = g->GetEntries();
	vector<double>* trapECal = 0;
	vector<double>* chan = 0;
	vector<double>* timestamp = 0;
	double startTime = 0;
	g->SetBranchAddress("trapECal",&trapECal);
	g->SetBranchAddress("channel",&chan);
	g->SetBranchAddress("startTime",&startTime);
	g->SetBranchAddress("timestamp",&timestamp);

	// rate spectra (channels go from 578 - 696)
	TH1D *rate = new TH1D("rate","Events over threshold",7200,1451857087,1451864290);
 	TH2D *rateByChan = new TH2D("rateByChan","rateByChan",7200,1451857087,1451864290,125,575,700); 

 	// Energy spectra (3 keV/bin)
 	TH1D *fullSpec = new TH1D("fullSpectrum","HG energy spectrum",1000,0,3000);

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
		
		// if (i % 10000 == 0) {
		// 	cout << tm/1E8 << " + " << (long)startTime << " = " << geTime << endl;
		// }

		// loop over HG channels
		for (int j = 0; j < (int)chan->size(); j++) 
		{ 
			ch = chan->at(j);
			en = trapECal->at(j);

			if (ch%2 == 0 && en > 2) // 2keV floor
			{
				rate->Fill(geTime);
				rateByChan->Fill(geTime,ch);	// x,y
				fullSpec->Fill(en);
			} 
		}
	}

	// output
	TFile *f = new TFile("GeRateTest.root","RECREATE");
	rate->Write();
	fullSpec->Write();
	rateByChan->Write();
	
	// TCanvas *c1 = new TCanvas("rateByChan","Bob Ross's Canvas",800,600);
	// c1->cd(0);
	// rateByChan->Draw("COLZ");
	// c1->Write();
	
	f->Close();

	return;
}
