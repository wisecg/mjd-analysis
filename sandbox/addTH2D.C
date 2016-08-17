// written to be implemented in channelSel/reportQuality.
// C. Wiseman

#include <tuple>

void addTH2D()
{
	TFile *f0 = new TFile("~/dev/channelSel/reportQuality/output/DS1_08.root");
	TFile *f1 = new TFile("~/dev/channelSel/reportQuality/output/DS1_09.root");

	TH2D *h0 = (TH2D*)f0->Get("664_EnergyVsRun");
	int xmin = h0->GetXaxis()->GetXmin();
	int h0Bins = h0->GetNbinsX();

	TCanvas *c0 = new TCanvas("c0","",800,600);
	h0->Draw("COLZ");

	TH2D *h1 = (TH2D*)f1->Get("664_EnergyVsRun");
	int xmax = h1->GetXaxis()->GetXmax();
	int h1Bins = h1->GetNbinsX();

	TCanvas *c1 = new TCanvas("c1","",800,600);
	h1->Draw("COLZ");

	int xBins = xmax-xmin;
	TH2D *test = new TH2D("","",xBins,xmin,xmax,300,0,3000);

	printf("reg method: bins: %i  xmin: %i  xmax: %i\n",xBins,xmin,xmax);

	// DS1_08 output:
	// th2d range: bins: 15  xmin: 9764  xmax: 9779
	// addTH2D output:
	// bins: 15  xmin: 9764  xmax: 9779

	// VERIFIED: X AND Y BINS MATCH.

	// loop over each x bin (run number)


	// create a vector tuple to hold the three numbers
	vector<tuple<double,double,double>> triplets;

	int numPushBacks = 0;

	TH1 *hbins0[xBins];
	for (int i = 0; i <= h0Bins; i++) {
		hbins0[i] = h0->ProjectionY(Form("bin%d",i),i,i); // name,firstxbin,lastxbin
		int size = hbins0[i]->GetNbinsX();
		int run = h0->GetXaxis()->GetBinCenter(i);
		for (int j = 0; j < size; j++){
			int w = hbins0[i]->GetBinContent(j);
			if (w > 0){
				double ene = hbins0[i]->GetBinCenter(j);
				if (ene < 0) ene = 1;
				test->Fill(run,ene,w);

				tuple<double,double,double> temp = make_tuple((double)run,(double)ene,(double)w);
				// cout << get<0>(temp) << " " << get<1>(temp) << " " << get<2>(temp) << endl;
				triplets.push_back(temp);
				numPushBacks++;
			}
		}
	}

	TH1 *hbins1[xBins];
	for (int i = 0; i <= h1Bins; i++) {
		hbins1[i] = h1->ProjectionY(Form("bin%d",i),i,i); // name,firstxbin,lastxbin
		int size = hbins1[i]->GetNbinsX();
		int run = h1->GetXaxis()->GetBinCenter(i);
		for (int j = 0; j < size; j++){
			int w = hbins1[i]->GetBinContent(j);
			if (w > 0){
				double ene = hbins1[i]->GetBinCenter(j);
				if (ene < 0) ene = 1;
				test->Fill(run,ene,w);

				tuple<double,double,double> temp = make_tuple((double)run,(double)ene,(double)w);
				// cout << get<0>(temp) << " " << get<1>(temp) << " " << get<2>(temp) << endl;
				triplets.push_back(temp);
				numPushBacks++;
			}
		}
	}

	TCanvas *c = new TCanvas("c","Bob Ross's Canvas",800,600);
	c->SetLogz();
	test->Draw("COLZ");

	// figure out the number of x bins (last run - first run) and the
	tuple<double,double,double> firstEntry = triplets.front();
	tuple<double,double,double> lastEntry = triplets.back();

	xmin = get<0>(firstEntry);
	xmax = 1 + get<0>(lastEntry);
	xBins = 1 + xmax - xmin;

	printf("tuples: h0bins: %i  h1bins: %i pushbacks: %i  bins: %i  xmin: %i  xmax: %i\n",h0Bins,h1Bins,numPushBacks,xBins,xmin,xmax);
	// reg method: bins: 65  xmin: 9764  xmax: 9829
	// tuple method: bins: 63  xmin: 9764  xmax: 9827
	// tuple method doesn't push back when there's no data.
	// maybe it's OK.

	TH2D *balls = new TH2D("","",xBins,xmin,xmax,300,0,3000);

	// fill the new histogram with the values
	for ( vector < tuple<double,double,double> >::const_iterator it = triplets.begin(); it != triplets.end();  it++)
	{
		tuple<double,double,double> temp = *it;
		int run = get<0>(temp);
		int ene = get<1>(temp);
		int cts = get<2>(temp);
		balls->Fill(run,ene,cts);
		// cout << get<0>(temp) << " " << get<1>(temp) << " " << get<2>(temp) << endl;
	}

	TCanvas *b = new TCanvas("b","the tuple one",800,600);
	balls->Draw("COLZ");

}