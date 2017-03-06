// Low-energy data cleaning plots.
// Clint Wiseman, USC/Majorana
// 9/19/2016
//
// v1: used for LEWG report, 9/27/16.
// v2: used for DNP conference, 10/11/16.

#include <iostream>

#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include "TEntryList.h"
#include "TGraph.h"
#include "GATDataSet.hh"

using namespace std;

// binning:
// kris's binning: 15 bins/kev.  Ex: 750,0,50
// clint's binning: 2 bins/kev.  Ex: 100,0,50
// kris's second one: 5 bins/kev.  Ex: 250,0,50

// channel references:
//
// found in current version of skim files
// 610 p3d2  600 p7d1  692 p2d1  598 p7d2
// 626 p6d3  648 p2d2  582 p1d2  640 p2d3
// 578 p1d4  580 p1d3  690 p6d4  592 p7d4
// 672 p5d3  608 p3d3  664 p3d4  632 p6d1
// (no channel 594 p7d3) (no channel 616 p3d1)

// int highGain[18] = {640, 648, 664, 672, 690, 692,
// 					578, 580, 582, 592, 594, 598,
// 					600, 608, 610, 616, 626, 632};

// corresponds to the highGain array
// string hgPos[18] = {"P2D3","P2D2","P3D4","P5D3","P6D4","P2D1",
// 					"P1D4","P1D3","P1D2","P7D4","P7D3","P7D2",
// 					"P7D1","P3D3","P3D2","P3D1","P6D3","P6D1"};

// detector active masses (Kris's DS-0 unidoc, 9 Sept 2016)
// P1D2 0.9788  P1D3 0.8113  P1D4 0.9679  P2D1 0.587  P2D2 0.7229
// P2D3 0.6591  P3D2 0.8863  P3D3 0.949  P3D4 1.0238  P5D3 0.632
// P6D1 0.7316  P6D3 0.7014  P6D4 0.5723  P7D1 0.588  P7D2 0.7099  P7D4 0.9643
// *Uncertainty is 0.01 kg for each detector

// natural detectors / beges: 600 p7d1, 692 p2d1
int chans[16] = {610, 600, 692, 598, 626, 648, 582, 640, 578, 580, 690, 592, 672, 608, 664, 632};

char theCut[1000];

char basicCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0";

char burstCut[1000] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";

char standardCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

// only difference from the standard cut is that this one goes out to 100 keV
char DNPCut[1000] = "trapENFCal < 100 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

char tuneTECut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

char tunedToeEnr[1000] = "(channel==578 && kvorrT/trapENFCal<1.8) || (channel==580 && kvorrT/trapENFCal<1.4) || (channel==582 && kvorrT/trapENFCal<2.0) || (channel==592 && kvorrT/trapENFCal<2.0) || (channel==598 && kvorrT/trapENFCal<2.0) || (channel==608 && kvorrT/trapENFCal<2.0) || (channel==610 && kvorrT/trapENFCal<2.0) || (channel==626 && kvorrT/trapENFCal<2.0) || (channel==632 && kvorrT/trapENFCal<2.0) || (channel==640 && kvorrT/trapENFCal<2.0) || (channel==648 && kvorrT/trapENFCal<2.0) || (channel==664 && kvorrT/trapENFCal<1.75) || (channel==672 && kvorrT/trapENFCal<2.0) || (channel==690 && kvorrT/trapENFCal<2.0)";

char tunedToeNat[1000] = "(channel==600 && kvorrT/trapENFCal<2.0) || (channel==692 && kvorrT/trapENFCal<2.0)";

char noisyRunsCut[1000] = "run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

char noisyRunsOnly[1000] = "(run==9648 || run==10663 || run==10745 || run==11175 || run==12445 || run==12723 || run==12735 || run==12745 || run==12746 || run==12765 || run==12766 || run==12767 || run==13004)";

char ds0_toeCut[1000] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";

char ds1_toeCut[1000] = "(((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";


void basicCutSpec(TChain *skim)
{
	TH1D *h1 = new TH1D("h1","h1",250,0,50);
	skim->Project("h1","trapENFCal",basicCut);
	h1->SetLineColor(kBlack);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");

	TAxis *ax2 = h1->GetXaxis();
	printf("Basic cut.  Total cts 0.8-50: %.0f  cts 0.8-5: %.0f  cts 5-50: %.0f\n", h1->Integral(ax2->FindBin(0.8),ax2->FindBin(50)),
	h1->Integral(ax2->FindBin(0.8),ax2->FindBin(5)),
	h1->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TFile *f = new TFile("./data/basicCutSpec.root","RECREATE");
	h1->Write();
	f->Close();
}


void TETMSpecRunChannel(TChain *skim)
{
	// all standard cuts EXCEPT trapETailMin
	sprintf(theCut,"trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && %s",noisyRunsCut);

	// speed up draw time
	TBenchmark b;
	b.Start("entrylist");
	skim->Draw(">>entList",theCut,"entrylist GOFF");
	TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
	skim->SetEntryList(entList);
	b.Show("entrylist");

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);

	// what is the range of trapETailMin?
	b.Start("tetm-1d");
	c1->SetLogy();
	TH1D *h1 = new TH1D("h1","h1",100,-20,10);
	skim->Project("h1","trapETailMin",theCut);
	h1->GetXaxis()->SetTitle("TrapETailMin");
	h1->Draw();
	c1->Print("./output/TETM_withStdCuts.pdf");
	b.Show("tetm-1d");

	c1->SetLogy(0);
	c1->SetLogz(1);

	// TETM vs energy
	b.Start("tetm-energy");
	TH2D *h2 = new TH2D("h2","h2",750,0,50,100,-20,10);
	skim->Project("h2","trapETailMin:trapENFCal",theCut);
	c1->SetLogy(0);
	c1->SetLogz(1);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->GetYaxis()->SetTitle("trapETailMin");
	h2->Draw("COLZ");
	c1->Print("output/TETMvsEnergy.pdf");
	b.Show("tetm-energy");

	// TETM vs run
	b.Start("tetm-run");
	TH2D *h3 = new TH2D("h3","h3",4970,9420,14390,100,-20,10);
	skim->Project("h3","trapETailMin:run",theCut);
	h3->GetXaxis()->SetTitle("Run");
	h3->GetXaxis()->SetNdivisions(011);
	h3->GetYaxis()->SetTitle("trapETailMin");
	h3->Draw("COLZ");
	c1->Print("output/TETMvsRun.pdf");
	b.Show("tetm-run");

	// TETM vs. channel
	b.Start("tetm-channel");
	TH2D *h4 = new TH2D("h4","h4",125,575,700,100,-20,10);
	skim->Project("h4","trapETailMin:channel",theCut);
	h4->GetXaxis()->SetTitle("Channel");
	h4->GetYaxis()->SetTitle("trapETailMin");
	h4->Draw("COLZ");
	c1->Print("output/TETMvsChannel.pdf");
	b.Show("tetm-channel");

	// TETM vs first few runs
	b.Start("tetm-firstruns");
	TH2D *h5 = new TH2D("h5","h5",2000,9500,11500,100,-20,10); // don't forget to change it here too.
	// sprintf(theCut,"%s && run<9500",theCut);
	// sprintf(theCut,"%s && run>12400 && run<12622",theCut);
	sprintf(theCut,"%s && run>9500 && run<11500",theCut);
	skim->Project("h5","trapETailMin:run",theCut);
	h5->GetXaxis()->SetTitle("Run");
	// h5->GetXaxis()->SetNdivisions(011);
	h5->GetXaxis()->SetNdivisions(509);
	h5->GetYaxis()->SetTitle("trapETailMin");
	h5->Draw("COLZ");
	c1->Print("output/TETMvsRun_pt5Band.pdf");
	b.Show("tetm-firstruns");
}


void compareCuts1D(TChain *skim)
{
	double tot=0,pt8to5=0,from5to50=0;
	double totPrev=4560752, pt8to5Prev=4556741, from5to50Prev=4117; // Basic cut
	// double totPrev = 33903, pt8to5Prev=30429, from5to50Prev = 3548; // Basic+TETM
	// double totPrev = 19493, pt8to5Prev=16842, from5to50Prev = 2717; // Basic+TETM+NoisyRuns
	// double totPrev = 13700, pt8to5Prev=11781, from5to50Prev = 1971; // Basic+TETM+NoisyRuns+mH
	// double totPrev = 10599, pt8to5Prev=8694, from5to50Prev = 1957; // Basic+TETM+NoisyRuns+mH+LN
	// double totPrev = 10417, pt8to5Prev=8694, from5to50Prev = 1775; // Basic+TETM+NoisyRuns+mH+LN+isGood
	// double totPrev = 10345, pt8to5Prev=8623, from5to50Prev = 1774; // Basic+TETM+NoisyRuns+mH+LN+isGood+muVeto
	// double totPrev = 10317, pt8to5Prev=8604, from5to50Prev = 1764; // Basic+TETM+NoisyRuns+mH+LN+isGood+muVeto+burst

	// 1 basic cut (from original skim files - large file size)
	TFile *f = new TFile("./data/basicCutSpec.root");
	TH1D *h1 = (TH1D*)f->Get("h1");
	TAxis *ax1 = h1->GetXaxis();
	tot = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(50));
	pt8to5 = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5));
	from5to50 = h1->Integral(ax1->FindBin(5),ax1->FindBin(50));
	printf("Basic Cut\t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 2 trapETailMin
	sprintf(theCut,"%s && trapETailMin < 0",basicCut);

	TH1D *h2 = new TH1D("h2","h2",250,0,50);
	skim->Project("h2","trapENFCal",theCut);
	TAxis *ax2 = h2->GetXaxis();
	tot = h2->Integral(ax2->FindBin(0.8),ax2->FindBin(50));
	pt8to5 = h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5));
	from5to50 = h2->Integral(ax2->FindBin(5),ax2->FindBin(50));
	printf("TETM      \t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracTETM = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 3 low-e noisy runs cut
	sprintf(theCut,"%s && trapETailMin < 0 && %s",basicCut,noisyRunsCut);

	TH1D *h3 = new TH1D("h3","h3",250,0,50);
	skim->Project("h3","trapENFCal",theCut);
	TAxis *ax3 = h3->GetXaxis();
	tot = h3->Integral(ax3->FindBin(0.8),ax3->FindBin(50));
	pt8to5 = h3->Integral(ax3->FindBin(0.8),ax3->FindBin(5));
	from5to50 = h3->Integral(ax3->FindBin(5),ax3->FindBin(50));
	printf("Noisy Runs\t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracNoisyRuns = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 4 single-detector (mH==1)
	sprintf(theCut,"%s && trapETailMin < 0 && mH==1 && %s",basicCut,noisyRunsCut);

	TH1D *h4 = new TH1D("h4","h4",250,0,50);
	skim->Project("h4","trapENFCal",theCut);
	TAxis *ax4 = h4->GetXaxis();
	tot = h4->Integral(ax4->FindBin(0.8),ax4->FindBin(50));
	pt8to5 = h4->Integral(ax4->FindBin(0.8),ax4->FindBin(5));
	from5to50 = h4->Integral(ax4->FindBin(5),ax4->FindBin(50));
	printf("mH==1    \t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracMH = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 5 !isLNFill
	sprintf(theCut,"%s && trapETailMin < 0 && !isLNFill && %s && mH==1",basicCut,noisyRunsCut);

	TH1D *h5 = new TH1D("h5","h5",250,0,50);
	skim->Project("h5","trapENFCal",theCut);
	TAxis *ax5 = h5->GetXaxis();
	tot = h5->Integral(ax5->FindBin(0.8),ax5->FindBin(50));
	pt8to5 = h5->Integral(ax5->FindBin(0.8),ax5->FindBin(5));
	from5to50 = h5->Integral(ax5->FindBin(5),ax5->FindBin(50));
	printf("!isLNFill\t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracLN = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 6 isGood
	sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill",basicCut,noisyRunsCut);

	TH1D *h6 = new TH1D("h6","h6",250,0,50);
	skim->Project("h6","trapENFCal",theCut);
	TAxis *ax6 = h6->GetXaxis();
	tot = h6->Integral(ax6->FindBin(0.8),ax6->FindBin(50));
	pt8to5 = h6->Integral(ax6->FindBin(0.8),ax6->FindBin(5));
	from5to50 = h6->Integral(ax6->FindBin(5),ax6->FindBin(50));
	printf("isGood  \t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracGood = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 7 !muVeto
	sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill && !muVeto",basicCut,noisyRunsCut);

	TH1D *h7 = new TH1D("h7","h7",250,0,50);
	skim->Project("h7","trapENFCal",theCut);
	TAxis *ax7 = h7->GetXaxis();
	tot = h7->Integral(ax7->FindBin(0.8),ax7->FindBin(50));
	pt8to5 = h7->Integral(ax7->FindBin(0.8),ax7->FindBin(5));
	from5to50 = h7->Integral(ax7->FindBin(5),ax7->FindBin(50));
	printf("!muVeto  \t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracVeto = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 8 burst cut
	sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill && !muVeto && %s",basicCut,noisyRunsCut,burstCut);

	TH1D *h8 = new TH1D("h8","h8",250,0,50);
	skim->Project("h8","trapENFCal",theCut);
	TAxis *ax8 = h8->GetXaxis();
	tot = h8->Integral(ax8->FindBin(0.8),ax8->FindBin(50));
	pt8to5 = h8->Integral(ax8->FindBin(0.8),ax8->FindBin(5));
	from5to50 = h8->Integral(ax8->FindBin(5),ax8->FindBin(50));
	printf("!burst  \t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracBurst = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.

	// 9 !wfDCBits
	sprintf(theCut,"%s && trapETailMin < 0 && isGood && %s && mH==1 && !isLNFill && !muVeto && %s && !wfDCBits",basicCut,noisyRunsCut,burstCut);

	TH1D *h9 = new TH1D("h9","h9",250,0,50);
	skim->Project("h9","trapENFCal",theCut);
	TAxis *ax9 = h9->GetXaxis();
	tot = h9->Integral(ax9->FindBin(0.8),ax9->FindBin(50));
	pt8to5 = h9->Integral(ax9->FindBin(0.8),ax9->FindBin(5));
	from5to50 = h9->Integral(ax9->FindBin(5),ax9->FindBin(50));
	printf("!wfDCBits\t[0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracWFDC = 100*(tot/totPrev);
	totPrev=tot; pt8to5Prev=pt8to5; from5to50Prev=from5to50;	// only when looking at reduction from prev cut.


	// Now that we've established the proper order of the cuts,
	// draw them all together to see the reduction.

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();
	h1->SetLineColor(1);
	h1->Draw();
	h2->SetLineColor(2);
	h2->Draw("same");
	h3->SetLineColor(3);
	h3->Draw("same");
	h4->SetLineColor(4);
	h4->Draw("same");
	h5->SetLineColor(6);
	h5->Draw("same");
	h6->SetLineColor(7);
	h6->Draw("same");
	h7->SetLineColor(kOrange);
	h7->Draw("same");
	h8->SetLineColor(kMagenta+3);
	h8->Draw("same");
	h9->SetLineColor(kMagenta+3);
	h9->Draw("same");
	TLegend* leg1 = new TLegend(0.5,0.3,0.87,0.92);
	leg1->AddEntry(h1,"HG Only (% reduction)","l");
	leg1->AddEntry(h2,TString::Format("+TETM (%.1f%%)",100-fracTETM),"l");
	leg1->AddEntry(h3,TString::Format("+!NoisyRuns (%.1f%%)",100-fracNoisyRuns),"l");
	leg1->AddEntry(h4,TString::Format("+mH==1 (%.1f%%)",100-fracMH),"l");
	leg1->AddEntry(h5,TString::Format("+!isLNFill (%.1f%%)",100-fracLN),"l");
	leg1->AddEntry(h6,TString::Format("+isGood (%.1f%%)",100-fracGood),"l");
	leg1->AddEntry(h7,TString::Format("+!muVeto (%.1f%%)",100-fracVeto),"l");
	leg1->AddEntry(h8,TString::Format("+!burst (%.1f%%)",100-fracBurst),"l");
	leg1->AddEntry(h9,TString::Format("+!wfDCBits (%.1f%%)",100-fracWFDC),"l");

	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/cutReduction.pdf");
}


void rateVrun(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	// plot rate under 5 keV, make a list of runs w/ rates over 100Hz
	// that are NOT accounted for by the standard cut.
	TH1D *h1 = new TH1D("h1","h1",4970,9420,14390);
	sprintf(theCut,"%s && trapENFCal < 5",standardCut);
	skim->Project("h1","run",theCut);
	h1->GetXaxis()->SetTitle("Run");
	h1->GetYaxis()->SetTitle("Rate < 5keV (cts/run)");
	h1->GetXaxis()->SetNdivisions(509);
	h1->Draw();
	for (int i = 0; i < 4970; i++)
	{
		double run = h1->GetXaxis()->GetBinCenter(i) - 0.5;
		double rate = h1->GetBinContent(i);
		if (rate > 50)
			cout << run << "  " << rate << endl;
	}
	c1->Print("./output/rate_vs_run.pdf");
}


void findLiveTime()
{
	// run this on PDSF
	vector<int> noisyRuns = {9648,10663,10745,11175,12445,12723,12735,12745,12746,12765,12766,12767,13004};

	GATDataSet bs;
	double totRunTime = 0;
	for (auto i : noisyRuns)
	{
		GATDataSet ds(i);
		bs.AddRunNumber(i);
		totRunTime += ds.GetRunTime()/1e9;
		cout << i << "  " << ds.GetRunTime()/1e9 << endl;
	}
	cout << "bs: " << bs.GetRunTime()/1e9 << "  ds: " << totRunTime << endl;

	// now need to add a channel-by-channel plot like in Kris' unidoc ...?
}


void checkRunCut(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	// standard cut
	TH1D *h1 = new TH1D("h1","h1",150,0,10);
	sprintf(theCut,"%s",standardCut);
	skim->Project("h1","trapENFCal",standardCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->SetLineColor(kBlack);
	h1->Draw();
	TAxis *ax1 = h1->GetXaxis();
	double tot = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(10));

	// remove noisy runs
	TH1D *h2 = new TH1D("h2","h2",150,0,10);
	sprintf(theCut,"%s && %s",standardCut,noisyRunsCut);
	skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kBlue);
	h2->Draw("same");
	TAxis *ax2 = h2->GetXaxis();
	double fracLoNoise = 100*((h2->Integral(ax2->FindBin(0.8),ax2->FindBin(10)))/tot);

	// noisy runs only
	TH1D *h3 = new TH1D("h3","h3",150,0,10);
	sprintf(theCut,"%s && %s",standardCut,noisyRunsOnly);
	skim->Project("h3","trapENFCal",theCut);
	h3->SetLineColor(kRed);
	h3->Draw("same");
	TAxis *ax3 = h3->GetXaxis();
	double fracNoisy = 100*((h3->Integral(ax3->FindBin(0.8),ax3->FindBin(10)))/tot);

	TLegend* leg1 = new TLegend(0.3,0.65,0.87,0.92);
	leg1->AddEntry(h1,"Standard DC Cut","l");
	leg1->AddEntry(h2,TString::Format("With Run Cut (%.1f%% persist)",100-fracLoNoise),"l");
	leg1->AddEntry(h3,TString::Format("Noisy Runs Only (%.1f%% removed)",100-fracNoisy),"l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("./output/noisyRunsCut.pdf");

	// channels without noisy runs
	c1->SetLogy(0);
	c1->SetLogz(1);
	TH2D *h4 = new TH2D("h4","h4",150,0,10,125,575,700);
	sprintf(theCut,"%s && %s && trapENFCal < 10",standardCut,noisyRunsCut);
	skim->Project("h4","channel:trapENFCal",theCut);
	h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h4->GetYaxis()->SetTitle("Channel");
	h4->Draw("COLZ");
	c1->Print("output/noisyRuns_EnergyVsChannel.pdf");

	// channels w/ noisy runs only.
	TH2D *h5 = new TH2D("h5","h5",150,0,10,125,575,700);
	sprintf(theCut,"%s && %s && trapENFCal < 10",standardCut,noisyRunsOnly);
	skim->Project("h5","channel:trapENFCal",theCut);
	h5->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h5->GetYaxis()->SetTitle("Channel");
	h5->Draw("COLZ");
	c1->Print("output/noisyRunsOnly_EnergyVsChannel.pdf");
}


void tuneTheTECut(TChain *skim)
{
	char hname[100];

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",2600,800);
	c1->Divide(3,1,0,0);
	TH1D *htoe[16];
	TH1D *htoe2[16];
	TH2D *htoe2D[16];
	for (int i=0; i<16; i++)
	{
		c1->cd(1);
		TVirtualPad *p1 = c1->GetPad(1);
		p1->SetTopMargin(0.05);
		p1->SetLogx(0);

		sprintf(hname,"h_%i",chans[i]);
		htoe[i] = new TH1D(hname,hname,200,0,10);

		sprintf(theCut,"%s && channel==%i",tuneTECut,chans[i]);
		skim->Project(hname,"kvorrT/trapENFCal",theCut);

		htoe[i]->GetXaxis()->SetTitle(TString::Format("kvorrT/trapENFCal (%i)",chans[i]));
		htoe[i]->GetXaxis()->SetTitleOffset(1.2);
		htoe[i]->Draw();

		c1->cd(2);

		TVirtualPad *p2 = c1->GetPad(2);
		p2->SetTopMargin(0.05);
		p2->SetLeftMargin(0.1);
		p2->SetRightMargin(0.1);
		p2->SetLogx(1);
		p2->SetLogz(1);

		sprintf(hname,"h_2d_%i",chans[i]);
		htoe2D[i] = new TH2D(hname,hname,1000,0,200,200,0,10);
		skim->Project(hname,"kvorrT/trapENFCal:trapENFCal",theCut);
		htoe2D[i]->GetXaxis()->SetTitle(TString::Format("trapENFCal (%i)",chans[i]));
		htoe2D[i]->GetXaxis()->SetTitleOffset(1.2);
		htoe2D[i]->GetYaxis()->SetTitle("T/E");
		htoe2D[i]->SetMinimum(1);
		if (htoe2D[i]->GetMaximum() < 100) htoe2D[i]->SetMaximum(100);
		htoe2D[i]->Draw("COLZ");

		c1->cd(3);
		TVirtualPad *p3 = c1->GetPad(3);
		p3->SetTopMargin(0.05);
		p3->SetLeftMargin(0.05);
		p3->SetRightMargin(0.05);
		p3->SetLogx(0);

		sprintf(hname,"h_%i",chans[i]);
		htoe2[i] = new TH1D(hname,hname,500,0.,2.5);

		sprintf(theCut,"%s && channel==%i",tuneTECut,chans[i]);
		skim->Project(hname,"kvorrT/trapENFCal",theCut);
		htoe2[i]->GetXaxis()->SetTitle(TString::Format("kvorrT/trapENFCal (%i)",chans[i]));
		htoe2[i]->GetXaxis()->SetTitleOffset(1.2);
		double ymax = htoe2[i]->GetMaximum();
		double binmax = htoe2[i]->GetMaximumBin();
		double xmax = htoe2[i]->GetBinCenter(binmax);
		htoe2[i]->GetXaxis()->SetRange(binmax-80,binmax+80);
		htoe2[i]->Draw();

		// htoe2[i]->Fit("gaus");
		cout << chans[i] << "  " << ymax << "  " << binmax << "  " << xmax << endl;

		sprintf(hname,"./output/toe_%i.pdf",chans[i]);
		c1->Print(hname);
	}
}


void checkTECut(TChain *skim)
{
	double tot=0, pt8to5=0, from5to50=0;
	double totPrev=10299, pt8to5Prev=8599, from5to50Prev=1751; // Standard Cut (no T/E)

	sprintf(theCut,"%s",standardCut);
	TH1D *h1 = new TH1D("h1","h1",250,0,50);
	skim->Project("h1","trapENFCal",theCut);
	TAxis *ax1 = h1->GetXaxis();
	tot = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(50));
	pt8to5 = h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5));
	from5to50 = h1->Integral(ax1->FindBin(5),ax1->FindBin(50));
	printf("Standard  [0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));

	// sprintf(theCut,"%s && %s",standardCut,ds0_toeCut);
	sprintf(theCut,"%s && %s",standardCut,ds1_toeCut);
	TH1D *h2 = new TH1D("h2","h2",250,0,50);
	skim->Project("h2","trapENFCal",theCut);
	TAxis *ax2 = h2->GetXaxis();
	tot = h2->Integral(ax2->FindBin(0.8),ax2->FindBin(50));
	pt8to5 = h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5));
	from5to50 = h2->Integral(ax2->FindBin(5),ax2->FindBin(50));
	printf("+T/E      [0.8-50]: %.0f (%.2f%%) [0.8-5]: %.0f (%.2f%%) [5-50]: %.0f (%.2f%%)\n", tot,100*(tot/totPrev),pt8to5,100*(pt8to5/pt8to5Prev),from5to50,100*(from5to50/from5to50Prev));
	double fracTE = 100*(tot/totPrev);


	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();
	h1->SetLineColor(kBlue);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->Draw();
	h2->SetLineColor(kRed);
	h2->Draw("same");
	TLegend* leg1 = new TLegend(0.4,0.7,0.87,0.92);
	leg1->AddEntry(h1,"Standard Cut","l");
	leg1->AddEntry(h2,TString::Format("+T/E (%.1f%% reduction)",100-fracTE),"l");

	leg1->Draw("SAME");
	c1->Update();
	// c1->Print("output/TECutReduction_DS0.pdf");
	c1->Print("output/TECutReduction_DS1.pdf");
}

void fullTESpec(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogz();
	c1->SetLogx();

	TH2D *h1 = new TH2D("h1","h1",1000,0,200,800,0,10);

	sprintf(theCut,"%s",tuneTECut);
	skim->Project("h1","kvorrT/trapENFCal:trapENFCal",theCut);

	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("T / E : enriched");
	h1->GetXaxis()->SetTitleOffset(1.2);
	h1->SetMinimum(1);
	// h1->SetMaximum(100);
	h1->Draw("COLZ");

	c1->Update();
	c1->Print("output/enrichedT-E.pdf");

}

// Clint, this is code you went with for DNP.  Best not to change it!
void enrVsNat(TChain* skim)
{
	double binsPerKeV = 5;

	double ds1livetime = 60.054;  // full DS-1 with noisy runs, LN, and veto times cut
	double natMass=0, enrMass=0;
	skim->SetBranchAddress("mAct_nat_kg",&natMass);
	skim->SetBranchAddress("mAct_enr_kg",&enrMass);
	skim->GetEntry(0);
	cout << "Active masses (kg).  Natural " << natMass << "  enr: " << enrMass << endl;

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	// c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",(int)binsPerKeV*25,5,30);
	sprintf(theCut,"%s && %s && !isEnr",standardCut,ds1_toeCut);
	skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (keV)");
	h1->GetYaxis()->SetTitle("Counts");
	h1->SetLineColor(kBlue);
	h1->Draw();

	TH1D *h2 = new TH1D("h2","h2",(int)binsPerKeV*25,5,30);
	sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
	skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");

	TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
	leg1->AddEntry(h1,TString::Format("Natural: %.1f kg-d",natMass*ds1livetime),"l");
	leg1->AddEntry(h2,TString::Format("Enriched: %.1f kg-d",enrMass*ds1livetime),"l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("./output/enrVsNat_DS1.pdf");

	TFile *f = new TFile("./output/enrVsNat_DS1.root","RECREATE");
	h1->Write("nat");
	h2->Write("enr");
	f->Close();
}

// Clint, this is code you went with for DNP.  Best not to change it!
void EnrVsNat_DS0_DS1(TChain *skim)
{
	double ds0enrExp = 478;  // taken from RS&DC DS-0 document
	double ds0natExp = 195;  //
	double ds1enrExp = 679;  // calculated in the above code, using LN fill, veto, and noisy run cuts.
	double ds1natExp = 71;   //
	double binsPerKeV = 5;	 // must match "enrVsNat" function above.

	// DS-0 spectra from Kris

	TFile *f1 = new TFile("./data/m1DiagnosticTree4Aug2016.root");
	TTree *t1 = (TTree*)f1->Get("diagTree");
	cout << t1->GetEntries() << endl;

	TH1D *ds0enr = new TH1D("ds0enr","ds0enr",(int)binsPerKeV*25,5,30);
	t1->Project("ds0enr","calENF","enr");

	TH1D *ds0nat = new TH1D("ds0nat","ds0nat",(int)binsPerKeV*25,5,30);
	t1->Project("ds0nat","calENF","!enr");

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	ds0nat->GetXaxis()->SetTitle("Energy (keV)");
	ds0nat->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
	// ds0nat->GetYaxis()->SetTitleOffset(1.1);
	ds0nat->SetLineColor(kBlue);
	ds0nat->Scale(binsPerKeV/ds0natExp);
	ds0nat->SetMinimum(0.01);
	ds0nat->SetMaximum(6);
	ds0nat->Draw();

	ds0enr->SetLineColor(kRed);
	ds0enr->Scale(binsPerKeV/ds0enrExp);
	ds0enr->Draw("same");

	TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
	leg1->AddEntry(ds0nat,TString::Format("Natural: %.0f kg-d",ds0natExp),"l");
	leg1->AddEntry(ds0enr,TString::Format("Enriched: %.0f kg-d",ds0enrExp),"l");
	leg1->Draw("SAME");
	c1->Update();

	c1->Print("./output/kris-ds0-enrvsnat.pdf");

	// DS-1 spectra from above code

	TH1D *ds1enr = new TH1D("ds1enr","ds1enr",(int)binsPerKeV*25,5,30);
	sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
	skim->Project("ds1enr","trapENFCal",theCut);

	TH1D *ds1nat = new TH1D("ds1nat","ds1nat",(int)binsPerKeV*25,5,30);
	sprintf(theCut,"%s && %s && !isEnr",standardCut,ds1_toeCut);
	skim->Project("ds1nat","trapENFCal",theCut);

	ds1nat->GetXaxis()->SetTitle("Energy (keV)");
	ds1nat->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
	// ds1nat->GetYaxis()->SetTitleOffset(1.1);
	ds1nat->SetLineColor(kBlue);
	ds1nat->Scale(binsPerKeV/ds1natExp);
	ds1nat->SetMinimum(0.01);
	ds1nat->SetMaximum(6);
	ds1nat->Draw();

	ds1enr->SetLineColor(kRed);
	ds1enr->Scale(binsPerKeV/ds1enrExp);
	ds1enr->Draw("same");

	TLegend* leg2 = new TLegend(0.45,0.75,0.87,0.92);
	leg2->AddEntry(ds1nat,TString::Format("Natural: %.0f kg-d",ds1natExp),"l");
	leg2->AddEntry(ds1enr,TString::Format("Enriched: %.0f kg-d",ds1enrExp),"l");
	leg2->Draw("SAME");
	c1->Update();

	c1->Print("./output/clint-ds1-enrvsnat.pdf");


	// Brandon asked for a natural detector comparison.

	c1->SetLogy(0);
	ds0nat->Draw();
	ds1nat->SetLineColor(kRed);
	ds1nat->Draw("same");
	c1->Print("./output/clint-ds0-ds1-natural.pdf");
}

// Clint, this is code you went with for DNP.  Best not to change it!
void EnrSpec_DS0_DS1(TChain *skim)
{
	double ds0enrExp = 478;  // taken from RS&DC DS-0 document
	double ds1enrExp = 679;  // calculated in the above code, using LN fill, veto, and noisy run cuts.
	double binsPerKeV = 2;	   // matches "enrVsNat" function above.

	// DS-0 spectra from Kris

	TFile *f1 = new TFile("./data/m1DiagnosticTree4Aug2016.root");
	TTree *t1 = (TTree*)f1->Get("diagTree");
	cout << t1->GetEntries() << endl;

	TH1D *ds0enr = new TH1D("ds0enr","ds0enr",(int)binsPerKeV*25,5,30);
	t1->Project("ds0enr","calENF","enr");

	// compare DS-0 and DS-1 enriched spectra

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);

	ds0enr->GetXaxis()->SetTitle("Energy (keV)");
	ds0enr->GetYaxis()->SetTitle("Counts kg^{-1} day^{-1} keV^{-1}");
	ds0enr->GetYaxis()->SetTitleOffset(1.1);
	ds0enr->SetLineColor(kBlue);
	ds0enr->Scale(binsPerKeV/ds0enrExp);

	ds0enr->SetMinimum(0);
	ds0enr->Draw();

	TH1D *ds1enr = new TH1D("ds1enr","ds1enr",(int)binsPerKeV*25,5,30);
	sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
	skim->Project("ds1enr","trapENFCal",theCut);

	ds1enr->SetLineColor(kRed);
	ds1enr->Scale(binsPerKeV/ds1enrExp);
	ds1enr->Draw("same");

	TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
	leg1->AddEntry(ds0enr,TString::Format("Enriched, DS-0: %.0f kg-d",ds0enrExp),"l");
	leg1->AddEntry(ds1enr,TString::Format("Enriched, DS-1: %.0f kg-d",ds1enrExp),"l");
	leg1->Draw("SAME");
	c1->Update();

	c1->Print("./output/kris-clint-enriched.pdf");

	TAxis *ax2 = ds0enr->GetXaxis();
	double ds0from15to30 = ds0enr->Integral(ax2->FindBin(15),ax2->FindBin(29),"width");

	printf("DS-0: From 15-30 keV: %.3f cts/kg-d,  %.3f cts/kg-d-keV\n",ds0from15to30, ds0from15to30/15.);

	TAxis *ax1 = ds1enr->GetXaxis();
	double ds1from15to30 = ds1enr->Integral(ax1->FindBin(15),ax1->FindBin(29),"width");

	printf("DS-1: From 15-30 keV: %.3f cts/kg-d,  %.3f cts/kg-d-keV\n",ds1from15to30, ds1from15to30/15.);

	cout << "Reduction Factor: " << ds0from15to30/ds1from15to30 << endl;

	TFile *f = new TFile("./output/ds0-ds1-enriched.root","RECREATE");
	ds0enr->Write("ds0enr");
	ds1enr->Write("ds1enr");
	f->Close();
}

// Clint, this is code you went with for DNP.  Best not to change it!
void IntegrateCounts()
{
	cout << "\nChecking integral of spectra ...\n";

	TFile *f = new TFile("./output/ds0-ds1-enriched.root");
	TH1D *ds0enr = (TH1D*)f->Get("ds0enr");
	TH1D *ds1enr = (TH1D*)f->Get("ds1enr");

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);

	ds0enr->GetXaxis()->SetRangeUser(5,30);
	ds0enr->SetMinimum(0);
	ds0enr->Draw();
	ds1enr->Draw("same");
	c1->Print("./output/check-ds0-ds1.pdf");

	// Clint, this computes the TOTAL counts between 15-30 keV,
	// so the result is in cts/kg-d
	// Then to get the average, divide by the number of bins

	TAxis *ax1 = ds0enr->GetXaxis();
	TAxis *ax2 = ds1enr->GetXaxis();
	double ds0tot = ds0enr->Integral(ax1->FindBin(15),ax1->FindBin(29),"width"); // don't include 30, it's overflow
	double ds1tot = ds1enr->Integral(ax2->FindBin(15),ax2->FindBin(29),"width");
	cout << "cts/kg-d                : ds0: " << ds0tot << "  ds1: " << ds1tot << "  reduction: " << ds0tot/ds1tot << endl;
	cout << "cts/kg-d-keV            : ds0: " << ds0tot/15. << "  ds1: " << ds1tot/15. << "  reduction: " << ds0tot/ds1tot << endl;

	// Clint, this computes the AVERAGE value between 15-30 keV
	// So the result is in cts/kg-d-keV

	double ds0avg=0, ds1avg=0;	// 15-30 keV
	for (int i = 20; i <= 49; i++)
	{
		// cout << i << "  " << (double)i/2 + 5 << "  \t" << ds0enr->GetBinContent(i) << "  " << ds1enr->GetBinContent(i) << endl;
		ds0avg += ds0enr->GetBinContent(i);
		ds1avg += ds1enr->GetBinContent(i);
	}
	ds0avg = ds0avg/30;
	ds1avg = ds1avg/30;
	cout << "cts/kg-d-keV, avg method: ds0: " << ds0avg << "  ds1: " << ds1avg << "  reduction: " << ds0avg/ds1avg << endl;

}


// Clint, this is code you went with for DNP.  Best not to change it!
void EnrSpec_DS0_DS1_Above15(TChain *skim)
{
	double lo = 20;
	double hi = 100;

	double ds0enrExp = 478;  // taken from RS&DC DS-0 document
	double ds1enrExp = 679;  // calculated in the above code, using LN fill, veto, and noisy run cuts.
	double binsPerKeV = 2;	   // matches "enrVsNat" function above.

	// DS-0 spectra from Kris

	TFile *f1 = new TFile("./data/m1DiagnosticTree4Aug2016.root");
	TTree *t1 = (TTree*)f1->Get("diagTree");
	cout << t1->GetEntries() << endl;

	TH1D *ds0enr = new TH1D("ds0enr","ds0enr",(int)binsPerKeV*(hi-lo),lo,hi);
	t1->Project("ds0enr","calENF","enr");

	// compare DS-0 and DS-1 enriched spectra

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	// c1->SetLogy();

	ds0enr->GetXaxis()->SetTitle("Energy (keV)");
	ds0enr->GetYaxis()->SetTitle("counts / kg / day / keV");
	ds0enr->GetYaxis()->SetTitleOffset(1.1);
	ds0enr->SetLineColor(kBlue);
	ds0enr->Scale(binsPerKeV/ds0enrExp);

	ds0enr->SetMinimum(0);
	ds0enr->SetMaximum(0.1);
	ds0enr->Draw();

	TH1D *ds1enr = new TH1D("ds1enr","ds1enr",(int)binsPerKeV*(hi-lo),lo,hi);
	sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
	skim->Project("ds1enr","trapENFCal",theCut);

	ds1enr->SetLineColor(kRed);
	ds1enr->Scale(binsPerKeV/ds1enrExp);
	ds1enr->Draw("same");

	TLegend* leg1 = new TLegend(0.45,0.75,0.87,0.92);
	leg1->AddEntry(ds0enr,TString::Format("Enriched, DS-0: %.0f kg-d",ds0enrExp),"l");
	leg1->AddEntry(ds1enr,TString::Format("Enriched, DS-1: %.0f kg-d",ds1enrExp),"l");
	leg1->Draw("SAME");
	c1->Update();

	char plotname[200];
	sprintf(plotname,"./output/kris-clint-enriched-lo%.0f-hi%.0f.pdf",lo,hi);
	c1->Print(plotname);

	TAxis *ax2 = ds0enr->GetXaxis();
	double ds0int = ds0enr->Integral(ax2->FindBin(lo),ax2->FindBin(hi-1),"width");

	printf("DS-0: From %.0f-%.0f keV: %.3f cts/kg-d,  %.3f cts/kg-d-keV\n",lo,hi,ds0int,ds0int/(hi-lo));

	TAxis *ax1 = ds1enr->GetXaxis();
	double ds1int = ds1enr->Integral(ax1->FindBin(lo),ax1->FindBin(hi-1),"width");

	printf("DS-1: From %.0f-%.0f keV: %.3f cts/kg-d,  %.3f cts/kg-d-keV\n",lo,hi,ds1int, ds1int/(hi-lo));

	cout << "Reduction Factor: " << ds0int/ds1int << endl;

	sprintf(plotname,"./output/ds0-ds1-enriched-lo%.0f-hi%.0f.root",lo,hi);
	TFile *f = new TFile(plotname,"RECREATE");
	ds0enr->Write("ds0enr");
	ds1enr->Write("ds1enr");
	f->Close();
}

// Clint, this is code you went with for DNP.  Best not to change it!
void IntegrateCounts_Above15()
{
	double lo = 20;
	double hi = 100;

	cout << "\nChecking integral of spectra ...\n";

	char plotname[200];
	sprintf(plotname,"./output/ds0-ds1-enriched-lo%.0f-hi%.0f.root",lo,hi);
	TFile *f = new TFile(plotname);
	TH1D *ds0enr = (TH1D*)f->Get("ds0enr");
	TH1D *ds1enr = (TH1D*)f->Get("ds1enr");

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);

	ds0enr->GetXaxis()->SetRangeUser(lo,hi);
	ds0enr->SetMinimum(0);
	ds0enr->Draw();
	ds1enr->Draw("same");
	c1->Print("./output/check-ds0-ds1-above15.pdf");

	// Clint, this computes the TOTAL counts between "lo" and "hi" keV,
	// so the result is in cts/kg-d
	// Then to get the average, divide by the number of bins

	TAxis *ax1 = ds0enr->GetXaxis();
	TAxis *ax2 = ds1enr->GetXaxis();
	double ds0tot = ds0enr->Integral(ax1->FindBin(lo),ax1->FindBin(hi-1),"width");
	double ds1tot = ds1enr->Integral(ax2->FindBin(lo),ax2->FindBin(hi-1),"width");
	cout << "cts/kg-d                : ds0: " << ds0tot << "  ds1: " << ds1tot << "  reduction: " << ds0tot/ds1tot << endl;
	cout << "cts/kg-d-keV            : ds0: " << ds0tot/(hi-lo) << "  ds1: " << ds1tot/(hi-lo) << "  reduction: " << ds0tot/ds1tot << endl;

	// Clint, this computes the AVERAGE value between 15-50 keV
	// So the result is in cts/kg-d-keV

	// double ds0avg=0, ds1avg=0;	// 15-50 keV
	// for (int i = 20; i <= 49; i++)
	// {
	// 	// cout << i << "  " << (double)i/2 + 5 << "  \t" << ds0enr->GetBinContent(i) << "  " << ds1enr->GetBinContent(i) << endl;
	// 	ds0avg += ds0enr->GetBinContent(i);
	// 	ds1avg += ds1enr->GetBinContent(i);
	// }
	// ds0avg = ds0avg/30;
	// ds1avg = ds1avg/30;
	// cout << "cts/kg-d-keV, avg method: ds0: " << ds0avg << "  ds1: " << ds1avg << "  reduction: " << ds0avg/ds1avg << endl;

}


// int main(int argc, char** argv)
int main()
{
	gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
	gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");

	// TChain *skim = new TChain("skimTree");
	// skim->Add("~/ds1-low/skimDS1_pt8*");
	// basicCutSpec(skim);	// done for new skim file
	// TETMSpecRunChannel(skim);  // done for new skim file
	// cout << "\007";	// beep when done

	// made using tcut-skimmer - basic cut + trapETailMinCut
	// TChain *cutSkim = new TChain("skimTree");
	// cutSkim->Add("./data/developCut_DS1.root");
	// printf("\n\nFound %lli entries.\n",cutSkim->GetEntries());
	// compareCuts1D(cutSkim);

	// made using tcut-skimmer - standard cut
	// TChain *cutSkim = new TChain("skimTree");
	// cutSkim->Add("./data/standardCut_DS1.root");
	// printf("\nFound %lli entries.\n",cutSkim->GetEntries());
	// rateVrun(cutSkim);
	// findLiveTime();
	// checkRunCut(cutSkim);
	// checkTECut(cutSkim);
	// enrVsNat(cutSkim);
	// EnrVsNat_DS0_DS1(cutSkim);
	// EnrSpec_DS0_DS1(cutSkim);
	// IntegrateCounts();

	TChain *cutSkim = new TChain("skimTree");
	cutSkim->Add("./data/DNPCut_DS1.root");
	printf("\nFound %lli entries.\n",cutSkim->GetEntries());
	EnrSpec_DS0_DS1_Above15(cutSkim);
	IntegrateCounts_Above15();

	// made using tcut-skimmer - standard cut but with a 200 keV limit
	// TChain *cutSkim = new TChain("skimTree");
	// cutSkim->Add("./data/tuneTECut_DS1.root");
	// printf("\n\nFound %lli entries.\n",cutSkim->GetEntries());
	// tuneTheTECut(cutSkim);
	// fullTESpec(cutSkim);

}
