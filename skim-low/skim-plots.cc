// Low-energy data cleaning plots.
// Clint Wiseman, USC/Majorana
// 9/19/2016

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

// natural detectors / beges: 600 p7d1, 692 p2d1
int chans[16] = {610, 600, 692, 598, 626, 648, 582, 640, 578, 580, 690, 592, 672, 608, 664, 632};

char theCut[1000];

char basicCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0";

char burstCut[1000] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";

char standardCut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && trapETailMin<0";

char tunedToeEnr[1000] = "(channel==578 && kvorrT/trapENFCal<1.8) || (channel==580 && kvorrT/trapENFCal<1.4) || (channel==582 && kvorrT/trapENFCal<2.0) || (channel==592 && kvorrT/trapENFCal<2.0) || (channel==598 && kvorrT/trapENFCal<2.0) || (channel==608 && kvorrT/trapENFCal<2.0) || (channel==610 && kvorrT/trapENFCal<2.0) || (channel==626 && kvorrT/trapENFCal<2.0) || (channel==632 && kvorrT/trapENFCal<2.0) || (channel==640 && kvorrT/trapENFCal<2.0) || (channel==648 && kvorrT/trapENFCal<2.0) || (channel==664 && kvorrT/trapENFCal<1.75) || (channel==672 && kvorrT/trapENFCal<2.0) || (channel==690 && kvorrT/trapENFCal<2.0)";

char tunedToeNat[1000] = "(channel==600 && kvorrT/trapENFCal<2.0) || (channel==692 && kvorrT/trapENFCal<2.0)";

void compareCuts1D(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("rawENF","rawENF",750,0,50);
	int cts = skim->Project("rawENF","trapENFCal",basicCut);
	h1->SetLineColor(kBlack);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->Draw();
	cout << "raw cts: " << cts << endl;

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s && isGood",basicCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kBlue);
	h2->Draw("SAME");
	cout << "isGood: " << cts << endl;

	TH1D *h3 = new TH1D("h3","h3",750,0,50);
	sprintf(theCut,"%s && mH==1",basicCut);
	cts = skim->Project("h3","trapENFCal",theCut);
	h3->SetLineColor(kViolet);
	h3->Draw("SAME");
	cout << "singlesite: " << cts << endl;

	TH1D *h4 = new TH1D("h4","h4",750,0,50);
	sprintf(theCut,"%s && !isLNFill",basicCut);
	cts = skim->Project("h4","trapENFCal",theCut);
	h4->SetLineColor(kOrange);
	h4->Draw("SAME");
	cout << "not lnFill: " << cts << endl;

	TH1D *h5 = new TH1D("h5","h5",750,0,50);
	sprintf(theCut,"%s && !muVeto",basicCut);
	cts = skim->Project("h5","trapENFCal",theCut);
	h5->SetLineColor(kGreen+2);
	h5->Draw("SAME");
	cout << "not muon veto: " << cts << endl;

	TH1D *h6 = new TH1D("h6","h6",750,0,50);
	sprintf(theCut,"%s && %s",basicCut,burstCut);
	cts = skim->Project("h6","trapENFCal",theCut);
	h6->SetLineColor(kCyan+1);
	h6->Draw("SAME");
	cout << "burst cut: " << cts << endl;

	TH1D *h7 = new TH1D("h7","h7",750,0,50);
	sprintf(theCut,"%s && !wfDCBits",basicCut);
	cts = skim->Project("h7","trapENFCal",theCut);
	h7->SetLineColor(kOrange+3);
	h7->Draw("SAME");
	cout << "dc cut: " << cts << endl;

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.94);
	leg1->AddEntry(h1,"basic","l");
	leg1->AddEntry(h2,"isGood","l");
	leg1->AddEntry(h3,"mH==1","l");
	leg1->AddEntry(h4,"!lnFill","l");
	leg1->AddEntry(h5,"!muVeto","l");
	leg1->AddEntry(h6,"!burst","l");
	leg1->AddEntry(h7,"!wfDCBits","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/makeStdCut.pdf");
}


void combineCuts1D(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",750,0,50);
	sprintf(theCut,"%s && mH==1",basicCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->SetLineColor(kBlack);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->Draw();
	cout << "raw cts with mH==1: " << cts << endl;

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s && mH==1 && isGood",basicCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kBlue);
	h2->Draw("SAME");
	cout << "isGood: " << cts << endl;

	TH1D *h3 = new TH1D("h3","h3",750,0,50);
	sprintf(theCut,"%s && mH==1 && !isLNFill",basicCut);
	cts = skim->Project("h3","trapENFCal",theCut);
	h3->SetLineColor(kOrange);
	h3->Draw("SAME");
	cout << "not lnFill: " << cts << endl;

	TH1D *h4 = new TH1D("h4","h4",750,0,50);
	sprintf(theCut,"%s && mH==1 && !muVeto",basicCut);
	cts = skim->Project("h4","trapENFCal",theCut);
	h4->SetLineColor(kGreen+2);
	h4->Draw("SAME");
	cout << "not muon veto: " << cts << endl;

	TH1D *h5 = new TH1D("h5","h5",750,0,50);
	sprintf(theCut,"%s && mH==1 && %s",basicCut,burstCut);
	cts = skim->Project("h5","trapENFCal",theCut);
	h5->SetLineColor(kCyan+1);
	h5->Draw("SAME");
	cout << "burst cut: " << cts << endl;

	TH1D *h6 = new TH1D("h6","h6",750,0,50);
	sprintf(theCut,"%s && mH==1 && !wfDCBits",basicCut);
	cts = skim->Project("h6","trapENFCal",theCut);
	h6->SetLineColor(kOrange+3);
	h6->Draw("SAME");
	cout << "dc cut: " << cts << endl;

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.94);
	leg1->AddEntry(h1,"basic","l");
	leg1->AddEntry(h2,"isGood","l");
	leg1->AddEntry(h3,"!lnFill","l");
	leg1->AddEntry(h4,"!muVeto","l");
	leg1->AddEntry(h5,"!burst","l");
	leg1->AddEntry(h6,"!wfDCBits","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/mHCutPlusOthers.pdf");
}


void checkLNFill(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",150,0,10);
	sprintf(theCut,"%s && mH==1",basicCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->SetLineColor(kBlack);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->Draw();
	cout << "raw cts with mH==1: " << cts << endl;

	TH1D *h2 = new TH1D("h2","h2",150,0,10);
	sprintf(theCut,"%s && mH==1 && !isLNFill",basicCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");
	cout << "add LN fill: " << cts << endl;

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.94);
	leg1->AddEntry(h1,"mH==1","l");
	leg1->AddEntry(h2,"mH==1 && !isLNFill","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/checkLNFill.pdf");

	TH2D *h3 = new TH2D("runENF","runENF",750,0,50,125,575,700);
	c1->SetLogy(0);
	c1->SetLogz(1);
	sprintf(theCut,"%s && !isLNFill",standardCut);
	skim->Project("runENF","channel:trapENFCal",theCut);
	h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h3->GetYaxis()->SetTitle("channel");
	h3->Draw("COLZ");
	c1->Print("output/lnFillByChannel.pdf");

	TH2D *h4 = new TH2D("runENF","runENF",750,0,50,4970,9420,14390);
	c1->SetLogz(1);
	skim->Project("runENF","run:trapENFCal",theCut);
	h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h4->GetYaxis()->SetTitle("run");
	h4->Draw("COLZ");
	c1->Print("output/lnFillByRun.pdf");
}


void checkBurstCut(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",150,0,10);
	sprintf(theCut,"%s",standardCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->SetLineColor(kBlack);
	h1->Draw();
	cout << "standard cut: " << cts << endl;

	TH1D *h2 = new TH1D("h2","h2",150,0,10);
	sprintf(theCut,"%s && %s",standardCut,burstCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");
	cout << "add burst cut: " << cts << endl;

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.94);
	leg1->AddEntry(h1,"standard","l");
	leg1->AddEntry(h2,"!burst","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/checkBurstCut.pdf");

	TH2D *h3 = new TH2D("runENF","runENF",750,0,50,125,575,700);
	c1->SetLogy(0);
	c1->SetLogz(1);
	sprintf(theCut,"%s && %s",standardCut,burstCut);
	skim->Project("runENF","channel:trapENFCal",theCut);
	h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h3->GetYaxis()->SetTitle("channel");
	h3->Draw("COLZ");
	c1->Print("output/burstByChannel.pdf");

	TH2D *h4 = new TH2D("runENF","runENF",750,0,50,4970,9420,14390);
	c1->SetLogz(1);
	skim->Project("runENF","run:trapENFCal",theCut);
	h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h4->GetYaxis()->SetTitle("run");
	h4->Draw("COLZ");
	c1->Print("output/burstByRun.pdf");
}


void checkTETMCut(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",750,0,50);
	sprintf(theCut,"%s",basicCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->SetLineColor(kBlack);
	h1->Draw();

	TAxis *ax1 = h1->GetXaxis();
	printf("basic cut (no TETM).  Total cts %i  cts 0.8-5 %.0f  cts 5-50 %.0f\n", cts, h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h1->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s && trapETailMin<0",basicCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");
	cout << "basic + TETM " << cts << endl;

	TAxis *ax2 = h2->GetXaxis();
	printf("basic cut w/ TETM.  Total cts %i  cts 0.8-5 %.0f  cts 5-50 %.0f\n", cts, h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5)), h2->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.92);
	leg1->AddEntry(h1,"basic,no TETM","l");
	leg1->AddEntry(h2,"standard,w/ TETM","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/checkTETMCut.pdf");
}


void standardCutReduction(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",750,0,50);
	sprintf(theCut,"%s",basicCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->SetLineColor(kBlack);
	h1->Draw();
	cout << "basic cut: " << cts << endl;

	TAxis *ax1 = h1->GetXaxis();
	printf("Basic cut.  Total cts %i  cts 0.8-5 %.0f  cts 5-50 %.0f\n", cts, h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h1->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s",standardCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");
	cout << "standard cut: " << cts << endl;

	TAxis *ax2 = h2->GetXaxis();
	printf("Standard Cut.  Total cts %i  cts 0.8-5 %.0f  cts 5-50 %.0f\n", cts, h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5)), h2->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.92);
	leg1->AddEntry(h1,"basic","l");
	leg1->AddEntry(h2,"standard","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/standardCutReduction.pdf");
}


void cutTreeStandard(TChain* skim)
{
	skim->Draw(">>entList",standardCut,"entrylist GOFF");
	TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
	skim->SetEntryList(entList);
	TFile *f2 = new TFile("standardCutLE_DS1.root","recreate");
	TTree *small = skim->CopyTree("");
	small->Write();
	// small->Print();	// prints branch names & info
	f2->Close();

	// you must always still use the same cuts!
	// even after making this small file!
}


void enrVsNat1D(TChain* skim)
{
	double natMass=0, enrMass=0;
	skim->SetBranchAddress("mAct_nat_kg",&natMass);
	skim->SetBranchAddress("mAct_enr_kg",&enrMass);
	skim->GetEntry(0);
	cout << "Active masses (kg).  Natural " << natMass << "  enr: " << enrMass << endl;

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	// c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",750,0,50);
	sprintf(theCut,"%s && !isEnr",standardCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("cts / kg");
	h1->SetLineColor(kCyan+4);
	h1->Scale(1/natMass);
	h1->Draw();
	h1->SetMaximum(25);

	TAxis *ax1 = h1->GetXaxis();
	printf("Natural det's.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h1->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s && isEnr",standardCut);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->SetLineColor(kRed);
	h2->Scale(1/enrMass);
	h2->Draw("SAME");

	TAxis *ax2 = h2->GetXaxis();
	printf("Enriched det's.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5)), h2->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TLegend* leg1 = new TLegend(0.6,0.6,0.87,0.92);
	leg1->AddEntry(h1,TString::Format("nat: %.2fkg",natMass),"l");
	leg1->AddEntry(h2,TString::Format("enr: %.1fkg",enrMass),"l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("output/enrVsNat1D.pdf");
}


void runVsEnergy(TChain* skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogz(1);

	TH2D *h1 = new TH2D("h1","h1",50,0,50,4970,9420,14390);
	skim->Project("h1","run:trapENFCal",standardCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("Run");
	h1->GetYaxis()->SetTitleOffset(1.4);
	h1->Draw("COLZ");
	c1->Print("output/runVsEnergy.pdf");

	TH2D *h2 = new TH2D("h2","h2",50,0,50,4970,9420,14390);
	sprintf(theCut,"%s && isEnr",standardCut);
	skim->Project("h2","run:trapENFCal",theCut);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->GetYaxis()->SetTitle("Run");
	h2->GetYaxis()->SetTitleOffset(1.4);
	h2->Draw("COLZ");
	c1->Print("output/runVsEnergyEnr.pdf");
}


void chanVsEnergy(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogz(1);

	TH2D *h1 = new TH2D("h1","h1",50,0,50,125,575,700);
	skim->Project("h1","channel:trapENFCal",standardCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("Channel");
	h1->GetYaxis()->SetTitleOffset(1.1);
	h1->Draw("COLZ");
	c1->Print("output/channelVsEnergy.pdf");

	TH2D *h2 = new TH2D("h2","h2",50,0,50,125,575,700);
	sprintf(theCut,"%s && isEnr",standardCut);
	skim->Project("h2","channel:trapENFCal",theCut);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->GetYaxis()->SetTitle("Channel");
	h2->GetYaxis()->SetTitleOffset(1.1);
	h2->Draw("COLZ");
	c1->Print("output/channelVsEnergyEnr.pdf");
}


void beforeAfter10500(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogz(1);

	TH2D *h1 = new TH2D("h1","h1",50,0,50,125,575,700);
	sprintf(theCut,"%s && run < 10500",standardCut);
	skim->Project("h1","channel:trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("Channel");
	h1->GetYaxis()->SetTitleOffset(1.1);
	h1->Draw("COLZ");
	c1->Print("output/channelVsEnergy_before10500.pdf");

	TH2D *h2 = new TH2D("h2","h2",50,0,50,125,575,700);
	sprintf(theCut,"%s && run > 10500",standardCut);
	skim->Project("h2","channel:trapENFCal",theCut);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->GetYaxis()->SetTitle("Channel");
	h2->GetYaxis()->SetTitleOffset(1.1);
	h2->Draw("COLZ");
	c1->Print("output/channelVsEnergy_after10500.pdf");

	c1->SetLogy();

	TH1D *h3 = new TH1D("h3","h3",150,0,10);
	sprintf(theCut,"%s && channel==610 && run < 10500",standardCut);
	int cts = skim->Project("h3","trapENFCal",theCut);
	h3->SetLineColor(kBlack);
	h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h3->Draw();

	TH1D *h4 = new TH1D("h4","h4",150,0,10);
	sprintf(theCut,"%s && channel==610 && run > 10500",standardCut);
	cts = skim->Project("h4","trapENFCal",theCut);
	h4->SetLineColor(kRed);
	h4->Draw("SAME");

	TLegend* leg1 = new TLegend(0.5,0.6,0.85,0.9);
	leg1->AddEntry(h3,"P3D2 before 10500","l");
	leg1->AddEntry(h4,"P3D2 after 10500","l");
	leg1->Draw("SAME");
	c1->Update();
	c1->Print("./output/P3D2_beforeAfter10500.pdf");

	// is the p3d2 noise from ln fills?  (it's not bursts ... it's too evenly spread across runs)

	TH1D *h5 = new TH1D("h5","h5",150,0,10);
	sprintf(theCut,"isLNFill && channel==610 && run < 10500 && trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && trapETailMin<0");
	cts = skim->Project("h5","trapENFCal",theCut);
	h5->SetLineColor(kOrange);
	h5->Draw("SAME");

	TH1D *h6 = new TH1D("h6","h6",150,0,10);
	sprintf(theCut,"isLNFill && channel==610 && run > 10500 && trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && trapETailMin<0");
	cts = skim->Project("h6","trapENFCal",theCut);
	h6->SetLineColor(kGreen);
	h6->Draw("SAME");

	TLegend* leg2 = new TLegend(0.5,0.6,0.85,0.9);
	leg2->AddEntry(h3,"P3D2 before 10500","l");
	leg2->AddEntry(h4,"P3D2 after 10500","l");
	leg2->AddEntry(h5,"P3D2 LN before 10500","l");
	leg2->AddEntry(h6,"P3D2 LN after 10500","l");
	leg2->Draw("SAME");
	c1->Update();
	c1->Print("./output/P3D2_LN_beforeAfter10500.pdf");
}


void enrichedSpec(TChain *skim)
{
	double natMass=0, enrMass=0;
	skim->SetBranchAddress("mAct_nat_kg",&natMass);
	skim->SetBranchAddress("mAct_enr_kg",&enrMass);
	skim->GetEntry(0);
	cout << "Active masses (kg).  Natural " << natMass << "  enr: " << enrMass << endl;

	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	// c1->SetLogy();

	TH1D *h1 = new TH1D("h1","h1",225,0,15);
	sprintf(theCut,"%s && isEnr",standardCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("cts / kg");
	h1->SetLineColor(kBlue);
	h1->Scale(1/enrMass);
	h1->Draw();
	h1->SetMaximum(2.5);

	TAxis *ax1 = h1->GetXaxis();
	printf("Enriched det's.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h1->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	c1->Print("output/enrichedSpec.pdf");
}


void PrintSkimChannels(TChain *skim)
{
	vector<int> chanList;
	sprintf(theCut,"%s && !isEnr",standardCut);
	int cts = skim->Draw("channel",theCut,"GOFF");
	cout << "found " << cts << " entries\n";
	for (int i = 0; i < cts; i++){
		int chan = skim->GetV1()[i];

		bool newChan = true;
		for (int j = 0; j<(int)chanList.size(); j++){
			if (chan == chanList[j]) {
				newChan = false;
				break;
			}
		}
		if (newChan) chanList.push_back(chan);
	}

	cout << "Final channel list: \n";
	for (int i = 0; i < (int)chanList.size(); i++) cout << chanList[i] << "  ";
	cout << endl;
}


void toeSpec(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",850,600);
	c1->SetLogz();
	c1->SetLogx();

	TH2D *h1 = new TH2D("h1","h1",1000,0,200,800,0,10);

	sprintf(theCut,"%s && isEnr",standardCut);
	skim->Project("h1","kvorrT/trapENFCal:trapENFCal",theCut);

	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->GetYaxis()->SetTitle("T / E : enriched");
	h1->GetXaxis()->SetTitleOffset(1.2);
	h1->SetMinimum(1);
	// h1->SetMaximum(100);
	h1->Draw("COLZ");

	c1->Update();
	c1->Print("output/enrichedT-E.pdf");

	// a dumb t/e cut, not tuned by detector

	c1->SetLogy();
	c1->SetLogx(0);

	TH1D *h2 = new TH1D("h2","h2",750,0,50);
	sprintf(theCut,"%s && isEnr",standardCut);
	int cts = skim->Project("h2","trapENFCal",theCut);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->GetYaxis()->SetTitle("cts / kg");
	h2->SetLineColor(kBlue);
	h2->Draw();
	// h2->SetMaximum(2.5);

	TAxis *ax1 = h2->GetXaxis();
	printf("Enriched det's.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h2->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h2->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	TH1D *h3  = new TH1D("h3","h3",750,0,50);
	sprintf(theCut,"%s && isEnr && (kvorrT/trapENFCal) >= 0 && (kvorrT/trapENFCal) <= 2.5",standardCut);
	cts = skim->Project("h3","trapENFCal",theCut);
	h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h3->GetYaxis()->SetTitle("cts / kg");
	h3->SetLineColor(kRed);
	h3->Draw("SAME");

	TAxis *ax2 = h3->GetXaxis();
	printf("With T/E.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h3->Integral(ax2->FindBin(0.8),ax2->FindBin(5)), h3->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TLegend* leg1 = new TLegend(0.5,0.6,0.85,0.9);
	leg1->AddEntry(h2,"standard","l");
	leg1->AddEntry(h3,"with simple T/E","l");
	leg1->Draw("SAME");

	c1->Print("output/enrichedDumbT-E.pdf");
}


void tuneDetToe(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",1800,800);
	c1->Divide(2,1,0,0);

	TH1D *htoe[16];
	TH2D *htoe2D[16];
	char hname[100];
	for (int i=0; i<16; i++){

		c1->cd(1);

	 	TVirtualPad *p1 = c1->GetPad(1);
		p1->SetTopMargin(0.05);
		p1->SetLogx(0);

		sprintf(hname,"h_%i",chans[i]);
		htoe[i] = new TH1D(hname,hname,200,0,10);

		sprintf(theCut,"%s && channel==%i",standardCut,chans[i]);
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

		sprintf(hname,"./output/toe_%i.pdf",chans[i]);
		c1->Print(hname);
	}
}

void enrTunedToe(TChain *skim)
{
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas", 850,600);
	c1->SetLogy();

	TH1D *h1  = new TH1D("h1","h1",100,0,20);
	sprintf(theCut,"%s && isEnr",standardCut);
	int cts = skim->Project("h1","trapENFCal",theCut);
	h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h1->SetLineColor(kBlack);
	h1->Draw();

	TAxis *ax1 = h1->GetXaxis();
	printf("With T/E.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h1->Integral(ax1->FindBin(0.8),ax1->FindBin(5)), h1->Integral(ax1->FindBin(5),ax1->FindBin(50)));

	TH1D *h2  = new TH1D("h2","h2",100,0,20);
	sprintf(theCut,"%s && %s && isEnr",standardCut,tunedToeEnr);
	cts = skim->Project("h2","trapENFCal",theCut);
	h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
	h2->SetLineColor(kRed);
	h2->Draw("SAME");

	TAxis *ax2 = h2->GetXaxis();
	printf("With tuned T/E.  Total cts %i  int 0.8-5 %.3f  int 5-50 %.3f\n", cts, h2->Integral(ax2->FindBin(0.8),ax2->FindBin(5)), h2->Integral(ax2->FindBin(5),ax2->FindBin(50)));

	TLegend* leg1 = new TLegend(0.5,0.6,0.85,0.9);
	leg1->AddEntry(h1,"standard","l");
	leg1->AddEntry(h2,"with tuned T/E","l");
	leg1->Draw("SAME");

	c1->Print("./output/tunedToeEnrCoarse.pdf");

}
// ============================================================================================


int main(int argc, char** argv)
{
	gROOT->ProcessLine(".x /Users/wisecg/env/MJDClintPlotStyle.C");

	// TChain *skim = new TChain("skimTree");
	// skim->Add("/Users/wisecg/datasets/ds1-low/*.root");
	// // skim->Add("/Users/wisecg/datasets/ds1-low/skimDS1_0.root");
	// long entries = skim->GetEntries();
	// printf("\n\nFound %li entries.\n",entries);
	// cout << "\a";	// beep when done

	// TBenchmark b;

	// b.Start("compareCuts1D");
	// compareCuts1D(skim);
	// b.Show("compareCuts1D");

	// combineCuts1D(skim);
	// checkLNFill(skim);
	// checkBurstCut(skim);
	// checkTETMCut(skim);
	// standardCutReduction(skim);
	// b.Reset();

	// This function cuts down the tree based on the "standard low-energy cut"
	// cutTreeStandard(skim);

	TChain *cutSkim = new TChain("skimTree");
	cutSkim->Add("standardCutLE_DS1.root");
	cout << "\nUsing cut tree with " << cutSkim->GetEntries() << " entries.\n";

	// enrVsNat1D(cutSkim);

	// gStyle->SetPadLeftMargin(0.16);
	// runVsEnergy(cutSkim);

	// gStyle->SetPadLeftMargin(0.12);	// back to normal
	// chanVsEnergy(cutSkim);

	// beforeAfter10500(cutSkim);

	// enrichedSpec(cutSkim);

	// PrintSkimChannels(cutSkim);

	// toeSpec(cutSkim);

	// tuneDetToe(cutSkim);

	enrTunedToe(cutSkim);
}
