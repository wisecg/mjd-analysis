#include "iostream"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TLine.h"

#include "GATDataSet.hh"

using namespace std;

int main(int argc, char** argv)
{
	TChain *skim = new TChain("skimTree");
	skim->Add("/Users/wisecg/dev/ds1/*.root");
	long entries = skim->GetEntries();
	printf("\n\nFound %li entries.\n",entries);

	// initialize all skim file branches, because what the hell
	unsigned int gatrev=0,EventDC1Bits=0;
	int run=0,iEvent=0,mH=0,mL=0;
	double startTime=0,stopTime=0,sumEH=0,sumEL=0;
	vector<int> *iHit=0,*channel=0,*P=0,*D=0,*gain=0,*mageID=0,*detID=0,*dateMT=0,*muType=0;
	vector<double> *mAct_g=0,*tloc_s=0,*time_s=0,*timeMT=0,*trapECal=0,*trapENFCal=0,*aenorm=0,*kvorrT=0,*toe=0,*dcrSlope=0,*trapETailMin=0,*dtmu_s=0;
	vector<bool> *isEnr=0,*isNat=0,*isGood=0,*isLNFill=0,*badScaler=0,*muVeto1ms=0;
	vector<string> *detName;
	vector<unsigned int> *wfDCBits=0;
	skim->SetBranchAddress("gatrev",&gatrev);
	skim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
	skim->SetBranchAddress("run",&run);
	skim->SetBranchAddress("iEvent",&iEvent);
	skim->SetBranchAddress("mH",&mH);
	skim->SetBranchAddress("mL",&mL);
	skim->SetBranchAddress("startTime",&startTime);
	skim->SetBranchAddress("stopTime",&stopTime);
	skim->SetBranchAddress("sumEH",&sumEH);
	skim->SetBranchAddress("sumEL",&sumEL);
	skim->SetBranchAddress("iHit",&iHit);
	skim->SetBranchAddress("channel",&channel);
	skim->SetBranchAddress("P",&P);
	skim->SetBranchAddress("D",&D);
	skim->SetBranchAddress("gain",&gain);
	skim->SetBranchAddress("mageID",&mageID);
	skim->SetBranchAddress("detID",&detID);
	skim->SetBranchAddress("dateMT",&dateMT);
	skim->SetBranchAddress("muType",&muType);
	skim->SetBranchAddress("mAct_g",&mAct_g);
	skim->SetBranchAddress("tloc_s",&tloc_s);
	skim->SetBranchAddress("time_s",&time_s);
	skim->SetBranchAddress("timeMT",&timeMT);
	skim->SetBranchAddress("trapECal",&trapECal);
	skim->SetBranchAddress("trapENFCal",&trapENFCal);
	skim->SetBranchAddress("aenorm",&aenorm);
	skim->SetBranchAddress("kvorrT",&kvorrT);  // in GAT: trirt100nsft10nsMax / trapECal 
	skim->SetBranchAddress("toe",&toe);
	skim->SetBranchAddress("dcrSlope",&dcrSlope);
	skim->SetBranchAddress("trapETailMin",&trapETailMin);	// events with positive trapETailMin are removed
	skim->SetBranchAddress("dtmu_s",&dtmu_s);
	skim->SetBranchAddress("isEnr",&isEnr);
	skim->SetBranchAddress("isNat",&isNat);
	skim->SetBranchAddress("isGood",&isGood);
	skim->SetBranchAddress("isLNFill",&isLNFill);
	skim->SetBranchAddress("badScaler",&badScaler);
	skim->SetBranchAddress("muVeto1ms",&muVeto1ms);
	skim->SetBranchAddress("detName",&detName);
	skim->SetBranchAddress("wfDCBits",&wfDCBits);

	// For reference, full event dump
	// for (long i = 0; i < 3; i++)
	// {
	// 	skim->GetEntry(i);
	// 	printf("Event: %li  Hits this event: %lu\n",i,channel->size());
	// 	printf("gatrev %u  EventDC1Bits %u  run %i  iEvent %i  mH %i  mL %i\n",gatrev,EventDC1Bits,run,iEvent,mH,mL);
	// 	printf("start %.0f  stop: %.0f  sumEH %.2f  sumEL %.2f\n\n",startTime,stopTime,sumEH,sumEL);
	// 	for (int j = 0; j < channel->size(); j++)
	// 	{
	// 		printf("\tiHit %i  chan %i  P %i  D %i  gain %i  mageID %i  detID %i  dateMT %i  muType %i\n"
	// 			,iHit->at(j),channel->at(j),P->at(j),D->at(j),gain->at(j),mageID->at(j),detID->at(j),dateMT->at(j),muType->at(j));
	// 		printf("\tmAct_g %.2f  tloc_s %.2f  time_s %.2f  timeMT %.2f trapECal %.2f  trapENFCal %.2f\n"
	// 			,mAct_g->at(j),tloc_s->at(j),time_s->at(j),timeMT->at(j),trapECal->at(j),trapENFCal->at(j));
	// 		printf("\taenorm %.2f  kvorrT %.2f  toe %.2f  dcrSlope %.2f  trapETailMin %.2f  dtmu_s %.2f\n"
	// 			,aenorm->at(j),kvorrT->at(j),toe->at(j),dcrSlope->at(j),trapETailMin->at(j),dtmu_s->at(j));
	// 		printf("\tisEnr %i  isNat	%i  isGood %i  isLNFill %i  badScaler %i  muVeto1ms %i\n\n"
	// 			,(bool)isEnr->at(j),(bool)isNat->at(j),(bool)isGood->at(j),(bool)isLNFill->at(j),(bool)badScaler->at(j),(bool)muVeto1ms->at(j));
	// 	}
	// 	cout << "\n";
	// }

	// For reference, an array of all channels (./py/skim.py)
	// channels range from 575 to 700
	//
	int enabled[36] = {640, 641, 648, 649, 664, 665, 672, 673, 690, 691, 692, 693, 
					   578, 579, 580, 581, 582, 583, 592, 593, 594, 595, 598, 599, 
					   600, 601, 608, 609, 610, 611, 616, 617, 626, 627, 632, 633};
	
	int highGain[18] = {640, 648, 664, 672, 690, 692, 
					   	578, 580, 582, 592, 594, 598,
					   	600, 608, 610, 616, 626, 632};

	int lowGain[18] = {641, 649, 665, 673, 691, 693, 
			    	   579, 581, 583, 593, 594, 599,
					   601, 609, 611, 617, 627, 633};

	gROOT->ProcessLine(".x /Users/wisecg/dev/MJDClintPlotStyle.C");

	// drawing
	char drawstr[500];
	char eDraw[50] = "trapENFCal";
	char t_eDraw[50] = "kvorrT/trapENFCal:trapENFCal";
	char t_eAccDraw[50] = "kvorrT/trapENFCal";
	char e_chanDraw[50] = "channel:trapENFCal";
	char chanVRunDraw[50] = "channel:run";

 	// cutting
	char cutstr[500];
	char energyCut[50] = "trapENFCal>0 && trapENFCal<50";
	char channelCut[50] = "channel==648";
	char HGCut[50] = "gain==0";
	char SSCut[50] = "mH==1";
	char IsEnr[50] = "isEnr";
	char DCBitCut[50] = "!wfDCBits";
	char LNFillCut[50] = "!isLNFill";
	char muonCut[50] = "(dtmu_s < -0.2e-3 || dtmu_s > 1)";
	char TECut[100] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";
	char TailMinCut[50] = "trapETailMin < 0";

	///////////////////////////////////////////////////////////
	//                                                       //
	// ===================== 1-D PLOTS ===================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	// %-------------------------------------% \\	

	// Comparison of high-gain, single-site, and T/E cuts
	if(0)
	{
		TCanvas *c1_1 = new TCanvas("c1_1", "Bob Ross's Canvas",800,600);
		c1_1->SetLogy();

		TH1D *h0 = new TH1D("h0","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s",energyCut,HGCut);
		skim->Project("h0", drawstr, cutstr);

		TH1D *h1 = new TH1D("h1","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SSCut);
		skim->Project("h1", drawstr, cutstr);

		TH1D *h2 = new TH1D("h2","",100,0,50);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,TECut);
		skim->Project("h2", drawstr, cutstr);

		h0->SetLineColor(kBlack);
		h0->GetXaxis()->SetTitle(drawstr);
		h0->GetYaxis()->SetTitle("Counts");
		h0->Draw();

		h1->SetLineColor(kBlue);
		h1->Draw("same");

		h2->SetLineColor(kRed);
		h2->Draw("same");

		TLegend* leg1 = new TLegend(0.6,0.7,0.87,0.94);
		leg1->AddEntry(h0,"HG","l");
		leg1->AddEntry(h1,"HG & SS","l");
		leg1->AddEntry(h2,"HG & SS & T/E","l");
		leg1->Draw();

		c1_1->Update();
		c1_1->Print("output/DS1_Spectrum.pdf");
	}

	// %-------------------------------------% \\

	// Adding LN fill and Data Cleaning bits to HG+SS+T/E cut
	if(0)
	{
		TCanvas *c1_2 = new TCanvas("c1_2", "Bob Ross's Canvas",800,600);
		c1_2->SetLogy();

		TH1D *h3 = new TH1D("h3","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,TECut);
		skim->Project("h3",drawstr,cutstr);

		TH1D *h4 = new TH1D("h4","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,TECut,LNFillCut);
		skim->Project("h4",drawstr,cutstr);	

		TH1D *h5 = new TH1D("h5","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,TECut,DCBitCut);
		skim->Project("h5",drawstr,cutstr);	

		h3->SetLineColor(kRed);
		h3->SetLineWidth(5);
		h3->GetXaxis()->SetTitle(drawstr);
		h3->GetYaxis()->SetTitle("Counts");
		h3->Draw();

		h4->SetLineColor(kGreen);
		h4->SetLineWidth(3);
		h4->Draw("same");

		h5->SetLineColor(kBlack);
		h5->SetLineWidth(1);
		h5->Draw("same");

		TLegend* leg2 = new TLegend(0.6,0.7,0.87,0.94);
		leg2->AddEntry(h3,"HG & SS & T/E","l");
		leg2->AddEntry(h4,"HG & SS & T/E & LN","l");
		leg2->AddEntry(h5,"HG & SS & T/E & DC","l");
		leg2->Draw();

		c1_2->Update();
		c1_2->Print("output/DS1_LNTag.pdf");
	}

	// %-------------------------------------% \\

	// Comparing HG+SS+T/E cut with adding TailMinCut
	if(0)
	{
		TCanvas *c1_3 = new TCanvas("c1_3", "Bob Ross's Canvas",800,600);
		c1_3->SetLogy();

		TH1D *h6 = new TH1D("h6","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,TECut);
		skim->Project("h6",drawstr,cutstr);

		TH1D *h7 = new TH1D("h7","",100,0,50);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,TECut,TailMinCut);
		skim->Project("h7",drawstr,cutstr);	

		h6->SetLineColor(kRed);
		h6->GetXaxis()->SetTitle(drawstr);
		h6->GetYaxis()->SetTitle("Counts");
		h6->Draw();

		h7->SetLineColor(kBlue);
		h7->Draw("same");

		TLegend* leg3 = new TLegend(0.6,0.7,0.87,0.94);
		leg3->AddEntry(h6,"HG & SS & T/E","l");
		leg3->AddEntry(h7,"HG & SS & T/E & TM","l");
		leg3->Draw();

		c1_3->Update();
		c1_3->Print("output/DS1_TE_TailMin.pdf");
	}

	// %-------------------------------------% \\

	// Comparing HG+SS+T/E+TM cut on enriched & natural detectors 
	if(0)
	{
		TCanvas *c1_4 = new TCanvas("c1_4", "Bob Ross's Canvas",800,600);
		// c1_4->SetLogy();

		TH1D *h8 = new TH1D("h8","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s && %s",energyCut,HGCut,SSCut,TECut,TailMinCut,IsEnr);
		skim->Project("h8",drawstr,cutstr);

		TH1D *h9 = new TH1D("h9","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s && !%s",energyCut,HGCut,SSCut,TECut,TailMinCut,IsEnr);
		skim->Project("h9",drawstr,cutstr);	

		h8->SetLineColor(kRed);
		h8->GetXaxis()->SetTitle(drawstr);
		h8->GetYaxis()->SetTitle("Counts");
		h8->Draw();

		h9->SetLineColor(kBlue);
		h9->Draw("same");

		TLegend* leg4 = new TLegend(0.6,0.7,0.87,0.94);
		leg4->SetHeader("HG+SS+TE+TM Cuts");
		leg4->AddEntry(h8,"Enriched","l");
		leg4->AddEntry(h9,"Natural","l");
		leg4->Draw();

		c1_4->Update();
		c1_4->Print("output/DS1_TETM_EnrVsNat.pdf");
	}

	// %-------------------------------------% \\

	// Acceptance of HG+SS+TM cut on T/E for enriched & natural detectors 
	if(0)
	{
		TCanvas *c1_5 = new TCanvas("c1_5", "Bob Ross's Canvas",800,600);
		// c1_5->SetLogy();

		TH1D *h10 = new TH1D("h10","",250,0,5);
		sprintf(drawstr,"%s",t_eAccDraw);
		sprintf(cutstr,"%s && %s && %s && %s && !%s",energyCut,HGCut,SSCut,TailMinCut,IsEnr);
		skim->Project("h10",drawstr,cutstr);

		TH1D *h11 = new TH1D("h11","",250,0,5);
		sprintf(drawstr,"%s",t_eAccDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,TailMinCut,IsEnr);
		skim->Project("h11",drawstr,cutstr);	

		h10->SetLineColor(kRed);
		h10->GetXaxis()->SetTitle(drawstr);
		h10->GetYaxis()->SetTitle("Counts");
		h10->Draw();

		h11->SetLineColor(kBlue);
		h11->Draw("same");

		TLegend* leg4 = new TLegend(0.6,0.7,0.87,0.94);
		leg4->AddEntry(h10,"Natural","l");
		leg4->AddEntry(h11,"Enriched","l");
		leg4->Draw();

		c1_5->Update();
		c1_5->Print("output/DS1_TETM_Acceptance.pdf");
	}

	///////////////////////////////////////////////////////////
	//                                                       //
	// ===================== 2-D PLOTS ===================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	// T/E plot for HG channels and SS events
	if(0)
	{
		TCanvas *c2_1 = new TCanvas("c2_1","Bob Ross's Canvas",850,600);
		c2_1->SetLogz();

		TH2D *h2D_1 = new TH2D("h2D_1","",100,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SSCut);
		skim->Project("h2D_1",drawstr,cutstr);

		h2D_1->GetXaxis()->SetTitle(eDraw);
		h2D_1->GetYaxis()->SetTitle("T / E");
		h2D_1->SetMinimum(1);
		h2D_1->SetMaximum(100);
		h2D_1->Draw("COLZ");

		TLine *line1 = new TLine(0,1.2,50,1.2);
	 	line1->SetLineColor(kRed);
	 	line1->SetLineWidth(3);
	 	line1->Draw();

	 	TLine *line2 = new TLine(0,2.1,50,2.1);
	 	line2->SetLineColor(kRed);
	 	line2->SetLineWidth(3);
	 	line2->Draw();
		
		c2_1->Update();
		c2_1->Print("output/DS1_TvsE.pdf");
	}

	// %-------------------------------------% \\

	// Comparing Enriched and Natural T/E (2 plots)
	if(0)
	{
		TCanvas *c2_2 = new TCanvas("c2_2","Bob Ross's Canvas",850,600);
		c2_2->SetLogz();

		TH2D *h2D_2 = new TH2D("h2D_2","",100,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,IsEnr);
		skim->Project("h2D_2",drawstr,cutstr);

		h2D_2->GetXaxis()->SetTitle(eDraw);
		h2D_2->GetYaxis()->SetTitle("T / E : enriched");
		h2D_2->SetMinimum(1);
		h2D_2->SetMaximum(100);
		h2D_2->Draw("COLZ");

		c2_2->Update();
		c2_2->Print("output/DS1_TE_Enriched.pdf");

		TCanvas *c2_3 = new TCanvas("c2_3","Bob Ross's Canvas",850,600);
		c2_3->SetLogz();

		TH2D *h2D_3 = new TH2D("h2D_3","",100,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s & !%s",energyCut,HGCut,SSCut,IsEnr);
		skim->Project("h2D_3",drawstr,cutstr);

		h2D_3->GetXaxis()->SetTitle(eDraw);
		h2D_3->GetYaxis()->SetTitle("T / E : natural");
		h2D_3->SetMinimum(1);
		h2D_3->SetMaximum(100);
		h2D_3->Draw("COLZ");

		c2_3->Update();
		c2_3->Print("output/DS1_TE_Natural.pdf");
	}

	// %-------------------------------------% \\

	// Comparing the effect of the tail min cut on enriched T/E
	if(0)
	{
		TCanvas *c2_4 = new TCanvas("c2_4","Bob Ross's Canvas",850,600);
		c2_4->SetLogz();

		TH2D *h2D_4 = new TH2D("h2D_4","",100,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,IsEnr,TailMinCut);
		skim->Project("h2D_4",drawstr,cutstr);

		h2D_4->GetXaxis()->SetTitle(eDraw);
		h2D_4->GetYaxis()->SetTitle("T / E : enriched w/ TM Cut");
		h2D_4->SetMinimum(1);
		h2D_4->SetMaximum(100);
		h2D_4->Draw("COLZ");

		c2_4->Update();
		c2_4->Print("output/DS1_TETM_Enriched.pdf");
	}

	// %-------------------------------------% \\

	// Energy vs. Channel, enriched & natural, with and without TM Cut (4 plots)
	if(0)
	{
		// all detectors, no TM cut
		TCanvas *c2_5 = new TCanvas("c2_5","Bob Ross's Canvas",850,600);
		c2_5->SetLogz();

		TH2D *h2D_5 = new TH2D("h2D_5","",40,0,20,125,575,700); 
		sprintf(drawstr,"%s",e_chanDraw);
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SSCut);
		skim->Project("h2D_5",drawstr,cutstr);

		h2D_5->GetXaxis()->SetTitle("trapENFCal");
		h2D_5->GetYaxis()->SetTitle("Channel");
		h2D_5->GetYaxis()->SetNdivisions(515);
		h2D_5->GetYaxis()->SetLabelSize(20);
		h2D_5->SetMinimum(1);
		h2D_5->SetMaximum(100);
		h2D_5->Draw("COLZ");

		c2_5->Update();
		c2_5->Print("output/DS1_ENFvsChannel.pdf");

		// all detectors, with TM cut
		TCanvas *c2_6 = new TCanvas("c2_6","Bob Ross's Canvas",850,600);
		c2_6->SetLogz();

		TH2D *h2D_6 = new TH2D("h2D_6","",40,0,20,125,575,700); 
		sprintf(drawstr,"%s",e_chanDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,TailMinCut);
		skim->Project("h2D_6",drawstr,cutstr);

		h2D_6->GetXaxis()->SetTitle("trapENFCal (TM Cut)");
		h2D_6->GetYaxis()->SetTitle("Channel");
		h2D_6->GetYaxis()->SetNdivisions(515);
		h2D_6->GetYaxis()->SetLabelSize(20);
		h2D_6->SetMinimum(1);
		h2D_6->SetMaximum(100);
		h2D_6->Draw("COLZ");

		c2_6->Update();
		c2_6->Print("output/DS1_TM_ENFvsChannel.pdf");

		// enriched detectors, no TM cut
		TCanvas *c2_7 = new TCanvas("c2_7","Bob Ross's Canvas",850,600);
		c2_7->SetLogz();

		TH2D *h2D_7 = new TH2D("h2D_7","",40,0,20,125,575,700); 
		sprintf(drawstr,"%s",e_chanDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SSCut,IsEnr);
		skim->Project("h2D_7",drawstr,cutstr);

		h2D_7->GetXaxis()->SetTitle("trapENFCal");
		h2D_7->GetYaxis()->SetTitle("Channel (enriched)");
		h2D_7->GetYaxis()->SetNdivisions(515);
		h2D_7->GetYaxis()->SetLabelSize(20);
		h2D_7->SetMinimum(1);
		h2D_7->SetMaximum(100);
		h2D_7->Draw("COLZ");

		c2_7->Update();
		c2_7->Print("output/DS1_Enr_ENFvsChannel.pdf");

		// enriched detectors, with TM cut
		TCanvas *c2_8 = new TCanvas("c2_8","Bob Ross's Canvas",850,600);
		c2_8->SetLogz();

		TH2D *h2D_8 = new TH2D("h2D_8","",40,0,20,125,575,700); 
		sprintf(drawstr,"%s",e_chanDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SSCut,IsEnr,TailMinCut);
		skim->Project("h2D_8",drawstr,cutstr);

		h2D_8->GetXaxis()->SetTitle("trapENFCal");
		h2D_8->GetYaxis()->SetTitle("Channel (enriched)");
		h2D_8->GetYaxis()->SetNdivisions(515);
		h2D_8->GetYaxis()->SetLabelSize(20);
		h2D_8->SetMinimum(1);
		h2D_8->SetMaximum(100);
		h2D_8->Draw("COLZ");

		c2_8->Update();
		c2_8->Print("output/DS1_Enr_TM_ENFvsChannel.pdf");
	}

	// Channel vs. Run, and Channels w/ 20 count minimum at low energy
	if(0)
	{
		TCanvas *c2_9 = new TCanvas("c2_9","Bob Ross's Canvas",850,600);
		c2_9->SetLogz();

		TH2D *h2D_9 = new TH2D("h2D_9","",3950,9420,13370,125,575,700); 
		sprintf(drawstr,"%s",chanVRunDraw);
		sprintf(cutstr,"%s && %s",energyCut,HGCut);
		skim->Project("h2D_9",drawstr,cutstr);

		h2D_9->GetXaxis()->SetTitle("Run");
		h2D_9->GetYaxis()->SetTitle("Channel");
		h2D_9->GetYaxis()->SetNdivisions(515);
		h2D_9->GetYaxis()->SetLabelSize(20);
		h2D_9->SetMinimum(1);
		h2D_9->Draw("COLZ");

		c2_9->Update();
		c2_9->Print("output/DS1_ChanVsRun.pdf");

		TCanvas *c2_10 = new TCanvas("c2_10","Bob Ross's Canvas",850,600);
		c2_10->SetLogz();

		TH2D *h2D_10 = new TH2D("h2D_10","",3950,9420,13370,125,575,700); 
		sprintf(drawstr,"%s",chanVRunDraw);
		sprintf(cutstr,"%s && %s",energyCut,HGCut);
		skim->Project("h2D_10",drawstr,cutstr);

		h2D_10->GetXaxis()->SetTitle("Run");
		h2D_10->GetYaxis()->SetTitle("Channel (>20 counts/run)");
		h2D_10->GetYaxis()->SetNdivisions(515);
		h2D_10->GetYaxis()->SetLabelSize(20);
		h2D_10->SetMinimum(20);
		h2D_10->Draw("COLZ");

		c2_10->Update();
		c2_10->Print("output/DS1_ChanVsRun_20CountMin.pdf");
		c2_10->Print("output/DS1_ChanVsRun_20CountMin.C");	// to make a run list
	}


}