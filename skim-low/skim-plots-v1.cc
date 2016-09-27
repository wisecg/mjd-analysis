#include "iostream"

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TLine.h"

#include "GATDataSet.hh"
#include "GATPeakShape.hh"

using namespace std;

// fitting function
// p0: const, p1: amp, p2: mu, p3: sigma
double gausFBG(Double_t* x, Double_t* p)
{
  Double_t dummy=(x[0]-p[2])/p[3];
  Double_t result=p[0] + p[1]*exp(-0.5*dummy*dummy)/ (2.506 * p[3]);
  return result;
}

int main(int argc, char** argv)
{
	TChain *skim = new TChain("skimTree");
	skim->Add("/Users/wisecg/dev/datasets/ds1/*.root");
	long entries = skim->GetEntries();
	printf("\n\nFound %li entries.\n",entries);

	// skim file branches.
	// May 2016 version.
	// get the list from a raw file with:
	// root[0] skimTree->GetListOfBranches()->Print()
	unsigned int gatrev=0,EventDC1Bits=0;
	int run=0,iEvent=0,mH=0,mL=0;
	double startTime=0,stopTime=0,sumEH=0,sumEL=0;
	vector<int> *iHit=0,*channel=0,*P=0,*D=0,*gain=0,*mageID=0,*detID=0,*dateMT=0,*muType=0;
	vector<double> *mAct_g=0,*tloc_s=0,*time_s=0,*timeMT=0,*trapECal=0,*trapENFCal=0,*aenorm=0,*kvorrT=0,*toe=0,*dcrSlope95=0,*dcrSlope99=0,*trapETailMin=0,*dtmu_s=0,*trap4usMax=0,*t150=0;
	vector<bool> *isEnr=0,*isNat=0,*isGood=0,*isLNFill=0,*badScaler=0,*muVeto1ms=0;
	// vector<string> *detName;
	vector<unsigned int> *wfDCBits=0;
	skim->SetBranchAddress("gatrev",&gatrev);
	skim->SetBranchAddress("EventDC1Bits",&EventDC1Bits);
	skim->SetBranchAddress("run",&run);
	skim->SetBranchAddress("iEvent",&iEvent);
	skim->SetBranchAddress("mH",&mH);
	skim->SetBranchAddress("mL",&mL);
	skim->SetBranchAddress("t150",&t150);
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
	skim->SetBranchAddress("trap4usMax",&trap4usMax);
	skim->SetBranchAddress("aenorm",&aenorm);
	skim->SetBranchAddress("kvorrT",&kvorrT);  // in GAT: trirt100nsft10nsMax / trapECal
	skim->SetBranchAddress("toe",&toe);
	skim->SetBranchAddress("dcrSlope95",&dcrSlope95);
	skim->SetBranchAddress("dcrSlope99",&dcrSlope99);
	skim->SetBranchAddress("trapETailMin",&trapETailMin);	// events with positive trapETailMin are removed
	skim->SetBranchAddress("dtmu_s",&dtmu_s);
	skim->SetBranchAddress("isEnr",&isEnr);
	skim->SetBranchAddress("isNat",&isNat);
	skim->SetBranchAddress("isGood",&isGood);
	skim->SetBranchAddress("isLNFill",&isLNFill);
	skim->SetBranchAddress("badScaler",&badScaler);
	skim->SetBranchAddress("muVeto1ms",&muVeto1ms);
	// skim->SetBranchAddress("detName",&detName);
	skim->SetBranchAddress("wfDCBits",&wfDCBits);

	// For reference, an array of all channels (./py/skim.py)
	// channels range from 575 to 700
	//
	int enabled[36] = {640, 641, 648, 649, 664, 665, 672, 673, 690, 691, 692, 693, 
					   578, 579, 580, 581, 582, 583, 592, 593, 594, 595, 598, 599, 
					   600, 601, 608, 609, 610, 611, 616, 617, 626, 627, 632, 633};
	
	int highGain[18] = {640, 648, 664, 672, 690, 692, 
					   	578, 580, 582, 592, 594, 598,
					   	600, 608, 610, 616, 626, 632};

   	// corresponds to the highGain array
	string hgPos[18] = {"P2D3","P2D2","P3D4","P5D3","P6D4","P2D1",
						"P1D4","P1D3","P1D2","P7D4","P7D3","P7D2",
						"P7D1","P3D3","P3D2","P3D1","P6D3","P6D1"};

    string detEnr[18] = {"enr","enr","enr","enr","enr","nat",
    					 "enr","enr","enr","enr","enr","enr",
			  			 "nat","enr","enr","nat","enr","enr"};

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
	char energyCut[50] = "trapENFCal>0 && trapENFCal<100";
	char channelCut[50] = "channel==648";
	char HGCut[50] = "gain==0";
	char SDCut[50] = "mH==1";
	char IsEnr[50] = "isEnr";
	char DCBitCut[50] = "!wfDCBits";
	char LNFillCut[50] = "!isLNFill";
	char muonCut[50] = "(dtmu_s < -0.2e-3 || dtmu_s > 1)";
	char TECut[100] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";
	char TMCut[50] = "trapETailMin < 0";
	// char burstCut[500] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";
	char burstCut[500] = "1";
	// char runCut[500] = "run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863";
	char runCut[500] = "1";


	// binning
	// kris's binning: 15 bins/kev.  Ex: 750,0,50
	// clint's binning: 2 bins/kev.  Ex: 100,0,50
	// kris's second one: 5 bins/kev.  Ex: 250,0,50

	///////////////////////////////////////////////////////////
	//                                                       //
	// ===================== 1-D PLOTS ===================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	// %-------------------------------------% \\

	// Comparison of high-gain, single-detector, and T/E cuts
	if(0)
	{
		TCanvas *c1_1 = new TCanvas("c1_1", "Bob Ross's Canvas",800,600);
		c1_1->SetLogy();

		TH1D *h0 = new TH1D("h0","",750,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s",energyCut,HGCut);
		skim->Project("h0", drawstr, cutstr);

		TH1D *h1 = new TH1D("h1","",750,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SDCut);
		skim->Project("h1", drawstr, cutstr);

		TH1D *h2 = new TH1D("h2","",750,0,50);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,TECut);
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
		leg1->AddEntry(h1,"HG & SD","l");
		leg1->AddEntry(h2,"HG & SD & T/E","l");
		leg1->Draw();

		c1_1->Update();
		c1_1->Print("output/DS1_Spectrum.pdf");
	}

	// %-------------------------------------% \\

	// Adding LN fill and Data Cleaning bits to HG+SD+T/E cut
	if(0)
	{
		TCanvas *c1_2 = new TCanvas("c1_2", "Bob Ross's Canvas",800,600);
		c1_2->SetLogy();

		TH1D *h3 = new TH1D("h3","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",HGCut,SDCut,TECut);
		skim->Project("h3",drawstr,cutstr);

		TH1D *h4 = new TH1D("h4","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",HGCut,SDCut,DCBitCut);
		skim->Project("h4",drawstr,cutstr);	

		TH1D *h5 = new TH1D("h5","",100,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",HGCut,SDCut,TECut,DCBitCut);
		skim->Project("h5",drawstr,cutstr);	

		h3->SetLineColor(kRed);
		h3->SetLineWidth(5);
		h3->GetXaxis()->SetTitle(drawstr);
		h3->GetYaxis()->SetTitle("Counts");
		h3->Draw();

		h4->SetLineColor(kBlue);
		h4->SetLineWidth(3);
		h4->Draw("same");

		h5->SetLineColor(kBlack);
		h5->SetLineWidth(1);
		h5->Draw("same");

		TLegend* leg2 = new TLegend(0.6,0.7,0.87,0.94);
		leg2->AddEntry(h3,"HG & SD & T/E","l");
		leg2->AddEntry(h4,"HG & SD & DC","l");
		leg2->AddEntry(h5,"HG & SD & T/E & DC","l");
		leg2->Draw();

		c1_2->Update();
		c1_2->Print("output/DS1_LNTag.pdf");
	}

	// %-------------------------------------% \\

	// Comparing DC cut to SD+T/E out to 500keV
	if(0)
	{
		TCanvas *c1_2a = new TCanvas("c1_2a", "Bob Ross's Canvas",800,600);
		// c1_2a->SetLogy();

		TH1D *h3a = new TH1D("h3a","",100,0,500);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",HGCut,SDCut,TECut);
		skim->Project("h3a",drawstr,cutstr);

		TH1D *h4a = new TH1D("h4a","",100,0,500);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s",HGCut,SDCut,DCBitCut);
		skim->Project("h4a",drawstr,cutstr);	

		TH1D *h5a = new TH1D("h5a","",100,0,500);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",HGCut,SDCut,TECut,DCBitCut);
		skim->Project("h5a",drawstr,cutstr);	

		TH1D *h5b = new TH1D("h5a","",100,0,500);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s",HGCut,SDCut);
		skim->Project("h5b",drawstr,cutstr);	

		h3a->SetLineColor(kRed);
		h3a->SetLineWidth(2);
		h3a->GetXaxis()->SetTitle(drawstr);
		h3a->GetYaxis()->SetTitle("Counts");
		h3a->GetYaxis()->SetRangeUser(0,200);
		h3a->Draw();

		h4a->SetLineColor(kGreen);
		h4a->SetLineWidth(2);
		h4a->Draw("same");

		h5a->SetLineColor(kBlack);
		h5a->SetLineWidth(1);
		h5a->Draw("same");

		h5b->SetLineColor(kBlue);
		h5b->SetLineWidth(1);
		h5b->Draw("same");

		TLegend* leg2 = new TLegend(0.6,0.7,0.87,0.94);
		leg2->AddEntry(h3a,"HG & SD & T/E","l");
		leg2->AddEntry(h4a,"HG & SD & DC","l");
		leg2->AddEntry(h5a,"HG & SD & T/E & DC","l");
		leg2->AddEntry(h5b,"HG & SD","l");
		leg2->Draw();

		c1_2a->Update();
		c1_2a->Print("output/DS1_WFDCBits.pdf");
	}

	// %-------------------------------------% \\

	// Comparing HG+SD+T/E cut with adding TMCut
	if(0)
	{
		TCanvas *c1_3 = new TCanvas("c1_3", "Bob Ross's Canvas",800,600);
		c1_3->SetLogy();

		TH1D *h6 = new TH1D("h6","",750,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,TECut);
		skim->Project("h6",drawstr,cutstr);

		TH1D *h7 = new TH1D("h7","",750,0,50);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SDCut,TECut,TMCut);
		skim->Project("h7",drawstr,cutstr);	

		h6->SetLineColor(kRed);
		h6->GetXaxis()->SetTitle(drawstr);
		h6->GetYaxis()->SetTitle("Counts");
		h6->Draw();

		h7->SetLineColor(kBlue);
		h7->Draw("same");

		TLegend* leg3 = new TLegend(0.6,0.7,0.87,0.94);
		leg3->AddEntry(h6,"HG & SD & T/E","l");
		leg3->AddEntry(h7,"HG & SD & T/E & TM","l");
		leg3->Draw();

		c1_3->Update();
		c1_3->Print("output/DS1_TE_TailMin.pdf");
	}

	// %-------------------------------------% \\

	// Comparing HG+SS+T/E+TM cut on enriched & natural detectors 
	if(1)
	{
		TCanvas *c1_4 = new TCanvas("c1_4", "Bob Ross's Canvas",800,600);
		// c1_4->SetLogy();

		// enriched
		TH1D *h8 = new TH1D("h8","",475,5,100);	
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && %s"
			,energyCut,muonCut,burstCut,runCut,HGCut,SDCut,TECut,TMCut,IsEnr);
		skim->Project("h8",drawstr,cutstr);

		// natural
		TH1D *h8a = new TH1D("nameOfh8a","",475,5,100);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s && %s && %s && %s && !%s"
			,energyCut,muonCut,burstCut,runCut,HGCut,SDCut,TECut,TMCut,IsEnr);
		skim->Project("h8a",drawstr,cutstr);	

		h8a->SetLineColor(kBlue);
		h8a->GetXaxis()->SetTitle("Energy (keV)");
		h8a->GetYaxis()->SetTitle("cts/exposure");
		h8a->Scale(1/61.);	// natural exposure is 61.27 pm 0.3 kg*days
		// h8a->Scale(1/5.);    // scale by bins/keV to get counts/keV.
		h8a->Draw();

		h8->SetLineColor(kRed);
		h8->Scale(1/651.); // enriched exposure is 651.04 pm 1.26 kg*days
		// h8->Scale(1/5.);	// scale by bins/keV to get counts/keV.s
		h8->Draw("same");

		TLegend* leg4 = new TLegend(0.6,0.8,0.85,0.94);
		// leg4->SetHeader("HG+SD+TE+TM Cuts");
		leg4->AddEntry(h8,"Enriched (651 kg-d)","l");
		leg4->AddEntry(h8a,"Natural (61 kg-d)","l");
		leg4->Draw();

		c1_4->Update();
		c1_4->Print("output/DS1_100_TETMLin_EnrVsNat.pdf");

		TFile *f = new TFile("output/DS1_100_TETMLin_EnrVsNat.root","RECREATE");
		h8->Write("enriched",TObject::kOverwrite);
		h8a->Write("natural",TObject::kOverwrite);
		c1_4->Write("both",TObject::kOverwrite);
		f->Close();

	}

	// %-------------------------------------% \\

	// Acceptance of HG+SD+TM cut on T/E for enriched & natural detectors 
	if(0)
	{
		TCanvas *c1_5 = new TCanvas("c1_5", "Bob Ross's Canvas",800,600);
		// c1_5->SetLogy();

		TH1D *h10 = new TH1D("h10","",250,0,5);
		sprintf(drawstr,"%s",t_eAccDraw);
		sprintf(cutstr,"%s && %s && %s && %s && !%s",energyCut,HGCut,SDCut,TMCut,IsEnr);
		skim->Project("h10",drawstr,cutstr);

		TH1D *h11 = new TH1D("h11","",250,0,5);
		sprintf(drawstr,"%s",t_eAccDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SDCut,TMCut,IsEnr);
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

	// %-------------------------------------% \\

	// Comparing HG+SD+T/E+TM cut, detector by detector (18 plots)
	if(0)
	{
		TCanvas *c1_dets[18]; 
		TH1D *h_det_Cuts[18]; 
		TH1D *h_det_NoCuts[18]; 
		char cname[50];
		char hname[50];
		for(int j = 0; j < 18; j++)  
		{
			int chan = highGain[j];	// analysis channel
			string det = hgPos[j];  // detector position

			sprintf(cname,"c1_%d",chan);
			c1_dets[j] = new TCanvas(cname,"Bob Ross's Canvas",800,600);
			c1_dets[j]->SetLogy();
			
			sprintf(hname,"HGSDCut_%d",chan);
			h_det_NoCuts[j] = new TH1D(hname, hname, 200,0,100);
			sprintf(drawstr,"%s",eDraw);
			sprintf(cutstr,"%s && %s && channel==%i",HGCut,SDCut,chan);
			skim->Project(hname,drawstr,cutstr);

			sprintf(hname,"TETMCut_%d",chan);
			h_det_Cuts[j]= new TH1D(hname, hname, 200,0,100);
			sprintf(cutstr,"%s && %s && %s && %s && channel==%i",HGCut,SDCut,TECut,TMCut,chan);
			skim->Project(hname,drawstr,cutstr);

			char xTitle[50];
			sprintf(xTitle,"trapENFCal (%d %s %s)",chan,det.c_str(),detEnr[j].c_str());
			h_det_NoCuts[j]->SetLineColor(kRed);
			h_det_NoCuts[j]->GetXaxis()->SetTitle(xTitle);
			h_det_NoCuts[j]->GetYaxis()->SetTitle("Counts");
			h_det_NoCuts[j]->Draw();

			h_det_Cuts[j]->SetLineColor(kBlue);
			h_det_Cuts[j]->Draw("same");

			TLegend* leg4 = new TLegend(0.6,0.8,0.87,0.94);
			leg4->AddEntry(h_det_NoCuts[j],"HG+SD","l");
			leg4->AddEntry(h_det_Cuts[j],"HG+SD+TE+TM","l");
			leg4->Draw();

			c1_dets[j]->Update();

			char file[100];
			sprintf(file,"output/DS1_Det%i_LESpectrum.pdf",chan);
			c1_dets[j]->Print(file);

			// may need to clear out memory
		}
	}

	// %-------------------------------------% \\

	// Acceptance of T/E with and without TM cut, detector by detector (18 plots)
	// And fitting the TM histo with a Gaussian
	if(0)
	{
		TCanvas *c1_dets2[18]; 
		TH1D *h_det_Cuts2[18];
		TH1D *h_det_Cuts2a[18];
		TF1 *fits[18];
		char cname2[50];
		char hname2[50];
		for(int j = 0; j < 1; j++)  
		{
			int chan = highGain[j];	// analysis channel
			string det = hgPos[j];	// detector position

			sprintf(cname2,"c1_%d",chan);
			c1_dets2[j] = new TCanvas(cname2,"Bob Ross's Canvas",800,600);

			sprintf(hname2,"TETMCut_%d",chan);
			h_det_Cuts2[j]= new TH1D(hname2, hname2, 100,0,5);
			sprintf(drawstr,"%s",t_eAccDraw);
			sprintf(cutstr,"%s && %s && %s && %s && channel==%i",HGCut,SDCut,TMCut,energyCut,chan);
			skim->Project(hname2,drawstr,cutstr);

			sprintf(hname2,"TEVals_%d",chan);
			h_det_Cuts2a[j]= new TH1D(hname2, hname2, 100,0,5);
			sprintf(drawstr,"%s",t_eAccDraw);
			sprintf(cutstr,"%s && %s && %s && channel==%i",HGCut,SDCut,energyCut,chan);
			skim->Project(hname2,drawstr,cutstr);

			char xTitle[50];
			sprintf(xTitle,"kvorrT/trapENFCal (%d %s %s, E<100)",chan,det.c_str(),detEnr[j].c_str());
			h_det_Cuts2a[j]->SetLineColor(kBlue);
			h_det_Cuts2a[j]->SetLineWidth(3);
			h_det_Cuts2a[j]->GetXaxis()->SetTitle(xTitle);
			h_det_Cuts2a[j]->GetYaxis()->SetTitle("Counts");
			h_det_Cuts2a[j]->Draw();

			h_det_Cuts2[j]->SetLineColor(kRed);
			h_det_Cuts2[j]->Draw("same");

			// fitting

			fits[j] = new TF1("f","gaus"); 
			// h_det_Cuts2[j]->Fit("gaus");
			// fits[j] = h_det_Cuts2[j]->GetFunction("gaus");

			int MaxBin = h_det_Cuts2[j]->GetMaximumBin();
			double mean = h_det_Cuts2[j]->GetXaxis()->GetBinCenter(MaxBin);

			fits[j]->FixParameter(1, mean);
			fits[j]->SetParLimits(2, 0.1, 1);
			h_det_Cuts2[j]->Fit(fits[j],"B");
			
			double fitResults[3][2] = {{0}}; // [rows][cols]
			for (int k = 0; k < 3; k++)	
			{
				fitResults[k][0]=fits[j]->GetParameter(k);
				fitResults[k][1]=fits[j]->GetParError(k);
			}
			
			TLegend* leg5 = new TLegend(0.6,0.8,0.87,0.94);
			char fitMean[200];
			sprintf(fitMean,"#mu = %.2f  #sigma = %.2f",fitResults[1][0],fitResults[2][0]);
			leg5->SetHeader(fitMean);
			leg5->AddEntry(h_det_Cuts2a[j],"Without TM","l");
			leg5->AddEntry(h_det_Cuts2[j],"With TM","l");
			leg5->Draw();


			c1_dets2[j]->Update();

			char file[100];
			sprintf(file,"output/DS1_Det%i_TETMAcceptance.pdf",chan);
			c1_dets2[j]->Print(file);

			// may need to clear out memory
		}
	}

	// %-------------------------------------% \\

	// Resolution of natural detectors at low-E
	// Trying to get fitting right
	if(0)
	{
		// spectral lines:
		// 10.37: 68Ge K-shell
		// 6.54: 55Fe
		//
		// full list:
		// 68-Ge (L) = 1.29 keV
		// 65-Zn (L) = 1.1 keV
		// 49-V = 4.97 keV
		// 51-Cr = 5.46 keV
		// 54-Mn = 5.99 keV
		// 55-Fe = 6.54  keV
		// 56,57,58-Co = 7.11 keV
		// 56-Ni = 7.71 keV
		// 65-Zn (K) = 8.98 keV 
		// 68-Ga = 9.66 keV
		// 68-Ge (K) = 10.37 keV
		// 73,74-As = 11.10 keV

		TCanvas *c1_6 = new TCanvas("c1_6", "Bob Ross's Canvas",800,600);
		// c1_6->SetLogy();

		// natural detectors: 
		// channel = 692 (p2d1)  highGain[5], hgPos[5]
		// channel = 600 (p7d1)  highGain[12], hgPos[12]
		// channel = 616 (p3d1)  highGain[15], hgPos[15]
		
		TH1D *h12 = new TH1D("h12","",250,0,50);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && channel==692",HGCut,SDCut,TECut,TMCut);
		skim->Project("h12",drawstr,cutstr);

		// MJD Peak Fitter
		// https://mjcalendar.npl.washington.edu/indico/event/2152/contribution/53/material/slides/0.pdf
		// GATPeakShape ps;
		// ps.EZFit(h12,10.37);
		// ps.EZFit(h12,6.54);
		// h12->GetXaxis()->SetRangeUser(0,20);
		// ps.Draw(h12);
		// ps.DrawComponents(h12);

		// Custom Gaus + Flat BG
		// p0: const, p1: amp, p2: mu, p3: sigma
		// https://root.cern.ch/root/roottalk/roottalk02/1393.html
    	TF1* func=new TF1("func",gausFBG,0,100,4);
    	func->SetParameters(1,100,10.3,0.1); 
    	h12->Fit("func");

		// EXT PARAMETER                                   STEP         FIRST
		// NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
		// 1  p0           1.94848e+00   1.47354e-01   3.11208e-04  -1.68438e-03
		// 2  p1           7.95154e+00   1.45590e+00  -1.08102e-02   1.60102e-05
		// 3  p2           1.00435e+01   4.39698e-02  -2.69875e-05  -7.33295e-03
		// 4  p3           2.16482e-01   4.04765e-02  -1.66817e-04  -3.49865e-03

		// GATPeakShape ps;
		// ps.GetInitialParGuess(h12,10.3);
		// // ps.DumpParameters();

		// // amplitude: 0.394704
		// // #mu: 10.1
		// // #sigma: 0.00874801
		// // f_{t}: 0
		// // #tau: 0.00874801
		// // f_{ht}: 0
		// // #tau_{ht}: 0.00874801
		// // h_{s}: 5.94059
		// // b_{BG}: 17.8218
		// // m_{BG}: 0
		// // q_{BG}: 0

		// ps.LimitPar(GATPeakShape::kFt,0,1);
		// ps.FitTo(h12,true,false,false);
		// ps.DumpParameters();
		// h12->GetXaxis()->SetRangeUser(0,20);
		// ps.Draw(h12);

		char xTitle[50];
		sprintf(xTitle,"trapENFCal (%i %s %s)",highGain[5],hgPos[5].c_str(),detEnr[5].c_str());
		h12->SetLineColor(kBlue);
		h12->GetXaxis()->SetTitle(xTitle);
		h12->GetYaxis()->SetTitle("Counts");
		h12->Draw("same");

		c1_6->Update();
		c1_6->Print("output/DS1_692_PeakFits.pdf");
	}

	// %-------------------------------------% \\

	// DS-1 Hot Boule detectors: P1D3 [580], P5D3 [672], P6D2 [630 - dead]
	if(0)
	{
		TCanvas *c1_7 = new TCanvas("c1_7", "Bob Ross's Canvas",800,600);
		// c1_7->SetLogy();

		TH1D *h13 = new TH1D("h13","",200,0,100);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && trapENFCal>5 && (channel==580||channel==672)",HGCut,SDCut);
		skim->Project("h13",drawstr,cutstr);

		char xTitle[50];
		sprintf(xTitle,"trapENFCal (%i %s & %i %s)",highGain[7],hgPos[7].c_str(),highGain[3],hgPos[3].c_str());
		h13->SetLineColor(kBlue);
		h13->GetXaxis()->SetTitle(xTitle);
		h13->GetYaxis()->SetTitle("Counts");
		h13->Draw();

		int bin1 = h13->GetXaxis()->FindBin(45);
		int bin2 = h13->GetXaxis()->FindBin(47);
		cout << " Counts: " << h13->Integral(bin1,bin2) << endl;

		c1_7->Update();
		c1_7->Print("output/DS1_HotBouleDetetectors.pdf");
	}

	// %-------------------------------------% \\

	// P7D3 : Spectra before and after run 12520 (2 plots)
	if(0)
	{
		// before 12520
		TCanvas *c1_8 = new TCanvas("c1_8", "Bob Ross's Canvas",800,600);
		c1_8->SetLogy();

		TH1D *h14 = new TH1D("h14","",100,0,20);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && (channel==594) && run < 12500",HGCut,SDCut);
		skim->Project("h14",drawstr,cutstr);

		TH1D *h15 = new TH1D("h15","",100,0,20);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && (channel==594) && run < 12500",HGCut,SDCut,TECut,TMCut);
		skim->Project("h15",drawstr,cutstr);

		char xTitle[50];
		sprintf(xTitle,"trapENFCal (%i %s %s)",highGain[10],hgPos[10].c_str(),detEnr[10].c_str());
		h14->SetLineColor(kBlue);
		h14->GetXaxis()->SetTitle(xTitle);
		h14->GetYaxis()->SetTitle("Counts before run 12520");
		h14->Draw();

		h15->SetLineColor(kRed);
		h15->Draw("same");

		TLegend* leg5 = new TLegend(0.6,0.8,0.87,0.94);
		leg5->AddEntry(h14,"HG+SD","l");
		leg5->AddEntry(h15,"HG+SD+TE+TM","l");
		leg5->Draw();

		c1_8->Update();
		c1_8->Print("output/DS1_P7D3Before12520.pdf");

		// after 12520
		TCanvas *c1_9 = new TCanvas("c1_9", "Bob Ross's Canvas",800,600);
		c1_9->SetLogy();

		TH1D *h16 = new TH1D("h16","",100,0,20);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && (channel==594) && run > 12500",HGCut,SDCut);
		skim->Project("h16",drawstr,cutstr);

		TH1D *h17 = new TH1D("h17","",100,0,20);
		sprintf(drawstr,"%s",eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && (channel==594) && run > 12500",HGCut,SDCut,TECut,TMCut);
		skim->Project("h17",drawstr,cutstr);

		sprintf(xTitle,"trapENFCal (%i %s %s)",highGain[10],hgPos[10].c_str(),detEnr[10].c_str());
		h16->SetLineColor(kRed);
		h16->GetXaxis()->SetTitle(xTitle);
		h16->GetYaxis()->SetTitle("Counts after run 12520");
		h16->Draw();

		h17->SetLineColor(kBlue);
		h17->Draw("same");

		TLegend* leg6 = new TLegend(0.6,0.8,0.87,0.94);
		leg6->AddEntry(h16,"HG+SD","l");
		leg6->AddEntry(h17,"HG+SD+TE+TM","l");
		leg6->Draw();

		c1_9->Update();
		c1_9->Print("output/DS1_P7D3After12520.pdf");
	}

	///////////////////////////////////////////////////////////
	//                                                       //
	// ===================== 2-D PLOTS ===================== //
	//                                                       //
	///////////////////////////////////////////////////////////

	// T/E plot for HG channels and SD events
	if(0)
	{
		TCanvas *c2_1 = new TCanvas("c2_1","Bob Ross's Canvas",850,600);
		c2_1->SetLogz();

		TH2D *h2D_1 = new TH2D("h2D_1","",750,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SDCut);
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

		TH2D *h2D_2 = new TH2D("h2D_2","",750,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,IsEnr);
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

		TH2D *h2D_3 = new TH2D("h2D_3","",750,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s & !%s",energyCut,HGCut,SDCut,IsEnr);
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

		TH2D *h2D_4 = new TH2D("h2D_4","",750,0,50,800,0,10); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SDCut,IsEnr,TMCut);
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
		sprintf(cutstr,"%s && %s && %s",energyCut,HGCut,SDCut);
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
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,TMCut);
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
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,IsEnr);
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
		sprintf(cutstr,"%s && %s && %s && %s && %s",energyCut,HGCut,SDCut,IsEnr,TMCut);
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

	// %-------------------------------------% \\

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

	// %-------------------------------------% \\

	// Comparing T/E spectra, detector by detector (18 plots)
	if(0)
	{
		TCanvas *c2_dets[18]; 
		TH2D *h2D_det[18]; 
		char c2name[50];
		char h2name[50];
		for(int j = 0; j < 18; j++)  
		{
			int chan = highGain[j];	// analysis channel
			string det = hgPos[j];	// detector position

			sprintf(c2name,"c1_%d",chan);
			c2_dets[j] = new TCanvas(c2name,"Bob Ross's Canvas",800,600);
			c2_dets[j]->SetLogy();
			c2_dets[j]->SetLogz();

			sprintf(h2name,"TE_chan%d",chan);
			h2D_det[j] = new TH2D(h2name,h2name,750,0,50,1000,0.1,100); 
			sprintf(drawstr,"%s",t_eDraw);
			sprintf(cutstr,"%s && %s && %s && channel==%i",energyCut,HGCut,SDCut,chan);
			skim->Project(h2name,drawstr,cutstr);

			char xTitle[50];
			sprintf(xTitle,"trapENFCal (%d %s %s)",chan,det.c_str(),detEnr[j].c_str());
			h2D_det[j]->GetXaxis()->SetTitle(xTitle);
			h2D_det[j]->GetYaxis()->SetTitle("T / E");
			h2D_det[j]->SetMinimum(1);
			h2D_det[j]->SetMaximum(100);
			h2D_det[j]->Draw("COLZ");
			c2_dets[j]->Update();

			char file2[100];
			sprintf(file2,"output/DS1_Det%i_TEvsE.pdf",chan);
			c2_dets[j]->Print(file2);

			// may need to clear out memory
		}
	}

	// %-------------------------------------% \\

	// T/E plots with TM cut added, detector by detector (18 plots)
	if(0)
	{
		TCanvas *c2_dets[18]; 
		TH2D *h2D_det[18]; 
		char c2name[50];
		char h2name[50];
		for(int j = 0; j < 18; j++)  
		{
			int chan = highGain[j];	// analysis channel
			string det = hgPos[j];	// detector position

			sprintf(c2name,"c1_%d",chan);
			c2_dets[j] = new TCanvas(c2name,"Bob Ross's Canvas",800,600);
			c2_dets[j]->SetLogy();
			c2_dets[j]->SetLogz();

			sprintf(h2name,"TE_chan%d",chan);
			h2D_det[j] = new TH2D(h2name,h2name,750,0,50,1000,0.1,100); 
			sprintf(drawstr,"%s",t_eDraw);
			sprintf(cutstr,"%s && %s && %s && %s && channel==%i",energyCut,HGCut,SDCut,TMCut,chan);
			skim->Project(h2name,drawstr,cutstr);

			char xTitle[50];
			sprintf(xTitle,"trapENFCal (%d %s %s)",chan,det.c_str(),detEnr[j].c_str());
			h2D_det[j]->GetXaxis()->SetTitle(xTitle);
			h2D_det[j]->GetYaxis()->SetTitle("T / E (TM Cut)");
			h2D_det[j]->SetMinimum(1);
			h2D_det[j]->SetMaximum(100);
			h2D_det[j]->Draw("COLZ");
			c2_dets[j]->Update();

			char file2[100];
			sprintf(file2,"output/DS1_Det%i_TETMvsE.pdf",chan);
			c2_dets[j]->Print(file2);

			// may need to clear out memory
		}
	}

	// %-------------------------------------% \\

	// P7D3 : T/E spectra before and after run 12520 (2 plots)
	if(0)
	{
		// before 12520

		TCanvas *c2_10 = new TCanvas("c2_10","Bob Ross's Canvas",850,600);
		c2_10->SetLogy();
		c2_10->SetLogz();

		TH2D *h2D_10 = new TH2D("h2D_10","",750,0,50,1000,0.1,100); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && (channel==594) && run < 12500",HGCut,SDCut);
		skim->Project("h2D_10",drawstr,cutstr);

		char xTitle[50];
		sprintf(xTitle,"trapENFCal (%i %s %s)",highGain[10],hgPos[10].c_str(),detEnr[10].c_str());
		h2D_10->GetXaxis()->SetTitle(xTitle);
		h2D_10->GetYaxis()->SetTitle("T / E : before 12520");
		h2D_10->SetMinimum(1);
		h2D_10->SetMaximum(100);
		h2D_10->Draw("COLZ");

		c2_10->Update();
		c2_10->Print("output/DS1_P7D3_TE_Before12520.pdf");

		// after 12520

		TCanvas *c2_11 = new TCanvas("c2_11","Bob Ross's Canvas",850,600);
		c2_11->SetLogy();
		c2_11->SetLogz();

		TH2D *h2D_11 = new TH2D("h2D_11","",750,0,50,1000,0.1,100); 
		sprintf(drawstr,"%s",t_eDraw);
		sprintf(cutstr,"%s && %s && (channel==594) && run > 12500",HGCut,SDCut);
		skim->Project("h2D_11",drawstr,cutstr);

		sprintf(xTitle,"trapENFCal (%i %s %s)",highGain[10],hgPos[10].c_str(),detEnr[10].c_str());
		h2D_11->GetXaxis()->SetTitle(xTitle);
		h2D_11->GetYaxis()->SetTitle("T / E : after 12520");
		h2D_11->SetMinimum(1);
		h2D_11->SetMaximum(100);
		h2D_11->Draw("COLZ");

		c2_11->Update();
		c2_11->Print("output/DS1_P7D3_TE_After12520.pdf");
	}

}