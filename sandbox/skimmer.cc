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

using namespace std;

void plotter();

int main(int argc, char** argv)
{
	plotter();	
}

void plotter()
{
	TChain *skim = new TChain("skimTree");
	skim->Add("/Users/wisecg/datasets/ds1/*.root");
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

	
	// binning
	// kris's binning: 15 bins/kev.  Ex: 750,0,50
	// clint's binning: 2 bins/kev.  Ex: 100,0,50
	// kris's second one: 5 bins/kev.  Ex: 250,0,50

	// "Method 1"
	if(0)
	{
		// Loop over the skim chain and make a low-E plot
		TH1D *enfHist = new TH1D("enfHist","name of histo",250,0,50);
		for (int i = 0; i < entries; i++)
		{
			skim->GetEntry(i);

			// event-level
			// printf("run %i  mH %i\n",run,mH);

			double enf;

			// channel-level
			for (int j = 0; j < (int)channel->size(); j++)
			{
				enf = trapENFCal->at(j);

				if (mH == 1 && channel->at(j)%2==0){
					// printf("channel %i  trapENFCal %.2f\n",channel->at(j),enf);
					enfHist->Fill(enf);
				}
			}
		}
		TCanvas *can = new TCanvas("can","Bob Ross's Canvas",800,600);
		can->SetLogy();
		enfHist->GetXaxis()->SetTitle("trapENFCal");
		enfHist->Draw();
		can->Print("SS_HG_0_50_Spec.pdf");	// can say .C, .jpg, but I like .pdf
	}

	
	// "Method 2"

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

	if(1)
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

		TH1D *h3 = new TH1D("h3","",750,0,50);
		sprintf(cutstr,"%s && %s && %s && %s",energyCut,HGCut,SDCut,IsEnr);
		skim->Project("h3", drawstr, cutstr);

		h0->SetLineColor(kBlack);
		h0->GetXaxis()->SetTitle(drawstr);
		h0->GetYaxis()->SetTitle("Counts");
		h0->Draw();

		h1->SetLineColor(kBlue);
		h1->Draw("same");

		h2->SetLineColor(kRed);
		h2->Draw("same");

		h3->SetLineColor(kGreen);
		h3->Draw("same");

		TLegend* leg1 = new TLegend(0.6,0.7,0.87,0.94);
		leg1->AddEntry(h0,"HG","l");
		leg1->AddEntry(h1,"HG & SD","l");
		leg1->AddEntry(h2,"HG & SD & T/E","l");
		leg1->AddEntry(h3,"HG & SD & enr","l");
		leg1->Draw();

		c1_1->Update();
		c1_1->Print("DS1_Spectrum.pdf");
	}


}