// TCut-based skim file skimmer.
// Clint Wiseman, USC/Majorana
// 9/19/2016

#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TNamed.h"

using namespace std;

char theCut[1000];

char developCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin<0";

char standardCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

// only difference from the standard cut is that this one goes out to 100 keV
char DNPCut[1000] = "trapENFCal < 100 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

char tuneTECut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

void PrintSkimChannels(TChain *skim)
{
	vector<int> chanList;
	sprintf(theCut,"%s && isGood",developCut);
	int cts = skim->Draw("channel",theCut,"GOFF");
	cout << "found " << cts << " entries\n";
	for (int i = 0; i < cts; i++)
	{
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


int main()
{
	TChain *skim = new TChain("skimTree");
	skim->Add("~/ds1-low/skimDS1_pt8*");
	long entries = skim->GetEntries();
	printf("\n\nFound %li entries.\n",entries);

	sprintf(theCut,"%s",DNPCut);

	cout << "\nUsing this cut: " << theCut << "\n\n";

	skim->Draw(">>entList",theCut,"entrylist GOFF");
	TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
	skim->SetEntryList(entList);
	TFile *f2 = new TFile("./data/DNPCut_DS1.root","recreate");
	TTree *small = skim->CopyTree("");
	small->Write();

	TNamed n("cutUsedHere",theCut);	// roll the cut used into the file.
	n.Write();

	f2->Close();

	// Print active detectors for this cut
	// TChain *cutSkim = new TChain("skimTree");
	// cutSkim->Add("./data/tuneTECut_DS1.root");
	// PrintSkimChannels(cutSkim);

	cout << "\007"; // beep when done

	// G@meC0cks
}