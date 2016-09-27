// TCut-based skim file skimmer.
// Clint Wiseman, USC/Majorana
// 9/19/2016

#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TEntryList.h"

using namespace std;

int main()
{
	char standardCut[1000] = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && trapETailMin<0";

	TChain *skim = new TChain("skimTree");
	skim->Add("/Users/wisecg/datasets/ds1-low/*.root");
	long entries = skim->GetEntries();
	printf("\n\nFound %li entries.\n",entries);
	cutTreeStandard(skim);

	skim->Draw(">>entList",standardCut,"entrylist GOFF");
	TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
	skim->SetEntryList(entList);
	TFile *f2 = new TFile("standardCutLE_DS1.root","recreate");
	TTree *small = skim->CopyTree("");
	small->Write();
	f2->Close();

	cout << "\a"; // give me a beep when done
}