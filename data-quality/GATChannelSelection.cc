#include "GATChannelSelection.hh"

using namespace std;

void GATChannelSelection::SlaveBegin() 
{
	fInputStream.open(fInputTable);
	if(!fInputStream.good()) {
		cout << "Couldn't open " << fInputTable << endl;
		return;
	}
	int chan = 0;
	uint32_t qual = 0;
	while(!fInputStream.eof())
	{
		fInputStream >> chan >> qual;
		channelTable.push_back(chan);
		qualityTable.push_back(qual);
		// printf("chan: %i  qual: %i\n",chan,qual);
	}
	cout << "Channel Quality vectors loaded\n";
	
	ReqEventBranch();
	// ReqRunBranch();
	PublishObj(&fQuality);
}

void GATChannelSelection::Process()
{
	LoadEventBranch();
  	const MGTEvent* event = GetEvent();
  	fQuality.resize(event->GetNDigitizerData());
	for(size_t iDD = 0; iDD < fQuality.size(); iDD++) 
	{
		size_t channel = event->GetDigitizerData(iDD)->GetID();
		
		// loop over channelTable and find match
		bool foundMatch = false;
		for (size_t iCQ=0; iCQ<channelTable.size(); iCQ++)
		{
			if ((int)channel == (int)channelTable[iCQ]) 
			{
				fQuality[iDD] = qualityTable[iCQ];
				foundMatch = true;
				continue;
			}
		}
		if (!foundMatch) fQuality[iDD] = 999;
	}
}

void GATChannelSelection::SlaveTerminate()
{
	channelTable.clear();
	qualityTable.clear();
	fInputStream.close();
}