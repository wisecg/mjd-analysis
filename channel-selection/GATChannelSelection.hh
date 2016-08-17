// To implement, add these three lines to process_mjd_data:
// 
// GATChannelSelection channelSel("quality","inputFile.txt");
// selector.AddInput(&channelSel);
// mjdDataWriter.AddPostedVector(channelSel.GetNameOfPostedVector());

#ifndef GATChannelSelection_hh
#define GATChannelSelection_hh

#include "GATMGTEventProcBase.hh"

class GATChannelQuality : public MGTDataObject, public std::vector<uint32_t>
{
  public:
    GATChannelQuality (const char* aName = "quality", const char* aTitle = "quality") : 
      MGTDataObject(aName, aTitle) {}
    virtual ~GATChannelQuality() {}
};

class GATChannelSelection : public GATMGTEventProcBase 
{
  public:
    GATChannelSelection(const char* outputName, const char* inputTable) 
    : GATMGTEventProcBase("GATChannelSelection","channelreader"), fQuality() {
      fOutputName = outputName;
      fInputTable = inputTable;
    }

    const char* GetNameOfPostedVector() { return fOutputName; }

    virtual ~GATChannelSelection() {}
   
  protected:
    virtual void SlaveBegin(); 
    virtual void Process();
    virtual void SlaveTerminate();
    
    ifstream fInputStream;
    GATChannelQuality fQuality;
    std::vector<double> channelTable;
    std::vector<uint32_t> qualityTable;
    const char* fOutputName;
    const char* fInputTable;
};

#endif   