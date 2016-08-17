#include "MJSlowControlsDoc.hh"
using namespace std;
using namespace MJDB;

int main(int argc, char** argv)
{
	string start = "2016/02/09 15:00:00"; 
	string stop = "2016/02/09 17:00:00";  
	string zone = "MDT";
	string variable = "";

	// Find runs matching this date
	MJSlowControlsDoc doc;
	doc.SetDatabase(kHDB_DAQ2,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc.GrabHistory();
	doc.CreateRunList("MatchingRuns.txt",false);	// if true, includes time info
	
	//P4D1	B8455
	MJSlowControlsDoc doc1;
	variable = "Baseline Voltage,B8455";
	doc1.SetDatabase(kHDB_DetHist,kFeresaValues,variable,start,stop,zone,kIncDocs);
	doc1.GrabHistory();	
	doc1.RootifyDocument("BaselineInfo.root",variable,"P4D1_B8455");
	
	//P4D5	B8469
	MJSlowControlsDoc doc2;
	variable = "Baseline Voltage,B8469";
	doc2.SetDatabase(kHDB_DetHist,kFeresaValues,variable,start,stop,zone,kIncDocs);
	doc2.GrabHistory();	
	doc2.RootifyDocument("BaselineInfo.root",variable,"P4D5_B8469",2);
}
