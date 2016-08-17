// temps.cc
// Clint Wiseman, USC/Majorana
// 2/22/16
//
/*
	A few useful environment monitor variables:
	"DR Diff Press,Davis SCM Environmental Monitoring Processes"
	"DR Radon,Davis SCM Environmental Monitoring Processes"		
	"DR Rel. Humidity,Davis SCM Environmental Monitoring Processes"
	"DR RH Rad7,Davis SCM Environmental Monitoring Processes"		
	"DR Temperature,Davis SCM Environmental Monitoring Processes"	
	"DR 0.3 um count,Davis SCM Environmental Monitoring Processes"	
	"DR 0.5 um count,Davis SCM Environmental Monitoring Processes"	
	"DR 0.7 um count,Davis SCM Environmental Monitoring Processes"	
	"DR 1.0 um count,Davis SCM Environmental Monitoring Processes"	
	"DR 2.0 um count,Davis SCM Environmental Monitoring Processes"	
	"DR 5.0 um count,Davis SCM Environmental Monitoring Processes"	
*/

#include "MJSlowControlsDoc.hh"
using namespace std;
using namespace MJDB;

int main(int argc, char** argv)
{
	string start = "2016/02/09 15:00:00"; 
	string stop = "2016/02/09 17:00:00";  
	string zone = "MDT";
	string variable = "DR Temperature,Davis SCM Environmental Monitoring Processes";

	// Create the ROOT file
	MJSlowControlsDoc doc;
	doc.SetDatabase(kHDB_SCM,kFeresaValues,variable,start,stop,zone,"");
	doc.GrabHistory(true);
	doc.RootifyDocument("DRTemperature.root",variable);

	// Store a list of runs (with unix time info) that matches this date range
	MJSlowControlsDoc doc2;
	doc2.SetDatabase(kHDB_DAQ2,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc2.GrabHistory();
	doc2.CreateRunList("RunList.txt",true);	// true: include unix time info
}
