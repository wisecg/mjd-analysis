#include <iostream>
#include <fstream>
#include "MJSlowControlsDoc.hh"

using namespace std;
using namespace MJDB;

void FindRunDateTimes();

int main(int argc, char* argv[])
{
	MJSlowControlsDoc doc;

	string start = "2016/08/25 00:00:00";
	string stop = "2016/09/15 00:00:00";
	string zone = "MDT";
	doc.SetDatabase(kHDB_DAQ1,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc.GrabHistory();

	ofstream output;
	output.open("DS3-m2runInfo.txt");
	doc.FindRunInfo("DS3-m2runs.txt",output);
	output.close();
}
