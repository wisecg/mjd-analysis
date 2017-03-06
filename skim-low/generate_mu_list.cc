#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TCut.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TEntryList.h"
#include "TROOT.h"

#include "GATDataSet.hh"
#include "MJTChannelMap.hh"
#include "TClonesArray.h"
#include "MGDOUtils.hh"

using namespace std;
using namespace CLHEP;

void LoadDS0(GATDataSet& ds);
void LoadDS1(GATDataSet& ds);
void LoadDS3(GATDataSet& ds);

int main(int argc, char** argv)
{
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [DS number]" << endl;
    return 1;
  }
  int dsNumber = atoi(argv[1]);
  cout << "generating mu list for DS" << dsNumber << endl;
  GATDataSet ds;
  string filename;
  string fileOutName;
  if(dsNumber == 0) {
    LoadDS0(ds);
    filename = "MuonList_DS0.txt";
    fileOutName = "MuonList_DS0_aug.txt";
  }
  else if(dsNumber == 1) {
    LoadDS1(ds);
    filename = "MuonList_DS1.txt";
    fileOutName = "MuonList_DS1_aug.txt";
  }
  else if(dsNumber == 3) {
    LoadDS3(ds);
    filename = "MuonList_DS3_Aug.txt";
    fileOutName = "MuonList_DS3_Aug_aug.txt";
  }
  else {
    cout << "Unknown dataset number " << dsNumber << endl;
    return 1;
  }
  TChain* gatChain = ds.GetGatifiedChain(false);

  cout << "Loading muon data..." << endl;
  vector<int> muRuns;
  vector<double> muTimes_s;
  vector<int> muTypes;
  vector<bool> badScalers;
  ifstream muFile(filename.c_str());
  if(!muFile.good()) {
    cout << "Couldn't open " << filename << endl;
    return 0;
  }
  int rundummy, typedummy;
  double timedummy;
  bool baddummy;
  double dummydummy;
  if(dsNumber == 0) {
    while(muFile >> rundummy >> timedummy >> typedummy >> baddummy) {
      muRuns.push_back(rundummy);
      muTimes_s.push_back(timedummy);
      muTypes.push_back(typedummy);
      badScalers.push_back(baddummy);
    }
  }
  else if(dsNumber == 1) {
    while(muFile >> rundummy >> dummydummy >> timedummy >> typedummy >> baddummy) {
      muRuns.push_back(rundummy);
      muTimes_s.push_back(timedummy);
      muTypes.push_back(typedummy);
      badScalers.push_back(baddummy);
    }
  }
  else if(dsNumber == 3) {
    while(muFile >> rundummy >> dummydummy >> timedummy >> typedummy >> baddummy) {
      muRuns.push_back(rundummy);
      muTimes_s.push_back(timedummy);
      muTypes.push_back(typedummy);
      badScalers.push_back(baddummy);
    }
  }
  else {
    cout << "Shouldn't get here! dsNumber = " << dsNumber << endl;
    return 1;
  }
  muFile.close();
  if(muTimes_s.size() == 0) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }

  cout << "Loading run start times..." << endl;
  int nRuns = gatChain->Draw("run:startTime:stopTime", "LocalEntry$ == 0", "GOFF");
  map<int, double> runStartTimes;
  vector<int> gapRuns;
  double lastStop = 0.0;
  for(int iRun = 0; iRun < nRuns; iRun++) {
    int run = gatChain->GetV1()[iRun];
    double tStart = gatChain->GetV2()[iRun];
    runStartTimes[run] = tStart;
    // mark >10 second intervals between runs as gaps (mu type = 3).
    if(tStart - lastStop > 10) gapRuns.push_back(run);
    lastStop = gatChain->GetV3()[iRun];
  }
  if(gapRuns.size() == 0) {
    cout << "should have gotten at least one gap, on first run" << endl;
    return 1;
  }

  ofstream muOutFile(fileOutName);
  if(!muOutFile.good()) {
    cout << "Couldn't open " << fileOutName << endl;
    return 0;
  }
  muOutFile.precision(12);
  size_t iGap = 0;
  for(size_t i=0; i<muTimes_s.size(); i++) {
    while(iGap < gapRuns.size() && gapRuns[iGap] <= muRuns[i]) {
      muOutFile << gapRuns[iGap] << ' ' << runStartTimes[gapRuns[iGap]]
                << " 0.0 3 0" << endl;
      iGap++;
    }
    muOutFile << muRuns[i] << ' ' << runStartTimes[muRuns[i]] << ' '
              << muTimes_s[i] << ' ' << muTypes[i] << ' '
              << badScalers[i] << endl;
  }
  muOutFile.close();

  return 0;
}

void LoadDS0(GATDataSet& ds)
{
  ds.AddRunRange(2580,2580); //1
  ds.AddRunRange(2582,2612); //31
  ds.AddRunRange(2614,2629); //16
  ds.AddRunRange(2644,2649);//6
  ds.AddRunRange(2658,2673);//16
  ds.AddRunRange(2689,2715);//28
  ds.AddRunRange(2717,2750);//33
  ds.AddRunRange(2751,2784);//33
  ds.AddRunRange(2785,2820);//35
  ds.AddRunRange(2821,2855);//34
  ds.AddRunRange(2856,2890);//34
  ds.AddRunRange(2891,2907);//17
  ds.AddRunRange(2909,2920);//12
  ds.AddRunRange(3137,3166);//30
  ds.AddRunRange(3167,3196);//30
  ds.AddRunRange(3197,3226);//30
  ds.AddRunRange(3227,3256);//30
  ds.AddRunRange(3257,3271);//14
  ds.AddRunRange(3293,3310);//18
  ds.AddRunRange(3311,3340);//30
  ds.AddRunRange(3341,3370);//30
  ds.AddRunRange(3371,3400);//30
  ds.AddRunRange(3401,3432);//32
  ds.AddRunRange(3461,3462);//2
  ds.AddRunRange(3464,3500);//36
  ds.AddRunRange(3501,3530);//30
  ds.AddRunRange(3531,3560);//30
  ds.AddRunRange(3561,3580);//20
  ds.AddRunRange(3596,3610);//15
  ds.AddRunRange(3611,3645);//50
  ds.AddRunRange(4034,4035);//2
  ds.AddRunRange(4038,4040);//3
  ds.AddRunRange(4045,4074);//30
  ds.AddRunRange(4075,4104);//30
  ds.AddRunRange(4105,4134);//30
  ds.AddRunRange(4239,4245);//7
  ds.AddRunRange(4248,4254);//7
  ds.AddRunRange(4256,4268);//13
  ds.AddRunRange(4270,4271);//2
  ds.AddRunRange(4273,4283);//11
  ds.AddRunRange(4285,4311);//27
  ds.AddRunRange(4313,4318);//6
  ds.AddRunRange(4320,4320);//1
  ds.AddRunRange(4322,4326);//5
  ds.AddRunRange(4328,4336);//9
  ds.AddRunRange(4338,4361);//24
  ds.AddRunRange(4363,4382);//20
  ds.AddRunRange(4384,4401);//18
  ds.AddRunRange(4403,4428);//26
  ds.AddRunRange(4436,4454);//19
  ds.AddRunRange(4457,4489);//33
  ds.AddRunRange(4491,4493);//3
  ds.AddRunRange(4497,4503);//7
  ds.AddRunRange(4505,4518);//14
  ds.AddRunRange(4549,4590);//32
  ds.AddRunRange(4591,4624);//35
  ds.AddRunRange(4625,4654);//30
  ds.AddRunRange(4655,4684);//30
  ds.AddRunRange(4685,4714);//30
  ds.AddRunRange(4715,4744);//30
  ds.AddRunRange(4745,4777);//33
  ds.AddRunRange(4789,4797);//10
  ds.AddRunRange(4800,4831);//32
  ds.AddRunRange(4854,4872);//19
  ds.AddRunRange(4874,4883);//10
  ds.AddRunRange(4885,4907);//23
  ds.AddRunRange(4938,4960);//23
  ds.AddRunRange(4962,4968);//7
  ds.AddRunRange(4970,4980);//11
  ds.AddRunRange(5007,5038);//32
  ds.AddRunRange(5040,5061);//22
  ds.AddRunRange(5090,5118);//29
  ds.AddRunRange(5125,5154);//30
  ds.AddRunRange(5155,5184);//30
  ds.AddRunRange(5185,5224);//40
  ds.AddRunRange(5225,5252);//28
  ds.AddRunRange(5277,5300);//34
  ds.AddRunRange(5301,5330);//30
  ds.AddRunRange(5372,5393);//22
  ds.AddRunRange(5405,5414);//10
  ds.AddRunRange(5449,5479);//30
  ds.AddRunRange(5480,5501);//23
  ds.AddRunRange(5525,5527);//3
  ds.AddRunRange(5531,5534);//4
  ds.AddRunRange(5555,5589);//35
  ds.AddRunRange(5591,5608);//18
  ds.AddRunRange(5610,5639);//30
  ds.AddRunRange(5640,5669);//30
  ds.AddRunRange(5670,5699);//30
  ds.AddRunRange(5700,5729);//30
  ds.AddRunRange(5730,5751);//22
  ds.AddRunRange(5753,5764);//12
  ds.AddRunRange(5766,5795);//30
  ds.AddRunRange(5796,5822);//27
  ds.AddRunRange(5826,5850);//24
  ds.AddRunRange(5889,5890);//2
  ds.AddRunRange(5894,5902);//9
  ds.AddRunRange(6553,6577);//25
  ds.AddRunRange(6775,6775);//1
  ds.AddRunRange(6776,6809);//34
  ds.AddRunRange(6811,6830);//20
  ds.AddRunRange(6834,6853);//20
  ds.AddRunRange(6887,6903); //17
  ds.AddRunRange(6957,6963);//7
}

void LoadDS1(GATDataSet& ds)
{
  ds.AddRunRange(9422, 9440); // 18
  ds.AddRunRange(9471, 9487); // 16
  ds.AddRunRange(9492, 9492); // 0
  ds.AddRunRange(9536, 9565); // 29
  ds.AddRunRange(9638, 9648); // 10
  ds.AddRunRange(9650, 9668); // 18
  ds.AddRunRange(9674, 9676); // 2
  ds.AddRunRange(9678, 9678); // 0
  ds.AddRunRange(9711, 9727); // 16
  ds.AddRunRange(9763, 9780); // 17
  ds.AddRunRange(9815, 9821); // 6
  ds.AddRunRange(9823, 9832); // 9
  ds.AddRunRange(9848, 9849); // 1
  ds.AddRunRange(9851, 9854); // 3
  ds.AddRunRange(9856, 9912); // 56
  ds.AddRunRange(9928, 9928); // 0
  ds.AddRunRange(9952, 9966); // 14
  ds.AddRunRange(10019, 10035); // 16
  ds.AddRunRange(10074, 10090); // 16
  ds.AddRunRange(10114, 10125); // 11
  ds.AddRunRange(10129, 10149); // 20
  ds.AddRunRange(10150, 10171); // 21
  ds.AddRunRange(10173, 10203); // 30
  ds.AddRunRange(10204, 10231); // 27
  ds.AddRunRange(10262, 10278); // 16
  ds.AddRunRange(10298, 10299); // 1
  ds.AddRunRange(10301, 10301); // 0
  ds.AddRunRange(10304, 10308); // 4
  ds.AddRunRange(10312, 10342); // 30
  ds.AddRunRange(10344, 10350); // 6
  ds.AddRunRange(10378, 10394); // 16
  ds.AddRunRange(10552, 10558); // 6
  ds.AddRunRange(10608, 10648); // 40
  ds.AddRunRange(10651, 10677); // 26
  ds.AddRunRange(10679, 10717); // 38
  ds.AddRunRange(10745, 10761); // 16
  ds.AddRunRange(10788, 10803); // 15
  ds.AddRunRange(10830, 10845); // 15
  ds.AddRunRange(10963, 10976); // 13
  ds.AddRunRange(11002, 11008); // 6
  ds.AddRunRange(11010, 11019); // 9
  ds.AddRunRange(11046, 11066); // 20
  ds.AddRunRange(11083, 11113); // 30
  ds.AddRunRange(11114, 11144); // 30
  ds.AddRunRange(11145, 11175); // 30
  ds.AddRunRange(11176, 11200); // 24
  ds.AddRunRange(11350, 11350); // 0
  ds.AddRunRange(11403, 11410); // 7
  ds.AddRunRange(11414, 11417); // 3
  ds.AddRunRange(11419, 11426); // 7
  ds.AddRunRange(11428, 11432); // 4
  ds.AddRunRange(11434, 11444); // 10
  ds.AddRunRange(11446, 11451); // 5
  ds.AddRunRange(11453, 11453); // 0
  ds.AddRunRange(11455, 11458); // 3
  ds.AddRunRange(11466, 11476); // 10
  ds.AddRunRange(11477, 11483); // 6
  ds.AddRunRange(12445, 12445); // 0
  ds.AddRunRange(12466, 12467); // 1
  ds.AddRunRange(12477, 12483); // 6
  ds.AddRunRange(12486, 12493); // 7
  ds.AddRunRange(12520, 12550); // 30
  ds.AddRunRange(12551, 12580); // 60
  ds.AddRunRange(12607, 12625); // 18
  ds.AddRunRange(12636, 12647); // 11
  ds.AddRunRange(12652, 12653); // 1
  ds.AddRunRange(12664, 12675); // 11
  ds.AddRunRange(12677, 12724); // 47
  ds.AddRunRange(12735, 12765); // 33
  ds.AddRunRange(12766, 12798); // 32
  ds.AddRunRange(12816, 12816); // 0
  ds.AddRunRange(12818, 12819); // 1
  ds.AddRunRange(12821, 12821); // 0
  ds.AddRunRange(12823, 12824); // 1
  ds.AddRunRange(12826, 12831); // 5
  ds.AddRunRange(12833, 12838); // 6
  ds.AddRunRange(12842, 12842); // 0
  ds.AddRunRange(12843, 12861); // 18
  ds.AddRunRange(12875, 12875); // 0
  ds.AddRunRange(13000, 13028); // 28
  ds.AddRunRange(13029, 13053); // 24
  ds.AddRunRange(13055, 13056); // 1
  ds.AddRunRange(13066, 13074); // 8
  ds.AddRunRange(13076, 13092); // 16
  ds.AddRunRange(13094, 13096); // 2
  ds.AddRunRange(13099, 13115); // 16
  ds.AddRunRange(13117, 13119); // 2
  ds.AddRunRange(13123, 13137); // 14
  ds.AddRunRange(13148, 13150); // 2
  ds.AddRunRange(13153, 13156); // 3
  ds.AddRunRange(13186, 13189); // 3
  ds.AddRunRange(13191, 13211); // 20
  ds.AddRunRange(13212, 13242); // 30
  ds.AddRunRange(13243, 13275); // 32
  ds.AddRunRange(13276, 13287); // 11
  ds.AddRunRange(13304, 13304); // 0
  ds.AddRunRange(13306, 13325); // 19
  ds.AddRunRange(13326, 13350); // 24
  ds.AddRunRange(13362, 13368); // 6
  ds.AddRunRange(13369,13383); // 15
  ds.AddRunRange(13395,13411); // 17
  ds.AddRunRange(13519,13548); // 30
  ds.AddRunRange(13572,13573); // 2
  ds.AddRunRange(13667,13688); // 22
  ds.AddRunRange(13699,13704); // 6
  ds.AddRunRange(13715,13719); // 5
  ds.AddRunRange(14010,14040); // 31
  ds.AddRunRange(14041,14041); // 1
  ds.AddRunRange(14342,14372); // 31
  ds.AddRunRange(14386,14387); // 2
}

void LoadDS3(GATDataSet& ds)
{
	ds.AddRunRange(9422, 9440); // 18
	ds.AddRunRange(16797,16826); // 31
	ds.AddRunRange(16827,16835); // 9
	ds.AddRunRange(16857,16886); // 31
	ds.AddRunRange(16887,16910); // 24
	ds.AddRunRange(16931,16936); // 7
	ds.AddRunRange(16947,16952); // 6
	ds.AddRunRange(16957,16959); // 3
	ds.AddRunRange(16970,16999); // 31
	ds.AddRunRange(17000,17009); // 10
	ds.AddRunRange(17035,17057); // 24
	ds.AddRunRange(17060,17090); // 31
	ds.AddRunRange(17091,17121); // 31
	ds.AddRunRange(17122,17127); // 6
	ds.AddRunRange(17129,17131); // 3
	ds.AddRunRange(17138,17156); // 20
	ds.AddRunRange(17158,17181); // 24
	ds.AddRunRange(17304,17318); // 15
	ds.AddRunRange(17322,17343); // 22
	ds.AddRunRange(17351,17381); // 31
	ds.AddRunRange(17382,17412); // 31
	ds.AddRunRange(17413,17422); // 10
	ds.AddRunRange(17448,17477); // 31
	ds.AddRunRange(17478,17493); // 16
	ds.AddRunRange(17500,17519); // 21
	ds.AddRunRange(17531,17553); // 23
	ds.AddRunRange(17555,17559); // 5
	ds.AddRunRange(17567,17597); // 31
	ds.AddRunRange(17598,17628); // 31
	ds.AddRunRange(17629,17659); // 31
	ds.AddRunRange(17660,17686); // 27
	ds.AddRunRange(17703,17717); // 15
	ds.AddRunRange(17720,17721); // 2
	ds.AddRunRange(17852,17882); // 31
	ds.AddRunRange(17883,17913); // 31
	ds.AddRunRange(17914,17944); // 31
	ds.AddRunRange(17945,17948); // 4
	ds.AddRunRange(17967,17980); // 20
}