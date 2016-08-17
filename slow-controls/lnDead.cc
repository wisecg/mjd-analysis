#include <iostream>
#include <vector>
#include "MJSlowControlsDoc.hh"

using namespace std;
using namespace MJDB;

void GenerateRunDateList();
void CalculateLivetimeReduction();

int main()
{
	// GenerateRunDateList();
	CalculateLivetimeReduction();
}

void CalculateLivetimeReduction()
{
	ifstream in("DS0_RunsAndDates.txt");
	string line;
	int run;
	float duration;
	string bDate;
	string bTime;
	long bUnixtime;
	string eDate;
	string eTime;
	long eUnixtime;
	vector<int> runs;
	vector<long> starts;
	vector<long> stops;
	while(!in.eof())
	{
		in >> run >> duration >> bDate >> bTime >> bUnixtime >> eDate >> eTime >> eUnixtime;
		// cout << run << " " << duration << " " << bDate << " " << bTime << " " << bUnixtime<< " " << eDate << " " << eTime << " " << eUnixtime << endl;
		runs.push_back(run);
		starts.push_back(bUnixtime);
		stops.push_back(eUnixtime);
	}

	vector<double> lnFillTimes;
	int dsNumber = 0;
	if(dsNumber == 0) {
		lnFillTimes.push_back(1435870160);
		lnFillTimes.push_back(1436020533);
		lnFillTimes.push_back(1436168622);
		lnFillTimes.push_back(1436316362);
		lnFillTimes.push_back(1436463827);
		lnFillTimes.push_back(1436614832);
		lnFillTimes.push_back(1436759402);
		lnFillTimes.push_back(1436908966);
		lnFillTimes.push_back(1437058884);
		lnFillTimes.push_back(1437060687);
		lnFillTimes.push_back(1437204895);
		lnFillTimes.push_back(1437350621);
		lnFillTimes.push_back(1437638022);
		lnFillTimes.push_back(1437781652);
		lnFillTimes.push_back(1437926383);
		lnFillTimes.push_back(1438070819);
		lnFillTimes.push_back(1438121812);
		lnFillTimes.push_back(1438212611);
		lnFillTimes.push_back(1438352837);
		lnFillTimes.push_back(1438496004);
		lnFillTimes.push_back(1438639835);
		lnFillTimes.push_back(1438782282);
		lnFillTimes.push_back(1438932169);
		lnFillTimes.push_back(1439080850);
		lnFillTimes.push_back(1439226624);
		lnFillTimes.push_back(1439371591);
		lnFillTimes.push_back(1439514291);
		lnFillTimes.push_back(1439657247);
		lnFillTimes.push_back(1439801649);
		lnFillTimes.push_back(1439944622);
		lnFillTimes.push_back(1440080091);
		lnFillTimes.push_back(1440192295);
		lnFillTimes.push_back(1440335863);
		lnFillTimes.push_back(1440477294);
		lnFillTimes.push_back(1440618411);
		lnFillTimes.push_back(1440759782);
		lnFillTimes.push_back(1440902658);
		lnFillTimes.push_back(1441046730);
		lnFillTimes.push_back(1441187019);
		lnFillTimes.push_back(1441325878);
		lnFillTimes.push_back(1441462000);
		lnFillTimes.push_back(1441600116);
		lnFillTimes.push_back(1441741779);
		lnFillTimes.push_back(1441883940);
		lnFillTimes.push_back(1442027368);
		lnFillTimes.push_back(1442169713);
		lnFillTimes.push_back(1442312599);
		lnFillTimes.push_back(1442453920);
		lnFillTimes.push_back(1442595578);
		lnFillTimes.push_back(1442737259);
		lnFillTimes.push_back(1442879000);
		lnFillTimes.push_back(1443021647);
	}
	if(dsNumber == 1) {
		lnFillTimes.push_back(1452519860);
		lnFillTimes.push_back(1452655867);
		lnFillTimes.push_back(1452791415);
		lnFillTimes.push_back(1453541032);
		lnFillTimes.push_back(1453670314);
		lnFillTimes.push_back(1453800137);
		lnFillTimes.push_back(1453929779);
		lnFillTimes.push_back(1453937178);
		lnFillTimes.push_back(1453998125);
		lnFillTimes.push_back(1454000591);
		lnFillTimes.push_back(1454002456);
		lnFillTimes.push_back(1454014971);
		lnFillTimes.push_back(1454107981);
		lnFillTimes.push_back(1454219392);
		lnFillTimes.push_back(1454332307);
		lnFillTimes.push_back(1454447212);
		lnFillTimes.push_back(1454559953);
		lnFillTimes.push_back(1454679496);
		lnFillTimes.push_back(1454769079);
		lnFillTimes.push_back(1454882301);
		lnFillTimes.push_back(1454946492);
		lnFillTimes.push_back(1454951907);
		lnFillTimes.push_back(1454954012);
		lnFillTimes.push_back(1454955958);
		lnFillTimes.push_back(1455225161);
		lnFillTimes.push_back(1455228289);
		lnFillTimes.push_back(1455440690);
		lnFillTimes.push_back(1455568585);
		lnFillTimes.push_back(1455696789);
		lnFillTimes.push_back(1455822149);
		lnFillTimes.push_back(1455952573);
		lnFillTimes.push_back(1456082166);
		lnFillTimes.push_back(1456206057);
		lnFillTimes.push_back(1456333236);
		lnFillTimes.push_back(1456460237);
		lnFillTimes.push_back(1456588495);
		lnFillTimes.push_back(1456717776);
		lnFillTimes.push_back(1456846882);
		lnFillTimes.push_back(1456943983);
		lnFillTimes.push_back(1456981934);
		lnFillTimes.push_back(1457110918);
		lnFillTimes.push_back(1457238098);
		lnFillTimes.push_back(1457365179);
		lnFillTimes.push_back(1457491997);
		lnFillTimes.push_back(1457619662);
		lnFillTimes.push_back(1457747884);
		lnFillTimes.push_back(1457874684);
		lnFillTimes.push_back(1458001143);
		lnFillTimes.push_back(1458130495);
		lnFillTimes.push_back(1458259402);
		lnFillTimes.push_back(1458387930);
		lnFillTimes.push_back(1458515657);
		lnFillTimes.push_back(1458639685);
		lnFillTimes.push_back(1458767518);
		lnFillTimes.push_back(1458896260);
		lnFillTimes.push_back(1459023862);
		lnFillTimes.push_back(1459150863);
		lnFillTimes.push_back(1459275989);
		lnFillTimes.push_back(1459402548);
		lnFillTimes.push_back(1459529663);
		lnFillTimes.push_back(1459654930);
		lnFillTimes.push_back(1459785212);
		lnFillTimes.push_back(1459912507);
		lnFillTimes.push_back(1460042708);
		lnFillTimes.push_back(1460169429);
		lnFillTimes.push_back(1460297779);
		lnFillTimes.push_back(1460426527);
		lnFillTimes.push_back(1460552742);
		lnFillTimes.push_back(1460678641);
		lnFillTimes.push_back(1460808916);
		lnFillTimes.push_back(1460939150);
		lnFillTimes.push_back(1461064619);
	}

	// loop over runs
	double TotalMinutesFilling = 0;

	for (int i = 0; i < (int)runs.size(); i++)
	{
		int run = runs[i];
		long start = starts[i];
		long stop = stops[i];

		for (int j = 0; j < (int)lnFillTimes.size(); j++)
		{
			long fill = lnFillTimes[j];
			long loFill = fill - 15 * 60;
			long hiFill = fill + 5 * 60;

			// Fill time is during run
			if (fill >= start && fill <= stop)
			{
				// upper and lower are both within run
				if (loFill >= start && hiFill <= stop)	{
					printf("Run %i  Fill %i  1: hiFill - loFill = %.2f\n",run,j,(double)(hiFill-loFill)/60);
					TotalMinutesFilling += (double)(hiFill-loFill)/60;
				}
				// upper is within run
				else if (hiFill >= start && loFill <= start)  {
					printf("Run %i  Fill %i  2: hiFill - start = %.2f\n",run,j,(double)(hiFill-start)/60);
					TotalMinutesFilling += (double)(hiFill-start)/60;
				}
				// lower is within run
				else if (loFill >= start && hiFill >= stop)	{
					printf("Run %i  Fill %i  3: stop - loFill = %.2f\n",run,j,(double)(stop-loFill)/60);
					TotalMinutesFilling += (double)(stop-loFill)/60;
				}
			}
			// Only hiFill boundary is during run
			else if (hiFill >= start && hiFill <= stop)
			{
				printf("Run %i  Fill %i  4: hiFill - start = %.2f\n",run,j,(double)(hiFill-start)/60);
				TotalMinutesFilling += (double)(hiFill-start)/60;
			}
			// Only loFill boundary is during run
			else if (loFill >= start && loFill <= stop)
			{
				printf("Run %i  Fill %i  5: stop - loFill = %.2f\n",run,j,(double)(stop-loFill)/60);
				TotalMinutesFilling += (double)(stop-loFill)/60;
			}
		}
	}

	cout << "MAKE SURE THERE ARE NO REPEAT RUNS IN THE OUTPUT!\n";
	cout << "Total minutes filling : " << TotalMinutesFilling << endl;
}

void GenerateRunDateList()
{
	vector<int> runList;
	// DS1 run list.
	if(0){
		for(int i = 9422; i <=9440; i++) runList.push_back(i);
		for(int i = 9455; i <=9455; i++) runList.push_back(i);
		for(int i = 9471; i <=9487; i++) runList.push_back(i);
		for(int i = 9492; i <=9492; i++) runList.push_back(i);
		for(int i = 9536; i <=9565; i++) runList.push_back(i);
		for(int i = 9638; i <=9648; i++) runList.push_back(i);
		for(int i = 9650; i <=9668; i++) runList.push_back(i);
		for(int i = 9674; i <=9676; i++) runList.push_back(i);
		for(int i = 9678; i <=9678; i++) runList.push_back(i);
		for(int i = 9711; i <=9727; i++) runList.push_back(i);
		for(int i = 9763; i <=9780; i++) runList.push_back(i);
		for(int i = 9815; i <=9821; i++) runList.push_back(i);
		for(int i = 9823; i <=9832; i++) runList.push_back(i);
		for(int i = 9848; i <=9849; i++) runList.push_back(i);
		for(int i = 9851; i <=9854; i++) runList.push_back(i);
		for(int i = 9856; i <=9912; i++) runList.push_back(i);
		for(int i = 9928; i <=9928; i++) runList.push_back(i);
		for(int i = 9952; i <=9966; i++) runList.push_back(i);
		for(int i = 10019; i <=10035; i++) runList.push_back(i);
		for(int i = 10074; i <=10090; i++) runList.push_back(i);
		for(int i = 10114; i <=10125; i++) runList.push_back(i);
		for(int i = 10129; i <=10171; i++) runList.push_back(i);
		for(int i = 10173; i <=10231; i++) runList.push_back(i);
		for(int i = 10262; i <=10278; i++) runList.push_back(i);
		for(int i = 10298; i <=10299; i++) runList.push_back(i);
		for(int i = 10301; i <=10301; i++) runList.push_back(i);
		for(int i = 10304; i <=10308; i++) runList.push_back(i);
		for(int i = 10312; i <=10342; i++) runList.push_back(i);
		for(int i = 10344; i <=10350; i++) runList.push_back(i);
		for(int i = 10378; i <=10394; i++) runList.push_back(i);
		for(int i = 10552; i <=10558; i++) runList.push_back(i);
		for(int i = 10608; i <=10648; i++) runList.push_back(i);
		for(int i = 10651; i <=10677; i++) runList.push_back(i);
		for(int i = 10679; i <=10717; i++) runList.push_back(i);
		for(int i = 10745; i <=10761; i++) runList.push_back(i);
		for(int i = 10788; i <=10803; i++) runList.push_back(i);
		for(int i = 10830; i <=10845; i++) runList.push_back(i);
		for(int i = 10963; i <=10976; i++) runList.push_back(i);
		for(int i = 11002; i <=11008; i++) runList.push_back(i);
		for(int i = 11010; i <=11019; i++) runList.push_back(i);
		for(int i = 11046; i <=11066; i++) runList.push_back(i);
		for(int i = 11083; i <=11200; i++) runList.push_back(i);
		for(int i = 11350; i <=11350; i++) runList.push_back(i);
		for(int i = 11403; i <=11410; i++) runList.push_back(i);
		for(int i = 11414; i <=11417; i++) runList.push_back(i);
		for(int i = 11419; i <=11426; i++) runList.push_back(i);
		for(int i = 11428; i <=11432; i++) runList.push_back(i);
		for(int i = 11434; i <=11444; i++) runList.push_back(i);
		for(int i = 11446; i <=11451; i++) runList.push_back(i);
		for(int i = 11453; i <=11453; i++) runList.push_back(i);
		for(int i = 11455; i <=11458; i++) runList.push_back(i);
		for(int i = 11466; i <=11483; i++) runList.push_back(i);
		for(int i = 12445; i <=12445; i++) runList.push_back(i);
		for(int i = 12466; i <=12467; i++) runList.push_back(i);
		for(int i = 12477; i <=12483; i++) runList.push_back(i);
		for(int i = 12486; i <=12493; i++) runList.push_back(i);
		for(int i = 12520; i <=12580; i++) runList.push_back(i);
		for(int i = 12607; i <=12625; i++) runList.push_back(i);
		for(int i = 12636; i <=12647; i++) runList.push_back(i);
		for(int i = 12652; i <=12653; i++) runList.push_back(i);
		for(int i = 12664; i <=12675; i++) runList.push_back(i);
		for(int i = 12677; i <=12724; i++) runList.push_back(i);
		for(int i = 12735; i <=12798; i++) runList.push_back(i);
		for(int i = 12816; i <=12816; i++) runList.push_back(i);
		for(int i = 12818; i <=12819; i++) runList.push_back(i);
		for(int i = 12821; i <=12821; i++) runList.push_back(i);
		for(int i = 12823; i <=12824; i++) runList.push_back(i);
		for(int i = 12826; i <=12831; i++) runList.push_back(i);
		for(int i = 12833; i <=12838; i++) runList.push_back(i);
		for(int i = 12842; i <=12861; i++) runList.push_back(i);
		for(int i = 12875; i <=12875; i++) runList.push_back(i);
		for(int i = 13000; i <=13053; i++) runList.push_back(i);
		for(int i = 13055; i <=13056; i++) runList.push_back(i);
		for(int i = 13066; i <=13074; i++) runList.push_back(i);
		for(int i = 13076; i <=13092; i++) runList.push_back(i);
		for(int i = 13094; i <=13096; i++) runList.push_back(i);
		for(int i = 13099; i <=13115; i++) runList.push_back(i);
		for(int i = 13117; i <=13119; i++) runList.push_back(i);
		for(int i = 13123; i <=13137; i++) runList.push_back(i);
		for(int i = 13148; i <=13150; i++) runList.push_back(i);
		for(int i = 13153; i <=13156; i++) runList.push_back(i);
		for(int i = 13186; i <=13189; i++) runList.push_back(i);
		for(int i = 13191; i <=13287; i++) runList.push_back(i);
		for(int i = 13304; i <=13304; i++) runList.push_back(i);
		for(int i = 13306; i <=13311; i++) runList.push_back(i);
		for(int i = 13313; i <=13350; i++) runList.push_back(i);
		for(int i = 13362; i <=13368; i++) runList.push_back(i);
	}

	// DS0 run list.
	if(1){
		for(int i = 2580; i <=2580; i++)runList.push_back(i);
		for(int i = 2582; i <=2612; i++)runList.push_back(i);
		for(int i = 2614; i <=2629; i++)runList.push_back(i);
		for(int i = 2644; i <=2649; i++)runList.push_back(i);
		for(int i = 2658; i <=2673; i++)runList.push_back(i);
		for(int i = 2689; i <=2715; i++)runList.push_back(i);
		for(int i = 2717; i <=2907; i++)runList.push_back(i);
		for(int i = 2909; i <=2920; i++)runList.push_back(i);
		for(int i = 3137; i <=3271; i++)runList.push_back(i);
		for(int i = 3293; i <=3370; i++)runList.push_back(i);
		for(int i = 3372; i <=3462; i++)runList.push_back(i);
		for(int i = 3464; i <=3580; i++)runList.push_back(i);
		for(int i = 3596; i <=3645; i++)runList.push_back(i);
		for(int i = 4034; i <=4035; i++)runList.push_back(i);
		for(int i = 4038; i <=4040; i++)runList.push_back(i);
		for(int i = 4045; i <=4134; i++)runList.push_back(i);
		for(int i = 4239; i <=4245; i++)runList.push_back(i);
		for(int i = 4248; i <=4254; i++)runList.push_back(i);
		for(int i = 4256; i <=4268; i++)runList.push_back(i);
		for(int i = 4270; i <=4271; i++)runList.push_back(i);
		for(int i = 4273; i <=4283; i++)runList.push_back(i);
		for(int i = 4285; i <=4311; i++)runList.push_back(i);
		for(int i = 4313; i <=4318; i++)runList.push_back(i);
		for(int i = 4320; i <=4320; i++)runList.push_back(i);
		for(int i = 4322; i <=4326; i++)runList.push_back(i);
		for(int i = 4328; i <=4336; i++)runList.push_back(i);
		for(int i = 4338; i <=4361; i++)runList.push_back(i);
		for(int i = 4363; i <=4382; i++)runList.push_back(i);
		for(int i = 4384; i <=4401; i++)runList.push_back(i);
		for(int i = 4403; i <=4428; i++)runList.push_back(i);
		for(int i = 4436; i <=4454; i++)runList.push_back(i);
		for(int i = 4457; i <=4489; i++)runList.push_back(i);
		for(int i = 4491; i <=4493; i++)runList.push_back(i);
		for(int i = 4497; i <=4503; i++)runList.push_back(i);
		for(int i = 4505; i <=4518; i++)runList.push_back(i);
		for(int i = 4549; i <=4654; i++)runList.push_back(i);
		for(int i = 4655; i <=4777; i++)runList.push_back(i);
		for(int i = 4789; i <=4797; i++)runList.push_back(i);
		for(int i = 4800; i <=4831; i++)runList.push_back(i);
		for(int i = 4854; i <=4872; i++)runList.push_back(i);
		for(int i = 4874; i <=4883; i++)runList.push_back(i);
		for(int i = 4885; i <=4907; i++)runList.push_back(i);
		for(int i = 4938; i <=4960; i++)runList.push_back(i);
		for(int i = 4962; i <=4968; i++)runList.push_back(i);
		for(int i = 4970; i <=4980; i++)runList.push_back(i);
		for(int i = 5007; i <=5038; i++)runList.push_back(i);
		for(int i = 5040; i <=5061; i++)runList.push_back(i);
		for(int i = 5090; i <=5118; i++)runList.push_back(i);
		for(int i = 5125; i <=5252; i++)runList.push_back(i);
		for(int i = 5277; i <=5330; i++)runList.push_back(i);
		for(int i = 5372; i <=5393; i++)runList.push_back(i);
		for(int i = 5405; i <=5414; i++)runList.push_back(i);
		for(int i = 5449; i <=5501; i++)runList.push_back(i);
		for(int i = 5525; i <=5527; i++)runList.push_back(i);
		for(int i = 5531; i <=5534; i++)runList.push_back(i);
		for(int i = 5555; i <=5589; i++)runList.push_back(i);
		for(int i = 5591; i <=5608; i++)runList.push_back(i);
		for(int i = 5610; i <=5751; i++)runList.push_back(i);
		for(int i = 5753; i <=5764; i++)runList.push_back(i);
		for(int i = 5766; i <=5822; i++)runList.push_back(i);
		for(int i = 5826; i <=5850; i++)runList.push_back(i);
		for(int i = 5889; i <=5890; i++)runList.push_back(i);
		for(int i = 5894; i <=5902; i++)runList.push_back(i);
		for(int i = 6553; i <=6577; i++)runList.push_back(i);
		for(int i = 6775; i <=6775; i++)runList.push_back(i);
		for(int i = 6776; i <=6809; i++)runList.push_back(i);
		for(int i = 6811; i <=6830; i++)runList.push_back(i);
		for(int i = 6834; i <=6853; i++)runList.push_back(i);
		for(int i = 6887; i <=6903; i++)runList.push_back(i);
		for(int i = 6957; i <=6963; i++)runList.push_back(i);
	}

	// ofstream out("DS1_RunList.txt");
	ofstream out("DS0_RunList.txt");
	for (int j = 0; j < (int)runList.size(); j++)
		out << runList[j] << endl;
	out.close();

	MJSlowControlsDoc doc;
	string start = "2015/06/29 00:00:00"; 	// 2580 - DS0 starts (2015-6-30-P3JDY_Run2580)
	string stop = "2015/09/24 00:00:00";	// 6963 - DS0 stops (2015-9-23-P3JDY_Run6963)
	// string start = "2016/01/10 00:00:00";  // 9422 - DS1 starts (2016-1-12-P3KJR_Run9422)
	// string stop = "2016/04/15 00:00:00";	// 13368 - DS1 stops (2016-4-14-P3KJR_Run13368)
	string zone = "MDT";
	doc.SetDatabase(kHDB_DAQ2,kFeresaRunInfo,"",start,stop,zone,kIncDocs);
	doc.GrabHistory();

	// Store the start/stop dates of every run.
	// ofstream out2("DS1_RunsAndDates.txt");
	// doc.FindRunInfo("DS1_RunList.txt",out2);

	ofstream out2("DS0_RunsAndDates.txt");
	doc.FindRunInfo("DS0_RunList.txt",out2);
}