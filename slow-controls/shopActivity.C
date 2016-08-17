void shopActivity() { 

	ifstream InputList;
	InputList.open("RunList.txt");
	int run = 0, duration = 0, unixtime1 = 0, unixtime2 = 0;
	string runStartDate = "", runStartTime = "", runStopDate = "", runStopTime = "";	// times are in GMT


	TFile *f = new TFile("MS3umCount_P3JDY.root");
	TTree *t = (TTree*)f->Get("MS 0.3 um count,Davis SCM Environmental Monitoring Processes");
	Long64_t nentries = t->GetEntries();
	int unixtime = 0;
	double value = 0; 
	bool gapInData = 0;
	t->SetBranchAddress("unixtime",&unixtime);
	t->SetBranchAddress("value",&value);
	t->SetBranchAddress("gapInData",&gapInData);

	
	double MSCount = 0;
	double MSAvg = 0;
	double SCEntries = 0;

	// store values for plotting
	vector<double> runs;
	vector<double> avgCount;

	int counter = 0;
	while(!InputList.eof()){
	//while(counter < 100){

		InputList >> run >> duration >> runStartDate >> runStartTime >> unixtime1 >> runStopDate >> runStopTime >> unixtime2;
		printf("\n%i %i %i %i %s %s %s %s \n",run,duration,unixtime1,unixtime2,runStartDate.c_str(),runStartTime.c_str(),runStopDate.c_str(),runStopTime.c_str());

		for (int i = 0; i < nentries; i++) {
			t->GetEntry(i);

			if (unixtime >= unixtime1 && unixtime <= unixtime2) { // scan up to a certain entry
				//printf("run: %i  time: %i  value: %.0f \n",run,unixtime,value);
				MSCount += value;
				SCEntries++;
			}
			if (unixtime > unixtime2) { // go to next run after you find the entries you need for this one
				
				if (SCEntries > 0) MSAvg = MSCount / SCEntries;
				else MSAvg = 0;
				runs.push_back(run);
				avgCount.push_back(MSAvg);
				printf("Average counts this run: %.0f \n",MSAvg);
				MSCount = 0;
				MSAvg = 0;
				SCEntries = 0;
				break;
			}
		}
		counter++;
	}

	// Plot the average count for each run
	TGraph *g = new TGraph(runs.size(), &(runs[0]), &(avgCount[0]));
	g->SetTitle("");
	g->Draw("ap");
	g->GetXaxis()->SetTitle("P3JDY Run");
	g->GetYaxis()->SetTitle("0.3 um counts");

}