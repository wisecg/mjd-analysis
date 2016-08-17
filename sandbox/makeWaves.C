void makeWaves()
{
	/*
	// This is the easiest way to view the waveforms of a run (12828) with a cut applied.
	// Type these into your ROOT prompt:
	GATDataSet ds(12828);
	TChain *built = ds.GetBuiltChain();
	GATWaveformBrowser wb(built,"trapENFCal < 100 && trapENFCal > 10");
	int i = 0;

	// hitting up-arrow / enter lets you page through waveforms
	wb.DrawWaveform(++i);
	*/

	// ------------------------------------------------------------------------------
	// This way lets you pull out some more info and print to pdf.  Run it as:
	// root[0] .X makeWaves.C

	GATDataSet ds(12828); // 15365
	TChain *waves = ds.GetBuiltChain();

	char cutstr[200] = "trapENFCal > 10 && trapENFCal < 100";	// optional cut, set to "" to bypass
	GATWaveformBrowser wb(waves, cutstr);
	int numWFs = (int)wb.GetNWaveforms();

	// Loop over waveforms and pull out the histos
	// Be careful, it prints a lot
	for (int i = 0; i < (int)numWFs; i++)
	{
		MGTWaveform *twf = wb.GetWaveform(i);
		if(twf==NULL) { cout << "couldn't get waveform!\n"; return; }
		MGTWaveform wf = *(twf);
		TH1D *hist = wf.GimmeHist();

		// Pull some variables out from the waveform object.
		// If you want to know run number, etc, put in the gatified data branch name.
		//
		double ch=0, pos=0, det=0, enf=0;
		const char* wfInfo = "P:D:channel:trapENFCal";
		waves->Draw(wfInfo,"","GOFF",1,i);
		if(waves->GetV1()==NULL||waves->GetV2()==NULL||waves->GetV3()==NULL||waves->GetV4()==NULL){
			cout << "Warning!  Variable doesn't exist!  Can't make waveform title!\n";
			break;
		}
		// V1, V2 etc change according to the variables you specify in "wfInfo"
		pos = waves->GetV1()[0];
		det = waves->GetV2()[0];
		ch = waves->GetV3()[0];
		enf = waves->GetV4()[0];
		TString title = TString::Format("P%gD%g(%g),%s=%g",pos,det,ch,waves->GetVar4()->GetTitle(),enf);
		hist->SetTitle(title.Data());

		// PDF Output
		TCanvas *c = new TCanvas("c","Bob Ross's Canvas",800,600);
		char wfName[500];
		sprintf(wfName,"./wave_%i_%.0f_%.0f.pdf",i,ch);
		// hist->SetMaximum(300);	// ADC y-axis limits if you want
		// hist->SetMinimum(0);
		hist->Draw("L");
		c->Print(wfName);
		delete c;
	}

}