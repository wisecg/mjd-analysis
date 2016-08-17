{
	// 6000323
	// 6000325
	// 6000327
	// 6000329
	// 6000331
	// 6000333

	int run = 6000370;
	char name[200];
	sprintf(name,"/global/project/projectdirs/majorana/data/mjd/surftest/data/built/NonLinearity/OR_run%i.root",run);
	TFile *f = new TFile(name);
	TTree *t = (TTree*)f->Get("MGTree");
	t->GetEntries();
	MGTEvent *event = new MGTEvent();
	t->SetBranchAddress("event",&event);
	t->GetEntry(0);
	event->GetWaveform(0);
	MGTWaveform *w = event->GetWaveform(0);
	TCanvas *c = new TCanvas("","",800,600);
	w.Draw();
	char cname[200];
	sprintf(cname,"wf_%i.png",run);
	c->Print(cname);
}