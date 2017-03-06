#include <MJTChannelMap.hh>
#include <MGTEvent.hh>
#include <MGTRun.hh>
#include <MGTWaveform.hh>
#include <MJTChannelData.hh>
#include <TFile.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TStyle.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TROOT.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <getopt.h>
#include <dirent.h>
#include <sys/types.h>

using namespace std;

TGraph* graph(string name, vector<double>& x, vector<double>& y, 
	      string sx, string sy, int color){
  if(x.size() == 0 || x.size() != y.size())
    return NULL;
  TGraph* g = new TGraph(x.size(), &x.front(), &y.front());
  g->SetName(name.c_str());
  g->SetTitle(("#font[132]{" + name + "}").c_str());
  g->GetXaxis()->SetTitle(("#font[132]{" + sx + "}").c_str());
  g->GetYaxis()->SetTitle(("#font[132]{" + sy + "}").c_str());
  g->GetXaxis()->SetLabelFont(132);
  g->GetYaxis()->SetTitleFont(132);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.2);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  return g;
}

TH1D* hist(string name, int nbins, vector<double>& x, 
	   string sx, string sy, int color){
  if(x.size() == 0)
    return NULL;
  double mean = 0.0;
  double xmax = -1.0e9;
  double xmin = 1.0e9;
  for(unsigned i=0; i<x.size(); i++){
    mean += x[i];
    xmax = std::max(xmax, x[i]);
    xmin = std::min(xmin, x[i]);
  }
  mean /= x.size();
  double rms = 0.0;
  for(unsigned i=0; i<x.size(); i++)
    rms += pow(x[i], 2);
  rms = sqrt(rms) / x.size();
  TH1D* h = new TH1D(name.c_str(), name.c_str(), nbins, 
		     std::max(mean - 8 * rms, xmin),
		     std::min(mean + 8 * rms, xmax));
  h->GetXaxis()->SetTitle(("#font[132]{" + sx + "}").c_str());
  h->GetYaxis()->SetTitle(("#font[132]{" + sy + "}").c_str());
  h->GetXaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelFont(132);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  for(unsigned i=0; i<x.size(); i++)
    h->Fill(x[i]);
  return h;
}

void halve_array(vector<double>& v, bool ave){
  for(unsigned i=0; i<v.size(); i+=2){
    if(ave)
      v[i/2] = (v[i] + v[i+1]) / 2;
    else
      v[i/2] = v[i];
  }
  v.resize(v.size() / 2);
}

int main(int argc, char* argv[]){

  string pn = "";
  int frun = 0;
  int lrun = 0;
  int fevent = 0;
  int levent = 0;

  static struct option opts[] = {
    {"partnumber", required_argument, NULL, 'p'},
    {"start", required_argument, NULL, 's'},
    {"end", required_argument, NULL, 'e'},
    {"firstevent", required_argument, NULL, 'f'},
    {"lastevent", required_argument, NULL, 'l'}
  };

  int opt = getopt_long(argc, argv, "p:s:e:f:l", opts, NULL);
  while(opt != -1){
    switch(opt){
    case 'p':
      pn = string(getenv("MJDDATADIR"))+"/surfmjd/data/built/"+string(optarg);
      break;
    case 's':
      frun = atoi(optarg); break;
    case 'e':
      lrun = atoi(optarg); break;
    case 'f':
      fevent = atoi(optarg); break;
    case 'l':
      levent = atoi(optarg); break;
    default:
      printf("unrecognized otption %s\n", optarg); break;
    }
    opt = getopt_long(argc, argv, "p:s:e:f:l", opts, NULL);
  }

  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");

  const int nbase = 200;
  const unsigned max_array_size = 2 << 16;

  vector<int> run_list;
  if(lrun <= frun)
    lrun = frun;
  if(frun == lrun){
    DIR* dir;
    struct dirent* infile;
    if((dir = opendir(pn.c_str())) == NULL)
      printf("unable to open directory %s\n", pn.c_str());
    else{
      while((infile = readdir(dir)) != NULL){
	string name(infile->d_name);
	if(name.find("run") + 3 < name.size() && name.find(".") < name.size()){
	  int run = atoi(name.substr(name.find("run") + 3,
				     name.find(".") - 2).c_str());
	  if(run >= frun && run <= lrun)
	    run_list.push_back(run);
	}
      }
    }
  }
  else{
    for(int i=frun; i<=lrun; i++)
      run_list.push_back(i);
  }
  sort(run_list.begin(), run_list.end());
  printf("test\n");

  stringstream ssfrun;
  ssfrun << run_list.front();
  stringstream sslrun;
  sslrun << run_list.back();

  TDirectory* dir = gROOT->CurrentDirectory();
  TFile* tmpfile = TFile::Open((pn+"/OR_run"+ssfrun.str()+".root").c_str());
  MJTChannelMap* chan_map =
    (MJTChannelMap*) ((MJTChannelMap*) tmpfile->Get("ChannelMap"))->Clone();
  tmpfile->Close();
  gROOT->cd(dir->GetPath());

  TChain* tree = new TChain("MGTree");
  for(unsigned i=0; i<run_list.size(); i++){
    stringstream ss;
    ss << run_list[i];
    tree->Add((pn + "/OR_run" + ss.str() + ".root").c_str());
    printf("Adding %s\n", (pn + "/OR_run" + ss.str() + ".root").c_str());
  }

  MGTEvent* event = new MGTEvent();
  tree->SetBranchAddress("event", &event);
  TClonesArray* channelData = new TClonesArray("MJTChannelData");
  tree->SetBranchAddress("channelData", &channelData);
  MGTRun* run = new MGTRun();
  tree->SetBranchAddress("run", &run);
  printf("test\n");
  // vector<uint32_t> enab = set->GetEnabledIDList();

  const int ndet = chan_map->NumberOfRows();
  map<int, int> index_map;
  vector<int> logain_chan, higain_chan;
  vector<string> logain_pos, higain_pos;
  vector<int> ids_chan, pos_chan;
  vector<pair<int, int> > phigain_chan, ppos_chan;
  for(int i=0; i<2000; i++){
    int id = chan_map->GetInt(i, "kIDLo");
    if(id == i){
      index_map[i] = (int) logain_chan.size();
      logain_chan.push_back(i);
      string postmp = chan_map->GetDetectorPos(id);
      string pos(4, '0');
      for(int j=2; j<6; j++) pos[j-2] = postmp[j];
      //int p = 10 * (int) (postmp[3] - '0');
      logain_pos.push_back(pos);
    }
    id = chan_map->GetInt(i, "kIDHi");
    if(id == i){
      index_map[i] = (int) higain_chan.size();
      higain_chan.push_back(i);
      ids_chan.push_back(id);
      string postmp = chan_map->GetDetectorPos(id);
      string pos(4, '0');
      for(int j=2; j<6; j++) pos[j-2] = postmp[j];
      int p = 10 * (int) (postmp[3] - '0');
      p += (int) (postmp[5] - '0');
      higain_pos.push_back(pos);
      ppos_chan.push_back(make_pair(p, i));
      phigain_chan.push_back(make_pair(id, i));
    }
  }

  vector<vector<double> > t(ndet);
  vector<vector<double> > v(ndet);
  vector<vector<double> > vbase(ndet);
  vector<vector<double> > vrms(ndet);
  vector<vector<double> > tevent(ndet);
  vector<vector<double> > trail(ndet);
  vector<string> detname(ndet, "");
  vector<TH1D*> spectra(ndet, NULL);  // plots for max/min value of waveform.

  vector<int> wf_prescale(ndet, 1);
  vector<int> ev_prescale(ndet, 1);
  vector<int> ev_count(ndet, 1);

  if(levent == 0)
    levent = (int) tree->GetEntries();

  tree->Show(0);

  // assumes the first event/run happens at t=0
  tree->GetEvent(0);
  int irun = run->GetRunNumber();
  printf("%i\n", irun);
  double frun_offset = run->GetStartTime();
  double run_offset = frun_offset;
  tree->GetEvent(fevent);
  printf("%e\n", frun_offset);
  double event_offset = event->GetTime();
  printf("test\n");

  //tree->GetEvent(tree->GetEntries() - 1);
  // double lrun_offset = run->GetEndTime();

  // loop over the entries in this run
  for(int iev=fevent; iev<levent; iev++){
    if(iev % 1000 == 0)
      printf("%i\t%i\n", iev, (int) tree->GetEntries());
    if(iev >= tree->GetEntries()) {
      cout << "user's event range is beyond the end of the run!\n";
      break;
    }
    tree->GetEvent(iev);
    stringstream ssiev;
    ssiev << iev; // make a string from the event number
    if(run->GetRunNumber() != irun){
      irun = run->GetRunNumber();
      //run_offset = run->GetStartTime();
      //printf("%e\n", event->GetTime()/1e9);
    }
    double time = (run_offset-frun_offset)*1e9-event_offset + event->GetTime();
    for(int ich=0; ich<(int)event->GetNWaveforms(); ich++){
      MJTChannelData* chan = (MJTChannelData*) channelData->At(ich);
      int ch_id = chan->GetID();
      int ch_index = index_map[ch_id];
      bool logain = chan_map->IsLoGain(ch_id);
      if(logain)
	continue;
    // get string and detector pos's
      int str = chan_map->GetDetectorPos(ch_id)[3] - '0' - 1;
      stringstream ssstr;
      ssstr << str;
      int det = chan_map->GetDetectorPos(ch_id)[5] - '0' - 1;
      stringstream ssdet;
      ssdet << det;
      string detpos = "P" + ssstr.str() + "D" + ssdet.str();

      // throw away pulser tagging channels
      if(str < 0 || det < 0 || ch_index < 0)
	continue;
      if(str < 0 || str >=7 || det < 0|| det >= 5)
	continue;
      if(detname[ch_index] == "")
	detname[ch_index] = detpos;

      MGTWaveform* wf = event->GetWaveform(ich);
      string hname = "h_" + ssiev.str();
      if(logain)
	hname += "_" + detpos + "_logain";
      else
	hname += "_" + detpos + "_higain";
      TH1D* h = wf->GimmeHist(hname.c_str());
      for(int i=0; i<4; i++)
	h->SetBinContent(h->GetNbinsX()-i, h->GetBinContent(h->GetNbinsX()-4));

      // time offset of waveforms that are grouped into the same event
      double offset = wf->GetTOffset();

      double base = 0.0;
      double rms = 0.0;
      for(int i=4; i<nbase; i++)
	base += h->GetBinContent(i);
      base /= nbase;
      for(int i=4; i<nbase; i++)
	rms += pow(h->GetBinContent(i) - base, 2);
      rms = sqrt(rms) / base;

      // max and min values of waveform
      // (useful to have both when we trigger on pos/neg pulses)
      double wmax = -1.0e6;
      double wmin =  1.0e6;
      bool rail = false;

      // looping over the waveform
      for(int i=1; i<=h->GetNbinsX(); i+=wf_prescale[ch_index]){
	t[ch_index].push_back((time + h->GetBinCenter(i) + offset) / 1e3);
	v[ch_index].push_back(h->GetBinContent(i));
	if(t[ch_index].size() == max_array_size){
	  halve_array(t[ch_index], false);
	  halve_array(v[ch_index], false);
	  wf_prescale[ch_index] *= 2;
	}

	// if the sample is a rail, save the time of the rail, and the bool
	if(abs(v[ch_index].back()) > 8000 && !rail)
          rail = true;
	
	wmax = std::max(h->GetBinContent(i) - base, wmax);
	wmin = std::min(h->GetBinContent(i) - base, wmin);
      }
      delete h;
      if(rail)
	trail[ch_index].push_back((time + offset) / 1e9);
      if(spectra[ch_index] == NULL){
	int dpos = atoi(higain_pos[ch_index].substr(3, 3).c_str());
	int color = (dpos + 1) * 2;
	if(dpos > 4)
	  color = 1;
	string title = "#font[132]{" + higain_pos[ch_index] + "}";
	h = new TH1D(("hspec_" + detname[ch_index]).c_str(), "",
		     2 << 13, -1.0 * (2 << 12), 2 << 12);
	h->GetXaxis()->SetTitle("Energy (ADC)");
	h->SetTitle(title.c_str());
	h->SetLineColor(color);
	spectra[ch_index] = h;
      }
      if(TMath::Abs(wmax) > TMath::Abs(wmin))
	spectra[ch_index]->Fill(wmax);
      else
	spectra[ch_index]->Fill(wmin);
      if(ev_count[ch_index] % ev_prescale[ch_index] == 0){
	tevent[ch_index].push_back(time / 1000 / ev_prescale[ch_index]);
	vbase[ch_index].push_back(base / ev_prescale[ch_index]);
	vrms[ch_index].push_back(rms / ev_prescale[ch_index]);
	ev_count[ch_index] = 0;
      }
      else{
	tevent[ch_index].back() += time / 1000 / ev_prescale[ch_index];
	vbase[ch_index].back() += base / ev_prescale[ch_index];
	vrms[ch_index].back() += rms / ev_prescale[ch_index];
      }
      ev_count[ch_index] ++;
      if(tevent.size() == max_array_size && 
	 ev_count[ch_index] == ev_prescale[ch_index]){
	printf("halving baseline arrays\n");
	halve_array(tevent[ch_index], false);
	halve_array(vbase[ch_index], true);
	halve_array(vrms[ch_index], true);
	ev_prescale[ch_index] *= 2;
	ev_count[ch_index] = 0;
      }
    }
  }

  string fname = "splice_waveforms_" + ssfrun.str();
  if(lrun > frun)
    fname += "_" + sslrun.str();
  TFile* outfile = new TFile((fname + ".root").c_str(), "recreate");
  outfile->cd();

  printf("test\n");

  for(int i=0; i<ndet; i++){
    if(t[i].size() == 0)
      continue;
    int dpos = atoi(higain_pos[i].substr(3, 3).c_str());
    int color = (dpos + 1) * 2;
    if(dpos > 4)
      color = 1;
    graph(higain_pos[i] + "_gwf", t[i], v[i], 
	  "Time (#mu s)", "ADC", color)->Write();
    graph(higain_pos[i] + "_gbase", tevent[i], vbase[i], 
	  "Time (#mu s)", "Baseline (ADC)", color)->Write();
    graph(higain_pos[i] + "_grms", tevent[i], vrms[i],
	  "Time (#mu s)", "Baseline RMS (ADC)", color)->Write();
    hist(higain_pos[i] + "_hbase", 300, vbase[i], 
	 "Baseline (ADC)", "", color)->Write();
    hist(higain_pos[i] + "_hrms", 300, vrms[i], 
	 "Baseline RMS (ADC)", "",color)->Write();
  }
  
  printf("test\n");

  for(int i=0; i<ndet; i++)
    if(spectra[i])
      spectra[i]->Write();
  
  outfile->Close();

  return 1;
}
