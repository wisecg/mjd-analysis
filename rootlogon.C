{
	cout << "Loading the CGW-MJD environment ...\n";

	// Plot style macro
	// gROOT->ProcessLine(".x /Users/wisecg/dev/MJDWavePlotStyle.C");
	gROOT->ProcessLine(".x /Users/wisecg/dev/MJDClintPlotStyle.C");

	// Load classes from mgsw scripts only:
	//gROOT->ProcessLine(".x $MGDODIR/Root/LoadMGDOClasses.C"); // loads automatically from MGDOMJ
	//gROOT->ProcessLine(".x $MGDODIR/Majorana/LoadMGDOMJClasses.C");
	//gROOT->ProcessLine(".x $GATDIR/LoadGATClasses.C");

	// Load classes manually:

	// load gat classes
	if(!TClass::GetDict("GATAnalysisEvent")) {
    	gSystem->Load("$GATDIR/lib/libGATBaseClasses");
    	gSystem->Load("$GATDIR/lib/libGATMGTEventProcessing");
    	gSystem->Load("$GATDIR/lib/libGATMGOutputMCRunProcessing");
    	gSystem->Load("$GATDIR/lib/libGATAnalysis");
		gROOT->ProcessLine(".include \"$GATDIR/BaseClasses\"");
	}

	// load MGDO classes
	if(!TClass::GetDict("TImage")) gSystem->Load("libGui");
	if(!TClass::GetDict("MGTWaveform")) {
	    gSystem->Load("$MGDODIR/lib/libMGDOBase" );
	    gSystem->Load("$MGDODIR/lib/libMGDOTransforms");
	    gSystem->Load("$MGDODIR/lib/libMGDORoot");
	    // gROOT->ProcessLine("using namespace CLHEP;");
	    // gROOT->ProcessLine(".include \"$CLHEP_INCLUDE_DIR\"");
	    gROOT->ProcessLine(".include \"$MGDODIR/Base\"");
	    gROOT->ProcessLine(".include \"$MGDODIR/Root\"");
	    gROOT->ProcessLine(".include \"$MGDODIR/Transforms\"");

	    // added by clint
	    // gSystem->Load("$MGDODIR/tam/lib/libTAM");
	    // gROOT->ProcessLine(".include \"$MGDODIR/tam/include\""); // contains header files
	}

	// load MGDOMJClasses
	if(!TClass::GetDict("MJTGretina4DigitizerData")) {
	    gSystem->Load("$MGDODIR/lib/libMGDOTabree");
	    gSystem->Load("$MGDODIR/lib/libMGDOMajorana");
	    gROOT->ProcessLine(".include \"$MGDODIR/Tabree\"");
	    gROOT->ProcessLine(".include \"$MGDODIR/Majorana\"");
	    gROOT->ProcessLine("#include <complex>");
	    gROOT->ProcessLine("((TExec*)TRef::GetListOfExecs()->At(TRef::AddExec(\"ReadMap\")))->SetAction(\"gFile->Get(\\\"ChannelMap\\\");\")");
	    gROOT->ProcessLine("((TExec*)TRef::GetListOfExecs()->At(TRef::AddExec(\"ReadSettings\")))->SetAction(\"gFile->Get(\\\"ChannelSettings\\\");\")");
	}

	//cout << "Loading MJDB classes .." << endl;
  	//gSystem->Load("$MGDODIR/lib/libMGDOMJDB");
  	//gROOT->ProcessLine(".include \"$MGDODIR/MJDB\"");
  	//gROOT->ProcessLine("using namespace MJDB;");

	cout << "I am Groot." << endl;
}
