#!/usr/common/usg/software/python/2.7.6/bin/python
###!/usr/local/bin/python

from ROOT import *
import numpy as np
import branch_libs as lib
import sys

def main():
	"""
	Saves all waveforms passing a cut into an output file with two trees.
	"waveTree" - Contains single waveforms, some basic physics info,
	             and the entry numbers of the skim and gat files.
	"eventTree" - Contains the full entries of each event (all skim, gat, & built branches)
	The idea is to look at the waveforms in "waveTree" and be able to grab any/all
	skim or gatified parameters from the other tree as needed.
	"""
	gROOT.SetBatch(kTRUE)  # don't display waveform plots
	writeET = 0  # write (1) or don't write (0) eventTree

	"""=========== 1. Skim-based wave-skimmer ===========
	This identifies events of interest from cuts based on skim file parameters,
	and then pulls the corresponding gatified & built data."""

	skimLoc = "$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/*.root"
	# skimLoc = "/Users/wisecg/datasets/ds1/*.root"
	outFile = "./output/waveSkim-skimTest.root"
	burstCut = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116"
	bigCut = "channel%2==0 && mH>1 && sumEH>1500 && !wfDCBits && isGood && " + burstCut

	SeymourSkimmer(bigCut, skimLoc, outFile, writeET)

	"""=========== 2. GAT-based wave-skimmer (takes a run list) ===========
	Take a list of runs from a text file, and apply a cut based on gatified parameters.
	Since the gatified data is not as "cleaned" as the skim file (pulsers are usually present),
	options are provided to first check how many events a given cut will retain,
	and a maximum number of waveforms to save per run.
	"""

	runList = "Debug.txt"
	outFile = "./output/waveSkim-cal-12470.root"
	# gatCut = "channel%2==0 && trapENFCal < 4000 && trapENFCal > 1 && (trapENFCal > 800 || trapENFCal < 500) && !wfDCBits"
	gatCut = "trapENFCal > 1591 && trapENFCal < 1593 && !wfDCBits"

	# GatGrabber(runList, gatCut, outFile, writeET, "",13)  # "test", or ""


# ==========================================================================================
# ==========================================================================================


def SeymourSkimmer(bigCut, skimLoc, outFile, writeET):
	print "Starting wave-skim.  Cutting on: \n\n",bigCut,"\n"

	skim = TChain("skimTree")
	skim.Add(skimLoc)
	lib.SetTreeInputs(skim, lib.skimDict)

	skim.SetEntryList(0)
	skim.Draw(">>elist", bigCut, "entrylist")
	elist = gDirectory.Get("elist")
	print "Found",elist.GetN(),"skim file entries."
	skim.SetEntryList(elist)

	# Save a list of runs containing events which pass the cut.
	skim.Draw("run",bigCut,"GOFF")
	v1 = skim.GetV1() # a ROOT.PyDoubleBuffer, very irritating to work with
	v1_list = list(set(v1[n] for n in xrange(skim.GetEntries())))
	RunList = sorted(v1_list)
	del RunList[0] # because the first element is zero for some reason

	print "Scanning skim file ..."
	EventMap = []
	for iSkimList in xrange(elist.GetN()):

		iSkimEntry = skim.GetEntryNumber(iSkimList);
		iSkimLocTreeEntry = skim.LoadTree( iSkimEntry )
		skim.GetEntry( iSkimEntry )
		skimTime = lib.s_tloc_s.at(0)
		skimGatEntry = lib.s_iEvent[0]
		skimRun = lib.s_run[0]

		thisEventCut = bigCut + "&& Entry$==" + str(iSkimEntry)
		numPass = int(skim.Draw("channel",thisEventCut,"GOFF"))
		chans = skim.GetV1()
		chanList = list(set(int(chans[n]) for n in xrange(numPass)))

		# save several parameters so we can verify we have the right gatified entry
		EventMap.append([skimRun, skimGatEntry, skimTime, chanList])


	# set output file
	out = TFile(outFile,"RECREATE")

	waveTree = TTree("waveTree","wave-skim single waveforms")
	enf = np.zeros(1,dtype=float)
	t50 = np.zeros(1,dtype=float)
	theRun = np.zeros(1,dtype=int)
	gatEvt = np.zeros(1,dtype=long)
	runTime = np.zeros(1,dtype=float)
	iEventTree = np.zeros(1,dtype=long)
	chn = np.zeros(1,dtype=int)
	itr = np.zeros(1,dtype=int)
	eventOut = MGTEvent()
	waveTree.Branch("event",eventOut)
	waveTree.Branch("iEventTree",iEventTree,"iEventTree/L")
	waveTree.Branch("itr",itr,"itr/I")
	waveTree.Branch("run",theRun,"theRun/I")
	waveTree.Branch("gatEvt",gatEvt,"gatEvt/L")
	waveTree.Branch("channel",chn,"channel/I")
	waveTree.Branch("trapENFCal",enf,"trapENFCal/D")
	waveTree.Branch("blrwfFMR50",t50,"blrwfFMR50/D")
	waveTree.Branch("runTime",runTime,"runTime/D")

	# loop over runs with events that pass the cut
	print "Scanning gatified + built data ..."
	eventMismatchCount, wfMismatchCount = 0, 0
	iEventTree[0] = 0
	for i in range(len(RunList)):
		run = int(RunList[i])
		ds = GATDataSet(run)

		gat = ds.GetGatifiedChain()
		lib.SetTreeInputs(gat, lib.gatDict)

		built = ds.GetBuiltChain()

		# loop over events in this run
		# EventMap - 0 run, 1 skimGatEntry, 2 skimTime, 3 channel list
		for evt in EventMap:
			if evt[0]==run:
				built.GetEntry(long(evt[1]))
				gat.GetEntry(long(evt[1]))
				gatTime = lib.timestamp.at(0)/1E8
				print "skm -- run %-6d  skimGatEntry %-7d  time %-8.2f  chans" % (evt[0],evt[1],evt[2]),evt[3]
				print "gat -- run %-6d  skimGatEntry %-7d  time %-8.2f" % (lib.run[0],evt[1],gatTime)

				if abs(evt[2] - gatTime) > 0.01:
					print "Timestamps don't match!  Skipping event ..."
					eventMismatchCount += 1
					continue

				iEventTree[0] += 1

				# ian's magic trick: assign it twice and it sticks.
				eventOut = built.event
				eventOut = built.event
				# "To be honest, I have no clue what is happening when you call that line, since MGTEvent does not have an = operator or a copy constructor explicitly defined. What this means is that an automatically generated one is being used, which might be the problem."
				numWFs = int(built.event.GetNWaveforms())

				# loop over individual hits in this event
				for i in xrange(numWFs):
					if lib.channel.at(i) in evt[3]:
						try:
							runTime[0] = evt[2]
							chn[0] = lib.channel.at(i)
							enf[0] = lib.trapENFCal.at(i)
							t50[0] = lib.blrwfFMR50.at(i)
							itr[0] = i
							theRun[0] = run
							gatEvt[0] = long(evt[1])
							waveTree.Fill()
						except:
							print "Filling failed for some damn reason.",sys.exc_info()[0]
							continue

	waveTree.Write()

	if writeET == 1:
		print "Now writing full-event tree ..."
		eventTree = TTree("eventTree","wave-skim full output")
		outDict = lib.CreateOutputDict('skimGat')
		lib.SetTreeOutputs(eventTree, outDict)

		eventOut = MGTEvent()
		eventTree.Branch("event",eventOut)

		# loop over events in this run
		# EventMap - 0 run, 1 skimGatEntry, 2 skimTime, 3 channel list
		for i in range(len(RunList)):
			run = int(RunList[i])
			ds = GATDataSet(run)
			gat = ds.GetGatifiedChain()
			built = ds.GetBuiltChain()
			lib.SetTreeInputs(gat, lib.gatDict)
			lib.SetTreeInputs(built, lib.builtDict)

			for evt in EventMap:
				if evt[0]==run:
					built.GetEntry(long(evt[1]))
					gat.GetEntry(long(evt[1]))

					# ian's magic trick: assign it twice and it sticks.
					eventOut = built.event
					eventOut = built.event

					if abs(evt[2] - lib.timestamp.at(0)/1E8) > 0.01:
						print "Timestamps don't match!  Skipping event ..."
						continue
					eventTree.Fill()

		eventTree.Write()

	out.Close()


# ==========================================================================================
# ==========================================================================================


def GatGrabber(runList, gatCut, outFile, writeET, action="test cut", maxwfs=-1):
	print "Starting wave-skim.  Cutting on: \n\n",gatCut,"\n"

	file_handle = open(runList, 'r')
	file_lines = file_handle.readlines()
	RunList = []
	for val in file_lines:
		RunList.append(int(val))
	print RunList

	# set output file
	out = TFile(outFile,"RECREATE")

	waveTree = TTree("waveTree","wave-skim single waveforms")
	enf = np.zeros(1,dtype=float)
	t50 = np.zeros(1,dtype=float)
	theRun = np.zeros(1,dtype=int)
	gatEvt = np.zeros(1,dtype=long)
	runTime = np.zeros(1,dtype=float)
	iEventTree = np.zeros(1,dtype=long)
	chn = np.zeros(1,dtype=int)
	itr = np.zeros(1,dtype=int)
	eventOut = MGTEvent()
	waveTree.Branch("event",eventOut)
	waveTree.Branch("iEventTree",iEventTree,"iEventTree/L")
	waveTree.Branch("itr",itr,"itr/I")
	waveTree.Branch("run",theRun,"theRun/I")
	waveTree.Branch("gatEvt",gatEvt,"gatEvt/L")
	waveTree.Branch("channel",chn,"channel/I")
	waveTree.Branch("trapENFCal",enf,"trapENFCal/D")
	waveTree.Branch("blrwfFMR50",t50,"blrwfFMR50/D")
	waveTree.Branch("runTime",runTime,"runTime/D")

	# loop over runs with events that pass the cut
	print "Scanning gatified + built data ..."
	iEventTree[0] = 0
	for i in range(len(RunList)):
		run = int(RunList[i])
		ds = GATDataSet(run)

		gat = ds.GetGatifiedChain()
		lib.SetTreeInputs(gat, lib.gatDict)

		built = ds.GetBuiltChain()

		# a TEntryListArray, the magical class that makes this possible
		gat.Draw(">>elist", gatCut, "entrylistarray")
		elist = gDirectory.Get("elist")
		gat.SetEntryList(elist)

		print "\n\nRun %d  entries %d  entry list %d  runtime %.2f" % (run, gat.GetEntries(), elist.GetN(), ds.GetRunTime()/1E9)

		if action == "test":
			continue

		evtsToKeep = maxwfs
		if maxwfs == -1:
			evtsToKeep = elist.GetN()
		if maxwfs > elist.GetN():
			evtsToKeep = elist.GetN()

		# loop over entries
		for iGatList in xrange(evtsToKeep):

			iGatEntry = gat.GetEntryNumber(iGatList);
			gat.LoadTree(iGatEntry)
			gat.GetEntry(iGatEntry)
			subEntries = elist.GetSubListForEntry(iGatEntry)
			subList = subEntries.GetSubLists()

			gatTime = lib.timestamp.at(0)/1E8
			numHits = lib.trapENFCal.size()

			# ian's magic trick: assign it twice and it sticks.
			eventOut = built.event
			eventOut = built.event
			# "To be honest, I have no clue what is happening when you call that line, since MGTEvent does not have an = operator or a copy constructor explicitly defined. What this means is that an automatically generated one is being used, which might be the problem."
			numWFs = int(built.event.GetNWaveforms())

			print "%d : run %d  entry %d  time %.2f  hits %d  numWFs %d" % (iGatList, run, iGatEntry, gatTime, numHits, numWFs)

			for hit in xrange(numHits):
				print "   Hit %d  chan %d  enf %.2f" % (hit,lib.channel.at(hit),lib.trapENFCal.at(hit))

			for i in xrange(numHits):
				if subEntries.Contains(i):
					try:
						runTime[0] = gatTime
						gatEvt[0] = iGatEntry
						chn[0] = lib.channel.at(i)
						enf[0] = lib.trapENFCal.at(i)
						t50[0] = lib.blrwfFMR50.at(i)
						itr[0] = i
						theRun[0] = run
						waveTree.Fill()
						print "   Output: Hit %d  chan %d  enf %.2f" % (itr[0], chn[0], enf[0])
					except:
						print "   Failed to fill for some damn reason.",sys.exc_info()[0]
						continue
			print ""

			iEventTree[0] += 1

	waveTree.Write()

	if writeET == 1:
		print "Now writing full-event tree ..."
		eventTree = TTree("eventTree","wave-skim full output")
		outDict = lib.CreateOutputDict('gat')
		lib.SetTreeOutputs(eventTree, outDict)
		eventOut = MGTEvent()
		eventTree.Branch("event",eventOut)

		for i in range(len(RunList)):
			run = int(RunList[i])
			ds = GATDataSet(run)
			gat = ds.GetGatifiedChain()
			lib.SetTreeInputs(gat, lib.gatDict)
			built = ds.GetBuiltChain()
			lib.SetTreeInputs(built, lib.builtDict)

			gat.Draw(">>elist", gatCut, "entrylistarray")
			elist = gDirectory.Get("elist")
			gat.SetEntryList(elist)

			evtsToKeep = maxwfs
			if maxwfs == -1:
				evtsToKeep = elist.GetN()
			if maxwfs > elist.GetN():
				evtsToKeep = elist.GetN()

			for iGatList in xrange(maxwfs):
				iGatEntry = gat.GetEntryNumber(iGatList);
				gat.LoadTree(iGatEntry)
				gat.GetEntry(iGatEntry)
				eventOut = built.event
				eventOut = built.event
				eventTree.Fill()

		eventTree.Write()

	out.Close()


if __name__ == "__main__":
	main()
