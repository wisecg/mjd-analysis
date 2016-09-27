#!/usr/common/usg/software/python/2.7.6/bin/python
###!/usr/local/bin/python

from ROOT import *
import numpy as np
import branch_libs as lib
import cPickle as pickle
import sys,time,collections

def main(argv):
	"""
	Saves all waveforms passing a cut into an output file with two trees.
	"waveTree" - Contains single waveforms, some basic physics info,
	             and the entry numbers of the skim and gat files.
	"eventTree" - Contains the full entries of each event (all skim, gat, & built branches)
	The idea is to look at the waveforms in "waveTree" and be able to grab any/all
	skim or gatified parameters from the other tree as needed.
	"""
	gROOT.SetBatch(kTRUE)  # don't display waveform plots
	writeET = 1  # write (1) or don't write (0) eventTree
	entryLimit = -1  # events to save (-1 saves all entries)
	useFile = 0
	if (len(argv) > 1):
		useFile = int(argv[1])
	print "useFile:",useFile

	"""=========== 1. Skim-based wave-skimmer ===========
	This identifies events of interest from cuts based on skim file parameters,
	and then pulls the corresponding gatified & built data."""

	# skimLoc = "$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/*.root"
	# skimLoc = "/Users/wisecg/datasets/ds1/*.root"
	# skimLoc = "$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/*.root"
	skimLoc = "~/wave-skim/input/standardCutLE_DS1.root"

	outFile = "./input/waveSkim-test.root"

	theCut = "trapENFCal < 200 && trapENFCal > 1 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && trapETailMin < 0"

	# theCut = "trapENFCal < 200 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && mH>1 && isGood && !isLNFill && !wfDCBits && !muVeto"

	SeymourSkimmer(theCut, skimLoc, outFile, writeET, useFile, entryLimit)

	"""=========== 2. GAT-based wave-skimmer (takes a run list) ===========
	Take a list of runs from a text file, and apply a cut based on gatified parameters.
	Since the gatified data is not as "cleaned" as the skim file (pulsers are usually present),
	options are provided to first check how many events a given cut will retain,
	and a maximum number of waveforms to save per run.
	"""
	runList = "Debug.txt"
	outFile = "./input/waveSkim-cal-12470.root"
	# gatCut = "channel%2==0 && trapENFCal < 4000 && trapENFCal > 1 && (trapENFCal > 800 || trapENFCal < 500) && !wfDCBits"
	gatCut = "trapENFCal > 1591 && trapENFCal < 1593 && !wfDCBits"

	# GatGrabber(runList, gatCut, outFile, writeET, "",13)  # "test", or ""


# ==========================================================================================
# ==========================================================================================


def SeymourSkimmer(theCut, skimLoc, outFile, writeET, useFile, entryLimit=-1):
	print "Starting wave-skim.  Cutting on: \n\n",theCut,"\n"

	if (useFile != 1):

		skim = TChain("skimTree")
		skim.Add(skimLoc)
		# skim.Add("$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/skimDS1_1.root")
		# skim.Add("$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/skimDS1_2.root")
		lib.SetTreeInputs(skim, lib.skimDict)

		skim.SetEntryList(0)
		skim.Draw(">>elist", theCut, "entrylistarray")
		elist = gROOT.FindObject("elist")
		print "Found",elist.GetN(),"skim file entries."
		skim.SetEntryList(elist)

		if entryLimit == -1:
			entryLimit = elist.GetN()
		if entryLimit > elist.GetN():
			entryLimit = elist.GetN()

		# Save a list of runs containing events which pass the cut.
		skim.Draw("run",theCut,"GOFF")
		v1 = skim.GetV1() # a ROOT.PyDoubleBuffer, very irritating to work with
		v1_list = list(set(v1[n] for n in xrange(skim.GetEntries())))
		RunList = sorted(v1_list)
		del RunList[0] # because the first element is zero for some reason

		print "Scanning skim file ..."
		EventMap = []
		for iSkimList in xrange(entryLimit):

			iSkimEntry = skim.GetEntryNumber(iSkimList);  # global index
			a = skim.LoadTree( iSkimEntry )	# returns entry number from this tree
			b = skim.GetEntry( iSkimEntry )
			subEntries = elist.GetSubListForEntry(iSkimEntry)

			skimTime = lib.s_tloc_s.at(0)
			skimGatEntry = lib.s_iEvent[0]
			skimRun = lib.s_run[0]

			numPass = skim.Draw("channel:Entry$","","GOFF",1,iSkimList)
			chans = skim.GetV1()
			chanList = list(set(int(chans[n]) for n in xrange(numPass)))

			# TEntryListArray method -- fails when there is more than one tree in the chain
			# numHits = lib.s_trapENFCal.size()
			# for i in xrange(numHits):
			# 	# if subEntries.Contains(i):
			# 	if subEntries.Contains(iSkimEntry,skim,i):
			# 		print "   Output: Entry %d  Hit %d  chan %d  enf %.2f" % (iSkimEntry,i,lib.s_channel.at(i),lib.s_trapENFCal.at(i))

			# save several parameters so we can verify we have the right gatified entry
			EventMap.append([skimRun, skimGatEntry, skimTime, chanList])

			if (iSkimList%100==0):
				print iSkimList,skimRun,skimGatEntry,skimTime,chanList

		pickle.dump(EventMap,open("EventMap.p","wb"))
		pickle.dump(RunList,open("RunList.p","wb"))

	EventMap = sorted(pickle.load(open("EventMap.p","rb")))
	RunList = pickle.load(open("RunList.p","rb"))

	# make a lookup dict for EventMap, to facilitate the upcoming loop over run numbers
	print "\n\nCreating lookup table ..."
	RangeDict = {}
	run, prevrun = 0,-1
	xlo, xhi = -1,-1
	lastxlo = 0
	for evt in EventMap:
		i = EventMap.index(evt)
		run = evt[0]
		xhi = i
		# print i,run
		if run != prevrun:
			# print "run:",prevrun,"lo:",xlo,"hi:",xhi-1
			RangeDict[prevrun] = [xlo,xhi-1]
			xlo = i
			lastxlo = xhi-1
		prevrun = run
	# print "last run:",run,"lo:",lastxlo,"hi:",len(EventMap)
	RangeDict[run] = [lastxlo,len(EventMap)]
	del RangeDict[-1]
	RangeLookup = collections.OrderedDict(sorted(RangeDict.items()))

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

	eventTree = TTree("eventTree","wave-skim full output")
	outDict = lib.CreateOutputDict('gat')
	lib.SetTreeOutputs(eventTree, outDict)

	# loop over runs with events that pass the cut
	print "\n\nScanning gatified + built data ... %d entries\n\n" % len(EventMap)
	eventMismatchCount = 0
	iEventTree[0] = -1
	start = time.time()
	for run in RangeLookup:

		ds = GATDataSet(run)
		gat = ds.GetGatifiedChain()
		lib.SetTreeInputs(gat, lib.gatDict)
		built = ds.GetBuiltChain()

		# loop over events in this run
		# reference: EventMap - 0 run, 1 skimGatEntry, 2 skimTime, 3 channel list
		lo = RangeDict[run][0]
		hi = RangeDict[run][1]+1
		if hi > len(EventMap):
			hi = RangeDict[run][1]

		for iEM in xrange(lo,hi):
			evt = EventMap[iEM]
			built.GetEntry(long(evt[1]))
			gat.GetEntry(long(evt[1]))

			gatTime = lib.timestamp.at(0)/1E8

			if abs(evt[2] - gatTime) > 0.01:
				print "Timestamps don't match!  Skipping event at index",iEM
				print "skim: run %d  gatEntry %d  skimTime %.2f  gatTime %.2f  chans" % (run,evt[1],evt[2],gatTime),evt[3]
				eventMismatchCount += 1
				continue

			iEventTree[0] += 1

			# ian's magic trick: assign it twice and it sticks.
			eventOut = built.event
			eventOut = built.event
			# "To be honest, I have no clue what is happening when you call that line, since MGTEvent does not have an = operator or a copy constructor explicitly defined. What this means is that an automatically generated one is being used, which might be the problem."

			# loop over individual hits in this event
			numWFs = int(built.event.GetNWaveforms())
			for i in xrange(numWFs):
				try:
					chan = lib.channel.at(i)
					if chan in evt[3]:
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

			if writeET == 1:
				eventTree.Fill()

			if iEM%100==0:
				elapsed = time.time() - start
				progress = "%.0f%%" % ((float(iEM)/len(EventMap))*100)
				print "%-3s  %d WF's saved.\tCurrent run: %i  Wall-clock time: %.2f" % (progress,iEM,run,elapsed)
				start = time.time()

	print "End of Scan.  Event mismatch count: ",eventMismatchCount,"of ",len(EventMap)," entries."
	waveTree.Write()

	if writeET == 1:
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
	main(sys.argv[:])
