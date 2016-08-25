#!/usr/common/usg/software/python/2.7.6/bin/python
###!/usr/local/bin/python

from ROOT import *
import numpy as np
import branch_libs as lib

def main():
	"""Take events in the skim file and pull all corresponding built/gatified data.
	Future: alternately be able to use a GATDataSet directly without the skim files."""

	# =========== Skim file & output file ===========
	skimLoc = "$MJDDATADIR/surfmjd/analysis/skim/DS1/20160621_265313037/*.root"
	# skimLoc = "/Users/wisecg/datasets/ds1/*.root"
	# wsOut = "./output/waveSkim-1550-1650.root"
	wsOut = "./output/waveSkim-1500-2000-mH-2.root"

	# =========== Skim file cuts ===========
	burstCut = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116"

	# low-energy noisy runs cut - need to research & refine
	# runCut = "run!=13312 && run!=13121 && run!=13004 && run!=12766 && run!=12735 && run!=12445 && run!=11175 && run!=12723 && run!=12746 && run!=12767 && run!=13071 && run!=13073 && run!=13074 && run!=13120 && run!=13205 && run!=13306 && run!=13307 && run!=9857 && run!=9862 && run!=9863"

	# bigCut = "channel%2==0 && mH==1 && (trapENFCal>1550 && trapENFCal<1650) && !wfDCBits && !muVeto && !isLNFill &&" + burstCut

	bigCut = "channel%2==0 && mH>1 && sumEH>1500 && !wfDCBits && isGood && " + burstCut

	# =========== Ready? Go! ===========
	skimmer(bigCut, skimLoc, wsOut)
	# skimChecker(wsOut)


def skimmer(bigCut, skimLoc, wsOut):
	"""Creates a wave-skim file.  Must be run on PDSF.

	NOTE: GATWaveformBrowser does not save the entire event.
	It saves only the waveforms which pass your cut.
	If you save other parameters from skim/built data,
	you may be confused because things like "sumEH" don't match
	the waveforms you have saved.

	My solution: write two trees into the output file.

	"waveTree" - Contains single waveforms, some basic physics info,
	             and the entry numbers of the skim and gat files.

	"eventTree" - Contains the full entries of each event (all skim and gat branches)

	The idea is to look at the waveforms in "waveTree" and be able to grab any/all
	skim or gatified parameters from the other tree, as needed.
	"""
	print "Cutting on: \n",bigCut

	skim = TChain("skimTree")
	skim.Add(skimLoc)
	lib.SetTreeInputs(skim, lib.skimDict)

	skim.SetEntryList(0)
	skim.Draw(">>elist", bigCut, "entrylist")
	elist = gDirectory.Get("elist")
	print "Found ",elist.GetN()," skim file entries."
	skim.SetEntryList(elist)

	ds = GATDataSet()
	wb = GATWaveformBrowser()
	ds = wb.LoadSkimWaveforms(skim, bigCut)
	print "Found ",wb.GetNWaveforms()," waveforms."

	gat = ds.GetGatifiedChain()
	lib.SetTreeInputs(gat, lib.gatDict)

	built = ds.GetBuiltChain()
	lib.SetTreeInputs(built, lib.builtDict)

	out = TFile(wsOut,"RECREATE")

	waveTree = TTree("waveTree","wave-skim single waveforms")
	enf = np.zeros(1,dtype=float)
	t50 = np.zeros(1,dtype=float)
	runTime = np.zeros(1,dtype=float)
	gatEnt, gatHit, skimEnt, builtEnt, eventTreeEntry, chn = 0., 0., 0., 0., 0., 0.
	wf = MGTWaveform()
	waveTree.Branch("gatEnt",long(gatEnt),"gatEnt/L")
	waveTree.Branch("gatHit",long(gatHit),"gatHit/L")
	waveTree.Branch("builtEnt",long(builtEnt),"builtEnt/L")
	waveTree.Branch("skimEnt",long(skimEnt),"skimEnt/L")
	waveTree.Branch("eventTreeEntry",long(eventTreeEntry),"eventTreeEntry/L")
	waveTree.Branch("waveform",wf)
	waveTree.Branch("channel",int(chn),"channel/I")
	waveTree.Branch("trapENFCal",enf,"trapENFCal/D")
	waveTree.Branch("blrwfFMR50",t50,"blrwfFMR50/D")
	waveTree.Branch("runTime",runTime,"runTime/D")

	# save the entry numbers for the full event tree
	EntryList = []

	# fill waveTree
	lastEvent = 0
	eventTreeEntry = -1
	eventMismatchCount = 0
	wfMismatchCount = 0
	for waveNum in xrange(wb.GetNWaveforms()):

		gatEnt = wb.GetEntryNumber(waveNum)
		gatHit = wb.GetIterationNumber(waveNum)
		gat.GetEntry(gatEnt)
		built.GetEntry(gatEnt)
		builtEnt = built.GetEntryNumber(gatEnt)
		skimEnt = 0
		for ientry in xrange(elist.GetN()):
			entryNumber = skim.GetEntryNumber(ientry)
			skim.LoadTree( entryNumber )
			skim.GetEntry( entryNumber )
			# gat.LoadTree returns the entry number of the original tree
			if skim.iEvent==gat.LoadTree(gatEnt):
				skimEnt = entryNumber
				break
		skim.GetEntry(skimEnt)

		if abs(lib.timestamp.at(0)/1E8 - lib.s_tloc_s[0]) > 0.001:
			print "waveform",waveNum,": mismatched events!"
			eventMismatchCount += 1
			print "skim - run %d  enf.at(0) %.3f  enf.size %d  time %.2f" % (lib.s_run[0], lib.s_trapENFCal.at(0), lib.s_trapENFCal.size(), lib.s_tloc_s[0])
			print "gat  - run %d  enf.at(0) %.3f  enf.size %d  time %.2f\n" % (lib.run[0], lib.trapENFCal.at(0), lib.trapENFCal.size(), lib.timestamp.at(0)/1E8)
			continue

		# output some physics
		wf = wb.GetWaveform(waveNum)

		nullchk = str(wf)
		if "nil" in nullchk:
			print "waveform",waveNum,",iteration ",gatHit,": unexpected number of waveforms ..."
			wfMismatchCount +=1
			print "skim - run %d  enf.at(0) %.3f  enf.size %d  time %.2f" % (lib.s_run[0], lib.s_trapENFCal.at(0), lib.s_trapENFCal.size(), lib.s_tloc_s[0])
			print "gat  - run %d  enf.at(0) %.3f  enf.size %d  time %.2f\n" % (lib.run[0], lib.trapENFCal.at(0), lib.trapENFCal.size(), lib.timestamp.at(0)/1E8)
			continue

		chn = lib.channel.at(gatHit)
		enf[0] = lib.trapENFCal.at(gatHit)
		t50[0] = lib.blrwfFMR50.at(gatHit)
		runTime[0] = lib.timestamp.at(gatHit)/1E8

		# so you can match waveTree to eventTree
		if lastEvent != gatEnt:
			eventTreeEntry += 1
			EntryList.append([gatEnt,skimEnt])

		waveTree.Fill()

		lastEvent = gatEnt

	print "\nDone. For %d waveforms:" % wb.GetNWaveforms()
	print "\t%d had mismatched gat/skim events based on timestamp differences (and were skipped)" % eventMismatchCount
	print "\t%d had fewer wf's than expected in the built data (and were skipped)." % wfMismatchCount

	waveTree.Write()

	# now fill the full event tree by looping over EntryList
	print "\n filling full event tree ...\n"
	eventTree = TTree("eventTree","wave-skim full output")
	outDict = lib.CreateOutputDict('all')
	lib.SetTreeOutputs(eventTree, outDict)

	for i in EntryList:
		gatEntry = i[0]
		skimEntry = i[1]

		gat.GetEntry(gatEntry)
		built.GetEntry(gatEntry)
		skim.GetEntry(skimEntry)

		# verify we get the same output as before
		# print "skim - run %d  enf-0 %.3f  size %d  time %.2f" % (lib.s_run[0],lib.s_trapENFCal.at(0),lib.s_trapENFCal.size(),lib.s_tloc_s[0])
		# print "gat  - run %d  enf-0 %.3f  size %d  time %.2f\n" % (lib.run[0],lib.trapENFCal.at(0),lib.trapENFCal.size(),lib.timestamp.at(0)/1E8)

		eventTree.Fill()

	eventTree.Write()
	out.Close()


def skimChecker(wsOut):
	"""Incomplete ..."""
	wsFile = TFile(wsOut)
	waveTree = wsFile.waveTree
	eventTree = wsFile.eventTree

	print "waveTree: %d  eventTree: %d" % (waveTree.GetEntries(),eventTree.GetEntries())

	for i in xrange(waveTree.GetEntries()):
		waveTree.GetEntry(i)
		eventTreeEntry = waveTree.eventTreeEntry
		chn = waveTree.channel

if __name__ == "__main__":
	main()
