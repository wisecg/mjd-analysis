#!/usr/local/bin/python
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
# import wave_libs as lib

# High-gain Module 1 detectors
chanDict = {584: "P1D1", 582: "P1D2", 580: "P1D3", 578: "P1D4", 692: "P2D1", 648: "P2D2", 640: "P2D3", 642: "P2D4", 616: "P3D1", 610: "P3D2", 608: "P3D3", 664: "P3D4", 624: "P4D1", 628: "P4D2", 688: "P4D3", 694: "P4D4", 614: "P4D5", 680: "P5D1", 678: "P5D2", 672: "P5D3", 696: "P5D4", 632: "P6D1", 630: "P6D2", 626: "P6D3", 690: "P6D4", 600: "P7D1", 598: "P7D2", 594: "P7D3", 592: "P7D4"}

def main():
	"""Create a .npz library of single-site waveforms for each detector."""

	gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

	waveTree = TChain("waveTree")
	waveTree.Add("./input/waveSkim-cal-12470.root")
	entries = waveTree.GetEntries()

	event = MGTEvent()
	waveTree.SetBranchAddress("theEvent",event)

	can = TCanvas("can","Bob Ross's Canvas",800,600)

	for iEnt in xrange(entries):
		waveTree.GetEntry(iEnt)

		print "\n%d run %d  gat %d  chan %d  enf %.2f  t50 %.2f  time %.2f" % (iEnt, waveTree.run, waveTree.gatEvt, waveTree.channel, waveTree.trapENFCal, waveTree.blrwfFMR50, waveTree.runTime)

		print type(event)

		wf = event.GetWaveform(waveTree.itr)

		try:
			waveName = "./output/wave_%d_%d_%d.pdf" % (waveTree.run,waveTree.gatEvt,waveTree.itr)
			wf.Draw()
			can.Print(waveName)
		except:
			print "balls, it didn't work"


if __name__ == "__main__":
	main()
