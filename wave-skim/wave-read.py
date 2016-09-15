#!/usr/common/usg/software/python/2.7.6/bin/python
###!/usr/local/bin/python

from ROOT import *
import numpy as np

"""Perform the all-important check that wave-skim's
waveforms are actually readable."""

gROOT.SetBatch(kTRUE)
gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

# f = TFile("./output/waveSkim-cal-12470.root")
f = TFile("./output/waveSkim-skimTest.root")
t = f.Get("waveTree")

can = TCanvas("can","Bob Ross's Canvas",800,600)

for i in xrange(t.GetEntries()):
	t.GetEntry(i)

	print "iEventTree %d  itr %d  run %d  gatEvt %d  chan %d  enf %.2f  t50 %.2f  time %.2f" % (t.iEventTree, t.itr, t.run, t.gatEvt, t.channel, t.trapENFCal, t.blrwfFMR50, t.runTime)

	wf = t.event.GetWaveform(t.itr)
	waveName = "./output/wave_%d.pdf" % (i)
	wf.Draw()
	can.Print(waveName)