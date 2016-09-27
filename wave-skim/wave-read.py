#!/usr/local/bin/python
##!/usr/common/usg/software/python/2.7.6/bin/python
from ROOT import *
# import numpy as np

"""Perform the all-important check that wave-skim's
waveforms are actually readable."""

gROOT.SetBatch(kTRUE)
gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
gStyle.SetOptStat(0)

f = TFile("./output/waveSkim-test.root")

t = f.Get("waveTree")
t1 = f.Get("eventTree")

can = TCanvas("can","Bob Ross's Canvas",800,600)

# for i in xrange(t.GetEntries()):
for i in xrange(10):
	t.GetEntry(i)
	t1.GetEntry(t.iEventTree)

	print "iEventTree %d  itr %d  run %d  gatEvt %d  chan %d  enf %.2f  t50 %.2f  time %.2f" % (t.iEventTree, t.itr, t.run, t.gatEvt, t.channel, t.trapENFCal, t.blrwfFMR50, t.runTime)

	itr = -1
	for j in xrange(t1.channel.size()):
		if t1.channel.at(j) == t.channel:
			itr = j
			print "  ET: chan %.2f  enf %.2f  t50 %.2f" % (t.channel, t1.trapENFCal.at(j), t1.blrwfFMR50.at(j))

	wf = t.event.GetWaveform(t.itr)
	h = wf.GimmeHist()

	# still need to (understand and) adjust for t_offset.
	# check out "GetFunction" in MGTWaveform.hh to get the parameter.  Need to make a TF1.

	h0 = TH1D("h","h",2018,0,20180)
	h0.SetTitle("Ch %i  E %.2f  t50 %.2f  tsc50 %.2f" % (t.channel, t.trapENFCal, t.blrwfFMR50, t1.TSCurrent50nsMax.at(itr)))
	h0.GetXaxis().SetTitle("t (ns)")
	h0.GetYaxis().SetTitle("ADC")
	h0.SetMaximum(1000)
	h0.Draw()
	h.Draw("same")

	can.Print("./output/wave_%d.pdf" % i)
	can.Print("./output/wave_%d.C" % i)