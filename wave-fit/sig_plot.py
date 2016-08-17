#!/usr/local/bin/python
"""Illustrate a couple different methods of accessing and plotting waveforms.

C. Wiseman, 6/12/16.
"""
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import rc
# import toy_model as tm
# import pymc
# import ROOT
from ROOT import gApplication, TChain, TCanvas


def main():
	"""."""
	gApplication.ExecuteFile("pyrootlogon.C")
	simplePlots()


def simplePlots():
	"""Use the if(1) switch to activate a block."""

	# The "wave" branch in this file is a single MGTWaveform, not a full MGTEvent.
	waves = TChain("waves")
	waves.Add("WSR_WFOutput.root")
	print "Got %i entries." % waves.GetEntries()

	# View member variables
	if(0):
		for i in xrange(waves.GetEntries()):
			waves.GetEntry(i)
			print "run %i\t iEvent %i\t channel %i\t trapENFCal %f\t trapECal %f" % (waves.run, waves.iEvent, waves.channel, waves.trapENFCal, waves.trapECal)

	# Two ways to recover the MGTWaveform
	waves.GetEntry(0)
	waveform1 = waves.wave.GetVectorData()  # gives a vector of doubles. checked with type(waveform1)
	waveform2 = waves.wave.GimmeHist()  # gives a th1d

	# Look at values
	if(0):
		print "waveform size : %i" % waveform1.size()
		for i in xrange(waveform1.size()):
			print "i %i  adc %f" % (i, waveform1[i])

	# 1. Plot the vector<double> with matplotlib
	if(1):
		plt.style.use('mjWaveforms')  # using stylesheet in ~/.matplotlib/stylelib
		fig = plt.figure(figsize=(10, 7), facecolor='w')
		fig.add_subplot(111)
		timestamps = np.arange(0, len(waveform1) * 10, 10)  # in ns
		plt.plot(timestamps, waveform1)
		plt.xlim(0, 20160)  # remove last sample at 2017
		plt.ylabel('ADC')
		plt.xlabel('time [ns]')
		plt.show()

	# 2. Plot the TH1D
	if(0):
		c = TCanvas("c", "Bob Ross's Canvas", 800, 600)
		waveform2.GetXaxis().SetRangeUser(0, 20160)	 # chop off the last sample
		waveform2.Draw()
		c.Print("c1.pdf")


if __name__ == "__main__":
	main()
