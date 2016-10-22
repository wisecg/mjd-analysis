#!/usr/local/bin/python
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import wave_model as wm
import sys

def main(argv):
	"""Interactively or rapid-draw waveforms that pass a given TCut."""
	if len(argv) < 1:
		print "Usage: ./wave-view.py [-R (rapid draw)  -I (interactive draw)]"
		return
	opt = argv[0]

	waveTree = TChain("waveTree")
	waveTree.Add("./data/wave-standardCutDS1-newBranch.root")
	print "Found",waveTree.GetEntries(),"input waveforms."

	eventTree = TChain("eventTree")
	eventTree.Add("./data/wave-standardCutDS1-newBranch.root")

	# initial "diagnosis" cut
	theCut = "trapENFCal > 0.8 && trapENFCal < 5 && channel!=600 && channel!=692 && d2wf0MHzTo50MHzPower < 4000"

	theCut = "trapENFCal > 0.8 && trapENFCal < 15 && channel!=600 && channel!=692 && d2wf0MHzTo50MHzPower < 4000"

	# extended run cut
	theCut += " && run != 9648 && run != 9663 && run != 10185 && run != 10663 && run != 10745 && run != 11175 && run != 12445 && run != 12723 && run != 12735 && run != 13004 && channel != 580"

	waveTree.Draw(">>elist", theCut , "entrylist")
	elist = gDirectory.Get("elist")
	waveTree.SetEntryList(elist)
	print "Found",elist.GetN(),"entries passing cut."

	if (opt == "-R"):
		drawRapid(waveTree, elist)

	if (opt == "-I"):
		drawInteractive(waveTree, elist)

def drawInteractive(waveTree, elist):
	print "Press Enter to start ..."

	plt.ion()
	fig = plt.figure(figsize=(9,5), facecolor='w')

	iList = -1
	while(True):
		iList += 1
		if iList > elist.GetN():
			exit(1)
		value = raw_input()
		if value=='q':		# quit
			exit(1)
		if value=='p': 		# previous entry
			iList -= 2
		if value.isdigit():	# go to entry number
			iList = int(value)

		entry = waveTree.GetEntryNumber(iList);
		waveTree.LoadTree(entry)
		waveTree.GetEntry(entry)

		trapENFCal = waveTree.trapENFCal	# waveTree parameters
		riseTime = waveTree.blrwfFMR50
		run = waveTree.run
		chan = waveTree.channel
		runTime = waveTree.runTime

		iEvent = waveTree.iEventTree		# eventTree parameters
		eventTree.GetEntry(iEvent)
		power = eventTree.d2wf0MHzTo50MHzPower.at(waveTree.itr)

		print "%d (%.0f %%)  %d  %d  %.2fkev  %.2fns  %.2fs" % (iList, 100*(iList/float(elist.GetN())), run, chan, trapENFCal, riseTime, runTime)

		signal = wm.CGWave(waveTree.event.GetWaveform(waveTree.itr))
		waveBLSub = signal.GetBLSub()
		waveFilt = signal.GetFilt()
		waveTS = signal.GetTS()

		plt.clf()
		ax = plt.subplot(111)
		ax.set_xlabel("time (ns)")
		ax.set_ylabel("ADC")

		ax.plot(waveTS,waveBLSub)
		ax.plot(waveTS,waveFilt, color='red', linewidth=2., alpha=0.5)
		xmin, xmax = np.amin(waveTS), np.amax(waveTS)
		ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
		ax.set_xlim([xmin,xmax])
		ax.set_ylim([ymin-0.1*ymin,ymax+0.1*ymax])



def drawRapid(waveTree, elist):

	ts = np.arange(0,2016*10,10)
	adc = np.ones(len(ts))

	plt.ion()
	fig = plt.figure(figsize=(9,5), facecolor='w')
	ax1 = plt.subplot(111)
	ax2 = plt.subplot(111)
	ax1.set_xlabel("time (ns)")
	ax1.set_ylabel("ADC")
	wave, = ax1.plot(ts, adc)
	waveF, = ax2.plot(ts, adc, color='red', linewidth=2., alpha=0.5)

	for iList in xrange(elist.GetN()):

		plt.pause(0.05)

		entry = waveTree.GetEntryNumber(iList);
		waveTree.LoadTree(entry)
		waveTree.GetEntry(entry)

		trapENFCal = waveTree.trapENFCal	# waveTree parameters
		riseTime = waveTree.blrwfFMR50
		run = waveTree.run
		chan = waveTree.channel
		runTime = waveTree.runTime

		iEvent = waveTree.iEventTree		# eventTree parameters
		eventTree.GetEntry(iEvent)
		power = eventTree.d2wf0MHzTo50MHzPower.at(waveTree.itr)

		print "%d (%.0f %%)  %d  %d  %.2fkev  %.2fns  %.2fs" % (iList, 100*(iList/float(elist.GetN())), run, chan, trapENFCal, riseTime, runTime)

		signal = wm.CGWave(waveTree.event.GetWaveform(waveTree.itr))
		waveBLSub = signal.GetBLSub()
		waveFilt = signal.GetFilt()
		waveTS = signal.GetTS()

		xmin, xmax = np.amin(waveTS), np.amax(waveTS)
		ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
		wave.set_ydata(waveBLSub)
		wave.set_xdata(waveTS)
		ax1.set_xlim([xmin,xmax])
		ax1.set_ylim([ymin-0.1*ymin,ymax+0.1*ymax])

		waveF.set_ydata(waveFilt)
		wave.set_xdata(waveTS)
		ax2.set_xlim([xmin,xmax])
		ax2.set_ylim([ymin-0.1*ymin,ymax+0.1*ymax])

		fig.canvas.draw()


if __name__ == "__main__":
	main(sys.argv[1:])
