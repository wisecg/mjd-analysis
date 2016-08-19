#!/usr/local/bin/python
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
	"""Good practice in Python is to document your functions at the beginning like this.
	-> This tutorial has four sub-programs: skimFiles, gatData, builtData, and eventLoop.
	   Each function demonstrates a few ways of processing MJD data and making some
	   simple plots.
	-> Below those are some "utility functions" that are called repeatedly in the sub-programs.
	-> One nice thing about python is that you can "import" this file, and use the utilities
	   in other programs.
	-> At the bottom is the magic invocation that allows you to run this from the command line:
		python tutorial.py
	Clint Wiseman, USC/Majorana, 8/11/2016.
	"""

	# ========= run routines =========
	# skimFiles()
	# gatData()
	# eventLoop()
	# builtData()
	pyWaves()


def skimFiles():
	"""Load skim files and make a few very simple plots."""

	# skim files are found on PDSF:
	# /project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS1/20160621_265313037
	# /project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/DS0/20160621_265313037

	skim = TChain("skimTree")
	skim.Add("~/datasets/ds1/*.root")
	print "Got %i entries." % skim.GetEntries()

	# print all the branches to the screen
	# skim.GetListOfBranches().Print()

	# draw a simple plot with a cut and save
	can = TCanvas("c1","Bob Ross's Canvas",800,600)
	skim.Draw("trapENFCal","mH==1")
	can.Print("energy.pdf")

	# project the data onto a user-specified histogram
	# and edit the ROOT plot a bit
	gStyle.SetOptStat(0) # kill stats box
	gStyle.SetOptTitle(0) # kill title
	can2 = TCanvas("c2","Bob Ross's Canvas",800,600)
	hist = TH1D("h1","Energy Plot",2900,100,3000)
	skim.Project("h1","trapENFCal","mH==1")
	hist.Draw()
	can2.SetLogy()
	hist.GetXaxis().SetTitle("trapENFCal (keV)")
	hist.GetYaxis().SetTitle("Counts (arb)")
	can2.Print("energy2.pdf")


def gatData():
	"""Look at gatified data and make simple plots in a more 'python' way."""

	ds = GATDataSet(12826,12831)
	gatified = ds.GetGatifiedChain()

	# choose only single-site events, and high-gain channels (even numbered channels)
	drawString = "trapENFCal"
	cutString = "mH==1 && channel%2==0"
	lo, hi = 10, 3000
	hist = TH1D("h1","Energy Plot",2900,lo,hi)
	gatified.Project("h1",drawString,cutString)

	# convert a ROOT histogram to a NumPy array.
	# this is MUCH more common datatype (stack overflow overfloweth with support)
	arr = rootToArray(hist,lo,hi)
	print "It's a ",type(arr),"!!"

	# make one more plot
	hist2 = TH1D("h2","High-M Energy",2900,lo,hi)
	gatified.Project("h2",drawString,"mH>1 && channel%2==0")
	arr2 = rootToArray(hist2,lo,hi)

	plt.style.use('mjWaveforms')  # using stylesheet in ~/.matplotlib/stylelib
	fig = plt.figure(figsize=(10, 7), facecolor='w')
	plt.plot(arr,"b")
	plt.plot(arr2,"r")
	axes = plt.gca()
	axes.set_yscale('log')
	plt.xlabel('trapENFCal')
	plt.ylabel('Counts (arb)')
	plt.legend(('single-site','multi-site'), loc='best')
	plt.show()
	plt.savefig("SSMSEnergies.pdf")


def eventLoop():
	"""Use ROOT's GetEntry to scan over every entry in a file."""

	ds = GATDataSet(12826)
	gatified = ds.GetGatifiedChain()
	entries = gatified.GetEntries()

	# need to declare individual branch access when we loop over events.
	timestamp = ROOT.std.vector("double")()
	trapENFCal = ROOT.std.vector("double")()
	channel = ROOT.std.vector("double")()
	theRun = np.zeros(1, dtype=float)
	startTime = np.zeros(1, dtype=float)
	stopTime = np.zeros(1, dtype=float)
	gatified.SetBranchAddress("timestamp", timestamp)
	gatified.SetBranchAddress("trapENFCal", trapENFCal)
	gatified.SetBranchAddress("channel", channel)
	gatified.SetBranchAddress("run", theRun)
	gatified.SetBranchAddress("startTime", startTime)
	gatified.SetBranchAddress("stopTime", stopTime)

	for iEnt in xrange(entries):
	    gatified.GetEntry(iEnt)
	    for i in xrange(channel.size()):
	        print trapENFCal.at(i)


def builtData():
	"""Look at the waveforms in the built data."""
	ds = GATDataSet(12826)
	built = ds.GetBuiltChain()
	wb = GATWaveformBrowser(built,"trapENFCal < 100 && trapENFCal > 10")
	num = wb.GetNWaveforms()

	# print waveforms, recover GAT parameters, and convert to numpy arrays
	can = TCanvas("c1","Bob Ross's Canvas",900,600)
	gStyle.SetOptStat(0)
	gROOT.ProcessLine(".X MJDWavePlotStyle.C")
	for iEnt in xrange(3): 	# can go up to 'num'
		waveName = "wave_%i.pdf" % iEnt

		# 1. simple grab 'n draw
		# wb.DrawWaveform(iEnt)
		# can.Print(waveName) # we could stop here, this is pretty good already

		# 2. now pull out the actual histogram and some variables from the reconstructed (gatified) data.
		# if you want to know run number, etc, put in the gatified data branch name.
		twf = wb.GetWaveform(iEnt)
		hist = twf.GimmeHist()
		built.Draw("P:D:channel:trapENFCal","","GOFF",1,iEnt)
		pos = built.GetV1()[0]
		det = built.GetV2()[0]
		chn = built.GetV3()[0]
		enf = built.GetV4()[0]
		title = built.GetVar4().GetTitle()
		waveTitle = "P%iD%i(%i),%s=%.2f" % (pos,det,chn,title,enf)
		hist.SetTitle(waveTitle)
		hist.Draw()
		# can.Print(waveName)


def pyWaves():
	"""convert waveforms to numpy arrays, baseline subtract, and get estmation parameters from gat."""

	ds = GATDataSet(12826)
	built = ds.GetBuiltChain() # automatically 'friends' the gatified chain
	recon = ds.GetGatifiedChain()
	wb = GATWaveformBrowser(built,"trapENFCal < 100 && trapENFCal > 10")

	waveNum = 2
	waveRaw = getWaveform(wb,waveNum)

	baseAvg, noiseAvg = findBaseline(waveRaw)
	waveBLSub = waveRaw
	waveBLSub[:] = [x - baseAvg for x in waveRaw]  # a python 'list comprehension'

	# pull out metadata for this waveform
	gatEnt = wb.GetEntryNumber(waveNum)
	gatItr = wb.GetIterationNumber(waveNum)
	built.GetEntry(gatEnt)
	waveInfo = {
		'run': recon.run,
		'channel': built.channel.at(gatItr),
		'trapENFCal': built.trapENFCal.at(gatItr),
		't0': built.blrwfFMR50.at(gatItr),
		'noiseAvg': noiseAvg,
		'baseAvg': baseAvg
		}
	print waveInfo

	plt.style.use('mjWaveforms')  # using stylesheet in ~/.matplotlib/stylelib
	fig = plt.figure(figsize=(10, 7), facecolor='w')

	loWin, hiWin = 0, 2017 # samples every 10 ns
	ts = np.arange(loWin * 10, hiWin * 10, 10)  # (start,stop,step) in ns
	plt.xlim(loWin * 10, hiWin * 10)  # remove last sample at 2017

	plt.ylabel('ADC')
	plt.xlabel('time [ns]')
	plt.plot(ts, waveBLSub, 'b', alpha=0.75)

	# draw the t0 calculated by GAT
	axes = plt.gca()
	plt.plot((waveInfo['t0'], waveInfo['t0']), axes.get_ylim(), 'r-', linewidth=1) # draw line at t0
	plt.show()
	plt.savefig("wave_py.pdf")


# ========= utilities used in routines =========

def rootToArray(hist, xLow, xHigh):
	"""Take a ROOT TH1D histogram and a range, get a numpy array."""
	binLow = hist.FindBin(xLow)
	binHigh = hist.FindBin(xHigh)
	loopArray = range(binLow, binHigh )
	npArray = np.empty_like(loopArray, dtype=np.float)
	for (iArray, iBin) in enumerate(loopArray):
		npArray[iArray] = np.float(hist.GetBinContent(iBin))
	return npArray


def getWaveform(wb, i):
	"""Take a GATWaveformBrowser, entry number, and give a numpy array."""
	twf = wb.GetWaveform(i) # MGTWaveform
	waveform = twf.GetVectorData() # vector<double>
	x = np.fromiter(waveform, dtype=np.double) # numpy array
	x1 = np.delete(x, x.size - 1)
	return x1


def findBaseline(signalRaw):
	"""Find the average starting baseline from an MJD waveform."""
	(hist, bins) = np.histogram(signalRaw, bins=np.arange(-8000, 8000, 1))
	fitfunc = lambda p, x: p[0] * np.exp(-0.5 * ((x - p[1]) / p[2])**2) + p[3]
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	c, pcov = curve_fit(gauss_function, bins[:-1], hist, p0=[np.amax(hist), bins[np.argmax(hist)], 5])
	mu = c[1]
	std = c[2]
	return mu, std


def gauss_function(x, a, x0, sigma):
	"""A homebrewed Gaussian."""
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


if __name__ == "__main__":
	main()
