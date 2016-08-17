#!/usr/local/bin/python
"""The one where we fit big waveforms to little ones."""

import numpy as np
import matplotlib.pyplot as plt
import pymc
from pymc import deterministic, DiscreteUniform, Normal, HalfNormal
from ROOT import gApplication, TChain
import sig_model as sig
from scipy.optimize import curve_fit
# from scipy.ndimage.filters import gaussian_filter


# High-gain Module 1 detectors
chanDict = {584: "P1D1", 582: "P1D2", 580: "P1D3", 578: "P1D4", 692: "P2D1", 648: "P2D2", 640: "P2D3", 642: "P2D4", 616: "P3D1", 610: "P3D2", 608: "P3D3", 664: "P3D4", 624: "P4D1", 628: "P4D2", 688: "P4D3", 694: "P4D4", 614: "P4D5", 680: "P5D1", 678: "P5D2", 672: "P5D3", 696: "P5D4", 632: "P6D1", 630: "P6D2", 626: "P6D3", 690: "P6D4", 600: "P7D1", 598: "P7D2", 594: "P7D3", 592: "P7D4"}


# -----------------------------------------------------------------
def main():
	"""."""
	gApplication.ExecuteFile("pyrootlogon.C")

	# Window parameters
	loWin = 800
	hiWin = 1200

	# Ben sez: use last ZERO CROSSING as the lo window! Then 200-250 samples after.

	# Get a template waveform
	temp = TChain("waves")
	temp.Add("WSR_DEP_WFOutput.root")
	templateRaw = getWaveform(temp, 2)
	temp.GetEntry(2)
	baselineAvg, noiseAvg = findBaseline(templateRaw)
	templateSubt = templateRaw
	templateSubt[:] = [x - baselineAvg for x in templateRaw]  # a python 'list comprehension'
	window = np.concatenate((np.arange(0, loWin), np.arange(hiWin, 2017)), axis=0)
	templateWin = np.delete(templateSubt, window)
	templateInfo = {'t0': temp.blrwfFMR50, 'trapENFCal': temp.trapENFCal, 'noiseSigma': noiseAvg, 'baseline': baselineAvg}

	# Get low-energy waveforms (put this in a loop later)
	wave = TChain("waves")
	wave.Add("WSR_WFOutput.root")
	signalRaw = getWaveform(wave, 0)
	window = np.concatenate((np.arange(0, loWin), np.arange(hiWin, 2017)), axis=0)
	baselineAvg, noiseAvg = findBaseline(signalRaw)
	signalSubt = signalRaw
	signalSubt[:] = [x - baselineAvg for x in signalRaw]  # a python 'list comprehension'
	signalWin = np.delete(signalSubt, window)
	signalInfo = {'t0': wave.blrwfFMR50, 'trapENFCal': wave.trapENFCal, 'noiseSigma': noiseAvg, 'baseline': baselineAvg}

	print signalInfo

	en = signalInfo['trapENFCal']
	tg = signalInfo['t0']
	na = signalInfo['noiseSigma']

	# set up a pymc model
	modelPyMC = pymc.Model(sig.createSignalModel_DEPWF(signalWin, en, tg, na, templateWin))
	M = pymc.MCMC(modelPyMC)

	# Change the proposal distribution from default (prior) to a discrtete normal distribution for t0
	M.use_step_method(pymc.DiscreteMetropolis, M.switchpoint, proposal_sd=4., proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.slowness_sigma, proposal_sd=0.1, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.wfScale, proposal_sd=0.1, proposal_distribution='Normal')

	# Run the MCMC
	M.sample(iter=10000)

	# how long do you wait to converge before you decide you're really sampling your posterior?
	burnin = 4000

	# pull out the fit parameters after burnin
	t0 = np.median(M.trace('switchpoint')[burnin:])
	slowness_sigma = np.median(M.trace('slowness_sigma')[burnin:])
	wfScale = np.median(M.trace('wfScale')[burnin:])

	print "t0 : gat %f  mcmc %f  template %f" % (tg, t0, templateInfo['t0'])
	print "slowness : mcmc %f" % slowness_sigma
	print "energy : gat %f  mcmc %f" % (en, wfScale)

	# scale the template by the mcmc parameters
	# remember t_offset is in samples, not ns
	e_scale = wfScale / templateInfo['trapENFCal']
	t_offset = (templateInfo['t0'] - t0) / 10
	loWinDEP = loWin + int(t_offset)
	hiWinDEP = hiWin + int(t_offset)
	windowDEP = np.concatenate((np.arange(0, loWinDEP), np.arange(hiWinDEP, 2017)), axis=0)
	scaledDEP = np.delete(templateSubt, windowDEP) * e_scale

	print "e_scale : this pulse %f  template gat %f  ratio %f" % (wfScale, templateInfo['trapENFCal'], e_scale)
	print "t_offset : template-signal ", t_offset
	print "loWin %i  hiWin %i  loDEP %i  hiDEP %i  offset %f  int offset %i," % (loWin, hiWin, loWinDEP, hiWinDEP, t_offset, int(t_offset))
	print "size - signalWin: %i  scaledDEP %i" % (signalWin.size, scaledDEP.size)

	# plot that shit
	# plt.style.use('mjWaveforms')  # stylesheet saved in ~/.matplotlib/stylelib
	# plt.figure(figsize=(10, 7), facecolor='w')
	# plt.plot(scaledDEP, color='red')
	# plt.plot(signalWin, color='blue')
	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif')
	# plt.xlabel(r't (samples)', fontsize=16)
	# plt.ylabel(r'ADC', fontsize=16)
	# plt.show()

	# plot the whole waveform (minus the template part we have to chop off)
	# (make it more like ben's plots with t0 drawn onto the whole WF)

	# check out the traces
	# MCMC full trace
	f, axarr = plt.subplots(3, sharex=True)
	axarr[0].plot(M.trace('switchpoint')[:])
	axarr[2].set_xlabel('MCMC Step Number')
	axarr[0].set_ylabel('t_0')
	axarr[1].plot(M.trace('slowness_sigma')[:])
	axarr[1].set_ylabel('slo')
	axarr[2].plot(M.trace('wfScale')[:])
	axarr[2].set_ylabel('enf')

	plt.show()
	plt.savefig("MCMC_steps.pdf")



# -----------------------------------------------------------------


def getWaveform(waves, i):
	"""Grab an MJ waveform."""
	waves.GetEntry(i)
	waveform = waves.wave.GetVectorData()  # gives a vector<double>

	# Look at values
	if(0):
		print "waveform size : %i" % waveform.size()
		for j in xrange(waveform.size()):
			print "i %i  adc %f" % (j, waveform[j])

	# Convert the vector<double> to a numpy array
	x = np.fromiter(waveform, dtype=np.double)

	# Remove last entry
	x1 = np.delete(x, x.size - 1)

	return x1


def gauss_function(x, a, x0, sigma):
	"""."""
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def findBaseline(signalRaw):
	"""Find the average starting baseline from an MJD waveform."""

	(hist, bins) = np.histogram(signalRaw, bins=np.arange(-8000, 8000, 1))
	fitfunc = lambda p, x: p[0] * np.exp(-0.5 * ((x - p[1]) / p[2])**2) + p[3]
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	c, pcov = curve_fit(gauss_function, bins[:-1], hist, p0=[np.amax(hist), bins[np.argmax(hist)], 5])
	mu = c[1]
	std = c[2]
	return mu, std


if __name__ == "__main__":
	main()
