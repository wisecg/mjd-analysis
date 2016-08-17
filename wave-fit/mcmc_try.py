#!/usr/local/bin/python
"""The one where Clint MCMC's waveforms."""

import numpy as np
import matplotlib.pyplot as plt
import pymc
from pymc import deterministic, DiscreteUniform, Normal, HalfNormal
from ROOT import gApplication, TChain
import toy_model as tm
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
# plt.style.use('mjWaveforms')  # stylesheet saved in ~/.matplotlib/stylelib


# High-gain Module 1 detectors
chanDict = {584: "P1D1", 582: "P1D2", 580: "P1D3", 578: "P1D4", 692: "P2D1", 648: "P2D2", 640: "P2D3", 642: "P2D4", 616: "P3D1", 610: "P3D2", 608: "P3D3", 664: "P3D4", 624: "P4D1", 628: "P4D2", 688: "P4D3", 694: "P4D4", 614: "P4D5", 680: "P5D1", 678: "P5D2", 672: "P5D3", 696: "P5D4", 632: "P6D1", 630: "P6D2", 626: "P6D3", 690: "P6D4", 600: "P7D1", 598: "P7D2", 594: "P7D3", 592: "P7D4"}


# -----------------------------------------------------------------
def main():
	"""."""
	gApplication.ExecuteFile("pyrootlogon.C")

	# The "wave" branch in this file is a single MGTWaveform, not a full MGTEvent.
	waves = TChain("waves")
	waves.Add("WSR_WFOutput.root")
	print "Got %i entries." % waves.GetEntries()

	# View member variables
	if(0):
		for i in xrange(waves.GetEntries()):
			waves.GetEntry(i)
			print "run %i\t iEvent %i\t channel %i\t trapENFCal %f\t trapECal %f\t blrwfFMR50 %f" % (waves.run, waves.iEvent, waves.channel, waves.trapENFCal, waves.trapECal, waves.blrwfFMR50)

	# Grab MJD waveform rising edge with windowing (~200 samples around rise time.)
	# and subtract out the baseline.
	loWin = 800
	hiWin = 1400
	window = np.concatenate((np.arange(0, loWin), np.arange(hiWin, 2017)), axis=0)

	signalRaw = getWaveform(waves, 0)  # 2017 samples (from signal.size)

	baselineAvg = subtractBaseline(signalRaw)
	signalSubt = signalRaw
	signalSubt[:] = [x - baselineAvg for x in signalRaw]  # a python 'list comprehension'

	signalWin = np.delete(signalSubt, window)

	# plot the waveform
	# plt.figure(0, figsize=(10, 7), facecolor='w')
	# plt.plot(signalWin)
	# plt.savefig("waveform.pdf")

	# the prior
	# plt.plot(np.arange(0, 100), 0.5 * np.ones(100), color="blue")
	# plt.plot(np.arange(100, 200), 0.5 * np.ones(100), color="red")
	# plt.axvline(100, color="black", linestyle=":")
	# plt.savefig("prior.pdf")

	# Generate the signal
	# signalWin = generateWaveform()

	# set up a pymc model
	siggen_model = pymc.Model(tm.createSignalModel(signalWin))

	M = pymc.MCMC(siggen_model)

	# Change the proposal distribution from default (prior) to a normal distribution for t0
	# And make t_0 discrete
	M.use_step_method(pymc.DiscreteMetropolis, M.switchpoint, proposal_sd=4., proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.early_mu, proposal_sd=0.1, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.late_mu, proposal_sd=0.1, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.early_sigma, proposal_sd=0.1, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.late_sigma, proposal_sd=0.1, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=0.1, proposal_distribution='Normal')

	# Run the MCMC
	M.sample(iter=10000)

	# how long do you wait to converge before you decide you're really sampling your posterior?
	burnin = 4000

	# pull out the fit parameters after burnin
	t0 = np.median(M.trace('switchpoint')[burnin:])
	early_mu = np.median(M.trace('early_mu')[burnin:])
	late_mu = np.median(M.trace('late_mu')[burnin:])
	slowness = np.median(M.trace('slowness')[burnin:])
	print "t0 = %f pm %f" % (t0, np.std(M.trace('switchpoint')[burnin:]))
	print "early_mu = %f pm %f" % (early_mu, np.std(M.trace('early_mu')[burnin:]))
	print "early_sig = %f pm %f" % (np.mean(M.trace('early_sigma')[burnin:]), np.std(M.trace('early_sigma')[burnin:]))
	print "late_mu = %f" % late_mu
	print "slowness = %f" % slowness

	# Recreate the best-fit sample to plot it

	t0_us = float(t0 * 10)

	tse = np.arange(0, int(t0_us), 10)
	tsl = np.arange(int(t0_us), len(signalWin) * 10, 10)
	timestamps = np.append(tse, tsl)

	e = np.ones(int(t0)) * early_mu  # early
	l = np.ones(signalWin.size - int(t0)) * late_mu  # late

	decay_late = exp_decay(l, tsl, t0_us, 72000.0)  # exp decay from rc time const (72000 ns = 72 us)
	s = np.append(e, decay_late)
	sigSmooth = gaussian_filter(s, sigma=slowness)

	plt.figure(1, figsize=(10, 7), facecolor='w')
	plt.plot(timestamps, sigSmooth, color="blue")
	plt.plot(timestamps, signalWin, color="red")
	plt.savefig("toy_fit.pdf")

	# What's the fit look like against the whole waveform?

	true_t0 = float(t0_us + loWin * 10.)
	print "true t0 : %i ns" % true_t0

	# timestamps
	tsfull_e = np.arange(0, int(true_t0), 10)
	tsfull_l = np.arange(int(true_t0), len(signalSubt) * 10, 10)
	tsfull = np.append(tsfull_e, tsfull_l)

	# fit signal
	full_e = np.ones(int(t0) + loWin) * early_mu
	# full_l = np.ones(signalSubt.size - (int(t0) + loWin)) * late_mu
	full_l = exp_decay(late_mu, tsfull_l, true_t0, 72000.0)
	sim = gaussian_filter(np.append(full_e, full_l), sigma=slowness)

	plt.figure(2, figsize=(10, 7), facecolor='w')
	plt.plot(tsfull, sim, color="red")
	plt.plot(tsfull, signalSubt, color="blue")

	# plt.clf()
	plt.show()

	# Histograms of posteriors

	# early_mu_arr = M.trace('early_mu')[burnin:]
	# early_sigma_arr = M.trace('early_sigma')[burnin:]
	#
	# n_bins = 50
	#
	#
	# plt.figure(1)
	# weights = np.ones_like(early_mu_arr) / float(len(early_mu_arr))
	# n, bins, patches = plt.hist(early_mu_arr, n_bins, histtype='step', linewidth=5, weights=weights)
	# plt.xlabel("mu_0 value")
	# plt.ylabel("probability")
	# plt.savefig("early_mu_pdf.pdf")

	# plt.figure(2)
	# weights = np.ones_like(early_sigma_arr) / float(len(early_sigma_arr))
	# n, bins, patches = plt.hist(early_sigma_arr, n_bins, histtype='step', linewidth=5, weights=weights)
	# plt.xlabel("sigma_0 value")
	# plt.ylabel("probability")
	# plt.savefig("early_sigma_pdf.pdf")

	# MCMC full trace
	# f, axarr = plt.subplots(5, sharex=True)
	# axarr[0].plot(M.trace('switchpoint')[:])
	# axarr[4].set_xlabel('MCMC Step Number')
	# axarr[0].set_ylabel('t_0')
	# axarr[1].plot(M.trace('early_mu')[:])
	# axarr[1].set_ylabel('mu_0')
	# axarr[2].plot(M.trace('late_mu')[:])
	# axarr[2].set_ylabel('mu_1')
	# axarr[3].plot(M.trace('early_sigma')[:])
	# axarr[3].set_ylabel('sigma_0')
	# axarr[4].plot(M.trace('late_sigma')[:])
	# axarr[4].set_ylabel('sigma_1')
	# plt.savefig("MCMC_steps.pdf")
	#
	# # Zoom in to see the convergence
	# f, axarr = plt.subplots(5, sharex=True)
	# axarr[0].plot(M.trace('switchpoint')[:])
	# axarr[4].set_xlabel('MCMC Step Number')
	# axarr[0].set_ylabel('t_0')
	# axarr[1].plot(M.trace('early_mu')[:])
	# axarr[1].set_ylabel('mu_0')
	# axarr[2].plot(M.trace('late_mu')[:])
	# axarr[2].set_ylabel('mu_1')
	# axarr[3].plot(M.trace('early_sigma')[:])
	# axarr[3].set_ylabel('sigma_0')
	# axarr[4].plot(M.trace('late_sigma')[:])
	# axarr[4].set_ylabel('sigma_1')
	# plt.xlim(0, 1000)
	# plt.savefig("MCMC_steps_zoom.pdf")
	#
	# plt.show()

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


def generateWaveform():
	"""Randomly generate a signal to fit to."""
	true_early_mean = 0
	true_late_mean = 2
	true_early_sigma = 1
	true_late_sigma = 0.75
	switchpoint = 133
	early_signal = np.random.normal(true_early_mean, true_early_sigma, switchpoint)
	late_signal = np.random.normal(true_late_mean, true_late_sigma, 200 - switchpoint)
	signal = np.append(early_signal, late_signal)
	return signal


def gauss_function(x, a, x0, sigma):
	"""."""
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def exp_decay(a, t, t0, tau):
	"""."""
	return a * np.exp(-(t - t0) / tau)


def subtractBaseline(signalRaw):
	"""Find the average starting baseline from an MJD waveform."""

	(hist, bins) = np.histogram(signalRaw, bins=np.arange(-8000, 8000, 1))
	fitfunc = lambda p, x: p[0] * np.exp(-0.5 * ((x - p[1]) / p[2])**2) + p[3]
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	c, pcov = curve_fit(gauss_function, bins[:-1], hist, p0=[np.amax(hist), bins[np.argmax(hist)], 5])

	print "Fit Coefficients:"
	print c[0], c[1], c[2]
	mu = c[1]
	std = c[2]
	print "mu is %f, std is %f" % (mu, std)

	return mu


if __name__ == "__main__":
	main()
