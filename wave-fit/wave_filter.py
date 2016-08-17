#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy import signal

def main():
	"""Look at the effects of different filters on an MJD waveform."""
	waves = TChain("waves")
	waves.Add("WSR_WFOutput.root")

	# waveEdge: windowed and baseline-subtracted, with metadata from GAT
	waveRaw = getWaveform(waves,0)
	baseAvg, noiseAvg = findBaseline(waveRaw)
	waveBLSub = waveRaw
	waveBLSub[:] = [x - baseAvg for x in waveRaw]  # a python 'list comprehension'
	zeros = np.asarray(np.where(abs(waveBLSub) < 0.1))
	lastZero = zeros[0,-1] # find last zero
	loWin, hiWin = lastZero-50, lastZero+200 # 50 and 200
	waveEnds = np.concatenate((np.arange(0, loWin), np.arange(hiWin, 2017)), axis=0)
	waveEdge = np.delete(waveBLSub, waveEnds)
	waveInfo = {
		'run': waves.run,
		'iEvent': waves.iEvent,
		'channel': waves.channel,
		'trapENFCal': waves.trapENFCal,
		'trapECal': waves.trapECal,
		't0': waves.blrwfFMR50,
		'noiseAvg': noiseAvg,
		'baseAvg': baseAvg
		}

	# ========= filtering =========

	# scipy.ndimage gaussian filter
	waveGaus = gaussian_filter(waveEdge, 3.)

	# 'filtfilt': an IIR filter with zero phase lag.
	# http://scipy-cookbook.readthedocs.io/items/FiltFilt.html
	b, a = signal.butter(3, 0.05)
	waveFilt = signal.filtfilt(b, a, waveEdge)


	# ========= plotting =========

	plt.style.use('mjWaveforms')  # using stylesheet in ~/.matplotlib/stylelib
	fig = plt.figure(figsize=(10, 7), facecolor='w')
	fig.add_subplot(111)
	ts = np.arange(loWin*10, hiWin*10, 10)  # (start,stop,step) in ns
	plt.xlim(loWin*10, hiWin*10)  # remove last sample at 2017
	plt.ylabel('ADC')
	plt.xlabel('time [ns]')
	plt.plot(ts, waveEdge, 'b', alpha=0.75)

	# gaussian_filter
	plt.plot(ts, waveGaus,'y')

	# filtfilt
	plt.plot(ts, waveFilt,'g')

	# plt.legend(('noisy signal', 'filtfilt'), loc='best')

	# draw the t0 calculated by GAT
	axes = plt.gca()
	plt.plot((waveInfo['t0'], waveInfo['t0']), axes.get_ylim(), 'r-', linewidth=1) # draw line at t0
	# plt.plot((lastZero*10, lastZero*10), axes.get_ylim(), 'r-', linewidth=1) # draw line at lastZero

	plt.show()

	# plt.savefig("fig.pdf")

# ===========================================================================================

def getWaveform(waves, i):
	"""Grab an MJ waveform. Convert a vector<double>
	to a numpy array, and remove the last entry."""
	waves.GetEntry(i)
	waveform = waves.wave.GetVectorData()  # vector<double>
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
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


if __name__ == "__main__":
	main()
