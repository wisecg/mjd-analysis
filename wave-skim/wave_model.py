#!/usr/local/bin/python
import numpy as np
from scipy.optimize import curve_fit
from pymc import TruncatedNormal, Normal, HalfNormal, deterministic, Uniform
from scipy.ndimage.filters import gaussian_filter
from scipy import signal as sg
import sys

def WaveformModel(waveEdge, waveEdgeTS, tempBLSub, tempTS, diff, amp, noiseAvg):

    # Old model:
    # switchpoint = Normal('switchpoint', mu=diff, tau=1.0, value=diff)
    # scale = Normal('scale', mu=amp, tau=sigToTau(1.0 * amp), value=amp)
    # noise = np.power(noiseAvg, -2)
    # slowness = HalfNormal('slowness', tau=0.001)
    # slowness = Normal('slowness', tau=0.001)

    # NOTE:
    # - tau needs to be whatever the sigma of the baseline noise in ADC counts is.
    # - tau is 1/sigma**2 for the normal distribution you're modeling
    #   ex. for a 10ns (1 sample) uncertainty in time (sigma=1), tau = 1/(1)^2 = 1
    #   ex. for a 1us (100 ns = 10 sample) uncertainty (sigma=10), tau = 0.01
    # - noiseAvg needs to be the sigma of the baseline RMS.
    # - we use TruncatedNormal so we're not anchored to 0

    # Current model:
    switchpoint = Normal('switchpoint', mu=diff, tau=1.0, value=diff)
    scale = Normal('scale', mu=amp, tau=sigToTau(1.0 * amp), value=amp)
    # noise = np.power(noiseAvg, -2)
    slowness = TruncatedNormal('slowness', mu=1., tau=0.01, a=0., b=np.inf)

    @deterministic(plot=False, name="depModel")
    def SignalModel(s=switchpoint, e=scale, sig=slowness):

        lo = waveEdgeTS[0]
        hi = waveEdgeTS[-1]+10
        idx = np.where((tempTS + s >= lo) & (tempTS + s <= hi))

        if (len(tempTS[idx]) != len(waveEdgeTS)):
            print "lengths don't match!",len(tempTS[idx]),len(waveEdgeTS)
            exit(1)

        out = e * tempBLSub[idx]
        blur = gaussian_filter(out, sigma=float(sig))

        return blur

    try:
        baseline_observed = Normal("baseline_observed", mu=SignalModel, tau=noiseAvg, value=waveEdge, observed=True)
    except:
        print "damn thing failed for some reason:",sys.exc_info()[0]
        return locals(),False

    return locals(),True

def sigToTau(sig):
    tau = np.power(sig, -2)
    return tau

def gauss_function(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def getWFArray(waveRaw):

    return npArr_copy

def rootToArray(hist, xLow, xHigh):
    """Take a ROOT TH1D histogram and a range, get a numpy array.
    Note on plotting:
    Can't use np.histogram, the hist has already been done.
    So make a plot that fakes it with :
        xaxis = np.arange(xlow, xhi, step)
        plt.bar(xaxis,hist,width=step)  # 1: fake it with bar
        plt.plot(xaxis,hist,ls='steps-post')  # 2: fake it with plot
    """
    binLow = hist.FindBin(xLow)
    binHigh = hist.FindBin(xHigh)
    loopArray = range(binLow, binHigh )
    npArray = np.empty_like(loopArray, dtype=np.float)
    for (iArray, iBin) in enumerate(loopArray):
        npArray[iArray] = np.float(hist.GetBinContent(iBin))
    return npArray


def findBaseline(signalRaw):
    """Find the average starting baseline from an MJD waveform."""
    (hist, bins) = np.histogram(signalRaw[:200], bins=np.arange(-8000, 8000, 1))
    fitfunc = lambda p, x: p[0] * np.exp(-0.5 * ((x - p[1]) / p[2])**2) + p[3]
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    c, pcov = curve_fit(gauss_function, bins[:-1], hist, p0=[np.amax(hist), bins[np.argmax(hist)], 5])
    mu = c[1]
    std = c[2]
    return mu, std

class processWaveform:
    def __init__(self, wave, removeNSamples=2, N=1, Wn=0.08, lo=50, hi=400):

        self.waveMGT = wave                            # input an MGTWaveform object
        self.offset = self.waveMGT.GetTOffset()        # time offset [ns]
        vec = wave.GetVectorData()
        npArr = np.fromiter(vec, dtype=np.double)      # raw numpy array

        hist = wave.GimmeUniqueHist()                  # get timestamp limits and make an array
        self.start = hist.GetXaxis().GetXmin() + 5     # add 5 ns to make it start at 0 (default is -5.)
        self.stop = hist.GetXaxis().GetXmax() + 5.

        self.binsPerNS = (self.stop - self.start) / hist.GetNbinsX()

        ts = np.arange(self.start,self.stop,self.binsPerNS)
        removeSamples = [npArr.size - i for i in xrange(1,removeNSamples+1)]
        self.ts = np.delete(ts,removeSamples)
        self.waveRaw = np.delete(npArr,removeSamples)   # force the size of the arrays to match

        self.baseAvg, self.noiseAvg = findBaseline(self.waveRaw)     # get baseline and baseline RMS
        self.waveBLSub = self.waveRaw
        self.waveBLSub[:] = [x - self.baseAvg for x in self.waveRaw] # subtract the baseline value

        self.b, self.a = sg.butter(N, Wn)
        self.waveFilt = sg.filtfilt(self.b, self.a, self.waveBLSub)  # do a smoothing

        # self.zeros = np.asarray(np.where(self.waveBLSub < 0.1))
        zeros = np.asarray(np.where(self.waveFilt < 0.01))      # make a list of negative points ("zeros")
        self.lastZero = 0
        if (zeros.size > 0):
            self.lastZero = zeros[0,-1]                         # find the last "zero crossing"
        self.lastZeroTime = ts[self.lastZero]

        self.loWin, self.hiWin = self.lastZero-lo, self.lastZero+hi  # indexes, not times
        waveEnds = np.concatenate((np.arange(0, self.loWin), np.arange(self.hiWin, self.waveRaw.size)), axis=0)

        self.waveEdge = np.delete(self.waveBLSub, waveEnds)     # rising edge of the waveform
        self.tsEdge = np.delete(self.ts, waveEnds)         # rising edge timestamps

    # constants
    def GetOffset(self): return self.offset
    def GetStartTime(self): return self.start
    def GetStopTime(self): return self.stop
    def GetBins(self): return self.binsPerNS
    def GetBaseNoise(self): return self.baseAvg, self.noiseAvg
    def GetWindowIndex(self): return self.loWin, self.hiWin

    # arrays
    def GetWaveEdge(self): return self.waveEdge
    def GetTSEdge(self): return self.tsEdge
    def GetTS(self): return self.ts
    def GetWaveRaw(self): return self.waveRaw
    def GetWaveBLSub(self): return self.waveBLSub
    def GetWaveFilt(self): return self.waveFilt
    def GetLastZeroIndex(self): return self.lastZero
    def GetLastZeroTime(self): return self.lastZeroTime