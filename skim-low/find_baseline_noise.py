#!/usr/local/bin/python
from ROOT import *
TROOT.gApplication.ExecuteFile("$MGDODIR/Root/LoadMGDOClasses.C")
TROOT.gApplication.ExecuteFile("$MGDODIR/Majorana/LoadMGDOMJClasses.C")

import sys, os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

gatDataName = "mjd_run"
gatTreeName = "mjdTree"
builtDataName = "OR_run"
builtTreeName = "MGTree"
dataSetName = "surfmjd"
detectorName = "P3KJR"


flatTimeSamples = 1500

####################################################################################################################################################################

chanDict = {584:"P1D1",582:"P1D2",580:"P1D3",578:"P1D4",692:"P2D1",648:"P2D2",640:"P2D3",642:"P2D4",616:"P3D1",610:"P3D2",608:"P3D3",664:"P3D4",624:"P4D1",628:"P4D2",688:"P4D3",694:"P4D4",614:"P4D5",680:"P5D1",678:"P5D2",672:"P5D3",696:"P5D4",632:"P6D1",630:"P6D2",626:"P6D3",690:"P6D4",600:"P7D1",598:"P7D2",594:"P7D3",592:"P7D4"}

def main(argv):
  runRange = (9914,9914)

  plt.ion()
  fig = plt.figure(1)

  #i do this one run at a time, instead of in a chain, because it makes it easier when we want to run on large data sets and create skim files for each run
  for iRun in range( runRange[0],  runRange[1]+1):
    print 'processing run', iRun
    gatFilePath =  os.path.expandvars("$MJDDATADIR/%s/data/gatified/%s/%s%d.root" % (dataSetName, detectorName, gatDataName, iRun  ) )
    builtFilePath =  os.path.expandvars("$MJDDATADIR/%s/data/built/%s/%s%d.root" % (dataSetName, detectorName, builtDataName, iRun  ) )

    if not os.path.isfile(gatFilePath):
      print ">>>Skipping file " + gatFilePath
      continue

    gat_file = TFile.Open(gatFilePath)
    gatTree = gat_file.Get(gatTreeName)
    built_file = TFile.Open(builtFilePath)
    builtTree = built_file.Get(builtTreeName)

    builtTree.AddFriend(gatTree)

    for channelNumber in chanDict.keys():
      chanHist = np.zeros_like(np.arange(-8000, 8000, 1)[:-1])

      chanCut =  "channel == %d" % channelNumber
      eCut = " energy > %d" % 1E5
      cut = chanCut + " && " + eCut
      #print "The cuts will be: " + cut

      gatTree.SetEntryList(0)
      gatTree.Draw(">>elist%d" % channelNumber, cut, "entrylist%d" % channelNumber)
      elist = gDirectory.Get("elist%d" % channelNumber)
      print "Number of entries in the entryList is " + str(elist.GetN())
      gatTree.SetEntryList(elist);
      builtTree.SetEntryList(elist);

      maxEntries = 1000
      numEntries = elist.GetN()

      eNum = np.amin((maxEntries, numEntries))

      if eNum == 0:
        print "No events for Detector %s (chan %d)" % (chanDict[channelNumber], channelNumber)
        continue

      for ientry in xrange( eNum ):
        update_progress(ientry/float(eNum))
        entryNumber = gatTree.GetEntryNumber(ientry);
        waveform = getWaveform(gatTree, builtTree, entryNumber, channelNumber)

        (hist, bins) = np.histogram(waveform, bins=np.arange(-8000, 8000, 1))
        chanHist += hist
  #      plt.clf()
  #
  #      value = raw_input('  --> Press q to quit, any other key to continue\n')
  #      if value == 'q':
  #        exit(1)

      plt.clf()

      sum = np.sum(chanHist)

      fitfunc  = lambda p, x: p[0]*exp(-0.5*((x-p[1])/p[2])**2)+p[3]
      errfunc  = lambda p, x, y: (y - fitfunc(p, x))

      c, pcov = curve_fit(gauss_function, bins[:-1], chanHist, p0 = [np.amax(chanHist), bins[np.argmax(chanHist)], 5])

      print "Fit Coefficients:"
      print c[0],c[1],c[2]

      mu = c[1]
      std = c[2]
      print "mu is %f, std is %f" % (mu, std)

      nonzero = np.argwhere(chanHist) #index number of bin
      nonzero_idx_low = np.amin(nonzero)
      nonzero_idx_hi = np.amax(nonzero)

      plt.xlim(bins[nonzero_idx_low], bins[nonzero_idx_hi])
      plt.plot(bins[:-1], chanHist, color="b")
      plt.plot(np.linspace(bins[nonzero_idx_low], bins[nonzero_idx_hi], 1000), gauss_function(np.linspace(bins[nonzero_idx_low], bins[nonzero_idx_hi], 1000), c[0], c[1], c[2]), color="r")

      plt.title("Baseline of %s (mean %0.2f, std %0.2f)" %   (chanDict[channelNumber], mu, std))
      plt.xlabel("Baseline ADC value")
      plt.ylabel("Counts (%d events, %d samples each" % (eNum, flatTimeSamples ))
      plt.savefig("baselinehist_det%s_chan%d_std%f.pdf" %(chanDict[channelNumber], channelNumber, std))

#      value = raw_input('  --> Press q to quit, any other key to continue\n')
#      if value == 'q':
#        exit(1)

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
####################################################################################################################################################################

def getWaveform(gatTree, builtTree, entryNumber, channelNumber):

    builtTree.GetEntry( entryNumber )
    gatTree.GetEntry( entryNumber )

    event = builtTree.event
    channelVec   = gatTree.channel
    energyVec   = gatTree.energy
    numWaveforms = event.GetNWaveforms()

    for i_wfm in xrange( numWaveforms ):
        channel = channelVec[i_wfm]
        e = energyVec[i_wfm]
        #print "---> channel is %d, energy %f" % (channel, e)
        if (channel != channelNumber): continue
        wf = event.GetWaveform(i_wfm)
        np_data = wf.GetVectorData()
        return np_data[:flatTimeSamples]

####################################################################################################################################################################

#I got impatient with loading longer runs, so this prints the progress.
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,2) , status)
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__=="__main__":
    main(sys.argv[1:])

