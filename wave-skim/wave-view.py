#!/usr/local/bin/python
import sys
from ROOT import TFile,TTree,TEntryList,gDirectory,gROOT
import wave_model as wm
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    """Interactive-draw or rapid-draw waveforms that pass a given TCut."""
    opt1, opt2 = "", ""
    setInteractive = False
    printWaveforms = False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-I" in (opt1, opt2):
        setInteractive = True
        print "Interactive mode selected."
    if "-S" in (opt1, opt2):
        printWaveforms = True
        print "Saving WF plots to current directory."
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    fileName = "./data/waveSkimDS1_test.root"
    inputFile = TFile(fileName)
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = ""
    theCut = inputFile.Get("theCut").GetTitle()
    # theCut = "channel != 580"
    print "Using cut:\n",theCut,"\n"

    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    print "Found",elist.GetN(),"entries passing cuts."

    scanSpeed = 0.5 # sec between plots
    plt.ion()
    fig = plt.figure(figsize=(9,5), facecolor='w')
    ax1 = plt.subplot(111)
    ax2 = plt.subplot(111)
    ax1.set_xlabel("time (ns)")
    ax1.set_ylabel("ADC")
    dummy = np.arange(0,2016*10,10)
    wave, = ax1.plot(dummy, np.ones(len(dummy)))
    waveF, = ax2.plot(dummy, np.ones(len(dummy)), color='red', linewidth=2., alpha=0.5)

    iList = -1
    while(True):
        plt.pause(scanSpeed/2)
        if iList >= elist.GetN()-1: exit(1)
        if setInteractive==True:
            iList += 1
            value = raw_input()
            if value=='q': exit(1)     # exit
            if value=='p': iList -= 2  # previous entry
            if (value.isdigit()):
               iList = int(value)      # go to entry number
        else:
            iList += 1

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWaves = waveTree.MGTWaveforms.size()

        # only get channels that pass the cut
        numPass = waveTree.Draw("channel","","GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        print "Entry",iList,"nCh",nChans,"nWF",nWaves,chanList

        for iHit in xrange(nChans):

            if waveTree.channel.at(iHit) in chanList:

                print "  chan %d  trapENF %.1f" % (waveTree.channel.at(iHit), waveTree.trapENFCal.at(iHit))

                signal = wm.processWaveform(waveTree.MGTWaveforms.at(iHit))
                waveBLSub = signal.GetWaveBLSub()
                waveFilt = signal.GetWaveFilt()
                waveTS = signal.GetTS()

                xmin, xmax = np.amin(waveTS), np.amax(waveTS)
                ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
                wave.set_ydata(waveBLSub)
                wave.set_xdata(waveTS)
                ax1.set_xlim([xmin,xmax])
                ax1.set_ylim([ymin-0.1*ymin,ymax+0.1*ymax])

                waveF.set_ydata(waveFilt)
                waveF.set_xdata(waveTS)
                ax2.set_xlim([xmin,xmax])
                ax2.set_ylim([ymin-0.1*ymin,ymax+0.1*ymax])

                plt.pause(scanSpeed/2)
                if (setInteractive): plt.pause(1.0)

                title = "Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (waveTree.run,waveTree.channel.at(iHit),iList,waveTree.trapENFCal.at(iHit))
                plt.title(title,loc='left')
                fig.canvas.draw()

                if (printWaveforms):
                    plt.savefig("waveforms/wave-%d-%d-%d.pdf" % (waveTree.run,iList,waveTree.channel.at(iHit)))


if __name__ == "__main__":
    main(sys.argv[1:])
