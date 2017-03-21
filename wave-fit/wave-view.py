#!/usr/local/bin/python
import sys
from ROOT import TFile,TTree,TEntryList,gDirectory,gROOT
import waveLibs as wl
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    """Interactive-draw or rapid-draw waveforms that pass a given TCut."""

    plt.ion()
    opt1, opt2 = "", ""
    intMode = False
    printWaveforms = False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-s" in (opt1, opt2):
        printWaveforms = True
        print "Saving WF plots to current directory."
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    # Set input file and cuts

    fileName = "./data/waveletSkimDS4.root"
    inputFile = TFile(fileName)
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    # theCut = inputFile.Get("cutUsedHere").GetTitle()
    theCut = "trapENFCal > 0 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    # S5
    # theCut += " && waveS5 > 4500 && waveS5 < 6000 && trapENFCal > 30 && trapENFCal < 50"  # DS4 blob #1
    # theCut += " && waveS5 > 1500 && waveS5 < 2500 && trapENFCal > 20 && trapENFCal < 30"  # DS4 blob #2
    # theCut += " && waveS5 > 1500 && waveS5 < 2500 && trapENFCal > 30 && trapENFCal < 40"  # DS4 blob #3
    # theCut += " && waveS5 > 1500 && waveS5 < 7000 && trapENFCal > 0 && trapENFCal < 5"    # DS4 blob #4
    # theCut += " && waveS5 < 1500 && trapENFCal > 0 && trapENFCal < 5" # what waveS5 would accept

    # S32/E
    # theCut += " && (waveS3-waveS2)/trapENFCal > 400 && (waveS3-waveS2)/trapENFCal < 600 && trapENFCal > 20 && trapENFCal < 30"
    # theCut += " && (waveS3-waveS2)/trapENFCal > 150 && (waveS3-waveS2)/trapENFCal < 320 && trapENFCal > 1.5 && trapENFCal < 5"

    # S34/E
    # theCut += " && (waveS3-waveS4)/trapENFCal > 25 && (waveS3-waveS4)/trapENFCal < 60 && trapENFCal > 0.8 && trapENFCal < 10"
    theCut += " && (waveS3-waveS4)/trapENFCal < 0 && trapENFCal > 0.8 && trapENFCal < 5"

    print "Using cut:\n",theCut,"\n"

    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."


    # Make a figure (only setting data in the loop is faster)
    # (With some more work, it could go much faster --
    # http://bastibe.de/2013-05-30-speeding-up-matplotlib.html )
    # See also clint's 'speed-plotter.py'

    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a2 = plt.subplot(111)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue')
    p2, = a2.plot(np.ones(1), np.ones(1), color='red', linewidth=2., alpha=0.5)


    # Loop over events
    iList = -1
    while(True):
        iList += 1
        plt.pause(0.00001)
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWaves = waveTree.MGTWaveforms.size()
        numPass = waveTree.Draw("channel","","GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH))
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

            # fill the figure
            p1.set_ydata(waveBLSub)
            p1.set_xdata(waveTS)
            p2.set_ydata(waveFilt)
            p2.set_xdata(waveTS)
            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
            a1.set_xlim([xmin,xmax])
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            plt.title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,energy))
            # plt.title(title,loc='left')
            # fig.canvas.draw()

            if (printWaveforms):
                plt.savefig("./plots/wave-%d-%d-%d.pdf" % (run,iList,chan))

if __name__ == "__main__":
    main(sys.argv[1:])
