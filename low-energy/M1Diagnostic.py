#!/usr/bin/python
#$ -V
#$ -j y
#$ -l gscratchio=1
#$ -l  projectio=1
#$ -o /project/projectdirs/majorana/users/krvorren/rateOutput
########################
# This script is designed to create some diagnostics for M1 to
# help determine good runs to use analysis.
#
# Created by Kris Vorren, krisvorren@unc.edu
#
# -Version History/Notes-
# Jul 31, 2015: Created.
# Sep 16, 2015: Diagnostic over range of runs used for the ring remover and more C1P4 only, also removing pulser tails
# Sep 17, 2015: All the detectors, no cuts initially, some histograms with cuts, new En v Run with cuts
# Sep 21, 2015: Less runs, actually include BG runs AFTER the external pulsing.
# Sep 22, 2015: There was an error with the trapMinCut. It was fixed and rerun
# Oct 3, 2015: No Timestamps, only cut T/E and cut Energy, all channels
# Oct 4, 2015: No triangle mins, look for coincidence with P6D1
# Oct 7, 2015: Plot singles
# Oct 8, 2015: Used smooth T/E
# Oct 11, 2015: Improved Granularity cut
# Oct 12, 2015: Only one channel events pass (noise)
# Oct 13, 2015: Only two channel events, remove bad runs from T/E (runs up to 5501)
# Oct 16, 2015: Only two channel events, fixed exposure
# Oct 19, 2015: Redid exposure (it was undercounting due to loop continue statement)
# Oct 21, 2015: All BG runs (couldn't get to work memory problem)
# Nov 3, 2015: Added the offset dictionary, runs to 5501
# Jan 22, 2016: fixed cal ext pul, removed singles condition and trap Min*****
# Jan 26, 2016: run 2368 - 5501: Compare toe with and w/o min cut
# Jan 27, 2016: same runs as yesterday, use trapEtailMin, only well behaved detectors
# Jan 28, 2016: test cut with min val
# Jan 29, 2016: Test cut seemed to work, now try to see what smooth T/E looks like (CHANGED T/E HIST LIMITS******)...
# Feb 2, 2016: All BG runs AFTER threshold lowered...
# Feb 3, 2016: Neglect range of funs above 5800 with bad noise issues...
# Feb 4, 2016: Need to break into two runs to get all the files... part 1 to run 4831, part 2 to 6963
# Feb 5, 2016: Had to add the min run as well
########################
import sys, os, re, glob, ROOT, array, numpy
ROOT.gROOT.SetBatch(True)

#----INITIALIZATION----
#Load some MGDO classes
ROOT.gROOT.Reset()
ROOT.gApplication.ExecuteFile("%s/Root/LoadMGDOClasses.C" %os.environ['MGDODIR'])
from ROOT import CLHEP

#Simple wait for return fuction
def WaitForReturn():
    val = raw_input("Hit enter or 'q' to quit:   ")
    if val == 'q':  sys.exit()

#Status bar output
def StatusBar(i, total):
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %.2f%%" % ('='*(100*(i+1)/total/5), 100.*(i+1)/total))
    sys.stdout.flush()
    if(i==total):print ""

##~~~~~~~~~~~~~~~MODULE 1~~~~~~~~~~~~~~~~
channelDict = {692.0:"P1D1", 690.0:"P1D2", 688.0:"P1D3", 640.0:"P1D4",
               584.0:"P2D1", 674.0:"P2D2", 576.0:"P2D3", 680.0:"P2D4",
               676.0:"P3D1", 616.0:"P3D2", 614.0:"P3D3", 610.0:"P3D4",
               608.0:"P4D1", 598.0:"P4D2", 600.0:"P4D3", 594.0:"P4D4", 592.0:"P4D5",
               664.0:"P5D1", 662.0:"P5D2", 656.0:"P5D3", 696.0:"P5D4",
               628.0:"P6D1", 626.0:"P6D2", 624.0:"P6D3", 632.0:"P6D4",
               646.0:"P7D1", 644.0:"P7D2", 642.0:"P7D3", 630.0:"P7D4"
               }
channelList = [692.0, 690.0, 688.0, 640.0,
               584.0, 674.0, 576.0, 680.0,
               676.0, 616.0, 614.0, 610.0,
               608.0, 598.0, 600.0, 594.0, 592.0,
               664.0, 662.0, 656.0, 696.0,
               628.0, 626.0, 624.0, 632.0,
               646.0, 644.0, 642.0, 630.0]

channelList = [692.0, 690.0, 688.0, 640.0,
                      674.0, 576.0,
                                    610.0,
               608.0, 598.0, 600.0, 594.0, 592.0,
               664.0, 662.0,        696.0,
                      626.0, 624.0, 
               646.0, 644.0, 642.0        ]

offsetDict = { #These were generated from the M1zeroCal.py script. Basically finds the mean of a gaussian peak from forced-acq data
    692.0: 0.318, 690.0: 0.124, 688.0: 0.469, 640.0: 0.432,
    584.0: 0.000, 674.0: -.022, 576.0: 0.075, 680.0: 0.000,
    676.0: 0.000, 616.0: 0.000, 614.0: 0.448, 610.0: -0.013,
    608.0: 0.381, 598.0: 0.322, 600.0: -.209, 594.0: 0.617, 592.0: 0.288,
    664.0: 0.474, 662.0: 0.159, 656.0: 0.196, 696.0: 0.449,
    628.0: 0.587, 626.0: -.181, 624.0: -.153, 632.0: 0.000,
    646.0: 0.022, 644.0: 0.091, 642.0: 0.009, 630.0: 0.000
}

#channelList = [608.0, 598.0, 600.0, 594.0, 592.0]
#channelList = [674.0, 624.0, 688.0, 662.0, 614.0, 608.0] #channels for external pulsing

#********Directory Setup*******    
#FOR PDSF
#serialNumberList = ["P3GKF/", "P3HUA/", "P3JCJ/", "P3JDY/"]
serialNumberList = ["P3JDY/"]

builtDirectory = "%s/surfmjd/data/built/" %os.environ['MJDDATADIR']
builtPrefix = "OR_run"
gatDirectory = "%s/surfmjd/data/gatified/" %os.environ['MJDDATADIR']
gatPrefix = "mjd_run"
outputFilePath = "%s/ENPA/project/m1Diagnostic4Feb2016_1.root" %os.environ['HOME']

#***runNumber is now a dictionary***
runNumber = {}
runNumber["P3GKF/"] = (
    range(1, 111 + 1) #For all the runs
    )
runNumber["P3HUA/"] = (
    #First run in here starts at 202...
    range(202, 1041 + 1) #For all the runs
    )
runNumber["P3JCJ/"] = (
    #First run is...
    range(1042, 2334 + 1) #For all the runs
    )
runNumber["P3JDY/"] = (
    #***ALL THE RUNS***
    ## range(2337, 4415 + 1) #For all the runs
    #***BACKGROUND RUNS***
    #range(2339,2360+1) + #medium thresholds
    #range(2368, 2440+1) + 
    #range(2547, 2570+1) +
    [2579] +
    range(2580, 2581) + #Run 2581 is a short garbage run
    range(2582, 2629+1) + #with pulsers on and raised thresholds on channels 6,7
    range(2644, 2649+1) + #with pulsers on and raised thresholds on channels 6,7
    range(2658, 2673+1) +
    range(2688, 2920 + 1) + #Wenqin says the Rn purge and seal is optimzed starting this run range,
    #Starting with 2924, the cryopump was turned off (running on the turbo). 
    range(3125, 3129+1) +
    range(3137,3271+1) + #Background runs on M1. BB-decay run bit was not set for first part
    #range(3293,3432+1) + #shadow shield is installed
    range(3293, 3432) +
    #range(3342, 3432) +
    range(3461,3462+1) + #Stopped for setting up the pulser positions
    range(3464,3556+1) + #"Really First background data after pulser adjustments. All pulsers should be at 640keV
    #- Matt opened valve to cryopump on Module 1 vacuum system, which caused it to start making a knocking noise. "
    range(3557,3580+1) + #Pulsers on/off, alternating between runs.
    range(3596,3645+1) +  #Pulsers on/off, alternating between runs. Restart of alternating pulsers on/off runs after rebooting SBC to fix bus error problems. 
    range(4034,4134+1) + #Pulsers on/off, alternating between runs. 
    range(4239, 4493+1) + #Re-enabled Card 11 with the correct VME firmware
    range(4497, 4518+1) + #After some tests, no veto
    range(4549, 4572+1) +#After some tests, no veto
    range(4573, 4831+1) + #veto is up and running
    range(4854, 4907+1) +
    range(4938, 4981+1) +
    range(5007, 5061+1) +
    range(5090, 5252+1) +
    range(5277, 5331+1) +
    range(5372, 5414+1) +
    range(5449, 5501+1) +
    
    range(5525, 5534+1) +
    range(5555, 5850+1) +
    ##range(5888, 5902+1) +
    ##range(5922, 5939+1) +
    ##range(6220, 6317+1) +
    ##range(6331, 6352+1) +
    range(6553, 6577+1) +
    range(6776, 6853+1) +
    range(6887, 6903+1) +
    range(6957, 6963+1) +
    [6965]
    #***QUESTIONABLE RUNS***
    #range(2335,2338+1) + #questionable
    #range(2682,2687+1) + #test/force acquision
    #range(2921,2930+1) + #The veto panel intallation may introduced some Rn
    #range(2975,3056+1) + #Wenqin says it is questionable since there is some Rn in the beginning of this period, due to veto panel installation. See 1.09 Elog #367
    #range(3078, 3124+1)  #Testing with all veto panels
    #range(3463,3463+1) + #First background runs after pulsers have been adjusted to all be at 640 keV.
    #range(3581,3595+1) + #"Starting with run 3581 (10am), the SBC started generating bus errors. At 1pm the even rates dropped to zero. These data are likely corrupted. 
    #Rebooting the SBC fixed the problem. "
    #range(3665,3687+1) + #detector biased down. E-box 2 opened (BAD)
    #range(4004,4033+1) + #bkg and test runs BAD
    #range(4171,4233+1) + #test/force acquision
    #range(4234,4238+1) + #Running 15-minute background runs with Card 11 (one with the incorrect VME firmware) removed.  This card's object has also been removed from ORCA.
    #range(4429,4429+1) + #missing runs from above
    #range(4433,4435+1) #missing runs from above
    #range(5873, 5886+1) #This is P2D2 calibration with high precision attenuated pulsers (200 ns)
    ## range(5947,5960+1) + #P2D2 calibration with 190 ns rt 300 s run len
    ## range(5964,5977+1) + #P6D3 external pulse
    ## range(5979,5992+1) + #P1D3 external pulse
    ## range(6191,6204+1) + #P5D2 external pulse
    ## range(6206,6219+1) #P3D3 external pulse
    ## range(7219, 7223+1) +
    ## range(7225, 7233+1) + #This and above are P2D2 ext pul 155 ns rt
    ## range(7234, 7246+1) + #P6D3 164ns rt
    ## range(7247, 7259+1) + #P1D3
    ## range(7260, 7272+1) #P5D2 138 ns
    )

badRuns = [2579, 2610, 3272, 3341, 3523, 3524, 3529, 3530, 4057, 4125, 4473, 4554, 5248, 5249, 5331, 6965]
#NOTE: The bad runs are only cut from the T/E

#FOR PERSONAL COMPUTER...
## serialNumberList = ["input/"]

## builtDirectory = "/Users/krvorren/Project/M1background/MJORData/"
## builtPrefix = "OR_run"
## gatDirectory = "/Users/krvorren/Project/M1background/GATifiedData/"
## gatPrefix = "mjd_run"
## outputFilePath = "/Users/krvorren/Project/M1background/output/DiagnosticTest.root"
## runNumber = {}
## runNumber["input/"] = [3000]

#************Chain Setup*************
builtFilesAdded = 0
gatFilesAdded = 0
chain = ROOT.TChain("mjdTree")
for sn in serialNumberList:
    for rn in runNumber[sn]:
        if(chain.Add(gatDirectory + sn + gatPrefix + str(rn) + ".root")): gatFilesAdded += 1

totalEntries = chain.GetEntries()
print "~~Added " + str(gatFilesAdded) + " built file(s) to the chain."
print "~~Total number of Entries is %i." %totalEntries

#*********Parameters of Interest*******
timestamp = ROOT.std.vector("double")()
calE = ROOT.std.vector("double")()
smoothT = ROOT.std.vector("double")()
unSmoothT = ROOT.std.vector("double")()
channel = ROOT.std.vector("double")()
theRun = numpy.zeros(1, dtype=float)
startTime = numpy.zeros(1, dtype=float)
stopTime = numpy.zeros(1, dtype=float)
trapMin = ROOT.std.vector("double")()
triMin = ROOT.std.vector("double")()
chi2 = ROOT.std.vector("double")()

chain.SetBranchStatus("*",0)
chain.SetBranchStatus("timestamp", 1)
chain.SetBranchStatus("trapECal", 1)
#chain.SetBranchStatus("energyCal", 1)
chain.SetBranchStatus("smoothTrirt100nsft10nsMax", 1)
chain.SetBranchStatus("trirt100nsft10nsMax", 1)
chain.SetBranchStatus("channel", 1)
chain.SetBranchStatus("run", 1)
chain.SetBranchStatus("startTime", 1)
chain.SetBranchStatus("stopTime", 1)
chain.SetBranchStatus("trapETailMin", 1)
chain.SetBranchStatus("triFilMin", 1)
chain.SetBranchStatus("RawWFblChi2", 1)

chain.SetBranchAddress("timestamp", timestamp)
chain.SetBranchAddress("trapECal", calE)
#chain.SetBranchAddress("energyCal", calE)
chain.SetBranchAddress("smoothTrirt100nsft10nsMax", smoothT)
chain.SetBranchAddress("trirt100nsft10nsMax", unSmoothT)
chain.SetBranchAddress("channel", channel)
chain.SetBranchAddress("run", theRun)
chain.SetBranchAddress("startTime", startTime)
chain.SetBranchAddress("stopTime", stopTime)
chain.SetBranchAddress("trapETailMin", trapMin)
chain.SetBranchAddress("triFilMin", triMin)
chain.SetBranchAddress("RawWFblChi2", chi2)

#*********Histogram Dictionary*********
histDict = {}
runMin = min(runNumber["P3JDY/"]) - 10
runMax = max(runNumber["P3JDY/"]) + 10
#runMin = 2500
#runMax = 4600
numBins = runMax - runMin
for iChan in channelList:
    #histDict["runTS"+str(iChan)] = ROOT.TH2F("runTS"+channelDict[iChan], "Run vs. Timestamp ("+channelDict[iChan]+");sec;Run", 3600, 0, 3600, numBins, runMin, runMax)
    #histDict["runEn"+str(iChan)] = ROOT.TH2F("energy"+channelDict[iChan], "Run vs. Energy ("+channelDict[iChan]+");Energy (keV);Run", 3000, 0, 100, numBins, runMin, runMax)
    histDict["runEnCut"+str(iChan)] = ROOT.TH2F("cutEnergy"+channelDict[iChan], "Run vs. Energy with Cuts ("+channelDict[iChan]+");Energy (keV);Run", 3000, 0, 100, numBins, runMin, runMax)
    #histDict["ToE"+str(iChan)] = ROOT.TH2F("ToE"+channelDict[iChan], "ToE vs. E ("+channelDict[iChan]+");Energy (keV);ToE", 3000, 0, 100, 1000, 0, 10)
    #Below are the new (maybe temporary) tools
    histDict["cutToE"+str(iChan)] = ROOT.TH2F("cutToE"+channelDict[iChan], "Clean ToE vs. E ("+channelDict[iChan]+");Energy (keV);ToE", 3000, 0, 100, 1000, 0, 10)
    #histDict["smoothToE"+str(iChan)] = ROOT.TH2F("smoothToE"+channelDict[iChan], "Cleaned Smoothed ToE vs. E ("+channelDict[iChan]+");Energy (keV);ToE", 3000, 0, 100, 2000, 0, 2)
    #histDict["chi2"+str(iChan)] = ROOT.TH2F("blRMS"+channelDict[iChan], "Baseline RMS vs E ("+channelDict[iChan]+");Energy (keV);RMS", 3000, 0, 100, 100, 0, 20)
histDict["exp"] = ROOT.TH2F("exposure", "Exposure by Channel; Channel; Run", 260, 570, 700, numBins, runMin, runMax)

#*********THE WORKING LOOP***********
for iEnt in xrange(chain.GetEntries()):
    #StatusBar(iEnt, chain.GetEntries())
    chain.GetEntry(iEnt)
    for i in xrange(channel.size()):
        axisXbin = histDict["exp"].GetXaxis().FindBin(channel.at(i))
        axisYbin = histDict["exp"].GetYaxis().FindBin(theRun[0])
        if stopTime[0] - startTime[0] > 0:
            histDict["exp"].SetBinContent(axisXbin, axisYbin, stopTime[0]-startTime[0])
    if channel.size() > 2: #noise triggers hi gain but not low gain, also want singles
        continue
    if channel.size() == 2:
        bigChannel = max([channel.at(0), channel.at(1)])
        if (bigChannel % 2 == 0): #lo gain channels are +1 hi gain channel, hi-gain channels are even numbers
            continue
        smallChannel = min([channel.at(0), channel.at(1)])
        if not (bigChannel - smallChannel == 1):
            continue
    for i in xrange(channel.size()):
        if(channel.at(i) in channelList):
            calibEn = calE.at(i) - offsetDict[channel.at(i)]
            if(calibEn > 0 and calibEn < 101):
                #histDict["runTS"+str(channel.at(i))].Fill(timestamp.at(i)/100000000, theRun[0])
                #histDict["runEn"+str(channel.at(i))].Fill(calE.at(i), theRun[0])
                #if not (theRun[0] in badRuns):
                #    histDict["ToE"+str(channel.at(i))].Fill(calibEn, unSmoothT.at(i)/calibEn)
                #if(trapMin.at(i) < 0):
                #blRMScut = -6.6 #NOTE, NOT ALL GATIFIED DATA HAS THIS CUT PARAMETER
                #histDict["chi2"+str(channel.at(i))].Fill(calE.at(i), blRMScut/-3)
                if(trapMin.at(i) < 0):
                    lowGoodFlag = True
                    if calibEn < 20:
                        if trapMin.at(i) < (-4.0/19.0*calibEn - 15.0/19.0):
                            lowGoodFlag = False
                    if lowGoodFlag:
                        histDict["runEnCut"+str(channel.at(i))].Fill(calibEn, theRun[0])
                        if not (theRun[0] in badRuns):
                            histDict["cutToE"+str(channel.at(i))].Fill(calibEn, unSmoothT.at(i)/calibEn)
                #histDict["smoothToE"+str(channel.at(i))].Fill(calE.at(i), smoothT.at(i)/calE.at(i))

#*********SAVE RESULTS***************
print ""
print "~Saving histograms to " + outputFilePath + "."
outFile = ROOT.TFile(outputFilePath, "RECREATE")
for key in histDict.keys():
    print "---'%s' histogram has %i entries, writing to file!---" % (histDict[key].GetTitle(), histDict[key].GetEntries())
    outFile.WriteTObject(histDict[key],histDict[key].GetName())
outFile.Close()
print "~~Saved"
print "*****END******\n"
