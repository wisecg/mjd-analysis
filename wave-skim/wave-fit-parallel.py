#!/usr/local/bin/python
import numpy as np
from ROOT import *
import sys

def main(argv):
    # Usage: ./wave-fit.py [-B (batch mode)]
    print "heeeey"
    run_mcmc = True
    opt = ""
    batch = False
    if (len(argv) > 0):
        opt = argv[0]
    if (opt == "-B"):
        batch = True

    theCut = "trapENFCal > 0.8 && trapENFCal < 5 && d2wf0MHzTo50MHzPower < 10000"
    theCut += "&& run!=9436 && run!=9648 && run!=9663 && run!=10185 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=13004"  # channel!=580
    theCut += "&& !(trapENFCal > 0.8 && trapENFCal < 1.5 && d2wf0MHzTo50MHzPower > 4000 && d2wf0MHzTo50MHzPower < 5500)"
    print "\nUsing cut: ",theCut

    calFile = "./data/wave-m1cal.root"
    inFile = "./data/waveSkimDS1_standardCut.root"
    outFile = "./data/waveFitDS1_extendedCut_2.root"
    waveTemplates = TChain("waveTree")
    waveTemplates.Add("./data/wave-m1cal.root")
    # PlotTemplateWaveforms(waveTemplates)
    # exit(1)
    f1 = TFile(inFile)
    waveTree = f1.Get("waveTree")
    print "Found",waveTree.GetEntries(),"input waveforms."

    waveTreeNew, mcmc_params = mcmc_wfs(waveTemplates,waveTree,eventTree,theCut,batch,run_mcmc)
    if (batch == False):
        exit(1)

def mcmc_wfs(waveTemplates, waveTree, eventTree, theCut, batch=False, run_mcmc=True):

    # 16 High-gain Module 1 detectors active in DS1.  BEGEs: 600 p7d1, 692 p2d1.
    chanDict = {610:"P3D2", 600:"P7D1", 692:"P2D1", 598:"P7D2", 626:"P6D3", 648:"P2D2",
    582:"P1D2", 640:"P2D3", 578:"P1D4", 580:"P1D3", 690:"P6D4", 592:"P7D4", 672:"P5D3",
    608:"P3D3", 664:"P3D4", 632:"P6D1"}
    calDict = {640:28, 672:18, 610:16, 580:19, 582:34, 648:38, 600:21, 578:39, 592:27,
    664:55, 626:62, 692:8, 598:22, 690:52, 632:9, 608:7}
    calList = [[key,calDict[key]] for key in calDict]

    waveTree.Draw(">>eList", theCut , "entrylist")
    eList = gDirectory.Get("eList")
    waveTree.SetEntryList(eList)
    print "Found",eList.GetN(),"entries passing cut."

if __name__ == "__main__":
    main(sys.argv[1:])
