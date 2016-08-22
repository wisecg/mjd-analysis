#!/usr/local/bin/python
###!/usr/common/usg/software/python/2.7.6/bin/python

from ROOT import *
import branch_libs as lib

def main():

	skim = TChain("skimTree")
	skim.Add("/Users/wisecg/datasets/ds1/*.root")
	lib.SetTreeInputs(skim, lib.skimDict)
	skim.GetEntry(0)
	print lib.s_trapENFCal.at(0)	# vector syntax
	print lib.s_run[0]				# scalar syntax

	ds = GATDataSet(12486)
	gat = ds.GetGatifiedChain()
	lib.SetTreeInputs(gat, lib.gatDict)
	gat.GetEntry(0)
	print lib.trapENFCal.at(0)
	print lib.run[0]

	out = TFile("./output/branch-check.root","RECREATE")
	tree = TTree("outputTree","tree of output")
	lib.SetTreeOutputs(tree, lib.outDict)

	tree.Fill()

	tree.Write()
	out.Close()



if __name__ == "__main__":
	main()
