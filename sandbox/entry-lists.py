"""
A snippet to remember how making entry lists works,
for three different cases. TChains in pyroot are stupid,
because you can't cast, i.e. there's no equivalent to:
TEntryList *elist = (TEntryList*)gDirectory->Get("elist")
"""
import ROOT
from ROOT import *
import branch_libs as lib

def main():
	# case1()
	# case2()
	# case3()
	case4()

def case1():
	"""TTree case"""
	gatFile = TFile("$MJDDATADIR/gatified/mjd_run12486.root")
	gatTree = gatFile.Get("mjdTree")
	gatTree.Draw(">>elist", "trapENFCal>2500", "entrylist")
	elist = gDirectory.Get("elist")
	print type(elist), elist.GetN()

def case2():
	"""GATDataSet case"""
	ds = GATDataSet(12826,12827)
	gatSkim = ds.GetGatifiedChain()
	gatSkim.Draw(">>elist2","trapENFCal>1000","entrylist")
	elist2 = gDirectory.Get("elist2")
	print type(elist2), elist2.GetN()

def case3():
	"""TChain case 1:
	Works iteratively, but can't handle arbitrary TCuts
	and requires we initialize branches"""
	skim = TChain("skimTree")
	skim.Add("/Users/wisecg/datasets/ds1/*.root")
	lib.SetTreeInputs(skim, lib.skimDict)
	elist3 = TEntryList()
	for i in xrange(skim.GetEntries()):
		skim.GetEntry(i)
		if (lib.s_run[0] == 9422):
			elist3.Enter(i,skim)
	print type(elist3), elist3.GetN()

def case4():
	""" TChain case 2: arbitrary TCuts."""
	skim = TChain("skimTree")
	skim.Add("/Users/wisecg/datasets/ds1/*.root")

	# attempt to play some python list tricks.
	# skim.Draw("Entry$","trapENFCal>1000","GOFF")
	# v1 = skim.GetV1() # a ROOT.PyDoubleBuffer, very irritating to work with
	# v1_list = list(set(v1[n] for n in xrange(skim.GetEntries())))
	# v1_sorted = sorted(v1_list)
	# elist = TEntryList()
	# for i in v1_sorted:
	# 	elist.Enter(long(i),skim)
	# print type(elist), elist.GetN()

	skim.SetEntryList(0)
	skim.Draw(">>elist", "trapENFCal>1000", "entrylist")
	elist = gDirectory.Get("elist")
	print "Number of entries in the entryList is " + str(elist.GetN())
	skim.SetEntryList(elist);

	# loop over entry list
	for ientry in xrange(elist.GetN()):
		entryNumber = skim.GetEntryNumber(ientry)
		skim.LoadTree( entryNumber )
		skim.GetEntry( entryNumber )
		print entryNumber,skim.run, skim.trapENFCal.at(0)

if __name__ == "__main__":
	main()
