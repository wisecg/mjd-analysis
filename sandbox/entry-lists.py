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
	""" TChain case 2:
	Handles arbitrary TCuts after playing some python list tricks.
	Gave same results as case 1 in a few test cases."""
	skim2 = TChain("skimTree")
	skim2.Add("/Users/wisecg/datasets/ds1/*.root")
	skim2.Draw("Entry$","run==9425","GOFF")
	v1 = skim2.GetV1() # a ROOT.PyDoubleBuffer, very irritating to work with
	v1_list = list(set(v1[n] for n in xrange(skim2.GetEntries())))
	v1_sorted = sorted(v1_list)
	elist4 = TEntryList()
	for i in v1_sorted:
		elist4.Enter(long(i),skim2)
		print v1_sorted[i]
	print type(elist4), elist4.GetN()

	# now show how to loop over a chain using the entry list
	# elist4.Print("all")
	skim2.SetEntryList(elist4);
	numEntries = elist4.GetN()
	for iEntry in xrange(numEntries):
		entryNumber = skimTree.GetEntryNumber(iEntry);
		# print entryNumber
		skim.GetEntry(entryNumber)

	# note that we could also loop directly over the 'v1_sorted' list:


if __name__ == "__main__":
	main()
