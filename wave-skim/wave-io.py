#!/usr/common/usg/software/python/2.7.6/bin/python
###!/usr/local/bin/python

from ROOT import *

def main():

	grabEvents()
	readEvents()

def grabEvents():

	ds = GATDataSet(12470)

	gat = ds.GetGatifiedChain()
	blt = ds.GetBuiltChain()

	output = TFile("./output/wave-io.root","RECREATE")
	tree = TTree("tree","tree")

	event = MGTEvent()
	wave = MGTWaveform()
	tree.Branch("event",event)
	tree.Branch("wave",wave)

	print "output:"
	for i in xrange(10):
		gat.GetEntry(i)
		blt.GetEntry(i)

		# ian's magic trick: assign it twice and it sticks.
		event = blt.event
		event = blt.event
		# "To be honest, I have no clue what is happening when you call that line, since MGTEvent does not have an = operator or a copy constructor explicitly defined. What this means is that an automatically generated one is being used, which might be the problem."

		# this doesn't work.
		wave = event.GetWaveform(0)
		wave = event.GetWaveform(0)

		print event.GetNWaveforms(),wf.GetLength()

		tree.Fill()

	tree.Write()
	output.Close()

def readEvents():

	infile = TFile("./output/wave-io.root")
	tree = infile.Get("tree")

	event = MGTEvent()
	wave = MGTWaveform()
	tree.SetBranchAddress("event",event)
	tree.SetBranchAddress("wave",wave)

	print "input:"
	for i in xrange(tree.GetEntries()):
		tree.GetEntry(i)

		print event.GetNWaveforms(),wave.GetLength()

		# this works because "event" is actually full
		# for i in xrange(event.GetNWaveforms()):
			# print event.GetWaveform(i).GetLength()

if __name__ == "__main__":
	main()
