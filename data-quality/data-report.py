#!/usr/local/bin/python
"""
data-report.py
C. Wiseman 9/10/2016

Wraps the program data-quality.cc, takes its output ROOT
file and generates a run summary PDF for this run.
Intended for near-term analysis.

Two plots are included for each detector:
	"gmax", the maximum or minimum value of each waveform (whichever has a greater absolute value)
	"gbase", the average value of the first four samples of each waveform.
"""

import sys, os
from pylatex import *
from ROOT import *
import numpy as np

M1Channels = {"P1D1":584, "P1D2":582, "P1D3":580, "P1D4":578, "P2D1":692, "P2D2":648, "P2D3":640, "P2D4":642, "P3D1":616, "P3D2":610, "P3D3":608, "P3D4":664, "P4D1":624, "P4D2":628, "P4D3":688, "P4D4":694, "P4D5":614, "P5D1":680, "P5D2":678, "P5D3":672, "P5D4":696, "P6D1":632, "P6D2":630, "P6D3":626, "P6D4":690, "P7D1":600, "P7D2":598, "P7D3":594, "P7D4":592}

M2Channels = {"P1D1":1140, "P1D2":1142, "P1D3":1110, "P1D4":1204, "P2D1":1174, "P2D2":1144, "P2D3":1106, "P2D4":1108, "P2D5":1138, "P3D1":1176, "P3D2":1172, "P3D3":1202, "P4D1":1170, "P4D2":1208, "P4D3":1206, "P4D4":1136, "P4D5":1168, "P5D1":1330, "P5D2":1304, "P5D3":1332, "P5D4":1302, "P6D1":1296, "P6D2":1298, "P6D3":1328, "P6D4":1234, "P7D1":1268, "P7D2":1238, "P7D3":1236, "P7D4":1232}

def main(argv):

	if len(argv) < 1:
		print "Usage: ./data-report.py [run number]"
		return

	run = int(argv[0])

	theFile = "./data-quality_%d.root" % run

	if os.path.isfile(theFile):
		infile = TFile(theFile)
	else:
		print "runnning data-quality (generates a ROOT file)"
		os.system("make && ./data-quality %d" % run)
		infile = TFile(theFile)
		return

	# plotTypes = ["gwf","gbase","grms","hbase","hrms","spec"]
	plotTypes = ["gmax","gbase"]

	print "Creating plots ..."
	plotList = createPlots(infile)

	print "Populating document ..."
	doc = Document('data-report')
	fillDocument(doc, plotList, M2Channels, plotTypes)

	print "Compiling document ..."
	doc.generate_pdf('data-report', clean_tex=True)

	print "Cleaning up plots..."
	for plot in plotList:
		os.remove(plot)


def createPlots(infile):

	# don't display "created canvas messages"
	gROOT.SetBatch(kTRUE)
	gROOT.ProcessLine( "gErrorIgnoreLevel = 2001;")

	plotList = []
	for key in gDirectory.GetListOfKeys():
		name = key.GetName()
		obj = key.ReadObj()
		plotname = "%s.png" % name
		plotList.append(plotname)

		if "gwf" in name or "gbase" in name or "gmax" in name:
			can = TCanvas("can","Bob Ross's Canvas",1600,550)
			obj.SetMarkerStyle(kFullDotMedium)
			obj.Draw("AP")
			can.Print(plotname)
		elif "spec" in name:
			can = TCanvas("can","Bob Ross's Canvas",750,450)
			can.SetLogy(1)
			obj.Draw()
			can.Print(plotname)
		elif "grms" in name:
			can = TCanvas("can","Bob Ross's Canvas",750,450)
			obj.SetMarkerStyle(kFullDotMedium)
			obj.Draw("AP")
			can.Print(plotname)
		else :
			can = TCanvas("can","Bob Ross's Canvas",750,450)
			obj.Draw()
			can.Print(plotname)

	return plotList


def fillDocument(doc, plotList, detList, plotTypes):

	# ======= preamble =======
	doc.preamble.append(NoEscape("""
	\usepackage[margin=1in]{geometry}
	\usepackage{fancyhdr}
	\usepackage{bold-extra}
	\usepackage{datetime}
	\usepackage{hyperref}
	\usepackage{graphicx}
	\usepackage[]{algorithm2e}
	\usepackage{framed}
	\usepackage{enumitem}
	\usepackage{multicol} % 2-col table of contents
	\extrafloats{100}
	\\newcommand{\\tty}{\\texttt}
	\\newcommand{\\ita}{\\textit}
	\\newcommand{\\bol}{\\textbf}
	\\newcommand{\\smc}{\\textsc}
	\\newenvironment{description_nospace}
	{ \\begin{description}
	    \setlength{\itemsep}{0pt}
	    \setlength{\parskip}{0pt}
	    \setlength{\parsep}{0pt}     }
	{ \end{description}                  }
	\\newenvironment{enumerate_nospace}
	{ \\begin{enumerate}
	    \setlength{\itemsep}{0pt}
	    \setlength{\parskip}{0pt}
	    \setlength{\parsep}{0pt}     }
	{ \end{enumerate}                  }
	"""))

	# ======= begin docment body =======
	doc.append(NoEscape("""
	\\pagestyle{fancy}
	\\fancyhf{}
	\\rfoot{\\thepage}
	\\begin{center}
	\Large\\bol{Majorana Data Quality Report}\\\\\\
	\large\ita{Auto-generated: \\today, \\hhmmsstime}
	\end{center}
	\\begin{centering}
	  \\tableofcontents
	\\end{centering}
	\\newpage
	\\lhead{\\ita{MJD Data Quality Report}}
	\\lfoot{\\ita{RS \& DC Group}}
	\\renewcommand{\\footrulewidth}{0.4pt}
	\\rfoot{\\thepage}
	\\renewenvironment{framed}[1][\hsize]
	  {\MakeFramed{\hsize#1\advance\hsize-\width \FrameRestore}}%
	  {\endMakeFramed}
	"""))

	dets = sorted(detList.keys())
	for det in dets:
		doc.append(NoEscape("\\newpage"))
		with doc.create(Section(det)):
			doc.append('')

		doc.append(NoEscape("\\begin{figure}[!ht]"))
		doc.append(NoEscape("\\centering"))

		# do the plots in the given order
		# plotTypes = ["gwf","gbase","grms","hbase","hrms","spec"]
		for plType in plotTypes:
			for filename in plotList:
				plot = ""
				if filename.endswith('.png'):
					plot = filename[:-4]
				if plType in plot and det in plot:
					if plType == "gmax":
						doc.append(NoEscape("\hspace*{-0.7cm}\includegraphics[width=510pt]{%s}" % plot))
					# if plType == "gwf":
						# doc.append(NoEscape("\\includegraphics[width=500pt]{%s}" % plot))
						# doc.append(NoEscape("\hspace*{-0.7cm}\includegraphics[width=510pt]{%s}" % plot))
					if plType == "gbase":
						doc.append(NoEscape("\hspace*{-0.7cm}\includegraphics[width=510pt]{%s}" % plot))
					if plType == "grms":
						doc.append(NoEscape("\\includegraphics[width=225pt]{%s}" % plot))
					if plType == "hbase":
						doc.append(NoEscape("\\includegraphics[width=225pt]{%s}" % plot))
					if plType == "hrms":
						doc.append(NoEscape("\\includegraphics[width=225pt]{%s}" % plot))
					if plType == "spec":
						doc.append(NoEscape("\\includegraphics[width=225pt]{%s}" % plot))

		doc.append(NoEscape("\\end{figure}"))

if __name__ == '__main__':
	main(sys.argv[1:])