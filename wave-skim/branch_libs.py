"""
Reference for skim file and gat file variables.
Import this library and use the "SetTree" functions
for one-line initialization of trees.

Skim-file VARIABLES have an "s_" prefix, to avoid repeated branch names.
Skim-file OUTPUT BRANCHES have an "s_" prefix as well.

Example: 'branchname', data_member
	skim:	'trapENFCal',s_trapENFCal	<- key must match input file branch name
	gat: 	'trapENFCal',trapENFCal

	output: 's_trapENFCal',s_trapENFCal	<- key can change, but data member must be the same
     	 	'trapENFCal',trapENFCal

Clint Wiseman, USC/Majorana
8/22/2016
"""

from ROOT import *
import numpy as np

def SetTreeInputs(tree, MJDict):
	for key in sorted(MJDict):
		if isinstance(MJDict[key], list):
			# print "key: %-10s  val: %-10s" % (key, MJDict[key][0])
			tree.SetBranchAddress(key, MJDict[key][0])
		else:
			# print "key: %-10s  val: %-10s" % (key, MJDict[key])
			tree.SetBranchAddress(key, MJDict[key])

def SetTreeOutputs(tree, MJDict):
	for key in sorted(MJDict):
		if isinstance(MJDict[key], list):
			# print "key: %-10s  val: %-10s" % (key, MJDict[key][0])
			tree.Branch(key, MJDict[key][0], MJDict[key][1])
		else:
			# print "key: %-10s  val: %-10s" % (key, MJDict[key])
			tree.Branch(key, MJDict[key])

# ================= Output dict (auto-generated) =================
def CreateOutputDict(whichDict='all'):
	"""."""
	outDict = {}
	if whichDict == 'skim':
		outDict = skimDict.copy()	# doesn't append the s_
	if whichDict == 'built':
		outDict = builtDict.copy()
	if whichDict == 'gat':
		outDict = gatDict.copy()
	if whichDict == 'gatBlt':
		outDict = gatDict.copy()
		outDict.update(builtDict)
	if whichDict == 'all':
		for key in skimDict:
			newKey = "s_"+key
			outDict[newKey] = skimDict[key]
		outDict.update(gatDict)
		outDict.update(builtDict)
	return outDict

# ================= Skim file var's (Aug. 2016) ================
s_skimgatrev = np.zeros(1,dtype=int)
s_gatrev = np.zeros(1,dtype=int)
s_run = np.zeros(1,dtype=int)
s_iEvent = np.zeros(1,dtype=int)
s_iHit = ROOT.std.vector("int")()
s_channel = ROOT.std.vector("int")()
s_P = ROOT.std.vector("int")()
s_D = ROOT.std.vector("int")()
s_gain = ROOT.std.vector("int")()
s_mageID = ROOT.std.vector("int")()
s_detID = ROOT.std.vector("int")()
s_detName = ROOT.std.vector("string")()
s_isEnr = ROOT.std.vector("bool")()
s_isNat = ROOT.std.vector("bool")()
s_mAct_g = ROOT.std.vector("double")()
s_isGood = ROOT.std.vector("bool")()
s_startTime = np.zeros(1,dtype=float)
s_stopTime = np.zeros(1,dtype=float)
s_tloc_s = ROOT.std.vector("double")()
s_time_s = ROOT.std.vector("double")()
s_timeMT = ROOT.std.vector("double")()
s_dateMT = ROOT.std.vector("int")()
s_trapECal = ROOT.std.vector("double")()
s_trap4usMax = ROOT.std.vector("double")()
s_trapENFCal = ROOT.std.vector("double")()
s_sumEH = np.zeros(1,dtype=float)
s_sumEL = np.zeros(1,dtype=float)
s_mH = np.zeros(1,dtype=int)
s_mL = np.zeros(1,dtype=int)
s_aenorm = ROOT.std.vector("double")()
s_t150 = ROOT.std.vector("double")()
s_kvorrT = ROOT.std.vector("double")()
s_toe = ROOT.std.vector("double")()
s_aenorm85 = ROOT.std.vector("double")()
s_dcrSlope85 = ROOT.std.vector("double")()
s_dcrSlope95 = ROOT.std.vector("double")()
s_dcrSlope98 = ROOT.std.vector("double")()
s_dcrSlope90 = ROOT.std.vector("double")()
s_dcrSlope99 = ROOT.std.vector("double")()
s_EventDC1Bits = np.zeros(1,dtype=int)
s_wfDCBits = ROOT.std.vector("unsigned int")()
s_isLNFill = ROOT.std.vector("bool")()
s_trapETailMin = ROOT.std.vector("double")()
s_dtmu_s = ROOT.std.vector("double")()
s_muType = ROOT.std.vector("int")()
s_badScaler = ROOT.std.vector("bool")()
s_muVeto = ROOT.std.vector("bool")()

skimDict = {
	'skimgatrev': [s_skimgatrev,'skimgatrev/i'],
	'gatrev': [s_gatrev,'gatrev/i'],
	'run': [s_run,'run/I'],
	'iEvent': [s_iEvent,'iEvent/I'],
	'startTime': [s_startTime,'startTime/D'],
	'stopTime': [s_stopTime,'stopTime/D'],
	'sumEH': [s_sumEH,'sumEH/D'],
	'sumEL': [s_sumEL,'sumEL/D'],
	'mH': [s_mH,'mH/I'],
	'mL': [s_mL,'mL/I'],
	'EventDC1Bits': [s_EventDC1Bits,'EventDC1Bits/i'],
	'iHit':s_iHit,
	'channel':s_channel,
	'P':s_P,
	'D':s_D,
	'gain':s_gain,
	'mageID':s_mageID,
	'detID':s_detID,
	'detName':s_detName,
	'isEnr':s_isEnr,
	'isNat':s_isNat,
	'mAct_g':s_mAct_g,
	'isGood':s_isGood,
	'tloc_s':s_tloc_s,
	'time_s':s_time_s,
	'timeMT':s_timeMT,
	'dateMT':s_dateMT,
	'trapECal':s_trapECal,
	'trap4usMax':s_trap4usMax,
	'trapENFCal':s_trapENFCal,
	'aenorm':s_aenorm,
	't150':s_t150,
	'kvorrT':s_kvorrT,
	'toe':s_toe,
	'aenorm85':s_aenorm85,
	'dcrSlope85':s_dcrSlope85,
	'dcrSlope95':s_dcrSlope95,
	'dcrSlope98':s_dcrSlope98,
	'dcrSlope90':s_dcrSlope90,
	'dcrSlope99':s_dcrSlope99,
	'wfDCBits':s_wfDCBits,
	'isLNFill':s_isLNFill,
	'trapETailMin':s_trapETailMin,
	'dtmu_s':s_dtmu_s,
	'muType':s_muType,
	'badScaler':s_badScaler,
	'muVeto':s_muVeto
}

# ================= GAT file var's (Aug. 2016) =================
timeMT    = ROOT.std.vector("double")()
dateMT    = ROOT.std.vector("int")()
detName   = ROOT.std.vector("string")()
detID     = ROOT.std.vector("int")()
C         = ROOT.std.vector("int")()
P         = ROOT.std.vector("int")()
D         = ROOT.std.vector("int")()
mageID    = ROOT.std.vector("int")()
isEnr     = ROOT.std.vector("bool")()
isNat     = ROOT.std.vector("bool")()
rawWFMax  = ROOT.std.vector("double")()
rawWFMin  = ROOT.std.vector("double")()
blrwfFMR0p1 = ROOT.std.vector("double")()
sgswfFMR0p1 = ROOT.std.vector("double")()
blrwfFMR1 = ROOT.std.vector("double")()
sgswfFMR1 = ROOT.std.vector("double")()
blrwfFMR3 = ROOT.std.vector("double")()
sgswfFMR3 = ROOT.std.vector("double")()
blrwfFMR10 = ROOT.std.vector("double")()
sgswfFMR10 = ROOT.std.vector("double")()
blrwfFMR20 = ROOT.std.vector("double")()
sgswfFMR20 = ROOT.std.vector("double")()
blrwfFMR50 = ROOT.std.vector("double")()
sgswfFMR50 = ROOT.std.vector("double")()
blrwfFMR80 = ROOT.std.vector("double")()
sgswfFMR80 = ROOT.std.vector("double")()
blrwfFMR90 = ROOT.std.vector("double")()
sgswfFMR90 = ROOT.std.vector("double")()
blrwfFMR97 = ROOT.std.vector("double")()
sgswfFMR97 = ROOT.std.vector("double")()
blrwfFMR99 = ROOT.std.vector("double")()
sgswfFMR99 = ROOT.std.vector("double")()
sgswft0   = ROOT.std.vector("double")()
TSCurrent50nsMax = ROOT.std.vector("double")()
TSCurrent100nsMax = ROOT.std.vector("double")()
TSCurrent200nsMax = ROOT.std.vector("double")()
RawWFblSlope = ROOT.std.vector("double")()
RawWFblOffset = ROOT.std.vector("double")()
RawWFblSlopeUnc = ROOT.std.vector("double")()
RawWFblOffsetUnc = ROOT.std.vector("double")()
RawWFblChi2 = ROOT.std.vector("double")()
RawWFftSlope = ROOT.std.vector("double")()
RawWFftSlopeUnc = ROOT.std.vector("double")()
RawWFftOffset = ROOT.std.vector("double")()
RawWFftOffsetUnc = ROOT.std.vector("double")()
RawWFftChi2 = ROOT.std.vector("double")()
RawWFwholeWFLinFitSlope = ROOT.std.vector("double")()
RawWFwholeWFLinFitSlopeUnc = ROOT.std.vector("double")()
RawWFwholeWFLinFitOffset = ROOT.std.vector("double")()
RawWFwholeWFLinFitOffsetUnc = ROOT.std.vector("double")()
RawWFwholeWFLinFitChi2 = ROOT.std.vector("double")()
trapENF55us = ROOT.std.vector("double")()
trapENF70us = ROOT.std.vector("double")()
trapENF85us = ROOT.std.vector("double")()
trapENF100us = ROOT.std.vector("double")()
trapENF115us = ROOT.std.vector("double")()
trapENF130us = ROOT.std.vector("double")()
trapBL    = ROOT.std.vector("double")()
trap500nsMax = ROOT.std.vector("double")()
trapENF500nsrt = ROOT.std.vector("double")()
trap1usMax = ROOT.std.vector("double")()
trapENF1usrt = ROOT.std.vector("double")()
trap2usMax = ROOT.std.vector("double")()
trapENF2usrt = ROOT.std.vector("double")()
trap4usMax = ROOT.std.vector("double")()
trapENF4usrt = ROOT.std.vector("double")()
trap6usMax = ROOT.std.vector("double")()
trapENF6usrt = ROOT.std.vector("double")()
trap8usMax = ROOT.std.vector("double")()
trapENF8usrt = ROOT.std.vector("double")()
trapE     = ROOT.std.vector("double")()
trapEMin  = ROOT.std.vector("double")()
longGapTrapMax = ROOT.std.vector("double")()
trap1ust0 = ROOT.std.vector("double")()
trapENF   = ROOT.std.vector("double")()
trapBL1us = ROOT.std.vector("double")()
trapETailMin = ROOT.std.vector("double")()
trirt50nsft10nsMax = ROOT.std.vector("double")()
trirt100nsft10nsMax = ROOT.std.vector("double")()
trirt200nsft10nsMax = ROOT.std.vector("double")()
triFilMin = ROOT.std.vector("double")()
trirt100nsft10nsIntegralW = ROOT.std.vector("double")()
blrwfIntegralW = ROOT.std.vector("double")()
smoothTrirt100nsft10nsMax = ROOT.std.vector("double")()
energyCal = ROOT.std.vector("double")()
trapECal  = ROOT.std.vector("double")()
trapEMinCal = ROOT.std.vector("double")()
trapENFCal = ROOT.std.vector("double")()
nlcblrwfSlope = ROOT.std.vector("double")()
RawWFdcblSlope = ROOT.std.vector("double")()
RawWFdcblOffset = ROOT.std.vector("double")()
RawWFdcblSlopeUnc = ROOT.std.vector("double")()
RawWFdcblOffsetUnc = ROOT.std.vector("double")()
RawWFdcblChi2 = ROOT.std.vector("double")()
RawWFdcftSlope = ROOT.std.vector("double")()
RawWFdcftOffset = ROOT.std.vector("double")()
RawWFdcftSlopeUnc = ROOT.std.vector("double")()
RawWFdcftOffsetUnc = ROOT.std.vector("double")()
blrwf2ndDiffPeakValue = ROOT.std.vector("double")()
blrwfNorm2ndDiffPeakValue = ROOT.std.vector("double")()
blrwfSSF  = ROOT.std.vector("double")()
BLMoverMerr = ROOT.std.vector("double")()
FTMoverMerr = ROOT.std.vector("double")()
delOffset = ROOT.std.vector("double")()
toe = ROOT.std.vector("double")()
d2wfnoiseTagNorm = ROOT.std.vector("double")()
d2wfDCPower = ROOT.std.vector("double")()
d2wf0MHzTo2MHzPower = ROOT.std.vector("double")()
d2wf2MHzTo10MHzPower = ROOT.std.vector("double")()
d2wf10MHzTo20MHzPower = ROOT.std.vector("double")()
d2wf30MHzTo35MHzPower = ROOT.std.vector("double")()
d2wf48MHzTo50MHzPower = ROOT.std.vector("double")()
d2wf0MHzTo50MHzPower = ROOT.std.vector("double")()
wfDCBits  = ROOT.std.vector("unsigned int")()
wfDC_Bit_0_SSSpikeBL = ROOT.std.vector("double")()
wfDC_Bit_1_SSSpikePhys = ROOT.std.vector("double")()
wfDC_Bit_2_EarlyTrigger = ROOT.std.vector("double")()
wfDC_Bit_3_LateTrigger = ROOT.std.vector("double")()
wfDC_Bit_4_PosSaturatedWFs = ROOT.std.vector("double")()
wfDC_Bit_5_NegSaturatedWFs = ROOT.std.vector("double")()
run       = np.zeros(1,dtype=float)
startTime = np.zeros(1,dtype=float)
stopTime  = np.zeros(1,dtype=float)
energy    = ROOT.std.vector("double")()
channel   = ROOT.std.vector("double")()
timestamp = ROOT.std.vector("double")()
gatrev    = np.zeros(1,dtype=int)
EventDC1Bits = np.zeros(1,dtype=int)

gatDict = {
	'timeMT'    : timeMT,
	'dateMT'    : dateMT,
	'detName'   : detName,
	'detID'     : detID,
	'C'         : C,
	'P'         : P,
	'D'         : D,
	'mageID'    : mageID,
	'isEnr'     : isEnr,
	'isNat'     : isNat,
	'rawWFMax'  : rawWFMax,
	'rawWFMin'  : rawWFMin,
	'blrwfFMR0p1' : blrwfFMR0p1,
	'sgswfFMR0p1' : sgswfFMR0p1,
	'blrwfFMR1' : blrwfFMR1,
	'sgswfFMR1' : sgswfFMR1,
	'blrwfFMR3' : blrwfFMR3,
	'sgswfFMR3' : sgswfFMR3,
	'blrwfFMR10' : blrwfFMR10,
	'sgswfFMR10' : sgswfFMR10,
	'blrwfFMR20' : blrwfFMR20,
	'sgswfFMR20' : sgswfFMR20,
	'blrwfFMR50' : blrwfFMR50,
	'sgswfFMR50' : sgswfFMR50,
	'blrwfFMR80' : blrwfFMR80,
	'sgswfFMR80' : sgswfFMR80,
	'blrwfFMR90' : blrwfFMR90,
	'sgswfFMR90' : sgswfFMR90,
	'blrwfFMR97' : blrwfFMR97,
	'sgswfFMR97' : sgswfFMR97,
	'blrwfFMR99' : blrwfFMR99,
	'sgswfFMR99' : sgswfFMR99,
	'sgswft0'   : sgswft0,
	'TSCurrent50nsMax' : TSCurrent50nsMax,
	'TSCurrent100nsMax' : TSCurrent100nsMax,
	'TSCurrent200nsMax' : TSCurrent200nsMax,
	'RawWFblSlope' : RawWFblSlope,
	'RawWFblOffset' : RawWFblOffset,
	'RawWFblSlopeUnc' : RawWFblSlopeUnc,
	'RawWFblOffsetUnc' : RawWFblOffsetUnc,
	'RawWFblChi2' : RawWFblChi2,
	'RawWFftSlope' : RawWFftSlope,
	'RawWFftSlopeUnc' : RawWFftSlopeUnc,
	'RawWFftOffset' : RawWFftOffset,
	'RawWFftOffsetUnc' : RawWFftOffsetUnc,
	'RawWFftChi2' : RawWFftChi2,
	'RawWFwholeWFLinFitSlope' : RawWFwholeWFLinFitSlope,
	'RawWFwholeWFLinFitSlopeUnc' : RawWFwholeWFLinFitSlopeUnc,
	'RawWFwholeWFLinFitOffset' : RawWFwholeWFLinFitOffset,
	'RawWFwholeWFLinFitOffsetUnc' : RawWFwholeWFLinFitOffsetUnc,
	'RawWFwholeWFLinFitChi2' : RawWFwholeWFLinFitChi2,
	# 'trapENF55us' : trapENF55us,
	'trapENF70us' : trapENF70us,
	# 'trapENF85us' : trapENF85us,
	'trapENF100us' : trapENF100us,
	# 'trapENF115us' : trapENF115us,
	'trapENF130us' : trapENF130us,
	# 'trapBL'    : trapBL,
	'trap500nsMax' : trap500nsMax,
	'trapENF500nsrt' : trapENF500nsrt,
	'trap1usMax' : trap1usMax,
	'trapENF1usrt' : trapENF1usrt,
	'trap2usMax' : trap2usMax,
	'trapENF2usrt' : trapENF2usrt,
	'trap4usMax' : trap4usMax,
	'trapENF4usrt' : trapENF4usrt,
	'trap6usMax' : trap6usMax,
	'trapENF6usrt' : trapENF6usrt,
	'trap8usMax' : trap8usMax,
	'trapENF8usrt' : trapENF8usrt,
	'trapE'     : trapE,
	'trapEMin'  : trapEMin,
	'longGapTrapMax' : longGapTrapMax,
	'trap1ust0' : trap1ust0,
	'trapENF'   : trapENF,
	'trapBL1us' : trapBL1us,
	'trapETailMin' : trapETailMin,
	'trirt50nsft10nsMax' : trirt50nsft10nsMax,
	'trirt100nsft10nsMax' : trirt100nsft10nsMax,
	'trirt200nsft10nsMax' : trirt200nsft10nsMax,
	'triFilMin' : triFilMin,
	'trirt100nsft10nsIntegralW' : trirt100nsft10nsIntegralW,
	'blrwfIntegralW' : blrwfIntegralW,
	'smoothTrirt100nsft10nsMax' : smoothTrirt100nsft10nsMax,
	'energyCal' : energyCal,
	'trapECal'  : trapECal,
	'trapEMinCal' : trapEMinCal,
	'trapENFCal' : trapENFCal,
	'nlcblrwfSlope' : nlcblrwfSlope,
	'RawWFdcblSlope' : RawWFdcblSlope,
	'RawWFdcblOffset' : RawWFdcblOffset,
	'RawWFdcblSlopeUnc' : RawWFdcblSlopeUnc,
	'RawWFdcblOffsetUnc' : RawWFdcblOffsetUnc,
	'RawWFdcblChi2' : RawWFdcblChi2,
	'RawWFdcftSlope' : RawWFdcftSlope,
	'RawWFdcftOffset' : RawWFdcftOffset,
	'RawWFdcftSlopeUnc' : RawWFdcftSlopeUnc,
	'RawWFdcftOffsetUnc' : RawWFdcftOffsetUnc,
	'blrwf2ndDiffPeakValue' : blrwf2ndDiffPeakValue,
	'blrwfNorm2ndDiffPeakValue' : blrwfNorm2ndDiffPeakValue,
	'blrwfSSF'  : blrwfSSF,
	'BLMoverMerr' : BLMoverMerr,
	'FTMoverMerr' : FTMoverMerr,
	'delOffset' : delOffset,
	'toe'       : toe,
	'd2wfnoiseTagNorm' : d2wfnoiseTagNorm,
	'd2wfDCPower' : d2wfDCPower,
	'd2wf0MHzTo2MHzPower' : d2wf0MHzTo2MHzPower,
	'd2wf2MHzTo10MHzPower' : d2wf2MHzTo10MHzPower,
	'd2wf10MHzTo20MHzPower' : d2wf10MHzTo20MHzPower,
	'd2wf30MHzTo35MHzPower' : d2wf30MHzTo35MHzPower,
	'd2wf48MHzTo50MHzPower' : d2wf48MHzTo50MHzPower,
	'd2wf0MHzTo50MHzPower' : d2wf0MHzTo50MHzPower,
	'wfDCBits'  :  wfDCBits,
	'wfDC_Bit_0_SSSpikeBL' : wfDC_Bit_0_SSSpikeBL,
	'wfDC_Bit_1_SSSpikePhys' : wfDC_Bit_1_SSSpikePhys,
	'wfDC_Bit_2_EarlyTrigger' : wfDC_Bit_2_EarlyTrigger,
	'wfDC_Bit_3_LateTrigger' : wfDC_Bit_3_LateTrigger,
	'wfDC_Bit_4_PosSaturatedWFs' : wfDC_Bit_4_PosSaturatedWFs,
	'wfDC_Bit_5_NegSaturatedWFs' : wfDC_Bit_5_NegSaturatedWFs,
	'run'       : [run,'run/D'],
	'startTime' : [startTime,'startTime/D'],
	'stopTime'  : [stopTime,'stopTime/D'],
	'energy'    : energy,
	'channel'   : channel,
	'timestamp' : timestamp,
	'gatrev'    : [gatrev,'gatrev/i'],
	'EventDC1Bits' : [EventDC1Bits,'EventDC1Bits/i']
}

# ================= Built file var's (Aug. 2016) =================
event = MGTEvent()

builtDict = {
	'event':event
}

