#!/usr/local/bin/python
##!/usr/common/usg/software/python/2.7.6/bin/python

def main():
	"""Create the input to skim_mjd_data: a list of M2 runs and muon times.
	Must be sure that the ranges of m2 runs defined below EXACTLY match what runs are used in skim_mjd_data.

	Output of muFinder:
	(m1) run, unix start, xtime, mutype, badscaler

	Output of generate_mu_list:
	run, unix start, xtime, mutype, badscaler

	The only difference is that generate_mu_list accounts for gaps between runs.

	1st M1 BG DS3 run: 16796 (2016-8-25)
	             Last: 17493 (2016-9-14)

	1st M2 BG DS3 run: 60000802 (2016-8-25)
	             Last: 60001506 (2016-9-14)

	In muFinder, I checked for muons in all veto runs, regardless of DAQ state (calibration, transition, etc.)
	Noted any gaps between stop/start of VETO (i.e. M1) runs and flagged the first entry time as
	"potential missed muon."  This means we skip what generate_mu_list does.
	"""

	m1runs, m1start, m1xTime, muType, badScaler = list(), list(), list(), list(), list()
	with open('MuonList_DS3AllVetoRuns.txt') as file:
		for line in file:
			line = line.strip()
			col = line.split()
			m1runs.append(int(col[0]))
			m1start.append(int(col[1]))
			m1xTime.append(float(col[2]))
			muType.append(int(col[3]))
			badScaler.append(int(col[4]))

	# create a list of muon times (unix time of m1 run start + muon time during run)
	muonTimes = list()
	for i in range(len(m1start)):
		muonTimes.append(m1start[i]+m1xTime[i])

	m2runs = list(range(60000802,60000823+1))
	m2runs += list(range(60000827,60000828+1))
	m2runs += list(range(60000830,60000830+1))
	m2runs += list(range(60000847,60000847+1))
	m2runs += list(range(60000850,60000851+1))
	m2runs += list(range(60000854,60000855+1))
	m2runs += list(range(60000858,60000862+1))
	m2runs += list(range(60000868,60000897+1))
	m2runs += list(range(60000898,60000922+1))
	m2runs += list(range(60000928,60000943+1))
	m2runs += list(range(60000953,60000953+1))
	m2runs += list(range(60000970,60001000+1))
	m2runs += list(range(60001001,60001011+1))
	m2runs += list(range(60001013,60001013+1))
	m2runs += list(range(60001033,60001062+1))
	m2runs += list(range(60001063,60001086+1))
	m2runs += list(range(60001088,60001093+1))
	m2runs += list(range(60001094,60001124+1))
	m2runs += list(range(60001125,60001129+1))
	m2runs += list(range(60001132,60001132+1))
	m2runs += list(range(60001135,60001135+1))
	m2runs += list(range(60001139,60001139+1))
	m2runs += list(range(60001144,60001144+1))
	m2runs += list(range(60001146,60001147+1))
	m2runs += list(range(60001163,60001181+1))
	m2runs += list(range(60001183,60001185+1))
	m2runs += list(range(60001187,60001205+1))
	m2runs += list(range(60001308,60001319+1))
	m2runs += list(range(60001330,60001350+1))
	m2runs += list(range(60001379,60001382+1))
	m2runs += list(range(60001384,60001414+1))
	m2runs += list(range(60001415,60001441+1))
	m2runs += list(range(60001463,60001489+1))
	m2runs += list(range(60001491,60001506+1))

	# now figure out the unix start and stop times of all m2 runs
	# using the program GetRunInfo (downloads start/stop times from the DAQ1 feresa DB)

	# This is input to GetRunInfo
	# thefile = open('DS3-m2runs.txt', 'w')
	# for item in m2runs:
	# 	thefile.write("%s\n" % item)

	# This is output of GetRunInfo: DS3-m2runInfo.txt
	# format: run elapsedTime runStart unixStart runStop unixStop
	# ex. 60000802	3600.05	2016/08/26 00:00:41	1472169641	2016/08/26 01:00:41	1472173241

	db_m2Runs, db_m2Duration, db_m2unixStart, db_m2unixStop = list(), list(), list(), list()
	with open('DS3-m2runInfo.txt') as file:
		header1 = file.readline()
		for line in file:
			line = line.strip()
			col = line.split()
			db_m2Runs.append(int(col[0]))
			db_m2Duration.append(float(col[1]))
			db_m2unixStart.append(int(col[4]))
			db_m2unixStop.append(int(col[7]))

	# For each run in DS3-M2, search the info and figure out if there was a muon in M1,
	# and calculate the time offset.

	muonList = open("MuonList_DS3_aug.txt", 'w')

	for m2run in m2runs:

		# get the unix start/stop times for this run
		if m2run in db_m2Runs:
			i = db_m2Runs.index(m2run)
			m2start = db_m2unixStart[i]
			m2stop = db_m2unixStop[i]
			m2length = db_m2Duration[i]

			# loop over the muon list and find entries that fall in
			# the time range of this M2 run
			for time in muonTimes:
				if (time >= m2start and time <= m2stop):
					j = muonTimes.index(time)
					xTime = time-m2start

					# check that this doesn't happen before generating the final list
					if xTime > m2length:
						print "hey, this doesn't make sense!"
						break

					# this is the final list for DS3-M2
					print m2run,m2start,xTime,muType[j],badScaler[j]
					muonList.write("%i %i %.8f %i %i\n" % (m2run,m2start,xTime,muType[j],badScaler[j]))

		else:
			print "couldn't find run",m2run


if __name__ == "__main__":
	main()
