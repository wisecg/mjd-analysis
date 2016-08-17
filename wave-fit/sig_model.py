#!/usr/local/bin/python
"""."""
import numpy as np
from pymc import Normal, HalfNormal, deterministic
from scipy.ndimage.filters import gaussian_filter

def createSignalModel_DEPWF(signal, en, tg, na, template):
	""" Info format:
		en = trapENFCal guess
		tg = t0 guess
		na = noise estimate
	"""
	verbose = 1

	switchpoint = Normal('switchpoint', mu=tg, tau=.01)
	slowness_sigma = HalfNormal('slowness_sigma', tau=.01)
	wfScale = Normal('wfScale', mu=en, tau=sigToTau(.25 * en))
	noise_tau = np.power(na, -2)

	if verbose:
		print "  t0 guess is %d" % tg
		print "  t0 init is %d" % switchpoint
		print "  wfScale guess is %f" % en
		print "  wfScale (energy) init is %f" % wfScale
		print "  noise sigma guess is %0.2f" % na
		print "  noise tau is %0.2e" % noise_tau
		print "  signal type is ", type(signal)

	@deterministic(plot=False, name="depModel")
	def depwf_model(s=switchpoint, e=wfScale, sig=slowness_sigma):
		out = np.zeros(len(signal))
		dep_data = template
		dep_data = out  # temp placeholder
		if e < 0:
			e = 0
		dep_data *= e
		if s < 0:
			s = 0
		if s > len(signal):
			s = len(signal)
		s = np.around(s)
		out[s:] += dep_data[0:(len(signal) - s)] * e
		out = gaussian_filter(out, float(sig))
		return out

	baseline_observed = Normal("baseline_observed", mu=depwf_model, tau=noise_tau, value=signal, observed=True)
	return locals()


def sigToTau(sig):
	"""."""
	tau = np.power(sig, -2)
	return tau
