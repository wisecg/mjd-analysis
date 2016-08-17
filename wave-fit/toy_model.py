#!/usr/local/bin/python
"""."""
import numpy as np
from pymc import DiscreteUniform, Normal, HalfNormal, deterministic
from scipy.ndimage.filters import gaussian_filter


# This is a 'model factory' function, good for being called over and over.
# See https://pymc-devs.github.io/pymc/modelfitting.html
def createSignalModel(data):
	"""."""
	# set up your model parameters

	# stochastic variables
	# (not compleletely detmined by parent values; they have a prob. distribution.)

	switchpoint = DiscreteUniform('switchpoint', lower=0, upper=len(data))
	early_sigma = HalfNormal('early_sigma', tau=sigToTau(1))
	late_sigma = HalfNormal('late_sigma', tau=sigToTau(1))
	early_mu = Normal('early_mu', mu=.5, tau=sigToTau(1))
	late_mu = Normal('late_mu', mu=.5, tau=sigToTau(1))
	slowness = Normal('slowness', mu=5., tau=sigToTau(1))

	# set up the model for uncertainty (ie, the noise) and the signal (ie, the step function)
	# (deterministic variables : given by values of parents)

	@deterministic(plot=False, name="test")
	def uncertainty_model(s=switchpoint, n=early_sigma, e=late_sigma):
		# Concatenate Uncertainty sigmas (or taus or whatever) around t0
		s = np.around(s)
		out = np.empty(len(data))
		out[:s] = n
		out[s:] = e
		return out

	@deterministic
	def tau(eps=uncertainty_model):
		# pymc uses this tau parameter instead of sigma to model a gaussian.  its annoying.
		return np.power(eps, -2)

	@deterministic(plot=False, name="siggenmodel")
	def signal_model(s=switchpoint, e=early_mu, l=late_mu, b=slowness):
		# makes the step function using the means
		out = np.zeros(len(data))
		out[:s] = e
		out[s:] = l
		blur = gaussian_filter(out, sigma=float(b))
		return blur

	# Full model: normally distributed noise around a step function
	baseline_observed = Normal("baseline_observed", mu=signal_model, tau=tau, value=data, observed=True)

	# Usage :
	# myModel = pymc.Model(tm.createSignalModel(myNumpyArray))
	# M = pymc.MCMC(myModel)
	return locals()


def sigToTau(sig):
	"""."""
	tau = np.power(sig, -2)
	return tau
