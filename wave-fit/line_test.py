#!/usr/local/env/python
"""."""
import numpy as np
import matplotlib.pyplot as plt


def main():
	size = 100
	signal = np.ones(size)

	t0 = 50
	m = 1.1
	b = 0

	for i in xrange(t0-10,t0+10):
		signal[i] = m * signal[i] + b
		print "i %i  sig %f" % (i, signal[i])

	fin = signal[t0+10]
	print "fin: ",fin
	signal[t0+10:] = fin

	plt.figure(figsize=(10, 7), facecolor='w')
	plt.plot(signal)
	plt.show()

if __name__ == "__main__":
	main()
