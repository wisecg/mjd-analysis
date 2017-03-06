"""
The idea is to take plots like "gmax" and "gbase" in data-quality.cc

and apply a Kalman filter to them, to look for "extrema"
http://scipy-cookbook.readthedocs.io/items/KalmanFiltering.html

Intro to Kalman filters:
http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf

Ben says:

everything is “basically linear,” which in electronics engineering speak means made of gaussians
so you model it with a bunch of kalman filters
and that gives you a statistically robust way to look for discontinuities or other jumps
its how they monitor parameters at like a gas turbine plant or jet engine or shit like that
its called fault detection and is a big component of controls engineering
its like wildly unexciting
but sort of cool math

but, like, say you want to monitor stability of a peak or whatever
you can make a bunch of plots of that peak position and look at them by eye
or you can have a filter that looks at the position vs time and says WOAH WTF BRO if it jumps
kalman filters are markov chain way to do that
and you know we roll markov style up in this bitch
same with rates or whatever
"""