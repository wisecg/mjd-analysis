#!/usr/local/bin/python
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation as animate
import numpy as np
import time

def main():
    figPlotter()
    # figAnimator()
    # animator()


def figAnimator():
    """ This one rips at about 60fps.  I guess that's my mac's top speed.
    Porting this to wave-view would require animate to control the loop ... """

    x = np.arange(0,2016*10,10)
    y = np.random.randn(2016)

    fig, ax = plt.subplots()
    line, = ax.plot(x,y)

    def update(i):
        print i
        ynew = np.random.randn(2016)
        line.set_ydata(ynew)
        ymin, ymax = np.amin(ynew), np.amax(ynew)

        # it's faster the fewer pixels change
        ax.set_ylim([ymin-abs(100*ymin),ymax+abs(100*ymax)])
        return line,

    start = time.time()

    a = animate(fig, update, frames=xrange(1, 100), interval=0.1, blit=True, repeat=False)
    plt.show(block=False)

    stop = time.time()
    print 100/(stop-start),"fps"


def animator():
    """ From StackOverflow.  Gives about 30-40 fps.
    http://stackoverflow.com/questions/8955869/why-is-plotting-with-matplotlib-so-slow/8956211
    """
    x = np.arange(0, 2*np.pi, 0.1)
    y = np.sin(x)

    fig, axes = plt.subplots(nrows=6)

    styles = ['r-', 'g-', 'y-', 'm-', 'k-', 'c-']
    def plot(ax, style):
        return ax.plot(x, y, style, animated=True)[0]
    lines = [plot(ax, style) for ax, style in zip(axes, styles)]

    def update(i):
        for j, line in enumerate(lines, start=1):
            line.set_ydata(np.sin(j*x + i/10.0))
        return lines

    start = time.time()
    a = animate(fig, update, xrange(1, 200), interval=0, blit=True, repeat=False)
    plt.show(block=False)
    stop = time.time()
    print 100/(stop-start),"fps"


def figPlotter():
    """ Seems to give about 60 fps on my mac. Requires TkAgg.
        http://bastibe.de/2013-05-30-speeding-up-matplotlib.html
    """
    fig, ax = plt.subplots()
    line, = ax.plot(np.random.randn(100))
    plt.show(block=False)

    tstart = time.time()
    num_plots = 0
    while time.time()-tstart < 30:
        line.set_ydata(np.random.randn(100))
        ax.draw_artist(ax.patch)
        ax.draw_artist(line)
        # fig.canvas.draw() # a dumb mac-specific speed bottleneck
        fig.canvas.blit(ax.bbox) # requires TkAgg backend
        fig.canvas.flush_events()
        num_plots += 1
    print(num_plots/30),"fps"



if __name__ == "__main__":
    main()
