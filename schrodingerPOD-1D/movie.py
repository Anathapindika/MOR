# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:09:53 2015

@author: TinoValentin
"""

import matplotlib.animation as animation
from numpy import *
from matplotlib.pylab import*
from pylab import *





def ani_frame(X,Array, path, title):
    
    fig = plt.figure()
    fig.suptitle(title)
    ax = plt.axes(xlim=(-4, 4), ylim=(0, 1.5))
    
    ax.set_aspect('equal')
    length = shape(Array)[1]
    im, = ax.plot(X,Array[:,0])
    dpi = 100
    
    tight_layout()

    print(length)

    def update_img(n):
        tmp = Array[:,n]
        im.set_data(X,tmp)
        return im

    ani = animation.FuncAnimation(fig,update_img,length,interval=20)
    plt.show()
    writer = animation.writers['ffmpeg'](fps=10)

    ani.save(path,writer=writer,dpi=dpi)
    print("done with movie")
    plt.close()


def sig_print(Arrays, path):
    fig_s = plt.figure()
    ax_s = fig_s.add_subplot(111)
    im_s = ax_s.plot(Arrays[:])
    fig_s.suptitle("sigma" + "-Decay", fontsize=20)
    plt.xlabel("#-sigma", fontsize=18)
    plt.ylabel("sigma", fontsize=16)
    savefig(path)
    plt.close()
    
def print_frame(Array, path):
    fig, ax = plt.subplots()
    im = plt.plot(Array)
    plt.savefig(path)
    plt.close()
    