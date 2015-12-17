# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:09:53 2015

@author: TinoValentin
"""

import matplotlib.animation as animation
from numpy import *
from matplotlib.pylab import*
from pylab import *





def ani_frame(Array, path, title):
    
    fig = plt.figure()
    fig.suptitle(title)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    length = shape(Array)[2]
    im = ax.imshow(Array[:,:,0], animated = True, cmap = "gist_ncar")
    fig.colorbar(im)
    dpi = 150
    
    tight_layout()



    def update_img(n):
        tmp = Array[:,:,n]
        im.set_data(tmp)
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
    im = plt.imshow(Array, cmap="gist_ncar")
    fig.colorbar(im, ax=ax)
    plt.savefig(path)
    plt.close()
    