import pmag
import numpy as np
import matplotlib.pyplot as plt

def subplot_net(title="", quadrant='all'):
    ax=plt.gca()
    # modified from plotNET to look nicer, play more nicely with subfigures; 
    # also modified to allow quadrants/halves only to be plotted
    """
    draws circle and tick marks for equal area projection
    can zoom in on quadrant/half of stereonet using quadrant='NE','E','S' etc.
    """
    quadrants= {'all': [-1.05,1.05,-1.05,1.05],'N': [-1.05,1.05,-0.1,1.05], 'E': [-0.1,1.05,-1.05,1.05], 
            'S': [-1.05,1.05,-1.05,0.1],'W': [-1.05,0.1,-1.05,1.05], 'NE': [-0.1,1.05,-0.1,1.05],
            'SE': [-0.1,1.05,-1.05,0.1],'SW': [-1.05,0.1,-1.05,0.1], 'NW': [-0.1,1.05,-1.05,0.1]}
     
    # make the perimeter
    ax.axis("off")
    Dcirc=np.arange(0,361.)
    Icirc=np.zeros(361,'f')
    Xcirc,Ycirc=[],[]
    for k in range(361):
        XY= pmag.dimap(Dcirc[k],Icirc[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    ax.plot(Xcirc,Ycirc,'k', linewidth=2, zorder=2)
    ax.set_xlim(quadrants[quadrant][:2])
    ax.set_ylim(quadrants[quadrant][-2:])
    # put on the tick marks
    Xsym,Ysym=[],[]
    for I in range(10,100,10):
        XY=pmag.dimap(0.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    ax.scatter(Xsym,Ysym, color='grey', marker='+', s=30)
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(90.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    ax.scatter(Xsym,Ysym,color='grey', marker='+', s=30)
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(180.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    ax.scatter(Xsym,Ysym,color='grey', marker='+', s=30)
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(270.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    ax.scatter(Xsym,Ysym,color='grey', marker='+', s=30)
    for D in range(0,360,10):
        Xtick,Ytick=[],[]
        for I in range(4):
            XY=pmag.dimap(D,I)
            Xtick.append(XY[0])
            Ytick.append(XY[1])
        ax.plot(Xtick,Ytick,color='grey', linewidth=1.5,zorder=1) 
    ax.set_title(title)
    ax.set(aspect=1)
