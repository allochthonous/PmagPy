import pmag
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

def plot_Jel(sdata,plotsymbol='o',symbolsize=20,symbolcolor='grey',symbolborder=None, z=4):
    """
    takes sdata, calculates Jelinek P' and T parameters, and plots 
    """
    #this may be a bit too specialised to really be worth it. We'll see
    JelPars=pd.DataFrame([pmag.s_JelPT(data) for data in sdata], columns=['Pdash','T'])
    ax=plt.gca()
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.scatter(JelPars.Pdash,JelPars['T'], marker=plotsymbol, color=symbolcolor,edgecolor=symbolborder, s=symbolsize, zorder=z)
    plt.axhline(0, color='grey',zorder=0)
    plt.xlim(1,1.2)
    plt.ylim(-1,1) 
    plt.xlabel('P \'', fontsize=14)
    plt.ylabel('T', fontsize=14)
    return JelPars


def plot_AMS(sdata,pointsize=50,errors='None',bedding=[],incolor='N',):
    """
    adapted from pmagplotlib.plotANIS
    set errors to 'h' for hext ellipses, 'b' for bootstrap, 'p' for parametric bootstrap
    bootstrap trials currently hard-coded as 1000
    for bedding, give as [strike,dip] ;
    """
    subplot_net() #set up stereonet
    if incolor=='N': colours=['0.4','0.6','0.5'] #specify greyscale colours
    else: colours=['lightcoral','lightskyblue','lightgreen']
    #if bedding plot bedding plane
    if bedding !=[]:
        # hard-coded version of pmag.plotC
        D_c,I_c=pmag.circ(bedding[0]-90.,90.-bedding[1],90.)
        X_c_d,Y_c_d=[],[]
        for k in range(len(D_c)):
            XY=pmag.dimap(D_c[k],I_c[k])
            if I_c[k]>0:
                X_c_d.append(XY[0])
                Y_c_d.append(XY[1])
        plt.plot(X_c_d,Y_c_d,color='grey',dashes=[10,2],linewidth=3,zorder=1)
    Vs=[]
    #plot individual sample data
    for s in sdata:
        tau,V=pmag.doseigs(s)
        Vs.append(V)
    plotEVEC(Vs,pointsize,colours)
    #plot mean eigenvectors
    nf,sigma,avs=pmag.sbar(sdata)
    Vs=[]
    mtau,mV=pmag.doseigs(avs)
    Vs.append(mV)
    plotEVEC(Vs,pointsize*4,['w','w','w'], 'black',7)
    #plot confidence limits
    hpars=pmag.dohext(nf,sigma,avs)
    if errors=='h':
        ellpars=[hpars["v1_dec"],hpars["v1_inc"],hpars["e12"],hpars["v2_dec"],hpars["v2_inc"],hpars["e13"],hpars["v3_dec"],hpars["v3_inc"]]
        plotELL(ellpars,'black',1,1)
        ellpars=[hpars["v2_dec"],hpars["v2_inc"],hpars["e23"],hpars["v3_dec"],hpars["v3_inc"],hpars["e12"],hpars["v1_dec"],hpars["v1_inc"]]
        plotELL(ellpars,'black',1,1)
        ellpars=[hpars["v3_dec"],hpars["v3_inc"],hpars["e13"],hpars["v1_dec"],hpars["v1_inc"],hpars["e23"],hpars["v2_dec"],hpars["v2_inc"]]
        plotELL(ellpars,'black',1,1)
    elif errors=='b' or errors=='p':
        if errors=='p': ipar=1
        else: ipar=0
        Tmean,Vmean,Taus,BVs=pmag.s_boot(sdata,ipar,1000) # get eigenvectors of mean tensor
        bpars=pmag.sbootpars(Taus,BVs)
        bpars['t1']=hpars['t1']
        bpars['t2']=hpars['t2']
        bpars['t3']=hpars['t3']
        #have to pair zeta value with eta dec/inc, and vice-versa, to align error ellipses correctly.
        ellpars=[bpars["v1_dec"],bpars["v1_inc"],bpars["v1_eta"],bpars["v1_eta_dec"],bpars["v1_eta_inc"],bpars["v1_zeta"],bpars["v1_zeta_dec"],bpars["v1_zeta_inc"]]
        plotELL(ellpars,'black',1,1)
        ellpars=[bpars["v2_dec"],bpars["v2_inc"],bpars["v2_eta"],bpars["v2_eta_dec"],bpars["v2_eta_inc"],bpars["v2_zeta"],bpars["v2_zeta_dec"],bpars["v2_zeta_inc"]]
        plotELL(ellpars,'black',1,1)
        ellpars=[bpars["v3_dec"],bpars["v3_inc"],bpars["v3_eta"],bpars["v3_eta_dec"],bpars["v3_eta_inc"],bpars["v3_zeta"],bpars["v3_zeta_dec"],bpars["v3_zeta_inc"]]
        plotELL(ellpars,'black',1,1) 
           
def plotEVEC(Vs,symsize=40,colours=['lightcoral','lightskyblue','lightgreen'],symboledgecolor='none',level=5):
    """
    plots eigenvector directions of S vectors
    adapted from pmagplotlib.plotEVEC
    """
    symb,symkey=['s','v','o'],0 # plot V1s as squares, V2s as triangles and V3s as circles
    for VEC in range(3):
        X,Y=[],[]
        for Vdirs in Vs:
    #  plot the V1 data  first
            XY=pmag.dimap(Vdirs[VEC][0],Vdirs[VEC][1])
            X.append(XY[0])
            Y.append(XY[1])
        plt.scatter(X,Y,s=symsize,marker=symb[VEC],c=colours[VEC],edgecolors=symboledgecolor, zorder=level)
        
def plotELL(pars,col,lower,plot):
    """
    function to calculate points on an ellipse about Pdec,Pdip with angle beta,gamma
    """
    rad=np.pi/180.
    Pdec,Pinc,beta,Bdec,Binc,gamma,Gdec,Ginc=pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7]
    if beta > 90. or gamma>90:
        beta=180.-beta
        gamma=180.-beta
        Pdec=Pdec-180.
        Pinc=-Pinc
    beta,gamma=beta*rad,gamma*rad # convert to radians
    X_ell,Y_ell,X_up,Y_up,PTS=[],[],[],[],[]
    nums=201
    xnum=float(nums-1.)/2.
# set up t matrix
    t=[[0,0,0],[0,0,0],[0,0,0]]
    X=pmag.dir2cart((Pdec,Pinc,1.0)) # convert to cartesian coordintes
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
# set up rotation matrix t
    t[0][2]=X[0]
    t[1][2]=X[1]
    t[2][2]=X[2]
    X=pmag.dir2cart((Bdec,Binc,1.0))
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
    t[0][0]=X[0]
    t[1][0]=X[1]
    t[2][0]=X[2]
    X=pmag.dir2cart((Gdec,Ginc,1.0))
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
    t[0][1]=X[0]
    t[1][1]=X[1]
    t[2][1]=X[2]
# set up v matrix
    v=[0,0,0]
    for i in range(nums):  # incremental point along ellipse
        psi=float(i)*np.pi/xnum
        v[0]=np.sin(beta)*np.cos(psi)
        v[1]=np.sin(gamma)*np.sin(psi)
        v[2]=np.sqrt(1.-v[0]**2 - v[1]**2)
        elli=[0,0,0]
# calculate points on the ellipse
        for j in range(3):
            for k in range(3):
                elli[j]=elli[j] + t[j][k]*v[k]  # cartesian coordinate j of ellipse
        PTS.append(pmag.cart2dir(elli))
        R=np.sqrt( 1.-abs(elli[2]))/(np.sqrt(elli[0]**2+elli[1]**2)) # put on an equal area projection
        if elli[2]<0:
#            for i in range(3): elli[i]=-elli[i]
            X_up.append(elli[1]*R)
            Y_up.append(elli[0]*R)
            # Adding None values stops plotting of an additional straight line 
            # between the points where the ellipse crosses the edge of the stereonet
            X_ell.append(None)
            Y_ell.append(None)            
        else:
            X_ell.append(elli[1]*R)
            Y_ell.append(elli[0]*R)
    if plot==1:
        if X_ell!=[]:plt.plot(X_ell,Y_ell,color=col, linewidth=2, zorder=6)
        if X_up!=[]:plt.plot(X_up,Y_up,color=col,linewidth=2,linestyle=':')
    else:
        return PTS



