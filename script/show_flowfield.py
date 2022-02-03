#!/usr/bin/python

from itertools import count
import sys, getopt, re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

colorarray = ['brown','green','fuchsia','orange','blue','lime']

def parseopt(argv):
    path=""
    try:
      opts, args = getopt.getopt(argv,"hp:",["path="])
    except getopt.GetoptError:
      print ('show_flowfield.py -p <input filepath>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print ('show_flowfield.py -p <input filepath>')
         sys.exit()
      elif opt in ("-p"):
         path = arg
         print ('input path is',path)
      else:
         print ('show_flowfield.py -p <input filepath>')
         sys.exit()
    return path

def readfile(path):
    lines = []
    #print(path)    
    with open(str(path)) as f:
        lines = f.readlines()
    result = []
    count = 1
    M = 0
    N = 0
    U = []
    V = []
    P = []
    for line in lines:
        if count == 1:
            SZ = line.split(' ')
            M = int(SZ[0])
            N = int(SZ[1])
        elif (count>=3)&(count<=(2+N+1)): # U
            Ui = line.split(', ')[0:-1]
            U.append([float(i) for i in Ui])
        elif (count>=2+N+4)&(count<=(2+2*N+3)): # V
            Vi = line.split(', ')[0:-1]
            V.append([float(i) for i in Vi])
        elif (count>=2+2*N+6)&(count<=(2+3*N+5)): # P
            Pi = line.split(', ')[0:-1]
            P.append([float(i) for i in Pi])
        count = count+1
    return M, N, U, V, P

def plotStream(M,N,U,V):
    X = np.arange(0.5/N,1,1/N)
    Y = np.arange(0.5/M,1,1/M)
    gY, gX = np.mgrid[0:(N+0.5)/N:1/N,0:(M+0.5)/M:1/M]
    theU = np.zeros((N+1,M+1))
    theV = np.zeros((N+1,M+1))
    for i in range(0,N+1):
        #print(theU[i,:])
        #print(np.size(gX[0,:]))
        #print(X)
        #rint(U[i][:])
        theU[i,:] = np.interp(gX[0,:],X,U[i][:])
    for i in range(0,M+1):
        theV[:,i] = np.interp(gY[:,0],Y,np.array(V)[:,i])
    #print(theU)
    #print(theV)
    mag = np.power((np.power(theU.transpose(),2)+np.power(theV.transpose(),2)),0.5)
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax1 = fig.add_subplot(gs[0, 0])
    strm = ax1.streamplot(gX, gY, theU.transpose(), theV.transpose(), color=mag, linewidth=2, cmap='RdYlBu')
    ax2 = fig.colorbar(strm.lines)
    ax1.set_title('Velocity Stream Line')
    ax1.set_xlabel('X position(%)')
    ax1.set_ylabel('Y Position(%)')
    ax2.ax.set_ylabel('Velocity Magnitude (m/s)')
    plt.tight_layout()
    plt.show()

def plotPressure(M,N,P):
    X = np.arange(0.5/N,1,1/N)
    Y = np.arange(0.5/M,1,1/M)
    gX, gY = np.meshgrid(X,Y)
    #print(np.shape(gX))
    #print(np.shape(P))
    #print(P)
    #print(np.array(P)[:,0:M])
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ctr = ax1.contourf(gX, gY, np.array(P)[:,0:M].transpose(), 10, alpha=1)
    ax2 = fig.colorbar(ctr)
    ax1.set_title('Pressure Contour')
    ax1.set_xlabel('X position(%)')
    ax1.set_ylabel('Y Position(%)')
    ax2.ax.set_ylabel('Pressure Magnitude (Pa)')
    plt.show()

if __name__ == "__main__":
    path = parseopt(sys.argv[1:])
    M, N, U, V, P = readfile(path)
    plotStream(M,N,U,V)
    plotPressure(M,N,P)
    #print(result)
    #plotrect(result)

