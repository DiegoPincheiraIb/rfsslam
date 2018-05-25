#!/usr/bin/python

 #
 # Software License Agreement (New BSD License)
 #
 # Copyright (c) 2018, Felipe Inostroza
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     * Redistributions of source code must retain the above copyright
 #       notice, this list of conditions and the following disclaimer.
 #     * Redistributions in binary form must reproduce the above copyright
 #       notice, this list of conditions and the following disclaimer in the
 #       documentation and/or other materials provided with the distribution.
 #     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 #       Universidad de Chile, nor the names of its contributors may be
 #       used to endorse or promote products derived from this software without
 #       specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT
 # HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 # GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 # LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 # THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #


import sys
import os.path
import numpy as np
import time

import matplotlib
matplotlib.use("TkAgg");

#print matplotlib.__version__

import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker

matplotlib.rcParams.update({'font.size': 20})

saveMovie = True;
saveFig = True


nTrDrawMax = 500
if len(sys.argv) < 2:
    print "Usage: animate1dSim DATA_DIR\n";
    sys.exit(0);

# Setting for file names

dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'

gtPoseFile = 'gtPose.dat';
gtPoseFile = dataDir + gtPoseFile;
if os.path.exists(gtPoseFile):
    print('Opening ' + gtPoseFile);
else:
    print(gtPoseFile + ' does not exist')
    sys.exit(0);

drPoseFile = 'deadReckoning.dat';
drPoseFile = dataDir + drPoseFile;
if os.path.exists(drPoseFile):
    print('Opening ' + drPoseFile);
else:
    print(drPoseFile + ' does not exist')
    sys.exit(0);
gtMapFile = 'gtLandmark.dat';
gtMapFile = dataDir + gtMapFile;
if os.path.exists(gtMapFile):
    print('Opening ' + gtMapFile);
else:
    print(gtMapFile + ' does not exist')
    sys.exit(0);

estPoseFile = 'particlePose.dat';
estPoseFile = dataDir + estPoseFile;
if os.path.exists(estPoseFile):
    print('Opening ' + estPoseFile);
else:
    print(estPoseFile + ' does not exist');
    sys.exit(0);
estPoseFileHandle = open(estPoseFile, "r");

estMapFile = 'landmarkEst.dat';
estMapFile = dataDir + estMapFile;
if os.path.exists(estMapFile):
    print('Opening ' + estMapFile);
else:
    print(estMapFile + ' does not exist')
    sys.exit(0);
estMapFileHandle = open(estMapFile, "r");

measurementFile = 'measurement.dat';
measurementFile = dataDir + measurementFile;
if os.path.exists(measurementFile):
    print('Opening ' + measurementFile);
else:
    print(measurementFile + ' does not exist')
    sys.exit(0);
measurementFileHandle = open(measurementFile, "r");

estimateImageFile = 'estimate.pdf';
estimateImageFile = dataDir + estimateImageFile;
errorImageFile = 'error.pdf';
errorImageFile = dataDir + errorImageFile;
estimateMovieFile = 'estimate.mp4';
estimateMovieFile = dataDir + estimateMovieFile;

# Reading files

print('Reading ' + gtPoseFile);
gtPose = np.genfromtxt(gtPoseFile);
gtPose_t = gtPose[:,0];
gtPose_x = gtPose[:,1];
gtPose_y = gtPose[:,2];
gtPose_th = gtPose[:,3];
print('Number of poses: ' + str(len(gtPose)));
print('Reading ' + drPoseFile);
drPose = np.genfromtxt(drPoseFile);
drPose_t = drPose[:,0];
drPose_x = drPose[:,1];
drPose_y = drPose[:,2];
drPose_th = drPose[:,3];


print('Reading ' + gtMapFile);
gtMap = np.atleast_2d(np.genfromtxt(gtMapFile) );
gtMap_x = gtMap[:,0];
gtMap_y = gtMap[:,1];

print('Reading ' + measurementFile);
measurements = np.genfromtxt(measurementFile);
measurements_i = measurements[:,0].astype(int);
print(measurements_i);
measurements_r = measurements[:,1];
measurements_b = measurements[:,2];


# Plotting

fig = plt.figure( figsize=(12,10), facecolor='w')

axTr = plt.subplot2grid((2, 1), (0, 0))
axMap = plt.subplot2grid((2, 1), (1, 0))

gtMapHandle, = axMap.plot(gtMap_x, gtMap_y, 'r*', markersize=15, zorder=5);
ax = plt.gca();



gtPoseHandle, = axTr.plot(gtPose_x, gtPose_y, 'r-', zorder=12);

drPoseHandle, = axTr.plot(drPose_x, drPose_y, 'r--', zorder=10);

bestPoseHandle, = axTr.plot([], [], 'g-',linewidth=3,zorder = 9);


trajectories = [];
angles = [None] * nTrDrawMax
for i in range(0, nTrDrawMax) :
    trajectories_line, = axTr.plot([], [], 'b-');
    trajectories.append( trajectories_line );
    
measurementHandle, = axMap.plot([], [], 'g*',zorder = 9);
bestLandmarks, = axMap.plot([],[], linestyle='', marker='o', color='r', zorder=9)
landmarks = [];
for i in range(0, nTrDrawMax) :
    landmark_particle, = axMap.plot([],[], linestyle='', marker='o', color='b', zorder=3)
    landmarks.append(landmark_particle);




xLim = axTr.get_xlim();
yLim = axTr.get_ylim();
txt = axTr.text(xLim[1]*0.5, yLim[1]*0.9, " ",zorder=20);

#axMap.set_title("Map")
axMap.set_xlabel("x [m]")
axMap.set_ylabel("y [m]")
#axTr.set_title("Trajectory")
axTr.set_xlabel("x [m]")
axTr.set_ylabel("y [m]")

def animateInit():

    txt.set_text("Iteration: ");
    global p;
    global m;
    poseLine = estPoseFileHandle.readline()
    p =np.fromstring(poseLine,dtype=float,sep=' ');
    mapline = estMapFileHandle.readline()
    m = np.fromstring(mapline,dtype=float,sep=' ');

    for i in range(0, nTrDrawMax):
        landmarks[i].set_data([],[])
    return [];

def animate(i):

    global p;
    global m;
    global z;

    txt.set_text("Iteration: "+ str(i));


    drawnObjects = [];
    drawnObjects.append(txt);

    timeStart = gtPose_t[0]

    p_idx_maxWeight = 0
    p_maxWeight = 0


    nparticle=0;
    while len(p)>0 and p[0] < i:
      poseLine = estPoseFileHandle.readline()
      p =np.fromstring(poseLine,dtype=float,sep=' ');
    bestparticle = 0;
    bestweight = -float("inf");
    while len(p)>0 and p[0] == i:

        if(bestweight<p[1]):
          bestparticle = nparticle
          bestweight = p[1]

        
        trajectories[nparticle].set_data(p[2::3], p[3::3]);
        angles[nparticle] = p[4::3]
        drawnObjects.append(trajectories[nparticle]);

        poseLine = estPoseFileHandle.readline()
        p =np.fromstring(poseLine,dtype=float,sep=' ');
        nparticle = nparticle + 1

    #print('traj ' + str(nparticle) + ' i ' + str(i) + '  p   '+ str(p))
    bestPoseHandle.set_data(trajectories[bestparticle].get_xdata() , trajectories[bestparticle].get_ydata())
    nparticle=0;
    
    measurement_x = np.multiply(measurements_r , np.cos(measurements_b + angles[bestparticle][measurements_i])) + trajectories[bestparticle].get_xdata()[measurements_i];
    measurement_y = np.multiply(measurements_r , np.sin(measurements_b + angles[bestparticle][measurements_i])) + trajectories[bestparticle].get_ydata()[measurements_i];

    measurementHandle.set_data(measurement_x , measurement_y)
    while len(m)>0 and m[0] < i:
      mapline = estMapFileHandle.readline()
      m = np.fromstring(mapline,dtype=float,sep=' ');
    while len(m)>0 and m[0] == i :

      landmarks[nparticle].set_data(m[1::2], m[2::2])


      drawnObjects.append(landmarks[nparticle]);
      mapline = estMapFileHandle.readline()
      m = np.fromstring(mapline,dtype=float,sep=' ');
      nparticle = nparticle + 1;
    bestLandmarks.set_data(landmarks[bestparticle].get_xdata() , landmarks[bestparticle].get_ydata())
    #print('maps ' + str(nparticle) + ' i ' + str(i) + '  m   '+ str(m) )
    #print("i: "+str(i))
    drawnObjects.append(gtPoseHandle);
    drawnObjects.append(drPoseHandle);
    drawnObjects.append(gtMapHandle);
    drawnObjects.append(bestLandmarks);
    drawnObjects.append(bestPoseHandle);
    drawnObjects.append(measurementHandle);




    return drawnObjects;

animation = anim.FuncAnimation(plt.figure(1), animate, np.arange(0, 100000 ,1000), interval=1,
                               init_func=animateInit, blit=True,  repeat=False);
if saveMovie:
    FFMpegWriter = matplotlib.animation.writers['ffmpeg']
    animation.save(estimateMovieFile, writer=FFMpegWriter(fps = 30))
else:
    plt.show(block=False)

if saveFig:
    # Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['ps.useafm'] = True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rcParams['text.usetex'] = True
    #plt.rc('text', usetex=True)

    axMap.autoscale()
    axMap.margins(0.05)

    plt.setp(gtPoseHandle, linewidth=2.0)
    txt.set_text(" ");
    axMap.legend([ gtMapHandle, bestLandmarks, landmarks[0]], [r"$\mathcal{M}$" , r"$\widehat{\mathcal{M}}$" , r"$\mathcal{M}^i$" ], loc='best');
    axTr.legend([gtPoseHandle , bestPoseHandle , trajectories[0] ], ["Ground-truth trajectory" , "Estimated trajectory" , "Particle trajectory"], loc='best')
    
    #scale = 10;
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    #plt.gca().xaxis.set_major_formatter(ticks)
    #ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*scale))
    #plt.gca().yaxis.set_major_formatter(ticks)
    plt.savefig(estimateImageFile, format='pdf', bbox_inches='tight')
    #plt.savefig('estimate.eps', format='eps', bbox_inches='tight')
    errorfig = plt.figure( figsize=(12,10), facecolor='w')
    axError = errorfig.gca()
    np.savetxt(dataDir+'error.txt', np.column_stack( ( np.abs(gtPose_x-bestPoseHandle.get_ydata()) ,  np.abs(gtPose_y-bestPoseHandle.get_ydata()) ) ) )
    print( np.abs(gtPose_x-bestPoseHandle.get_xdata()) )
    print( np.abs(gtPose_y-bestPoseHandle.get_ydata()) )
    xerrorHandle, = axError.plot(gtPose_t , np.abs(gtPose_x-bestPoseHandle.get_xdata())  ,'b-')
    yerrorHandle, = axError.plot(gtPose_t , np.abs(gtPose_y-bestPoseHandle.get_ydata())  ,'r-')
    axError.legend([xerrorHandle , yerrorHandle ], ["Absolute error in x" , "Absolute error in y" ], loc='best')
    axError.set_xlabel("time [s]")
    axError.set_ylabel("error [m]")
    plt.savefig(errorImageFile, format='pdf', bbox_inches='tight')
    
measurementFileHandle.close();
