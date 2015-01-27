#!/usr/bin/python

 #
 # Software License Agreement (New BSD License)
 #
 # Copyright (c) 2013, Keith Leung
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

import matplotlib
#matplotlib.use("TKAgg");
#print matplotlib.__version__

import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker   

saveMovie = False;

nLandmarksDrawMax = 1000;
nMeasurementsDrawMax = 100;

if len(sys.argv) < 2:
    print "Usage: animate2dSim DATA_DIR\n";
    sys.exit(0);

# Setting for file names

dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'

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

measurementFile = 'measurements.dat';
measurementFile = dataDir + measurementFile;
if os.path.exists(measurementFile):
    print('Opening ' + measurementFile);
else:
    print(measurementFile + ' does not exist')
    sys.exit(0);
measurementFileHandle = open(measurementFile, "r");

estimateImageFile = 'estimate.pdf';
estimateImageFile = dataDir + estimateImageFile;
estimateMovieFile = 'estimate.mp4';
estimateMovieFile = dataDir + estimateMovieFile;

# Reading files

p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
p_idx = 0;
p_idx_maxWeight = 0;
p_maxWeight = p[5];
p_x = [];
p_y = [];
p_r = [];
p_w = [];
p_t = p[0]
while p[0] == p_t:
    p_x.append(p[2]);
    p_y.append(p[3]);
    p_r.append(p[4]);
    p_w.append(p[5]);
    if p[5] > p_maxWeight :
        p_maxWeight = p[5];
        p_idx_maxWeight = p_idx;
    p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
    p_idx += 1;
px_best = []
py_best = []
pr_best = []

m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);

z = np.fromfile(measurementFileHandle, count=4, sep=" ", dtype=float);

# Plotting 

fig = plt.figure( figsize=(12,12), facecolor='w')

particles, = plt.plot(p_x, p_y, 'b.');
ax = plt.gca()
plt.axis('equal');
plt.grid(True);
plt.xlim([-125, 225])
plt.ylim([-75, 275])

measurements = [];
for i in range(0, nMeasurementsDrawMax) : 
    measurement_line, = plt.plot([], [], 'b-');
    measurements.append( measurement_line );

landmarks = [];
for i in range(0, nLandmarksDrawMax) : 
    landmark_ellipse = Ellipse(xy=(0,0), width=0, height=0, angle=0);
    landmarks.append(landmark_ellipse); 
    ax.add_patch(landmarks[i]);

landmarkCenters, = plt.plot([], [], '+')
landmarkCenters.set_color([0.2,0.2,0.8])

trajectory, = plt.plot(0, 0, 'b-')

xLim = plt.getp(ax, 'xlim');
yLim = plt.getp(ax, 'ylim');
txt = plt.text(150, -65, " ");

def animateInit():

    txt.set_text("Time: ");
    particles.set_data([],[]);
    for i in range(0, nMeasurementsDrawMax) :
        measurements[i].set_data([],[]);
        measurements[i].set_color([1.0, 0.2 ,0.2]);
    for i in range(0, nLandmarksDrawMax):
        landmarks[i].center = (0,0);
        landmarks[i].width = 0;
        landmarks[i].height = 0;
        landmarks[i].set_facecolor([0.2,0.2,0.8])
    return [];

def animate(i):
    
    global p;
    global m;
    global z;

    if not p.any():
        print i
        print "No more messages"
        return []

    currentTime = p[0];
    drawnObjects = [];

    # Time
    txt.set_text("Time: {0}".format(currentTime));
    drawnObjects.append(txt);

    # Particles
    p_idx = 0;
    p_idx_maxWeight = 0;
    p_maxWeight = p[5];
    px_best.append(p[2]);
    py_best.append(p[3]);
    pr_best.append(p[4]);
    p_x = [];
    p_y = [];
    p_w = [];
    while p.any() and abs(p[0] - currentTime) < 1e-4:
        p_x.append(p[2]);
        p_y.append(p[3]);
        p_r.append(p[4]);
        p_w.append(p[5]);
        if p[5] > p_maxWeight :
            p_maxWeight = p[5];
            p_idx_maxWeight = p_idx;
            px_best[i] = p[2];
            py_best[i] = p[3];
            pr_best[i] = p[4];
        p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
        p_idx += 1;
    particles.set_data(p_x, p_y)
    trajectory.set_data(px_best, py_best)
    
    # Landmarks
    m_idx = 0;
    m_x = []
    m_y = []
    m_x_min = 0;
    m_x_max = 0;
    m_y_min = 0;
    m_y_max = 0;
    while m.any() and abs(m[0] - currentTime) < 1e-12:
   

        cov = np.array([ [ m[4], m[5] ], [ m[5], m[6] ] ]);
        w = m[7];
        eVal, eVec = np.linalg.eig(cov);
        eVal = eVal.real;
        a1 = 4*np.sqrt(eVal[0]); # Assume this is semi-major axis first
        a2 = 4*np.sqrt(eVal[1]); # 3 dof, 4 stdev is roughly prob = 0.997
        semiMajorAxis = eVec[:,0];
        if a2 > a1:
            aTmp = a1
            a1 = a2
            a2 = aTmp
            semiMajorAxis = eVec[:,1];
        a1Angle = np.arctan2(semiMajorAxis[1], semiMajorAxis[0]);
        
        m_x.append(m[2])
        m_y.append(m[3])
        landmarks[m_idx].set_alpha(min(w, 0.75));
        landmarks[m_idx].center = (m[2], m[3]);
        landmarks[m_idx].height = a2;
        landmarks[m_idx].width = a1;
        t_start = ax.transData;
        t_rot = transforms.Affine2D().rotate_around(m[2], m[3], a1Angle);
        t_compound = t_rot + t_start;
        landmarks[m_idx].set_transform(t_compound);
        drawnObjects.append(landmarks[m_idx]);

        if m[2] < m_x_min:
            m_x_min = m[2]
        if m[2] > m_x_max:
            m_x_max = m[2]
        if m[3] < m_y_min:
            m_y_min = m[3]
        if m[3] > m_y_max:
            m_y_max = m[3]        

        m_idx += 1;

        m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);

    #plt.xlim([m_x_min - 10, m_x_max + 10])
    #plt.ylim([m_y_min - 10, m_y_max + 10])
    #plt.gca().set_aspect('equal')

    while landmarks[m_idx].height != 0:
        landmarks[m_idx].set_alpha(0);
        landmarks[m_idx].center = (0, 0);
        landmarks[m_idx].height = 0;
        landmarks[m_idx].width = 0;
        m_idx += 1;

    landmarkCenters.set_data(m_x, m_y)
    drawnObjects.append(landmarkCenters)

    # Measurements
    nZ = 0;
    while z.any() and abs(z[0] -  currentTime) < 1e-12:
        z_dir = pr_best[i] + z[2] - np.pi / 2;
        z_end = [px_best[i] + z[1]*np.cos(z_dir), py_best[i] + z[1]*np.sin(z_dir) ];
        measurements[nZ].set_data([px_best[i], z_end[0]], [py_best[i], z_end[1]]);
        drawnObjects.append(measurements[nZ]);
        z = np.fromfile(measurementFileHandle, count=4, sep=" ", dtype=float);
        nZ += 1;
    while measurements[nZ].get_xdata() != []:
        measurements[nZ].set_data([], []);
        nZ += 1;
    
    drawnObjects.append(particles)
    drawnObjects.append(trajectory)

    return drawnObjects;

animation = anim.FuncAnimation(plt.figure(1), animate, np.arange(0, 7230), interval=1, 
                               init_func=animateInit, blit=True, repeat=False);

if saveMovie:
    FFMpegWriter = matplotlib.animation.writers['ffmpeg']
    animation.save(estimateMovieFile, writer=FFMpegWriter(fps = 30)) #extra_args=['-loglevel','quiet','-vcodec','libx264']
    #estPoseHandle, = plt.plot(px_best, py_best, 'b-');
    for i in range(0, nMeasurementsDrawMax) : 
        measurements[i].remove();
    #plt.setp(gtPoseHandle, linewidth=2.0)
    #plt.legend([gtPoseHandle, estPoseHandle, gtMapHandle, landmarks[0]], ["Ground-truth trajectory", "Estimated trajectory", "Ground-truth landmark", "Estimated landmark" ], loc=4);
    plt.legend([estPoseHandle, landmarks[0]], ["Estimated trajectory", "Estimated landmark" ], loc=4);
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
    plt.savefig(estimateImageFile, format='pdf', bbox_inches='tight')
else:
    plt.show()

measurementFileHandle.close();


