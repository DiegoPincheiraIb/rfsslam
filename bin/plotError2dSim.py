#!/usr/bin/python

 #
 # Software License Agreement (New BSD License)
 #
 # Copyright (c) 2013, Keith Leung, Felipe Inostroza
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
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print "Usage: plotError2dSim DATA_DIR\n";
    sys.exit(0);

# Open files
dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'
poseEstFile = 'poseEstError.dat';
poseEstFile = dataDir + poseEstFile;
if os.path.exists(poseEstFile):
    print('Opening ' + poseEstFile);
else:
    print(poseEstFile + ' does not exist')
    sys.exit(0);
mapEstFile = 'landmarkEstError.dat';
mapEstFile = dataDir + mapEstFile;
if os.path.exists(mapEstFile):
    print('Opening ' + mapEstFile);
else:
    print(mapEstFile + ' does not exist')
    sys.exit(0);


# Read File
print('Reading ' + poseEstFile);
poseEst = np.genfromtxt(poseEstFile);
poseTimesteps = poseEst[:,0];
poseErr_x = poseEst[:,1];
poseErr_y = poseEst[:,2];
poseErr_r = poseEst[:,3] * 180.0 / np.pi ;
poseErr_d = poseEst[:,4];
print('Reading ' + mapEstFile);
mapEst = np.genfromtxt(mapEstFile);
mapTimesteps = mapEst[:,0];
landmarksMeasured = mapEst[:,1];
landmarksEstimated = mapEst[:,2];
errorOSPA = mapEst[:,3];

# Generate plots
plt.figure(1);
p1, = plt.plot(poseTimesteps, poseErr_x, 'r-');
p2, = plt.plot(poseTimesteps, poseErr_y, 'b-');
plt.legend([p1, p2], ["x", "y"], loc=4);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel('Timestep');
plt.ylabel('Position error [m]');
plt.grid(True);
plt.ylim(-0.3, 0.3);

plt.figure(2);
plt.plot(poseTimesteps, poseErr_r, 'r-');
plt.xlabel('Timestep');
plt.ylabel('Rotation error [deg]');
plt.grid(True);
plt.ylim(-10, 10);

plt.figure(3);
plt.plot(poseTimesteps, poseErr_d, 'r-');
plt.xlabel('Timestep');
plt.ylabel('Position error [m]');
plt.grid(True);
plt.ylim(ymax = 0.3);

plt.figure(4);
plt.plot(mapTimesteps, errorOSPA, 'r-');
plt.xlabel('Timestep');
plt.ylabel('OSPA error');
plt.grid(True);
plt.ylim(ymax = 3);

plt.figure(5);
p1, = plt.plot(mapTimesteps, landmarksMeasured, 'k-', linewidth=3.0);
p2, = plt.plot(mapTimesteps, landmarksEstimated, 'r-');
plt.legend([p1, p2], ["Actual", "Estimated"], loc=4);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel('Timestep');
plt.ylabel('Number of landmarks observed');
plt.grid(True);
plt.ylim(ymax = 70);
#plt.show();

# Save plots
errorPosePosFile = dataDir + 'errorPosePos.pdf';
errorPoseRotFile = dataDir + 'errorPoseRot.pdf';
errorPoseDstFile = dataDir + 'errorPoseDst.pdf';
errorOSPAFile = dataDir + 'errorOSPA.pdf';
errorCardinalityFile = dataDir + 'errorCardinality.pdf';
plt.figure(1);
print('Saving  ' + errorPosePosFile);
plt.savefig(errorPosePosFile, format='pdf', bbox_inches='tight');
plt.figure(2);
print('Saving  ' + errorPoseRotFile);
plt.savefig(errorPoseRotFile, format='pdf', bbox_inches='tight');
plt.figure(3);
print('Saving  ' + errorPoseDstFile);
plt.savefig(errorPoseDstFile, format='pdf', bbox_inches='tight');
plt.figure(4);
print('Saving  ' + errorOSPAFile); 
plt.savefig(errorOSPAFile, format='pdf', bbox_inches='tight');
plt.figure(5);
print('Saving  ' + errorCardinalityFile); 
plt.savefig(errorCardinalityFile, format='pdf', bbox_inches='tight');


