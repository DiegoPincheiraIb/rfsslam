#!/usr/bin/python

import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print "Usage: plotError2dSim DATA_DIR\n";
    sys.exit(0);

# Open file
dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'
dataFile = 'landmarkEstError.dat';
dataFile = dataDir + dataFile;
if os.path.exists(dataFile):
    print('Opening ' + dataFile);
else:
    print(dataFile + ' does not exist')
    sys.exit(0);

# Read File
print('Reading ' + dataFile);
data = np.genfromtxt(dataFile);
timesteps = data[:,0];
landmarksMeasured = data[:,1];
landmarksEstimated = data[:,2];
errorOSPA = data[:,3];

# Generate plots
plt.figure(1);
plt.plot(timesteps, errorOSPA, 'r-');
plt.xlabel('Timestep');
plt.ylabel('OSPA Error');
plt.ylim(ymax = 3);

plt.figure(2);
p1, = plt.plot(timesteps, landmarksMeasured, 'k-', linewidth=3.0);
p2, = plt.plot(timesteps, landmarksEstimated, 'r-');
plt.legend([p1, p2], ["Actual", "Estimated"], loc=4);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel('Timestep');
plt.ylabel('Number of landmarks observed');
plt.ylim(ymax = 70);
#plt.show();

# Save plots
errorOSPAFile = dataDir + 'errorOSPA.pdf';
errorCardinalityFile = dataDir + 'errorCardinality.pdf';
plt.figure(1);
print('Saving  ' + errorOSPAFile); 
plt.savefig(errorOSPAFile, format='pdf', bbox_inches='tight');
plt.figure(2);
print('Saving  ' + errorCardinalityFile); 
plt.savefig(errorCardinalityFile, format='pdf', bbox_inches='tight');


