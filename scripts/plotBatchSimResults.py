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
matplotlib.use("TkAgg");

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#print matplotlib.__version__

import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker   
import argparse

parser = argparse.ArgumentParser(description="Batch Simulation Results Plotter")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--saveFig", help="save last animation frame as a pdf file", action="store_true")
parser.add_argument("resultFile", help="path to result log file")
args = parser.parse_args()

# Open and read result log
if not os.path.exists(args.resultFile):
    print(args.resultFile + ' does not exist');
    sys.exit(0);
print('Reading ' + args.resultFile);
data = np.genfromtxt(args.resultFile)
np.set_printoptions(threshold=np.nan)
data = data[ np.lexsort((data[:,1], data[:,0])) ]

data_stat = np.empty((0,6), float)
Pd_current = data[0,0]
c_current = data[0,1]
trajErrors = [data[0,2]]
mapColas = [data[0,3]]
c_unique = []
for i in range(1, data.shape[0]):

    if data[i,0] == Pd_current and data[i,1] == c_current :
        trajErrors.append(data[i,2])
        mapColas.append(data[i,3])
    else:
        # Calculate statistics
        data_add = np.array([Pd_current, np.log10(c_current), np.mean(trajErrors), np.std(trajErrors), np.mean(mapColas), np.std(mapColas) ])
        data_stat = np.vstack( [data_stat, data_add])
        c_unique.append(np.log10(c_current))

        print trajErrors

        Pd_current = data[i,0]
        c_current = data[i,1]
        trajErrors = [data[i,2]]
        mapColas = [data[i,3]]

data_add = np.array([Pd_current, np.log10(c_current), np.mean(trajErrors), np.std(trajErrors), np.mean(mapColas), np.std(mapColas) ])
data_stat = np.vstack( [data_stat, data_add])

data_stat = data_stat[ np.lexsort((data_stat[:,0], data_stat[:,1])) ]

c_unique = np.unique(c_unique)

plt.figure(1);
for i in c_unique:
    data_ = data_stat[ (data_stat[:,1] == i) ]
    plt.plot(data_[:,0], data_[:,2], 'r-');
plt.grid(True)
plt.xlabel(r"Probability of detection")
plt.ylabel(r"Error [m]")

plt.figure(2);
for i in c_unique:
    data_ = data_stat[ (data_stat[:,1] == i) ]
    plt.plot(data_[:,0], data_[:,4], 'r-');
plt.grid(True)
plt.xlabel(r"Probability of detection")
plt.ylabel(r"COLA error")
plt.show()

#print data_stat.shape

#dtype = [('Pd',float), ('c',float), ('PoseErr',float), ('COLA',float)]
#data2 = np.array(data, dtype=dtype)
#np.sort(data2, axis=0, order=['Pd', 'c'])
#print data2
