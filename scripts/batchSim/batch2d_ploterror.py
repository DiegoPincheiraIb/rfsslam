#!/usr/bin/python3

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
import os
import argparse
import subprocess 
import errno
import time
import signal
import shutil
import matplotlib
import numpy as np
matplotlib.use("TkAgg");


import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker

matplotlib.rcParams.update({'font.size': 20})

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


parser = argparse.ArgumentParser(description="plot the error of batch 1d results")
parser.add_argument("-r","--resultDir", help="path to result data", default='batchResults')


args = parser.parse_args()





if not os.path.exists(args.resultDir):
    os.makedirs(args.resultDir)


psopath = os.path.join(args.resultDir, 'rfspso2d' )
cerespath = os.path.join(args.resultDir, 'rfsceres2d' )

    

psoerrors= []
cereserrors = []
    
for folder in os.listdir(psopath):

    if not os.path.isdir(os.path.join(psopath,folder)):
        continue
    
    print(folder)
    error = np.genfromtxt(os.path.join(psopath,folder,'error.txt'))
    psoerrors.append(error);

for folder in os.listdir(cerespath):

    if not os.path.isdir(os.path.join(cerespath,folder)):
        continue
    print(folder)

    error = np.genfromtxt(os.path.join(cerespath,folder,'error.txt'))
    cereserrors.append(error);
 
psoavgerror = np.zeros(psoerrors[0].shape)
ceresavgerror = np.zeros(cereserrors[0].shape)

for error in psoerrors:
    print(error)
    psoavgerror += error
psoavgerror= psoavgerror/len(psoerrors)
for error in cereserrors:
    print(error)
    ceresavgerror += error
ceresavgerror= ceresavgerror/len(cereserrors)


xerrorpso = psoavgerror[:,0]
yerrorpso = psoavgerror[:,1]

ceresxerror = ceresavgerror[:,0]
ceresyerror = ceresavgerror[:,1]

print(xerrorpso)
gtPose_t = range(len(xerrorpso))


        
        
        
errorfig = plt.figure( figsize=(12,10), facecolor='w')
axError = errorfig.gca()
xpsoHandle, = axError.plot(gtPose_t , xerrorpso  ,'b-')
xceresHandle, = axError.plot(gtPose_t , ceresxerror  ,'r-')
ypsoHandle, = axError.plot(gtPose_t , yerrorpso  ,'b--')
yceresHandle, = axError.plot(gtPose_t , ceresyerror  ,'r--')
axError.legend([xpsoHandle ,ypsoHandle , xceresHandle ,yceresHandle ], ["PSO x error" ,"PSO y error" ,"Known DA x error" , "Known DA y error"  ], loc='best')
axError.set_xlabel("time [s]")
axError.set_ylabel("error [m]")
plt.savefig("2derror.pdf", format='pdf', bbox_inches='tight')



