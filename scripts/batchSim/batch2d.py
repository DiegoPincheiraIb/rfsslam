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



# Should probably replace most the OS calls with native python commands


parser = argparse.ArgumentParser(description="Run rfsceresslam2dSim and rfspsoslam2dSim  with several trajectories")
parser.add_argument("-n", "--numberOfRuns", help="number of monte carlo runs to execute",type=int, default=5, metavar='NUM')
parser.add_argument("-p","--plot", help="plot the results of each simulation", action="store_true")
parser.add_argument("-o","--resultDir", help="path to result data", default='batchResults')


args = parser.parse_args()





if not os.path.exists(args.resultDir):
    os.makedirs(args.resultDir)


psopath = os.path.join(args.resultDir, 'rfspso2d' )
cerespath = os.path.join(args.resultDir, 'rfsceres2d' )
if not os.path.exists(cerespath):
    os.makedirs(cerespath )
if not os.path.exists(psopath):
    os.makedirs(psopath )
    

if not os.path.exists('data'):
    os.makedirs('data')
    
for i in range(args.numberOfRuns):


    if not os.path.exists('data/rfsceres2dslam'):
        os.makedirs('data/rfsceres2dslam')
    logfile =  open('data/rfsceres2dslam/output.txt', 'w')
    ceresslam = subprocess.run(['../rfs-slam-build/bin/rfsceresslam2dSim','-t',str(i)], stdout=logfile)
    logfile.close()
    runFolder =os.path.join(cerespath,'run_'+str(i))
    shutil.move('data/rfsceres2dslam',runFolder)
    subprocess.Popen(['./scripts/sim/animate2dceres.py',runFolder])
    
        
    if not os.path.exists('data/rfspso2dslam'):
        os.makedirs('data/rfspso2dslam')
    
    
    logfile =  open('data/rfspso2dslam/output.txt', 'w')
    psoslam = subprocess.run(['../rfs-slam-build/bin/rfspsoslam2dSim','-t',str(i)], stdout=logfile)
    logfile.close()
    runFolder =os.path.join(psopath,'run_'+str(i))
    shutil.move('data/rfspso2dslam',runFolder)
    subprocess.Popen(['./scripts/sim/animate2dpso.py',runFolder])
        
        

        
        
        
        
        
    



