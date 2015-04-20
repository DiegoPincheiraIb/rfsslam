#!/bin/bash

#script to run simulation with and without amplitude information 

mkdir data/runs
for i in {1..5}
do
  build/bin/rbphdslam2dSim_amplitude
  mv data/rbphdslam data/runs/amp$i
  bin/animate2dSim.py data/runs/amp$i &
  
  build/bin/rbphdslam2dSim
  mv data/rbphdslam data/runs/noamp$i
  bin/animate2dSim.py data/runs/noamp$i &
  
done

