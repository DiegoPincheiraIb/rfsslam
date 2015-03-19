#!/bin/bash

cd ~/Projects/phdFilter
mkdir -p data/batch/rbphdslam

for Pd in 0.99 0.9 #$(seq 1.0 -0.1 0.1)
do
    if (( $(bc <<< "$Pd == 1.0") == 1 ))
    then
	Pd=0.99  
    fi

    for c in 0.0001 #0.0001 0.001 0.01 0.1 1 10
    do
	
	for trial in {1..10}
	do

	    echo "Pd = $Pd   c = $c"

	    # Edit xml config files
	    sed -e "s/<probDetection>.*<\/probDetection>/<probDetection>$Pd<\/probDetection>/" -e "s/<clutterIntensity>.*<\/clutterIntensity>/<clutterIntensity>$c<\/clutterIntensity>/" -e "s/<logDirPrefix>.*<\/logDirPrefix>/<logDirPrefix>data\/batch\/rbphdslam\/<\/logDirPrefix>/" cfg/rbphdslam2dSim.xml > data/batch/rbphdslam/rbphdslam2dSim.xml

	    # Run the simulator
	    build/bin/rbphdslam2dSim 28 data/batch/rbphdslam/rbphdslam2dSim.xml

	    # Analyze results
	    build/bin/analysis2dSim data/batch/rbphdslam/

	    # Get the position and map errors at the very end
	    posError=$(tail -n 1 data/batch/rbphdslam/poseEstError.dat | tr -s ' ' | cut -d ' ' -f5)
	    mapError=$(tail -n 1 data/batch/rbphdslam/landmarkEstError.dat | tr -s ' ' | cut -d ' ' -f4)

	    # Write results to file
	    echo "$Pd   $c   $posError   $mapError" >> data/batch/batch_results_rbphdslam.dat

	done

    done

done