#!/bin/bash

LANG=en_US # so that seq doesn't use commas as decimal separator
LC_NUMERIC=en_US.UTF-8 # so that seq doesn't use commas as decimal separator

cd ~/Projects/phdFilter
mkdir -p data/batch/rbphdslam

for Pd in $(seq 1.0 -0.1 0.1)
do
    if (( $(bc <<< "$Pd == 1.0") == 1 ))
    then
	Pd=0.99  
    fi

    for c in 0.0001 0.001 0.01 0.1 1 10
    do
	
	for trialseed in 0 1 2 3 4 5 6 7 8 28
	do

	    echo "Pd = $Pd   c = $c   seed = $trialseed"

	    # Edit xml config files
	    sed -e "s/<probDetection>.*<\/probDetection>/<probDetection>$Pd<\/probDetection>/" -e "s/<clutterIntensity>.*<\/clutterIntensity>/<clutterIntensity>$c<\/clutterIntensity>/" -e "s/<logDirPrefix>.*<\/logDirPrefix>/<logDirPrefix>data\/batch\/rbphdslam\/<\/logDirPrefix>/" cfg/rbphdslam2dSim.xml > data/batch/rbphdslam/rbphdslam2dSim.xml

	    # Run the simulator
	    build/bin/rbphdslam2dSim $trialseed data/batch/rbphdslam/rbphdslam2dSim.xml

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
