#!/bin/bash
n_events=100
n_particles=100
moms=(0.1 1)
degrees=(0)

BUILD_DIR=../../BELLA-traccc_build

for deg in "${degrees[@]}"
do
    for p in "${moms[@]}"
    do
	command="
	${BUILD_DIR}/bin/do_telescope_simulation
	--gen-events=${n_events} 				      
	--gen-nparticles=${n_particles} 
	--gen-theta=${deg}:${deg}
	--gen-mom-gev=${p}:${p}
	--gen-phi-degree=0:0
	--output-directory=${PWD}/sim_data/${p}_GeV_${deg}_theta/
	--field-file=${PWD}/bfield/BELLA_field.cvf
    "
	${command}

	command="
	${BUILD_DIR}/bin/do_truth_fitting_momentum_residual 
	--detector-file=telescope_detector_geometry.json
	--material-file=telescope_detector_material.json
	--input-directory=${PWD}/sim_data/${p}_GeV_${deg}_theta/
	--input-events=${n_events}
	--use-detray-detector
	--field-file=${PWD}/bfield/BELLA_field.cvf
	"
	${command}

	# Move the CSV file to data directory
	command="
 	mv residual.csv 
	../data/residual_${p}_GeV_${deg}_theta.csv
	"
	${command}

    done
done

# Remove the simulation file
rm -rf ${PWD}/sim_data
rm -rf ${PWD}/telescope_detector_geometry.json