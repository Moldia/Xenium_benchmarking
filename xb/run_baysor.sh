#!/bin/bash

#------------------------------------------------------------
# How to run the script

# 1. Set the parameters below
# 2. Get the docker image louisk92/txsim_baysor:v0.6.2bin:
#    docker pull louisk92/txsim_baysor:v0.6.2bin
# 3. Activate a docker container and bind mount the input and output directories:
#    docker run -it --rm -v /path/to/data:/path/to/data -v /path/to/script:/path/to/script louisk92/txsim_baysor:v0.6.2bin
#    (take data dirs from the parameters below: XENIUM_DIR(, OUT_DIR), and script dir from the path to this script)
#    Note: if you can't run the docker container you can convert it to a singularity/apptainer container.
# 4. Run this script inside the docker container:
#    cd /path/to/run/baysor/bash/script
#    bash run_baysor.sh

#------------------------------------------------------------
# User-defined parameters

# Set scale parameter here instead of in the .toml file (leads to errors)




SCALE=20
TOML_FILE="baysor_params.toml" # Path to the .toml file, adjust parameters in the toml file as needed
PRIOR=True # Whether to use a prior segmentation
PRIOR_CONFIDENCE=0.5


XENIUM_DIR=$1 # Path to the input directory
MOLECULES=$XENIUM_DIR"/spots.csv" # Path to the spots csv with coordinates scaled to pixel space
PRIOR_SEGMENTATION=$XENIUM_DIR"/label_image.tif" # Path to the prior segmentation .tif file
OUT_DIR=$XENIUM_DIR # Path to the output directory

#------------------------------------------------------------

if [[ $SCALE > 0 ]]; then
    BAYSOR_CLI="run -s $SCALE -c $TOML_FILE -o $OUT_DIR/ $MOLECULES"
else
    BAYSOR_CLI="run -c $TOML_FILE -o $OUT_DIR/ $MOLECULES"
fi

if [[ $PRIOR == True ]]; then
    BAYSOR_CLI+=" --prior-segmentation-confidence $PRIOR_CONFIDENCE $PRIOR_SEGMENTATION"
fi

echo "Running Baysor with the following command:"
echo $BAYSOR_CLI
baysor $BAYSOR_CLI
