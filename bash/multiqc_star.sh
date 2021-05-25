#!/bin/bash

# Expects multiqc to be downloaded to a conda environment.
source activate py3-env

star_outdir=$1
name=$2

cd $star_outdir
mkdir -p $star_outdir/multiqc_out  # Create MultiQC output folder.
multiqc -n $name -o $star_outdir/multiqc_out *Log.final.out
