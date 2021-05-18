#!/bin/bash

# Expects multiqc to be downloaded to a conda environment.
source activate py3-env

fastqc_outdir=$1
name=$2

cd $fastqc_outdir
mkdir -p $fastqc_outdir/multiqc_out  # Create MultiQC output folder.
multiqc -n $name -o $fastqc_outdir/multiqc_out *fastqc.zip
