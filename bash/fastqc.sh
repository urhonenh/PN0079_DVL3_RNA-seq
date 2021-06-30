#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p normal
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=henna.urhonen@tuni.fi
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=1G
#SBATCH --output=slurm-%j_fastqc.out

trimmed_folder=$1

fastqc_outdir=$trimmed_folder/fastqc_out_nogroup
mkdir -p $fastqc_outdir

fastqc -t $SLURM_CPUS_PER_TASK --nogroup -o $fastqc_outdir $trimmed_folder/*trimmed.fastq.gz
