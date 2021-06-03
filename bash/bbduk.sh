#!/bin/bash
#SBATCH -J bbduk
#SBATCH -p normal
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=henna.urhonen@tuni.fi
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=slurm-%j_bbduk.out

start=`date +%s`
module load compbio/fastqc/0.11.7

# INPUT ARGUMENTS:
datafolder=$1
trimmed_folder=$2

bbmapdir=/bmt-data/genomics/apps/bbmap

mkdir -p $trimmed_folder  # -p flag: no error if existing, creates parent directories if needed.
mkdir -p $trimmed_folder/trim_summary

for i in $datafolder/*.fastq.gz
do
	samplename=$(basename $i _R1_001.fastq.gz);
	cat $i | $bbmapdir/bbduk.sh in=stdin.fastq.gz out=$trimmed_folder/${samplename}_trimmed.fastq.gz \
	int=f ref=/home/hu425279/PN0079_DVL3_rna-seq/data/illumina_pcr_primers+truseq_adapters.fa \
	k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlen=20 t=4 \
	2> $trimmed_folder/trim_summary/${samplename}_bbduk.log;
done

mkdir -p $trimmed_folder/fastqc_out
fastqc -t $SLURM_CPUS_PER_TASK -o $trimmed_folder/fastqc_out $trimmed_folder/*trimmed.fastq.gz

end=`date +%s`
runtime=$((end-start))
echo $runtime