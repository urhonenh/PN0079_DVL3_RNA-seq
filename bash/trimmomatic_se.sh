#!/bin/bash
#SBATCH -J trimmomatic_se
#SBATCH -p normal
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=henna.urhonen@tuni.fi
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-%j_trimmomatic_se.out

start=`date +%s`
module load compbio/fastqc/0.11.7

datafolder=$1
trimmed_folder=$2

mkdir -p $trimmed_folder  # -p flag: no error if existing, creates parent directories if needed.
mkdir -p $trimmed_folder/trim_summary

# Removed: LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
for i in $datafolder/*.fastq.gz
do
	samplename=$(basename $i _R1_001.fastq.gz);
	java -jar /bmt-data/genomics/apps/Trimmomatic-0.36/trimmomatic-0.36.jar SE \
		-threads ${SLURM_CPUS_PER_TASK} \
		$i $trimmed_folder/${samplename}_trimmed.fastq.gz \
		ILLUMINACLIP:/bmt-data/genomics/apps/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 \
		2> $trimmed_folder/trim_summary/${samplename}_trim_summary.log;
done

mkdir -p $trimmed_folder/fastqc_output
fastqc -t $SLURM_CPUS_PER_TASK -o $trimmed_folder/fastqc_output $trimmed_folder/*trimmed.fastq.gz

end=`date +%s`
runtime=$((end-start))
echo $runtime
