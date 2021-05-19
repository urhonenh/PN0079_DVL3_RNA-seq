#!/bin/bash
#SBATCH -p normal 	# partition name
#SBATCH -t 3:00:00 	# hours:minutes runlimit after which job will be killed
#SBATCH --cpus-per-task=5 	# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --job-name STAR-index 	# Job name
#SBATCH --output=slurm-%j_star_index.out	# File to which standard out and error log will be written
#SBATCH --mail-user=henna.urhonen@tuni.fi
#SBATCH --mail-type=FAIL

genomedir=/bmt-data/genomics/reference/Mus_musculus
mkdir -p $genomedir/STAR_GRCm39_index

echo "$(STAR --version)"

STAR --runMode genomeGenerate \
	--genomeDir $genomedir/STAR_GRCm39_index \
	--genomeFastaFiles $genomedir/gencode_release_M27_GRCm39/GRCm39.primary_assembly.genome.fa \
	--runThreadN $SLURM_CPUS_PER_TASK \
	--sjdbGTFfile $genomedir/gencode_release_M27_GRCm39/gencode.vM27.primary_assembly.annotation.gtf \
	--limitGenomeGenerateRAM 15000000000 \
	--limitIObufferSize 50000000
