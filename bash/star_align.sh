#!/bin/bash
#SBATCH --job-name=star-align-dvl3 # Job name
#SBATCH -p normal
#SBATCH -t 07:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5   # 5 threads per task 
#SBATCH --mem=40G
#SBATCH --array=1-33%5   # Creates 1-x array jobs and runs max. y jobs (%y) at a time.
#SBATCH --output=/home/hu425279/PN0079_DVL3_rna-seq/bash/slurm-outs/star_align_trimmed_%A_%a.out
#SBATCH --mail-user=henna.urhonen@tuni.fi
#SBATCH --mail-type=FAIL

start=`date +%s`

echo "$(STAR --version)"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" ~/PN0079_DVL3_rna-seq/data/samples.txt`

# SET THESE VARIABLES
#fastqdir=/bmt-data/genomics/projects/dvl3_mouse_rna-seq/01_raw_data
fastqdir=/bmt-data/genomics/projects/dvl3_mouse_rna-seq/02_qc/trimmomatic_out_pcr_primers+truseq
outdir=/bmt-data/genomics/projects/dvl3_mouse_rna-seq/03_alignments/trimmed
IDX=/bmt-data/genomics/reference/Mus_musculus/STAR_GRCm39_index
GTF=/bmt-data/genomics/reference/Mus_musculus/gencode_release_M27_GRCm39/gencode.vM27.primary_assembly.annotation.gtf

[[ -d ${outdir} ]] || mkdir ${outdir}

echo "SAMPLE: ${sample}"

STAR --runThreadN 5 \
--genomeDir $IDX \
--sjdbGTFtagExonParentGene gene_id \
--sjdbGTFfile $GTF \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix ${outdir}/${sample}_ \
--readFilesCommand zcat \
--readFilesIn ${fastqdir}/${sample}_trimmed.fastq.gz  # Or _R1_001.fastq.gz

end=`date +%s`
runtime=$((end-start))
echo $runtime
