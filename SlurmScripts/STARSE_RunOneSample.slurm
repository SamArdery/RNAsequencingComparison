#!/bin/bash
#SBATCH --job-name=star_aligner
#SBATCH -p compute # partition (queue)
#SBATCH -q batch 
#SBATCH -N 1 # number of nodes
#SBATCH -n 12 # number of cores
#SBATCH --mem 64GB # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
# example use: sbatch --export=ALL,SAMPLE="PB367_45_ES_B",SEX='F',GENERATION='16',INDIR="/projects/munger-lab/raw/DO_mESC/Read1",OUTDIR="/projects/munger-lab/ArderyProject/STAR/" /projects/munger-lab/ArderyProject/Scripts/STARSE_RunOneSample.slurm


# check that variables SAMPLE, INDIR and OUTDIR exist
if [ -z "$SAMPLE" ]; then echo "SAMPLE is unset"; exit 1; fi # e.g. DO-1290
if [ -z "$SEX" ]; then echo "SEX is unset"; exit 1; fi # M or F
if [ -z "$GENERATION" ]; then echo "GENERATION is unset"; exit 1; fi # e.g. G8
if [ -z "$INDIR" ]; then echo "INDIR is unset"; exit 1; fi # FASTQ.GZ files folder
if [ -z "$OUTDIR" ]; then echo "OUTDIR is unset"; exit 1; fi # output folder


module load singularity # load singularity
cd $OUTDIR
# export GBRS_DATA=/projects/munger-lab/raw/DO_mESC/Read1

echo "0) Merging FASTQ files and setting variables"
echo `date`

#TEMPDIR=`mktemp -d`
TEMPDIR="/fastscratch/--/temp" # make sure my fastscratch directory exists
echo "Files in Temporal Directory ${TEMPDIR}"
echo "Sample $SAMPLE"

# merge gzipped fastq files into one FASTQ file
SEQFILE="${TEMPDIR}/${SAMPLE}.fastaq"
zcat `ls ${INDIR}/${SAMPLE}*` > $SEQFILE

echo "STAR files ready"

#file for RSEM
OUTFILE2="$OUTDIR$SAMPLE*.bam"

echo "RSEM files ready"

echo "1) Starting STAR run"
echo `date` 

STARINDEX="/projects/munger-lab/STAR/R84-REL1505/C57BL6J/star_index"
gtfFILE="/projects/munger-lab/STAR/R84-REL1505/C57BL6J/C57BL6J.gtf"

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/star:v020201 --genomeDir $STARINDEX --readFilesIn $SEQFILE --runThreadN 12 --outFileNamePrefix $OUTDIR$SAMPLE --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile $gtfFILE

echo "done with STAR"

echo "2) Starting RSEM"
echo `date`

singularity run library://jaxreg.jax.org/rnaseq-modelad/rsem:v1.3.3 rsem-calculate-expression --alignments -p 8 --strandedness reverse $OUTFILE2 /projects/munger-lab/RSEM/R84-REL1505/C57BL6J/C57BL6J $SAMPLE 

echo "done with RSEM"
