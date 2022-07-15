#!/bin/bash
#SBATCH --job-name=gbrs_one_slurm
#SBATCH -p compute # partition (queue)
#SBATCH -q batch 
#SBATCH -N 1 # number of nodes
#SBATCH -n 12 # number of cores
#SBATCH --mem 64GB # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
# example use: sbatch --export=ALL,SAMPLE="PB367_45_ES_B",SEX='F',GENERATION='16',INDIR="/projects/munger-lab/raw/DO_mESC",OUTDIR="/projects/munger_lab/ArderyProject/PE_Data" /projects/munger-lab/ArderyProject/Scripts/BowtiePE_RunOneSample.slurm

# check that variables SAMPLE, INDIR and OUTDIR exist
if [ -z "$SAMPLE" ]; then echo "SAMPLE is unset"; exit 1; fi # e.g. DO-1290
if [ -z "$SEX" ]; then echo "SEX is unset"; exit 1; fi # M or F
if [ -z "$GENERATION" ]; then echo "GENERATION is unset"; exit 1; fi # e.g. G8
if [ -z "$INDIR" ]; then echo "INDIR is unset"; exit 1; fi # FASTQ.GZ files folder
if [ -z "$OUTDIR" ]; then echo "OUTDIR is unset"; exit 1; fi # output folder

module load singularity
cd $OUTDIR
export GBRS_DATA=/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/hmm

echo "PE Bowtie Intersect Sample Rerun"

echo "0) Merging FASTQ files and setting variables"
echo `date`

#TEMPDIR=`mktemp -d`
TEMPDIR="/fastscratch/--/temp"
echo "Files in Temporal Directory ${TEMPDIR}"
echo "Sample $SAMPLE"

# merge gzipped fastq files into one FASTQ file
SEQFILE_R1="${TEMPDIR}/${SAMPLE}_R1.fastaq"
zcat `ls ${INDIR}/Read1/${SAMPLE}*` > ${SEQFILE_R1}

SEQFILE_R2="${TEMPDIR}/${SAMPLE}_R2.fastaq"
zcat `ls ${INDIR}/Read2/${SAMPLE}*` > ${SEQFILE_R2}

# 1) bowtie options
INDEXBASE="/projects/munger-lab/Genomes/bowtie1/transcripts"
OUTFILE1a="$TEMPDIR/$SAMPLE.8-way.transcriptome.sam"
OUTFILE1b="$TEMPDIR/$SAMPLE.8-way.transcriptome.bam"

# 2) gbrs bam2emase options
LIDFILE="/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/ref/emase.transcripts.info"
SUFFICES=A:B:C:D:E:F:G:H
OUTFILE2="$TEMPDIR/$SAMPLE.8-way.transcriptome.h5"
INDEX_DTYPE=uint32
DATA_DTYPE=uint8

# 3) gbrs compress options
OUTFILE2c="$TEMPDIR/$SAMPLE.8-way.transcriptome.compressed.h5"

# 4) gbrs quantify_multiway options - quantifying multiway allele specificity
GRPFILE="/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/ref/emase.gene2transcripts.tsv"
LENFILE="/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/8-way/emase.pooled.transcripts.info"
MODEL=4
OUTFILEBASE="$TEMPDIR/$SAMPLE"
#READLEN=100
#PSEUDOCOUNT=0.0
#MAXITERS=999
#TOLERANCE=0.0001

# Note the output files are:
# $OUTFILEBASE.multiway.genes.expected_read_counts
# $OUTFILEBASE.multiway.genes.tpm
# and the two corresponding *.isoforms.* files

# 5) gbrs reconstruct options
EXPRFILE="$OUTFILEBASE.multiway.genes.tpm"
AVECFILE="/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/hmm/avecs.npz"
TPROBFILE="/projects/churchill-lab/resource/Custom_Genomes/R84-REL1505/hmm/tranprob.DO.G${GENERATION}.${SEX}.npz"
GPOSFILE="/projects/munger-lab/gbrs/ref.gene_pos.ordered.npz"
EXPRCUTOFF=1.5
SIGMA=0.12

# What are the output files?
# .genoprobs.npz ?
# .genotypes.tsv

# 6) gbrs quantify2 options
GENOFILE="$OUTFILEBASE.genotypes.tsv"
OUTFILEBASE2="$TEMPDIR/$SAMPLE"

echo "1) Bowtie alignment"
echo `date`

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/bowtie:v1.1.2 -q -a --threads 12 --sam -v 3 ${INDEXBASE} -1 ${SEQFILE_R1} -2 ${SEQFILE_R2} > ${OUTFILE1a}

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/samtools:v0.1.18 view -bS -F 4 ${OUTFILE1a} > ${OUTFILE1b}

echo "2) BAM to EMASE conversion"
echo `date`

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/gbrs:v0.1.6 bam2emase -i ${OUTFILE1b} -m ${LIDFILE} -s ${SUFFICES//\:/\,} -o ${OUTFILE2} 

echo "3) Compress"
echo `date`

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/gbrs:v0.1.6 compress -i ${OUTFILE2} -o ${OUTFILE2c}

echo "4) Quantify multiway"
echo `date`

#run-emase -i ${OUTFILE2} -g ${GRPFILE} -M ${MODEL} -o ${OUTFILE3} -r ${READLEN} -p ${PSEUDOCOUNT} -m ${MAXITERS} -t ${TOLERANCE} -L ${LENFILE}
singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/gbrs:v0.1.6 quantify -i ${OUTFILE2c} -g ${GRPFILE} -L ${LENFILE} -M ${MODEL} -o ${OUTFILEBASE} --report-alignment-counts

echo "5) Reconstruct"
echo `date`

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/gbrs:v0.1.6 reconstruct -e ${EXPRFILE} -t ${TPROBFILE} -x ${AVECFILE} -g ${GPOSFILE} -c ${EXPRCUTOFF} -s ${SIGMA} -o ${OUTFILEBASE}

echo "6) Quantify"
echo `date`

singularity run library://jaxreg.jax.org/rna-seq-do-gbrs/gbrs:v0.1.6 quantify -i ${OUTFILE2c} -G ${GENOFILE} -g ${GRPFILE} -L ${LENFILE} -M ${MODEL} -o ${OUTFILEBASE2} --report-alignment-counts




echo "3) MOVE RESULTS FROM TEMP DIRECTORY TO OUTPUT DIRECTORY"
echo `date`

cp "$TEMPDIR/$SAMPLE.8-way.transcriptome.bam" "$OUTDIR"
cp "$TEMPDIR/$SAMPLE.8-way.transcriptome.sam" "$OUTDIR" 
 
cp "$TEMPDIR/$SAMPLE.8-way.transcriptome.h5" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.multiway.genes.expected_read_counts" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.multiway.genes.tpm" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.multiway.isoforms.expected_read_counts" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.multiway.isoforms.tpm" "$OUTDIR" 

cp "$TEMPDIR/$SAMPLE.diploid.genes.expected_read_counts" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.diploid.genes.tpm" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.diploid.isoforms.expected_read_counts" "$OUTDIR" 
cp "$TEMPDIR/$SAMPLE.diploid.isoforms.tpm" "$OUTDIR" 

cp "$TEMPDIR/$SAMPLE.genotypes.tsv" "$OUTDIR"
cp "$TEMPDIR/$SAMPLE.genoprobs.npz" "$OUTDIR"


rm "${TEMPDIR}/${SAMPLE}_R1.fastaq"
rm "${TEMPDIR}/${SAMPLE}_R2.fastaq"