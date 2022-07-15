#!/bin/bash
#SBATCH --job-name=rsem_index
#SBATCH -p compute # partition (queue)
#SBATCH -q batch 
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 64GB # memory pool for all cores
#SBATCH -t 0-1:00 # time (D-HH:MM)

module load singularity

RSEMINDEX="/projects/munger-lab/RSEM/R84-REL1505/C57BL6J/star_index"
FASTAFILES="/projects/munger-lab/STAR/R84-REL1505/C57BL6J/C57BL6J.fa"
gtfFILE="/projects/munger-lab/STAR/R84-REL1505/C57BL6J/C57BL6J.gtf"

singularity run library://jaxreg.jax.org/rnaseq-modelad/rsem:v1.3.3 rsem-prepare-reference --gtf $gtfFILE $FASTAFILES C57BL6J

echo "Done! Yaaay!"