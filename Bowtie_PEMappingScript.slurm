#!/bin/bash
#SBATCH --job-name=map_eQTL_sex_int
#SBATCH -q batch 
#SBATCH --mem 64GB # memory pool for all cores
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH -t 1-0:00 # time (D-HH:MM)
# example use: sbatch /projects/munger-lab/ArderyProject/Scripts/Bowtie_PEMappingScript.slurm

module load singularity
echo "Starting the Rscript!"
echo "singularity exec docker://rocker/tidyverse:4.0.3 Rscript /projects/munger-lab/ArderyProject/src/RNAseq_comp/Bowtie_PEMappingScript.R"

singularity exec docker://rocker/tidyverse:4.0.3 Rscript /projects/munger-lab/ArderyProject/src/RNAseq_comp/Bowtie_MappingScript.R

echo "Done!"