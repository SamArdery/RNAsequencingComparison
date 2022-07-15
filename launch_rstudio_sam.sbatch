#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --job-name=rstudio-arders  
#SBATCH --output=/projects/munger-lab/ArderyProject/logs/rstudio-%j.log  
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=512GB
#SBATCH --export=ALL
### SLURM HEADER
localcores=${SLURM_CPUS_PER_TASK}
simg_path="/projects/munger-lab/ArderyProject/sif/rstudio_renv.simg" 
work_dir="/projects/munger-lab/ArderyProject/work" 
set -euo pipefail
#cd ${work_dir}

# start modifications
# adding these to fix the bind error
# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
mkdir -p ${work_dir}/run ${work_dir}/tmp 
mkdir -p -m 700 ${work_dir}/var/lib/rstudio-server
cat > ${work_dir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${work_dir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
export R_LIBS_USER=/projects/munger-lab/ArderyProject/R/rocker-rstudio/4.1
exec rsession "\${@}"
END

chmod +x ${work_dir}/rsession.sh

export SINGULARITY_BIND="${work_dir}/run:/run,${work_dir}/tmp:/tmp,${work_dir}/database.conf:/etc/rstudio/database.conf,${work_dir}/rsession.sh:/etc/rstudio/rsession.sh,${work_dir}/var/lib/rstudio-server:/var/lib/rstudio-server"

## end modifications


export PASSWORD=$(openssl rand -base64 15)
PORT=$(shuf -i8899-11999 -n1)
hostname_with_port=$(echo $(hostname -A):${PORT} | tr -d " ")
cat 1>&2 <<EOF
Login to ${hostname_with_port} with
  username: ${USER}
  password: ${PASSWORD}
EOF

cd ${work_dir}
module load singularity
singularity exec ${simg_path} rserver \
    --www-port ${PORT} \
    --auth-none=0 --auth-pam-helper-path=pam-helper \
    --auth-timeout-minutes=0 \
    --auth-stay-signed-in-days=30 
    
