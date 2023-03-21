#!/bin/bash

#SBATCH --job-name=500cells           # Assign an short unique name to your job
#SBATCH --requeue
#SBATCH --nodes=1                   # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2           # Cores per task (>1 if multithread tasks)
#SBATCH --mem=25600                  # Real memory (RAM) required (MB)
#SBATCH --time=08:00:00             # Total run time limit (HH:MM:SS) // Days: 1-3:0:0
#SBATCH --array=1-1000
#SBATCH --output=/dev/null #%A_%a.out          # STDOUT output file
#SBATCH --error=/dev/null #%A_%a.err           # STDERR output file (optional)
#SBATCH --mail-user=XXXX@scarletmail.rutgers.edu 
#SBATCH --mail-type=BEGIN,END,FAIL

sleep $((RANDOM%30+1))

# load modulde
module purge
module load singularity/3.8.3


# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment
cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
export R_LIBS=/home/XXXX/R/rocker-rstudio/4.0
exec rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export SINGULARITY_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export SINGULARITYENV_RSTUDIO_SESSION_TIMEOUT=0


# run script
singularity exec /home/XXXX/rstudio_4.0.4.sif Rscript simulation-2conditions.R \
$SLURM_ARRAY_TASK_ID \
--gene_perc=50 \
--factor_lower=1.5 \
--factor_upper=3.0 \
--ncells_per_type=500
