#!/bin/bash -l
#SBATCH --job-name="Al-12"
#SBATCH --time=1-00:00:00
#SBATCH --account=mr22
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#========================================

############################################################################
# Variables definition
############################################################################
LAMMPS_HOME=/users/piaggip/bin/DAINT/LAMMPS/lammps-git/src
LAMMPS_RUN=${LAMMPS_HOME}/lmp_gnu-xc30
coresPerPart=6
partPerNode=6
cycles=2
Partitions=6
PartitionsMinus1=5
############################################################################

export PLUMED_USE_LEPTON=yes

############################################################################
# Run
############################################################################
if [ -e runno ] ; then
   #########################################################################
   # Restart runs
   #########################################################################
   nn=`tail -n 1 runno | awk '{print $1}'`
   srun -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK ${LAMMPS_RUN} -p ${Partitions}x${coresPerPart} -in restart.lmp

   #########################################################################
else
   #########################################################################
   # First run
   #########################################################################
   nn=1
   # Number of partitions
   srun -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK ${LAMMPS_RUN} -p ${Partitions}x${coresPerPart} -in start.lmp
   #########################################################################
fi
############################################################################

############################################################################
# Prepare next run
############################################################################
# Back up
for j in $(seq 0 ${PartitionsMinus1})
do
        cp log.lammps.${j} log.lammps.${j}.${nn}
        cp Restart2.lmp.${j} Restart2.lmp.${j}.${nn}
        cp Restart.lmp.${j} Restart.lmp.${j}.${nn}
        cp data.final.${j} data.final.${j}.${nn}
done

############################################################################
# Check number of cycles
############################################################################
mm=$((nn+1))
echo ${mm} > runno
#cheking number of cycles
if [ ${nn} -ge ${cycles} ]; then
  exit
fi
############################################################################

############################################################################
# Resubmitting again
############################################################################
sbatch < job.sh
############################################################################

