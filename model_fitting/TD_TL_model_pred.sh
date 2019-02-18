#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=58:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=TD_TL_model_pred
#SBATCH --mail-type=END
#SBATCH --mail-user=alm652@nyu.edu
#SBATCH --output=slurm_%j.out

#NAME="ADHD"
#WORKDIR="${SCRATCH}/${NAME}"

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab/2015b

export MATLABPATH=${MATLABPATH}:/${HOME}/${NAME}/matlab:${HOME}/MATLAB
source ${HOME}/MATLAB/setpath.sh

#Check if running as an array job

if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
        IID=${SLURM_ARRAY_TASK_ID}
fi

# Run the program

#echo ${WORKDIR} ${NAME} ${IID}.job

cat<<EOF | matlab -nodisplay
addpath(genpath('${HOME}/MATLAB'));
#cd('${WORKDIR}');
TD_TL_model_pred($IID);
EOF




