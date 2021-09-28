#!/bin/bash
#SBATCH --job-name Some_Trial
#SBATCH --comment "Trial"
#SBATCH --time=0-00:01:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Mert.Kurttutan@physik.uni-muenchen.de
#SBATCH --chdir=/home/m/Mert.Kurttutan/Academia/Codes/Physics/Projects/qh_fm_01/codes/excs
#SBATCH --output=/project/th-scratch/m/Mert.Kurttutan/slurm/DMRG.%j.%N.out
#SBATCH --constraint=avx
#SBATCH --ntasks=4
#SBATCH --mem=4GB
#SBATCH --partition=th-cl,cluster,th-ws

__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# input: Lx Ly alpha U N S(spin) PBC job

source /project/theorie/s/Sam.Mardazad/Group/syten_inc.sh

#python ${__dir}/trial.py 
python ./trial.py 
