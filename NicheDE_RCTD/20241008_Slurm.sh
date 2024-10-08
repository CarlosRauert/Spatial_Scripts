#!/bin/bash
#
#SBATCH --mem=1500G
#SBATCH --ntasks=32
#SBATCH --time=2-00:00:00

conda init
conda activate NicheDE

Rscript $HOME/work/Scripts_Git_Repos/Spatial_Scripts/NicheDE_RCTD/20241008_NicheDE_all_Slurm.R