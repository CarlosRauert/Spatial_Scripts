#!/bin/bash
#
#SBATCH --mem=800G
#SBATCH --ntasks=12
#SBATCH --time=2-00:00:00

conda init
conda activate NicheDE

Rscript $HOME/work/Scripts_Git_Repos/Spatial_Scripts/CellCharter/20241007_NicheDE_Slurm.R