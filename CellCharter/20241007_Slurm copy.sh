#!/bin/bash
#
#SBATCH --mem=600G
#SBATCH --ntasks=16
#SBATCH --time=2-00:00:00

conda init
conda activate NicheDE

Rscript $HOME/work/Scripts_Git_Repos/Spatial_Scripts/CellCharter/20241007_Slurm_Region3.R