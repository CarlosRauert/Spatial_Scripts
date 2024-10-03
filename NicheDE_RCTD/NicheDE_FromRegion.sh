#!/bin/sh

Region=$1

conda init
conda activate /data/cephfs-1/home/users/rauertc_c/work/miniconda/envs/seurat

Rscript /data/cephfs-2/unmirrored/projects/liposarcoma-wgs/Niche_DE_Rang20/Scripts/20240924_NicheDE_xR.R ${Region}

echo "done"