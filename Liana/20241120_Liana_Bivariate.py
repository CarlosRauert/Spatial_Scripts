import pandas as pd
import scipy
import scanpy as sc
import decoupler as dc
import liana as li
import squidpy as sq

from mudata import MuData

# test vignette
adata = sc.read("kuppe_heart19.h5ad", backup_url='https://figshare.com/ndownloader/files/41501073?private_link=4744950f8768d5c8f68c')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
li.ut.spatial_neighbors(adata, bandwidth=200, cutoff=0.1, kernel='gaussian', set_diag=True)
tt=li.mt.bivariate(adata,
                resource_name='consensus', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True,
                inplace=False
                )

# load data and subset to R2
adata=sc.read("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241104_CC_Xenium_SCVIEpoch/20241106_800Epoch/adata_PostCC_800Epoch.h5ad")
adata_r2=adata[adata.obs['sample'] == "2"].copy()

#As part of the LIANA+ manuscript, we performed two distinct tasks to evaluate the ability of these metrics 
#to preserve biological information, and saw that on average when used to identify local ligand-receptor relationships, 
#spatially-weighted Cosine similarity did best. 

# find good bandwidth
plot, _ = li.ut.query_bandwidth(coordinates=adata_r2.obsm['spatial'], start=0, end=500, interval_n=20)
plot.save("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Bivariate/Plots/Bandwidth_Plot.pdf", format="pdf")
# 100 seems to include ~100 neighbors

li.ut.spatial_neighbors(adata_r2, bandwidth=100, cutoff=0.1, kernel='gaussian', set_diag=False)

tt=li.mt.bivariate(adata_r2,
                resource_name='cellinker', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True,
                inplace=False
                )

tt[1].var.sort_values("mean", ascending=False).head(3)
tt[1].var.sort_values("std", ascending=False).head(3)

tt[0]

InterPlots=[]

#region loop
for xR in ["2","3","4","5","6","7","9"]:
    adata_xR=adata[adata.obs['sample'] == xR].copy()
    li.ut.spatial_neighbors(adata_xR, bandwidth=100, cutoff=0.1, kernel='gaussian', set_diag=False)
    tt=li.mt.bivariate(adata_xR,
                resource_name='cellinker', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                #nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True,
                inplace=False
                )
    sq.pl.spatial_scatter(
        tt[1],
        color=tt[0],
        img=None,
        library_key='sample',
        spatial_key="spatial",
        library_id = xR,
        size=10,
        save=f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Bivariate/Plots/LR_Cosine_Plot_Region{xR}.pdf"
    )