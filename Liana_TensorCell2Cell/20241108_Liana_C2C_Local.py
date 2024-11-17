C:\ProgramData\miniforge3\Scripts\activate

import cell2cell as c2c
import decoupler as dc
import liana as li

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

import matplotlib.pyplot as plt
import seaborn as sns

import muon as mu
import mofax as mofa

import plotnine as p9

from scipy.spatial.distance import squareform
import warnings
warnings.filterwarnings('ignore')

#Finally, for comparing tissues from multiple patients, a 4D tensor could be built for the same tissue region across patients. So the 4th dimension would be a 
#tissue region aligned and present across patients.

output_folder = 'C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/20241108_TensorC2C_Out'
c2c.io.directories.create_directory(output_folder)

#Load data
#Similar to the tutorial of using LIANA and MISTy, we will use an ischemic 10X Visium spatial slide from Kuppe et al., 2022. 
#It is a tissue sample obtained from a patient with myocardial infarction, specifically focusing on the ischemic zone of the heart tissue.
#The slide provides spatially-resolved information about the cellular composition and gene expression patterns within the tissue.

adata = sc.read("C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/AnnData/adata_PostCC_PostMDE_800Epoch.h5ad")
adata

# use cellcharter cluster as spatial context

li.mt.rank_aggregate.by_sample(
    adata,
    groupby='cell_types',
    resource_name='consensus', # Note: uses HUMAN gene symbols!
    sample_key='cluster_cellcharter', # sample key by which we which to loop
    expr_prop = 0.1,
    use_raw=False,
    n_perms=100, # reduce permutations for speed
    return_all_lrs=False, # we don't return all LR values to utilize MOFA's flexible views
    verbose=True, # use 'full' to show all information
    )

spatial_liana = adata.uns['liana_res']

#Using LIANA results to build a 4D-communication tensor, with grid windows as the contexts in the 4th dimension.
tensor = li.multi.to_tensor_c2c(liana_res=spatial_liana, # LIANA's dataframe containing results
                                sample_key='cluster_cellcharter', # Column name of the samples
                                source_key='source', # Column name of the sender cells
                                target_key='target', # Column name of the receiver cells
                                ligand_key='ligand_complex', # Column name of the ligands
                                receptor_key='receptor_complex', # Column name of the receptors
                                score_key='lr_means', # Column name of the communication scores to use
                                inverse_fun=None, # Transformation function
                                how='outer', # What to include across all samples
                                outer_fraction=1/4., # Fraction of samples as threshold to include cells and LR pairs.
                               )
#100%|██████████████████████████████████████████████████████████████| 7/7 [00:00<00:00, 27.49it/s]#Metadata for coloring the elements in the tensor.

#Here we do not assign major groups, each element is colored separately.

dimensions_dict = [None, None, None, None]
meta_tensor = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=dimensions_dict,
                                                  fill_with_order_elements=True
                                                 )
#Tensor properties

tensor.shape
#(7, 12, 17, 17)
tensor.excluded_value_fraction()
#0.9188498928983359
tensor.sparsity_fraction()
#0.0
#Run Tensor-cell2cell pipeline
#For simplicity, we factorize the tensor into 8 factors or CCC patterns instead of using the elbow analysis to determine the number of factors.

c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                           meta_tensor,
                                           rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis
                                           tf_optimization='regular', # To define how robust we want the analysis to be.
                                           random_state=0, # Random seed for reproducibility
                                           device='cuda', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                           cmaps=['plasma', 'Dark2_r', 'Set1', 'Set1'],
                                           output_folder="C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/20241108_TensorC2C_Out", # Whether to save the figures in files. If so, a folder pathname must be passed"
                                          )

#Running Tensor Factorization

# Generate columns in the adata.obs dataframe with the loading values of each factor
factor_names = list(tensor.factors['Contexts'].columns)
for f in factor_names:
    adata.obs[f] = adata.obs['cluster_cellcharter'].apply(lambda x: float(tensor.factors['Contexts'][f].to_dict()[x]))
    adata.obs[f] = pd.to_numeric(adata.obs[f])
# These columns are used for coloring the tissue regions to associate them with each factor
sq.pl.spatial_scatter(adata, color=factor_names, size=4, vmin=0., img=None,cmap="coolwarm", library_key="sample", spatial_key="spatial", library_id='2')
plt.suptitle(f'interactions by Cellcharter Cluster', fontsize=32, ha='left')
plt.savefig("C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/TensorC2C_Out/20241022_R2_Tensors_Consensus_elbow.pdf")

sq.pl.spatial_scatter(
    adata,
    color=["Factor 8"],
    img=None,
    #palette=custom_cmap,
    library_key='sample',
    spatial_key="spatial",
    library_id = ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
    size=10,
    ncols=4,
    title=['Region_1', 'Region_2', 'Region_3', 'Region_4', 'Region_5', 'Region_6', 'Region_7', 'Region_8', 'Region_9'], 
    save=f"C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/20241108_TensorC2C_Out/Factor8.pdf"
)

_ = c2c.plotting.loading_clustermap(loadings=tensor.factors['Ligand-Receptor Pairs'],
                                    loading_threshold=0.1,
                                    use_zscore=False,
                                    figsize=(28, 8),
                                    filename="C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/20241108_TensorC2C_Out/20241108_Tensors_LRpairs_Consensus.pdf",
                                    row_cluster=False
                                   )