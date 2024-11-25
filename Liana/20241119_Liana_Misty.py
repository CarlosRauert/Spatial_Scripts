import scanpy as sc
import decoupler as dc
import plotnine as p9
from plotnine import ggplot, aes, geom_point, ggtitle
import liana as li
import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel
from matplotlib import gridspec
from plotnine import data
from PIL import Image
import os

# load data
adata=sc.read("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241104_CC_Xenium_SCVIEpoch/20241106_800Epoch/adata_PostCC_800Epoch.h5ad")
# restrict to Region 2
adata_r2=adata[adata.obs['sample'] == "2"].copy()

# check key cell types
sq.pl.spatial_scatter(
    adata_r2,
    color=["cell_types"],
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    save=f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/spatial_plot.pdf"
)

# decoupler
# Generate PCA features
sc.tl.pca(adata_r2, svd_solver='arpack')
# Compute distances in the PCA space, and find spot neighbors
sc.pp.neighbors(adata_r2)
# Run leiden clustering algorithm and plot
sc.tl.leiden(adata_r2)
adata_r2.obs['leiden'] = ['Clust. {0}'.format(i) for i in adata_r2.obs['leiden']]
sq.pl.spatial_scatter(
    adata_r2,
    color="leiden",
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    save=f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/leiden_Plot.pdf"
)
# Spatial connectivity
li.ut.spatial_neighbors(
    adata_r2,
    bandwidth=150,
    cutoff=0.1,
    kernel='gaussian',
    set_diag=True,
    standardize=True
)
adata_r2.obs['conn'] = adata_r2.obsp['spatial_connectivities'][0].toarray().ravel()
sq.pl.spatial_scatter(
    adata_r2,
    color="conn",
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    save=f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/conn_Plot.pdf"
)
# Update X with spatially weighted gene exression
adata_r2_spW = adata_r2.copy()
adata_r2_spW.X = adata_r2_spW.obsp['spatial_connectivities'].toarray().dot(adata_r2_spW.X.toarray())
# compare spatially weighted and unweighted for OGN
fig, ax = plt.subplots(1, 2, figsize=(12, 6))  # 1 row, 2 columns
sq.pl.spatial_scatter(
    adata_r2,
    color="OGN",
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    ax=ax[0])
sq.pl.spatial_scatter(
    adata_r2_spW,
    color="OGN",
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    ax=ax[1])
plt.savefig("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/spW_Plot.pdf")
plt.close()
#tf activity inference
net = dc.get_collectri(organism='human', split_complexes=False)
# later

#back to misty
#estimate pathway activities
progeny = dc.get_progeny(organism='human', top=500) # obtain genesets
dc.run_mlm( # use multivariate linear model to estimate activity
    mat=adata_r2,
    net=progeny,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw=False,
)
acts_progeny = li.ut.obsm_to_adata(adata_r2, 'mlm_estimate') # extract progeny activities as an AnnData object
sq.pl.spatial_scatter(
    acts_progeny,
    color=["Hypoxia","JAK-STAT"],
    img=None,
    library_key='sample',
    spatial_key="spatial",
    library_id = "2",
    size=10,
    save=f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/Progeny_Plot.pdf"
)
misty = genericMistyData(intra=adata_r2, extra=acts_progeny, cutoff=0.05, bandwidth=200, n_neighs=6)
misty(model=RandomForestModel, n_jobs=128, verbose = True)

# ligand receptor misty
sc.pp.highly_variable_genes(adata_r2)  #For the sake of computational speed, let’s identify the highly variable genes (not used bc too few features)
hvg = adata_r2.var[adata_r2.var['highly_variable']].index
combined_resources=pd.DataFrame(columns=["ligand", "receptor"]) # combine all liana resources
for resource in li.resource.show_resources():
    combined_resources=pd.concat([combined_resources, li.resource.select_resource(resource)])
all_values = pd.unique(combined_resources.values.ravel()).tolist()
adata_subset = adata_r2[:, adata_r2.var_names.isin(all_values)]
misty=li.method.genericMistyData(adata_subset)
#misty = lrMistyData(adata_r2, bandwidth=200, set_diag=False, cutoff=0.01, nz_threshold=0.1, resource_name=None, resource=combined_resources) # Build LR Misty object
misty = lrMistyData(adata_r2, bandwidth=200, set_diag=False, cutoff=0.01, nz_threshold=0.1, resource_name="cellinker") # Build LR Misty object
misty(bypass_intra=True, model=LinearModel, verbose=True) # run speedy linear model
(
    li.pl.interactions(misty, view="extra", return_fig=True, figure_size=(6, 5), key=abs) +
    p9.scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    p9.labs(y='Receptor', x='Ligand')
).save("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Plots/MistyLR_Plot.pdf", format="pdf")


cellinker = li.resource.select_resource("cellinker")
subset_df = lrdb[lrdb["ligand"] == "VCAN"]
subset_df

# region loop
adata.obs["subtype"] = np.where(
    adata.obs["cluster_cellcharter"].isin([1, 3, 4, 5]), "DDLS",  # Condition 1
    np.where(
        adata.obs["cluster_cellcharter"] == 6, "WDLS",            # Condition 2
        "NA"                                                     # Default value
    )
)

xR="2"
xType=0

Misty_All=pd.DataFrame()

# Loop through each region
for xR in ["1", "2", "3", "4", "5", "6", "7", "9"]:
    adata_xR = adata[adata.obs['sample'] == xR].copy()
    plotlist = []
    subs = ["WDLS", "DDLS"]
    for xType in (0, 1):
        adata_xR_sub = adata_xR[adata_xR.obs['subtype'] == subs[xType]].copy()
        misty = lrMistyData(adata_xR_sub, bandwidth=200, set_diag=False, cutoff=0.01, nz_threshold=0.1, resource_name="cellinker")
        misty(bypass_intra=True, model=LinearModel, verbose=True)
        # create dataframe and append to Misty_All
        Misty_xR_sub=misty.uns['interactions']
        Misty_xR_sub["Region"]=xR
        Misty_xR_sub["subtype"]=subs[xType]
        Misty_All=Misty_All.append(Misty_xR_sub)
        # Create the plot for the current subtype
        p = li.pl.interactions(misty, view="extra", return_fig=True, figure_size=(6, 5), key=abs) + \
            p9.scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) + \
            p9.labs(y='Receptor', x='Ligand') + \
            p9.ggtitle(subs[xType])
        # Save each individual subplot to a temporary file
        subplot_filename = f"{subplot_dir}MistyLR_{xR}_{subs[xType]}.png"
        p.save(subplot_filename, format="png", dpi=800)  # High DPI for subplots
        plotlist.append(subplot_filename)
        plt.close()  # Close the figure after saving it
    # Create a combined image (side by side)
    img1 = Image.open(plotlist[0])  # WDLS plot
    img2 = Image.open(plotlist[1])  # DDLS plot
    img_combined = Image.new('RGB', (img1.width + img2.width, img1.height))  # Combined width
    img_combined.paste(img1, (0, 0))  # Paste the first image at the left
    img_combined.paste(img2, (img1.width, 0))  # Paste the second image at the right
    figsize = (img_combined.width / 800, (img_combined.height + 100) / 800)  # Inches for high DPI
    # Create a matplotlib figure to add the title
    fig, ax = plt.subplots(figsize=figsize, dpi=800)  # High DPI for publication quality
    ax.text(0.5, 0.95, f"Region {xR}", ha='center', va='center', fontsize=16, fontweight='bold', color='black')
    ax.axis('off')
    ax.imshow(img_combined)
    title_filename = f"/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Misty_Plots/MistyLR_combined_{xR}.png"
    plt.savefig(title_filename, format="png", bbox_inches='tight', pad_inches=0, dpi=800)  # High DPI
    plt.close()  # Close the plot
    # Delete the individual subplot files after combining
    for subplot_file in plotlist:
        os.remove(subplot_file)

# Save Misty_All to csv
Misty_All.to_csv("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/20241119_Liana/Misty_Plots/MistyLR_AllRes.csv")