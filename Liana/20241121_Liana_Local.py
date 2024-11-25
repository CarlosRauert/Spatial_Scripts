C:\ProgramData\miniforge3\Scripts\activate
conda deactivate
conda activate Liana

import pandas as pd
import liana as li
import squidpy as sq
from mudata import MuData
import plotnine as p9
from plotnine import ggplot, aes, geom_point, ggtitle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel
from matplotlib import gridspec
from plotnine import data
from PIL import Image
import scanpy as sc
import os
import matplotlib.image as mpimg



# load data
adata=sc.read("C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/AnnData/adata_PostCC_PostMDE_800Epoch.h5ad")

#subtype
adata.obs["subtype"] = np.where(
    adata.obs["cluster_cellcharter"].isin([1, 3, 4, 5]), "DDLS",  # Condition 1
    np.where(
        adata.obs["cluster_cellcharter"] == 6, "WDLS",            # Condition 2
        "NA"                                                     # Default value
    )
)

# Misty

xR="2"
xType=0

Misty_All=pd.DataFrame()
subplot_dir="C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Misty_Out/Plots/"

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
        Misty_All=Misty_All._append(Misty_xR_sub)
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
    title_filename = f"C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Misty_Out/Plots/MistyLR_combined_{xR}.png"
    plt.savefig(title_filename, format="png", bbox_inches='tight', pad_inches=0, dpi=800)  # High DPI
    plt.close()  # Close the plot
    # Delete the individual subplot files after combining
    for subplot_file in plotlist:
        os.remove(subplot_file)

# Save Misty_All to csv
Misty_All.to_csv("C:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Misty_Out/MistyLR_AllRes.csv")


# Bivariate

xR="2"
tt_All=pd.DataFrame()
InterPlots=[]
InterDict={"init":0}
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
                nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True,
                inplace=False
                )
    tt_df=tt[0]
    tt_df["Region"]=xR
    tt_All=tt_All._append(tt_df)
    for interaction in tt[0]["interaction"].tolist():
        if interaction in InterDict.keys():
            pass
        else:
            ind=max(InterDict.values())+1
            InterDict.update({interaction:ind})
    #for interaction in tt[0]["interaction"].tolist():
    #    sq.pl.spatial_scatter(
    #        tt[1],
    #        color=interaction,
    #        img=None,
    ##        library_key='sample',
    #        spatial_key="spatial",
    #        library_id=xR,
    #        size=10,
    #    )
    # #   morans=value = tt[0].loc[tt[0]['interaction'] == interaction, 'morans'].iloc[0]
    #    morans=round(morans, 3)
    #    plt.suptitle(f"Region: {xR} morans: {morans:.3f}")
    #    indD=InterDict[interaction]
    #    plt.savefig(f"c:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Bivariate_Out/Plot_By_LR/Subplots/Subplot_{xR}_{indD}.png", dpi=2000)
    #    plt.close()
    #    print(f"Plot{interaction} done")
    #sq.pl.spatial_scatter(
    #    tt[1],
    #    color=tt[0]["interaction"].tolist(),
    #    img=None,
    #    library_key='sample',
    #    spatial_key="spatial",
    #    library_id = xR,
    #    size=10,
    #    save=f"c:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Bivariate_Out/Plots/LR_Cosine_Plot_Region{xR}.pdf"
    #)
    #plt.close()

# tt_all save
tt_All.to_csv("c:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Bivariate_Out/tt_All.csv")

# concatenate LR Subplots

# Define the directory where the PNG files are stored
base_dir = "c:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Bivariate_Out/Plot_By_LR/Subplots"
save_dir = "c:/Users/carlo/Documents/Charite/Promotion/PromotionSpatialLS/20241022_Liana/Bivariate_Out/Plot_By_LR"
# Iterate over the range of i (1 to 21)
for i in range(8, 22):
    images = []
    for j in [2, 3, 4, 5, 6, 7, 9]:
        file_path = os.path.join(base_dir, f"Subplot_{j}_{i}.png")
        if os.path.exists(file_path):
            images.append(file_path)  # Store file paths instead of loading all images
    if images:
        fig, axes = plt.subplots(1, len(images), figsize=(5 * len(images), 5))
        if len(images) == 1:
            axes = [axes]
        for ax, img_path in zip(axes, images):
            img = Image.open(img_path)
            img = img.resize((1280, 960))  # Resize to reduce memory usage
            ax.imshow(img)
            ax.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"Subplots_{i}.png"), dpi=1000)
        plt.close()