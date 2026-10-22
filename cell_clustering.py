from pathlib import Path
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Define paths
input_file = Path("data/pbmc3k_pca_umap.h5ad")
results_folder = Path("results")
plot_folder = results_folder / "plots"

results_folder.mkdir(exist_ok=True)
plot_folder.mkdir(parents=True, exist_ok=True)

# Load processed dataset
adata = sc.read_h5ad(input_file)

print("=" * 45)
print(" INSIGHTCELL: CELL CLUSTERING")
print("=" * 45)

print("\nCells:", adata.n_obs)
print("Genes:", adata.n_vars)

# Perform Leiden clustering
sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="leiden",
    flavor="igraph",
    n_iterations=2,
    directed=False
)

print("\nCells in each cluster:")
print(adata.obs["leiden"].value_counts().sort_index())

# Visualize clusters on UMAP
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    show=False
)

plt.savefig(
    plot_folder / "cell_clusters_umap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close("all")

# Identify marker genes for each cluster
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)

# Export marker gene results
markers = sc.get.rank_genes_groups_df(
    adata,
    group=None
)

markers.to_csv(
    results_folder / "cluster_marker_genes.csv",
    index=False
)

# Display top 5 markers per cluster
top_markers = (
    markers.groupby("group", sort=False)
    .head(5)
)

print("\nTop marker genes by cluster:")
print(
    top_markers[["group", "names", "logfoldchanges", "pvals_adj"]]
    .to_string(index=False)
)

# Save dataset with cluster assignments
output_file = Path("data/pbmc3k_clustered.h5ad")
adata.write_h5ad(output_file)

print("\nUMAP plot saved to:", plot_folder / "cell_clusters_umap.png")
print("Marker genes saved to:", results_folder / "cluster_marker_genes.csv")
print("Clustered dataset saved to:", output_file.resolve())

print("\nPhase 8 completed!")