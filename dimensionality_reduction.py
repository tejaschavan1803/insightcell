from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

# File paths
input_file = Path("data/pbmc3k_normalized.h5ad")
plot_folder = Path("results/plots")

plot_folder.mkdir(parents=True, exist_ok=True)

# Load normalized dataset
adata = sc.read_h5ad(input_file)

print("=" * 45)
print(" INSIGHTCELL: DIMENSIONALITY REDUCTION")
print("=" * 45)

print("\nCells:", adata.n_obs)
print("Genes:", adata.n_vars)

# Use highly variable genes for PCA
adata_pca = adata[:, adata.var["highly_variable"]].copy()

# Scale gene expression values
sc.pp.scale(adata_pca, max_value=10)

# Perform Principal Component Analysis
sc.tl.pca(adata_pca, svd_solver="arpack")

print("\nPCA completed successfully!")

# Save PCA visualization
sc.pl.pca(
    adata_pca,
    show=False
)

plt.savefig(
    plot_folder / "pca_plot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close("all")

# Build nearest-neighbor graph
sc.pp.neighbors(
    adata_pca,
    n_neighbors=15,
    n_pcs=30
)

print("Nearest-neighbor graph created!")

# Calculate UMAP coordinates
sc.tl.umap(adata_pca)

print("UMAP completed successfully!")

# Save UMAP visualization
sc.pl.umap(
    adata_pca,
    show=False
)

plt.savefig(
    plot_folder / "umap_plot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close("all")

# Save processed data
output_file = Path("data/pbmc3k_pca_umap.h5ad")
adata_pca.write_h5ad(output_file)

print("\nPCA plot: results/plots/pca_plot.png")
print("UMAP plot: results/plots/umap_plot.png")
print("Processed data saved to:", output_file.resolve())

print("\nPhase 7 completed!")