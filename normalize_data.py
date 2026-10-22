from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

# File paths
input_file = Path("data/pbmc3k_filtered.h5ad")
output_file = Path("data/pbmc3k_normalized.h5ad")
plot_folder = Path("results/plots")

plot_folder.mkdir(parents=True, exist_ok=True)

# Load filtered dataset
adata = sc.read_h5ad(input_file)

print("=" * 45)
print(" INSIGHTCELL: NORMALIZATION AND GENE SELECTION")
print("=" * 45)

print("\nCells:", adata.n_obs)
print("Genes:", adata.n_vars)

# Normalize each cell to 10,000 total counts
sc.pp.normalize_total(adata, target_sum=10000)

# Log-transform normalized counts
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

# Count highly variable genes
hvg_count = int(adata.var["highly_variable"].sum())

print("\nHighly variable genes identified:", hvg_count)

# Save highly variable gene list
hvg = adata.var.loc[
    adata.var["highly_variable"]
].copy()

hvg.to_csv("results/highly_variable_genes.csv")

# Visualize highly variable genes
sc.pl.highly_variable_genes(adata, show=False)

plt.savefig(
    plot_folder / "highly_variable_genes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close("all")

# Save normalized dataset
adata.write_h5ad(output_file)

print("\nNormalized dataset saved to:", output_file.resolve())
print("HVG list saved to: results/highly_variable_genes.csv")
print("HVG plot saved to:", (plot_folder / "highly_variable_genes.png").resolve())

print("\nPhase 6 completed successfully!")