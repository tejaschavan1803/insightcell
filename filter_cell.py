from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

# Define file paths
input_file = Path("data/pbmc3k_raw.h5ad")
output_folder = Path("data")
results_folder = Path("results/plots")

results_folder.mkdir(parents=True, exist_ok=True)

# Load original dataset
adata = sc.read_h5ad(input_file)

print("=" * 45)
print("       INSIGHTCELL: QUALITY CONTROL")
print("=" * 45)

print("\nOriginal cells:", adata.n_obs)
print("Original genes:", adata.n_vars)

# Identify mitochondrial genes
adata.var["mt"] = (
    adata.var_names.str.upper().str.startswith("MT-")
)

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Display QC summary
print("\nQC metrics summary:")
print(
    adata.obs[
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    ].describe()
)

# Save QC plots
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.3,
    multi_panel=True,
    show=False
)

plt.savefig(
    results_folder / "qc_before_filtering.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close("all")

# Initial filtering thresholds
min_genes = 200
max_genes = 3000
max_mito_percent = 10

# Filter cells using QC metrics
keep_cells = (
    (adata.obs["n_genes_by_counts"] >= min_genes)
    & (adata.obs["n_genes_by_counts"] < max_genes)
    & (adata.obs["pct_counts_mt"] < max_mito_percent)
)

adata_filtered = adata[keep_cells].copy()

print("\nCells before filtering:", adata.n_obs)
print("Cells after filtering:", adata_filtered.n_obs)
print("Cells removed:", adata.n_obs - adata_filtered.n_obs)

if adata_filtered.n_obs == 0:
    raise ValueError(
        "No cells passed filtering. Review the QC summary and thresholds."
    )

# Save filtered dataset
output_file = output_folder / "pbmc3k_filtered.h5ad"
adata_filtered.write_h5ad(output_file)

print("\nFiltered dataset saved to:", output_file.resolve())
print("\nQC plot saved to:", (results_folder / "qc_before_filtering.png").resolve())
print("\nQuality control completed!")