import scanpy as sc
from pathlib import Path

file_path = Path("data/pbmc3k_raw.h5ad")

adata = sc.read_h5ad(file_path)

print("=" * 45)
print("       INSIGHTCELL DATASET CHECK")
print("=" * 45)

print("\nOriginal dataset:")
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

print("\nExpression matrix:")
print("Minimum value:", adata.X.min())
print("Maximum value:", adata.X.max())

print("\nExisting cell metadata:")
print(adata.obs.columns.tolist())

print("\nFirst five cells:")
print(adata.obs.head())

print("\nDataset check completed.")


filtered_path = Path("data/pbmc3k_filtered.h5ad")

if filtered_path.exists():
    filtered = sc.read_h5ad(filtered_path)

    print("\nFiltered dataset:")
    print("Cells:", filtered.n_obs)
    print("Genes:", filtered.n_vars)