from pathlib import Path
import scanpy as sc

data_folder = Path("data")
data_folder.mkdir(parents=True,exist_ok=True)

print("InsightCell: Downloading real scRNA-seq data...")
adata = sc.datasets.pbmc3k()

output_file = data_folder / "pbmc3k_raw.h5ad"
adata.write_h5ad(output_file)

print("\nDataset downloaded successfully!") 
print("File location:", output_file.resolve()) 
print("Number of cells:", adata.n_obs) 
print("Number of genes:", adata.n_vars) 
print("\nDataset summary:")
print(adata)