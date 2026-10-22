import scanpy as sc
from pathlib import Path

file_path = Path("data/pbmc3k_raw.h5ad")
if not file_path.exists():
    print("Dataset not found!")
    print("Pleasse run download_data.py first.")

else:
    ## load the dataset
    adata = sc.read_h5ad("data/pbmc3k_raw.h5ad")
    print("="*45)
    print(" INSIGHTCELL DATASET EXPLORATION")
    print("="*45)

    print("\n1. DATASET SIZE")
    print("Number of cells:", adata.n_obs)
    print("number of genes:",adata.n_vars)

    #First five cell identifiers
    print("\n2. FIRST FIVE CELL IDs") 
    print(adata.obs_names[:5].tolist()) 
    
    # First ten gene names print
    ("\n3. FIRST TEN GENES") 
    print(adata.var_names[:10].tolist())
     # Available cell metadata 
    print("\n4. CELL METADATA COLUMNS") 
    print(adata.obs.columns.tolist()) 
    #
    #  Available gene metadata 
    print("\n5. GENE METADATA COLUMNS") 
    print(adata.var.columns.tolist()) 
    # Data object summary 
    print("\n6. DATASET SUMMARY") 
    print(adata) 
    print("\nDataset exploration completed!")

