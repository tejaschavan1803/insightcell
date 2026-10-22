from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

## dataset load
input_file = Path("data/pbmc3k_raw.h5ad")
output_folder = Path("results/plots")
output_folder.mkdir(parents=True,exist_ok=True)

if not input_file.exists():
    raise FileNotFoundError(
        "Dataset not found. Run download_data.py first"
    )

adata = sc.read_h5ad(input_file)
print("InsightCell: Quality Control")
print("-"*40)
print("cell before QC:",adata.n_obs)
print("cell before QC:",adata.n_vars)
print("Genes before QC:",adata.n_vars)

## mt = mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

##QC matrics

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("\nQC metrics calculated successfully!")
print("\nFirst five cells:")
print(
    adata.obs[
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt"
    ]

    ].head()
)

## QC plots

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
     output_folder / "qc_violin_plots.png",
     dpi=300, 
     bbox_inches="tight" ) 
plt.close("all")

##scatter
sc.pl.scatter(
    adata,
    x="total_counts",
    y="pct_counts_mt",
    show=False

)

plt.savefig(
    output_folder / "qc_counts_vs_mitochondrial.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close("all")


## save the qc metrics
metrics_file = Path("results/qc_metrics.csv")
metrics_file.parent.mkdir(parents=True, exist_ok=True)

adata.obs[
    [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt"
    ]

].to_csv(metrics_file)

print("\nQuality-control analysis completed!") 
print("Violin plots:", output_folder / "qc_violin_plots.png") 
print( "Scatter plot:", output_folder / "qc_counts_vs_mitochondrial.png" ) 
print("QC metrics:", metrics_file)