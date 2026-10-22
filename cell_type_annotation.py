from pathlib import Path
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Define paths

input_file = Path("data/pbmc3k_clustered.h5ad")
results_folder = Path("results")
plot_folder = results_folder / "plots"

results_folder.mkdir(exist_ok=True)
plot_folder.mkdir(parents=True, exist_ok=True)

# Load clustered dataset

adata = sc.read_h5ad(input_file)

print("=" * 45)
print(" INSIGHTCELL: CELL-TYPE ANNOTATION")
print("=" * 45)

print("\nCells:", adata.n_obs)
print("Genes:", adata.n_vars)

# Common marker genes for PBMC immune-cell populations

marker_genes = {
"T_cells": ["CD3D", "CD3E", "IL7R"],
"B_cells": ["MS4A1", "CD79A", "CD74"],
"Monocytes": ["LYZ", "S100A8", "S100A9"],
"NK_cells": ["NKG7", "GNLY", "PRF1"],
"Platelets": ["PPBP", "PF4"],
"Dendritic_cells": ["FCER1A", "CST3"]
}

# Keep marker genes that exist in the dataset

available_markers = {
cell_type: [
gene for gene in genes
if gene in adata.var_names
]
for cell_type, genes in marker_genes.items()
}

print("\nAvailable marker genes:")
for cell_type, genes in available_markers.items():
    print(f"{cell_type}: {genes}")

# Plot marker expression across clusters

genes_to_plot = sorted({
gene
for genes in available_markers.values()
for gene in genes
})

if genes_to_plot:
    sc.pl.dotplot(
adata,
var_names=genes_to_plot,
groupby="leiden",
standard_scale="var",
show=False
)


plt.savefig(
    plot_folder / "cell_type_marker_dotplot.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close("all")

sc.pl.umap(
    adata,
    color=genes_to_plot[:6],
    ncols=3,
    show=False
)

plt.savefig(
    plot_folder / "marker_gene_umap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close("all")


# Export marker genes for biological interpretation

marker_table = pd.DataFrame([
{"Candidate_Cell_Type": cell_type, "Marker_Gene": gene}
for cell_type, genes in available_markers.items()
for gene in genes
])

marker_table.to_csv(
results_folder / "known_cell_type_markers.csv",
index=False
)

print("\nMarker table saved to:")
print(results_folder / "known_cell_type_markers.csv")

print("\nReview the plots before assigning cell types.")
print("Phase 9 completed!")
