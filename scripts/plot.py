from gtf_utils import parse_gene_map

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import re

h5ad_dir = snakemake.input["h5ad"]
sample = snakemake.params["sample"]
knee_out = snakemake.output["knee"]
scatter_out = snakemake.output["scatter"]

adata = sc.read_h5ad(h5ad_dir)
print(f"Loaded h5ad file from {h5ad_dir}")

# First, mark out the mito genes
gene_map = parse_gene_map(snakemake.input["gtf"])
bare_ids = adata.var_names.str.replace(r"\.\d+$", "", regex=True)
bare_map = {re.sub(r"\.\d+$", "", k): v for k, v in gene_map.items()}

adata.var["gene_symbol"] = [bare_map.get(g, (g, ""))[0] for g in bare_ids]
adata.var["chromosome"] = [bare_map.get(g, (g, ""))[1] for g in bare_ids]
adata.var["mt"] = adata.var["chromosome"] == "chrM"
print(f"{adata.var['mt'].sum()} flagged as mito genes")

# Next, calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

# Knee Plot -> total counts per barcode
counts = np.sort(np.asarray(adata.obs["total_counts"], dtype=float))[::-1]
plt.figure(figsize=(6, 5))
plt.loglog(np.arange(1, counts.size + 1), np.maximum(counts, 1e-1))
plt.xlabel("Barcode rank")
plt.ylabel("Total UMI counts")
plt.title(f"{sample}: Barcode-Rank Knee Plot")
plt.tight_layout()
plt.savefig(knee_out, dpi=150)

# Scatter of counts vs. genes; colored by mito%
fig, ax = plt.subplots(figsize=(6.5, 5))
sc_plot = ax.scatter(
    adata.obs["total_counts"],
    adata.obs["n_genes_by_counts"],
    c=adata.obs["pct_counts_mt"],
    cmap="viridis",
    s=15,
    linewidths=0,
)
cbar = fig.colorbar(sc_plot, ax=ax, label="% mito counts")
ax.set_xlabel("total_genes")
ax.set_ylabel("n_genes_by_counts")
ax.set_title(f"{sample}: Counts vs. Genes")
fig.tight_layout()
fig.savefig(scatter_out, dpi=150, bbox_inches="tight")
plt.close(fig)
