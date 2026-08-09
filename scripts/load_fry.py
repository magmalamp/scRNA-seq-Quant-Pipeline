from pyroe import load_fry

count_dir = snakemake.input["counts"]
format = snakemake.params["format"]
sample = snakemake.params["sample"]

# Generate anndata object
adata = load_fry(count_dir, output_format=format)
adata.obs["sample"] = sample
print(f"[{sample}] {adata.shape[0]} cells by {adata.shape[1]} genes.")

# Download counts as h5ad
adata.write_h5ad(snakemake.output["h5ad"])
print(f"[{sample}] written to {snakemake.output['h5ad']}")
