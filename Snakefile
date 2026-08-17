configfile: "config.yaml"

# These first two rules are just to grab the needed data
rule get_ref:
  output:
    data=temp("data/ref/ref.tar.gz"),
    genome="data/ref/genome.fa",
    gtf="data/ref/genes.gtf"
  params:
    url="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"
  shell:
    r"""
      curl -L -o {output.data} {params.url}
      tar -xzf {output.data} -O refdata-gex-GRCh38-2020-A/fasta/genome.fa > {output.genome}
      tar -xzf {output.data} -O refdata-gex-GRCh38-2020-A/genes/genes.gtf > {output.gtf}
    """

rule fetch_and_combine:
  output:
    r1="data/reads/{sample}_R1.fastq.gz",
    r2="data/reads/{sample}_R2.fastq.gz"
  params:
    tar_url=lambda wc: config["samples"][wc.sample]["tar"],
    tmp_sample="data/reads/tmp_{sample}"
  shell:
    r"""
      mkdir -p {params.tmp_sample}
      curl -L -o {params.tmp_sample}/fastqs.tar {params.tar_url}
      tar -xf {params.tmp_sample}/fastqs.tar -C {params.tmp_sample}
      cat $(find {params.tmp_sample} -name '*_R1_*.fastq.gz' | sort) > {output.r1}
      cat $(find {params.tmp_sample} -name '*_R2_*.fastq.gz' | sort) > {output.r2}
      rm -rf {params.tmp_sample}
      """
    

rule make_splici:
  input:
    genome="data/ref/genome.fa",
    gtf="data/ref/genes.gtf"
  output:
    fasta="data/splici/transcriptome_splici_fl85.fa",
    t2g_3col="data/splici/transcriptome_splici_fl85_t2g_3col.tsv"
  params:
    read_length=90,
    flank_trim=5,
    outdir="data/splici"
  shell:
    r"""
      pyroe make-splici {input.genome} {input.gtf} {params.read_length} {params.outdir} \
        --flank-trim-length {params.flank_trim} \
        --filename-prefix transcriptome_splici --dedup-seqs
    """

# needed index file for salmon alevin
rule generate_index:
  input:
    "data/splici/transcriptome_splici_fl85.fa"
  output:
    directory("data/splici/idx")
  threads: 4
  shell:
    r"""
      salmon index -t {input} -i {output} -k 31 -p {threads}
    """

# grab first two columns of t2g_3col for t2g_2col
rule t2g_2col:
  input:
    "data/splici/transcriptome_splici_fl85_t2g_3col.tsv"
  output:
    "data/splici/transcriptome_splici_fl85_t2g_2col.tsv"
  shell:
    r"""
      cut -f1,2 {input} > {output}
    """

# finally run salmon alevin
rule seq_to_ref_map:
  input:
    idx_loc="data/splici/idx",
    fq_r1="data/reads/{sample}_R1.fastq.gz",
    fq_r2="data/reads/{sample}_R2.fastq.gz",
    tg_map="data/splici/transcriptome_splici_fl85_t2g_2col.tsv"
  output:
    directory("output/alevin/{sample}/{sample}_map")
  params:
    lib="A",
    chem="chromiumV3"   # "chromium", "chromiumV3", or "dropseq"
  threads: 4
  shell:
    r"""
      salmon alevin -i {input.idx_loc} -p {threads} -l {params.lib} \
        --sketch -1 {input.fq_r1} -2 {input.fq_r2} -o {output} \
        --tgMap {input.tg_map} --dumpFeatures --{params.chem}
    """

rule gen_permit_list:
  input:
    "output/alevin/{sample}/{sample}_map"
  output:
    directory("output/alevin/{sample}/{sample}_quant")
  params:
    ori="both"    # "fw", "rc", or "both"
  shell:
    r"""
      alevin-fry generate-permit-list --input {input} --output-dir {output} \
        --expected-ori {params.ori} --knee-distance
    """

rule collate_and_quant:
  input:
    quant="output/alevin/{sample}/{sample}_quant",
    map="output/alevin/{sample}/{sample}_map",
    t2g_3col="data/splici/transcriptome_splici_fl85_t2g_3col.tsv"
  output:
    directory("output/alevin/{sample}/{sample}_count")
  threads: 4
  shell:
    r"""
      alevin-fry collate -t {threads} -i {input.quant} -r {input.map} --compress
      alevin-fry quant -t {threads} -i {input.quant} -o {output} --tg-map {input.t2g_3col} \
        --resolution cr-like --use-mtx
    """

rule load_fry:
  input:
    counts="output/alevin/{sample}/{sample}_count"
  output:
    h5ad="output/anndata/{sample}.h5ad"
  params:
    format="scRNA",     # "scRNA", "S+A", "raw", "U+S+A", "snRNA", "velocity" or "all"
    sample="{sample}"
  script:
    "scripts/load_fry.py"

rule QC_and_plots:
  input:
    h5ad="output/anndata/{sample}.h5ad",
    gtf="data/ref/genes.gtf"
  output:
    knee="plots/{sample}/knee.png",
    scatter="plots/{sample}/counts_vs_genes.png"
  params:
    sample="{sample}"
  script:
    "scripts/plot.py"
