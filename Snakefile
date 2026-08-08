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

def get_r1_urls(wildcards):
  return config["samples"][wildcards.sample]["r1"]

def get_r2_urls(wildcards):
  return config["samples"][wildcards.sample]["r2"]

rule fetch_and_combine:
  output:
    r1="data/reads/{sample}_R1.fastq.gz",
    r2="data/reads/{sample}_R2.fastq.gz"
  params:
    r1_urls=get_r1_urls,
    r2_urls=get_r2_urls
  shell:
    r"""
    i=0
    for url in {params.r1_urls}; do
      curl -L -o {output.r1}.lane$i "$url"; i=$((i+1))
    done
    cat {output.r1}.lane* > {output.r1} && rm {output.r1}.lane*

    i=0
    for url in {params.r2_urls}; do
      curl -L -o {output.r2}.lane$i "$url"; i=$((i+1))
    done
    cat {output.r2}.lane* > {output.r2} && rm {output.r2}.lane*
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

def get_r2_urls(wildcards):
  return config["samples"][wildcards.sample]["r2"]

# finally run salmon alevin
rule seq_to_ref_map:
  input:
    idx_loc="data/splici/idx",
    fq_r1="data/reads/{sample}_R1.fastq.gz",
    fq_r2="data/reads/{sample}_R2.fastq.gz",
    tg_map="data/splici/transcriptome_splici_fl85_t2g_2col.tsv"
  output:
    directory("alevin-output/{sample}/{sample}_map")
  params:
    lib="A",
    chem="chromium"   # "chromium", "chromiumV3", or "dropseq"
  threads: 4
  shell:
    r"""
      salmon alevin -i {input.idx_loc} -p {threads} -l {params.lib} \
        --sketch -1 {input.fq_r1} -2 {input.fq_r2} -o {output} \
        --tgMap {input.tg_map} --dumpFeatures --{params.chem}
    """

rule gen_permit_list:
  input:
    "alevin-output/{sample}/{sample}_map"
  output:
    directory("alevin-output/{sample}/{sample}_quant")
  params:
    ori="both"    # "fw", "rc", or "both"
  shell:
    r"""
      alevin-fry generate-permit-list --input {input} --output-dir {output} \
        --expected-ori {params.ori} --knee-distance
    """

rule collate_and_quant:
  input:
    quant="alevin-output/{sample}/{sample}_quant",
    map="alevin-output/{sample}/{sample}_map",
    t2g_3col="data/splici/transcriptome_splici_fl85_t2g_3col.tsv"
  output:
    directory("alevin-output/{sample}/{sample}_count")
  threads: 4
  shell:
    r"""
      alevin-fry collate -t {threads} -i {input.quant} -r {input.map} --compress
      alevin-fry quant -t {threads} -i {input.quant} -o {output} --tg-map {input.t2g_3col} \
        --resolution cr-like --use-mtx
    """
