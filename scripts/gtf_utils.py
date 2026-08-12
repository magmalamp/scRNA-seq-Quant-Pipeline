import re


def parse_gene_map(gtf_path):
    gene_map = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            attrs = fields[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            if gid:
                gene_map[gid.group(1)] = (
                    gname.group(1) if gname else gid.group(1),
                    chrom,
                )
    return gene_map
