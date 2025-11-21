def samples(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            parts  = line.strip().split(': ')
            key = parts[0]
            lineage = parts[1]
            mapping[key] = lineage
    return mapping

SAMPLE_DICT = samples(config["sample_map"])
SAMPLE = list(SAMPLE_DICT.keys())

def ncbi(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            parts  = line.strip().split(': ')
            key = parts[0]
            accession = parts[1]
            lineage = parts[2]
            mapping[key] = (accession,lineage)
    return mapping

NCBI_DICT = ncbi(config["ncbi_map"])
NCBI = list(NCBI_DICT.keys())

def all_input(wildcards):
    genes = get_genes()
    return {
        "compleasm": expand("qc/{sample}_compleasm", sample=SAMPLE),
        "gene_names": "gene.txt",
        "species_tree": "species_tree.tre"
    }


def minig(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            parts  = line.strip().split(': ')
            key = parts[0]
            g = parts[1].strip().split(',')
            g1 = g[0]
            g2 = g[1]
            div = g[2]
            name1 = g[3].strip('"')
            name2 = g[4].strip('"')
            cdna_key1 = g[5]
            cdna_key2 = g[6]
            mapping[key] = (g1, g2, div, name1, name2, cdna_key1, cdna_key2)
    return mapping

MINIG_DICT = minig(config["minimap_genomes"])
MINIG = list(MINIG_DICT.keys())

def minicDNA(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            parts  = line.strip().split(': ')
            key = parts[0]
            g = parts[1].strip().split(',')
            cDNA = g[0]
            g2 = g[1]
            mapping[key] = (cDNA, g2)
    return mapping

MINIC_DICT = minicDNA(config["minimap_cDNA"])
MINIC = list(MINIC_DICT.keys())

