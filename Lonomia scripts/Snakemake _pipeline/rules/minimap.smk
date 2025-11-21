import glob

def find_genome(iid):
    sample_path = f"{iid}.fa"
    if os.path.exists(sample_path):
        return sample_path
    ncbi_pattern = f"ncbi/{iid}_phylogeny/ncbi_dataset/data/*/*fna"
    ncbi_matches = glob.glob(ncbi_pattern)   
    if ncbi_matches:
        return ncbi_matches[0]
    raise FileNotFoundError(f"No genome file found for {iid}")


rule minimap_genomes_faidx: 
    input:
        g1 = lambda wildcards: find_genome( MINIG_DICT[wildcards.minig][0]), 
        g2 = lambda wildcards: find_genome( MINIG_DICT[wildcards.minig][1]), 
    output: 
        "minimap/{minig}_faidx.done"
    conda:
        config["environments"]["samtools"]
    resources:
        resources=config["default_resources"]
    shell:
        """
        samtools faidx {input.g1}
        samtools faidx {input.g2}
        
        touch {output}
        
        """

rule minimap_genomes: 
    input:
        g1 = lambda wildcards: find_genome( MINIG_DICT[wildcards.minig][0]), 
        g2 = lambda wildcards: find_genome( MINIG_DICT[wildcards.minig][1]), 
    output: 
        "minimap/{minig}.sam"
    params: 
        div = lambda wildcards: MINIG_DICT[wildcards.minig][2]
    conda:
        config["environments"]["minimap"]
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        
        minimap2 -ax {params.div} -r 5000 -p 0.3 -N 1 {input.g1} {input.g2} > {output}
        
        """

rule minimap_genomes_converted: 
    input:
        "minimap/{minig}.sam" 
    output: 
        "minimap/{minig}.paf"
    conda:
        config["environments"]["minimap"]
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        paftools.js sam2paf {input} > {output}
        
        """

rule minimap_dotplot: 
    input:
        "minimap/{minig}.paf"
    output: 
        "minimap_dotplot/{minig}.png"
    params: 
        synteny = config["scripts"]["synteny_paf2"]
    conda:
        config["environments"]["r-base"]
    resources:
        resources=config["default_resources"]
    shell:
        """

        Rscript {params.synteny} --input {input} --output {output} -l 1000 -f 100 -q 20
        
        """

#####################


rule minimap_cDNA: 
    input:
        cDNA = lambda wildcards: MINIC_DICT[wildcards.minic][0] + ".fasta", 
        g2 = lambda wildcards: find_genome( MINIC_DICT[wildcards.minic][1]), 
    output: 
        "minimap_cDNA/{minic}.sam"
    conda:
        config["environments"]["minimap"]
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        
        minimap2 -ax splice -uf {input.g2} {input.cDNA} > {output}
        
        """

rule minimap_cDNA_converted: 
    input:
        "minimap_cDNA/{minic}.sam"
    output: 
        "minimap_cDNA/{minic}.paf"
    conda:
        config["environments"]["minimap"]
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        paftools.js sam2paf {input} > {output}
        
        """


rule minimap_dotplot_cDNA: 
    input:
        gg = "minimap/{minig}.paf",
        cdna_target = lambda wildcards: f"minimap_cDNA/{MINIG_DICT[wildcards.minig][5]}.paf",
        cdna_query = lambda wildcards: f"minimap_cDNA/{MINIG_DICT[wildcards.minig][6]}.paf"
    output: 
        "minimap_dotplot2/{minig}.png"
    params: 
        synteny = config["scripts"]["synteny_paf2"],
        cdna_id = config["cdna_id"],
        name1 = lambda wildcards: MINIG_DICT[wildcards.minig][3],
        name2 = lambda wildcards: MINIG_DICT[wildcards.minig][4]
    conda:
        config["environments"]["r-base"]
    resources:
        resources=config["default_resources"]
    shell:
        """
        Rscript {params.synteny} \
        --input {input.gg} \
        --cdna_paf {input.cdna_target} \
        --cdna_paf_query {input.cdna_query} \
        --cdna_id {params.cdna_id} \
        --x_label "{params.name1}" \
        --y_label "{params.name2}" \
        --output {output} \
        --verbose -l 1000 -f 100 -q 20
        """
