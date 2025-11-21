rule compleasm: 
    input: 
        "{sample}.fa"
    output: 
        directory("qc/{sample}_compleasm")
    params:
        MODE= "genome",
        lineage = lambda wildcards: SAMPLE_DICT[wildcards.sample]
    conda:
        config["environments"]["compleasm"]
    resources:
        resources=config["default_resources_8thread"]
    log: 
        "logs/busco_quickmerge/{sample}"
    benchmark:
        "benchmarks/busco_quickmerge/{sample}"
    shell:
        """

        compleasm run -a {input} -o {output} -l {params.lineage} -t 8 

        """

rule compleasm_ncbi: 
    input: 
        "ncbi/{ncbi}_phylogeny/ncbi_dataset"
    output: 
        directory("phylogeny/compleasm/{ncbi}")
    params:
        MODE= "genome",
        lineage = lambda wildcards: NCBI_DICT[wildcards.ncbi][1]
    conda:
        config["environments"]["compleasm"]
    resources:
        resources=config["default_resources_8thread"]
    log: 
        "logs/busco_ncbi/{ncbi}"
    benchmark:
        "benchmarks/busco_ncbi/{ncbi}"
    shell:
        """

        compleasm run -a {input}/data/*/*fna -o {output} -l {params.lineage} -t 8 

        """