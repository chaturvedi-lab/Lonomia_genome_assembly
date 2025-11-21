rule blast_cDNA: 
    input:
        query = lambda wildcards: MINIC_DICT[wildcards.minic][0] + ".fasta", 
        genome = lambda wildcards: find_genome( MINIC_DICT[wildcards.minic][1]), 
    output: 
        "blast2/{minic}.txt"
    params:
        script = config["scripts"]["blast"]
    conda:
        config["environments"]["blast"]
    resources:
        resources=config["default_resources_8thread"]
    shell:
        """
        
        cd blast2
        bash {params.script} {input.genome} ../{input.query} {wildcards.minic}.txt
        
        """
