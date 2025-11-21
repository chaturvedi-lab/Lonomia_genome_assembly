rule get_genomes_phylogeny:
    output:
        directory("ncbi/{ncbi}_phylogeny/ncbi_dataset")
    params: 
        dir0 = lambda wildcards: "ncbi/" + wildcards.ncbi + "_phylogeny/",
        accession = lambda wildcards: NCBI_DICT[wildcards.ncbi][0]
    conda:
        config["environments"]["ncbi_datasets"]
    resources:
        resources=config["default_resources"]
    benchmark:
        "benchmarks/get_genomes_phylogeny/{ncbi}"
    shell:
        """
        
        if [ ! -d ncbi ]; then
        	mkdir -p ncbi;
        fi
        
        if [ ! -d {params.dir0} ]; then
        	mkdir -p {params.dir0};
        fi
        
        cd {params.dir0}
        
        datasets download genome accession {params.accession} --include gff3,rna,cds,protein,genome,seq-report --filename ncbi_dataset.zip
        
        if [ -e README.md ]; then
        rm README.md
        fi
        
        unzip ncbi_dataset.zip
        rm ncbi_dataset.zip
  
        """