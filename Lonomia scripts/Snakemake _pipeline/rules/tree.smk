checkpoint gene_names:
    input:
        samples_compleasm = expand("qc/{sample}_compleasm", sample=SAMPLE),
        ncbi_compleasm = expand("phylogeny/compleasm/{ncbi}", ncbi=NCBI)
    output:
        txt = "gene.txt"
    params:
        threshold = 15
    resources:
        resources=config["default_resources"]
    shell:
        """
        # Process sample genomes
        for item in {input.samples_compleasm}; do
            grep ">" ${{item}}/*odb10/gene_marker.fasta | sed 's/^>//'
        done > temp_sample_genes.txt
        
        # Process NCBI genomes
        for item in {input.ncbi_compleasm}; do
            grep ">" ${{item}}/*odb10/gene_marker.fasta | sed 's/^>//'
        done > temp_ncbi_genes.txt
        
        # Combine and count occurrences
        cat temp_sample_genes.txt temp_ncbi_genes.txt | sort | uniq -c | sort -nr > tmp_{output}
        
        # Filter by threshold
        awk '$1 >= {params.threshold} {{print $2}}' tmp_{output} > {output}
        
        # Clean up temporary files
        rm temp_sample_genes.txt temp_ncbi_genes.txt
        """
# checkpoint gene_names:
#     input:
#         expand("qc/{sample}_compleasm", sample=SAMPLE)
#     output:
#         txt = "gene.txt"
#     params:
#         threshold = 20
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         for item in {input}; do
#             grep ">" ${{item}}/*odb10/gene_marker.fasta | sed 's/^>//'
#         done | sort | uniq -c | sort -nr > tmp_{output}
#         
#         awk '$1 >= {params.threshold} {{print $2}}' tmp_{output} > {output}
#         
#         """

def get_genes():
    out = checkpoints.gene_names.get().output.txt
    with open(out) as f:
        genes = [line.strip() for line in f]
    return genes

rule gene:
    input:
        sample_inputs = expand("qc/{sample}_compleasm", sample=SAMPLE),
        ncbi_inputs = expand("phylogeny/compleasm/{ncbi}", ncbi=NCBI),
        txt = "gene.txt"
    output:
        "gene/{gene}.fasta"
    resources:
        resources=config["default_resources"]
    shell:
        """
        # Extract genes from sample genomes
        for file in {input.sample_inputs}; do
            base_name=$(basename "$file")
            grep -A 1 "\\b{wildcards.gene}\\b" "$file"/*odb10/gene_marker.fasta | sed '/^--$/d' | sed "s/^>/>${{base_name}}|/" >> {output} || true
        done
            

        # Extract genes from NCBI genomes
        for file in {input.ncbi_inputs}; do
            base_name=$(basename "$file")
            grep -A 1 "\\b{wildcards.gene}\\b" "$file"/*odb10/gene_marker.fasta | sed '/^--$/d' | sed "s/^>/>${{base_name}}|/" >> {output} || true
        done
        """

# rule gene:
#     input:
#         inputs = expand("qc/{sample}_compleasm", sample=SAMPLE),
#         txt = "gene.txt"
#     output:
#         "gene/{gene}.fasta"
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
# 
#         for file in {input.inputs}; do
#             base_name=$(basename "$file")
#             grep -A 1 {wildcards.gene} "$file"/*odb10/gene_marker.fasta | sed '/^--$/d' | sed "s/^>/>${{base_name}}|/" >> {output} || true
#         done
# 
#         """


# def gather_gene(wildcards):
#     item = get_genes()
#     return [f"gene/{gene}.fasta" for gene in item]

rule align_gene:
    input:
        "gene/{gene}.fasta"
    output:
        "aligned/{gene}.faa"
    resources:
        resources=config["default_resources"]
    shell:
        "/global/scratch/users/diler/software/mafft-linux64/mafft.bat --auto {input} > {output}"



# rule edit_seq_names:
#     input:
#         "aligned/{gene}.faa"
#     output:
#         "edited/{gene}.faa"
#     shell:
#         "sed '/>/s/[:_].*//' {input} > {output}"


rule infer_gene_trees:
    input:
        "aligned/{gene}.faa"
    output:
        "aligned/{gene}.faa.treefile"
    resources:
        resources=config["default_resources"]
    shell:
        "/global/scratch/users/diler/software/iqtree-1.6.12-Linux/bin/iqtree -s {input} -nt AUTO"


def gather_trees(wildcards):
    item = get_genes()
    return [f"aligned/{gene}.faa.treefile" for gene in item]


rule species_tree:
    input:
        gather_trees
    output:
        "species_tree.tre"
    resources:
        resources=config["highmem_resources"]
    shell:
        """

        cat {input} > all_trees_tmp.tre
        
        python convert_gene_tree.py all_trees_tmp.tre all_trees.tre

        module load java
        java -jar /global/scratch/users/diler/software/Astral/astral.5.7.8.jar -i all_trees.tre -o {output} 2>{output}.log

        """
