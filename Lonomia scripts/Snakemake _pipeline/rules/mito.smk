rule find_mito_reference:
    output:
        ref_fasta="PV936473.1.fasta",
        ref_gb="PV936473.1.gb"
    resources:
        resources=config["default_resources"]
    shell:
        """
        singularity exec --bind /global/scratch/users/diler/lonomia mitohifi_master.sif \
            findMitoReference.py \
            --species "Hylesia metabus" \
            --outfolder . \
            --min_length 14000
        """

rule mitohifi:
    input:
        reads="/global/scratch/users/samridhichaturvedi/lonomia/combined_run1_run2_LINO.fasta",
        ref_fasta="PV936473.1.fasta",
        ref_gb="PV936473.1.gb"
    output:
        assembly="lino_mitogenome.fasta"
    resources:
        resources=config["default_resources_24cpus"]
    threads: 4
    shell:
        """
        singularity exec --bind /global/scratch/users/diler/lonomia mitohifi_master.sif \
            mitohifi.py \
            -r {input.reads} \
            -f {input.ref_fasta} \
            -g {input.ref_gb} \
            -t {threads} \
            --mitos
        """