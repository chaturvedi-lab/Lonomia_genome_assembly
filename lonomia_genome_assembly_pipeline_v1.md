# STEP 1: Adapter filtration using HiFiAdapterFilt v. 0.2.3

It is possible that the adapters have been already filtered by the UC Davis core before they handed over the data. Check this with Oahn or use the following code for filtering out adapters. 

Link: https://github.com/sheinasim-USDA/HiFiAdapterFilt


# STEP 2: Assembly using HiFiASM and purge duplicates
If you recieved the data in seperate files, combine all the files to create a master fasta file. This can be done using "cat" command from command line.

Then, create the assembly using HiFiASM. You might have to install this in your local directory following instructions here: https://github.com/chhylp123/hifiasm

Install HiFiASM using the following isntructions:

```bash
# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# Run on test data (use -f0 for small datasets)
wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz
./hifiasm -o test -t4 -f0 chr11-2M.fa.gz 2> test.log
awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa  # get primary contigs in FASTA
```

Below is the Bash script I used for assembling the genome.

```bash
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

hifiasm/hifiasm/hifiasm -o lonomia_assembly1 -t 44 fasta/lonomia_combined.fasta \
        --hom-cov 19 -s 0.25 --write-paf --write-ec /dev/null
```

I then ran Genomescope on this assembly and found that the duplicate rate was pretty high and the assembly was not a great quality. I then decided to use the purge duplicate option to do the assembly. This might be required for your plant genomes.

Install purge_dups as follows using conda:

```bash
#conda install
conda create -n purgedups_env -c conda-forge -c bioconda purge_dups
conda activate purgedups_env
```

Do the preliminary data prep to run the command using the assembly created in previous step:

```bash
#Step_1 create config file from command line
pd_config.py -l pg_dup /global/scratch/users/samridhichaturvedi/lonomia/hifiasm_data/lonomia_assembly1.bp.hap1.p_ctg.fa pb.fofn

#Step 2 run_purge_dups.py to run pipeline using SLURM
run_purge_dups.py config.json global/home/users/samridhichaturvedi/.conda/envs/purgedups_env/bin/ lonomia -p bash
```

Here is a bash script to run the code using run_purge_dups.py script:

```bash
#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

# I ran busco from its own conda environment
module load anaconda3
conda activate purgedups_env

run_purge_dups.py config.json /global/home/users/samridhichaturvedi/.conda/envs/purgedups_env/bin lonomia -p bash
```

# STEP 3: Meryl and Merqury to estimate K-mers and assess assembly completeness

I used the combination of these two programs to assess the assembly completeness and I followed the instructions from here: https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/QAQC

Here is the Bash script to run Meryl and Merqury. You will have to install them using Conda.

```bash
#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

# I ran busco from its own conda environment
module load anaconda3
conda activate merqury_env

#command to run busco
meryl count k=16.42 /global/scratch/users/samridhichaturvedi/lonomia/hifiasm_data/lonomia_assembly1.bp.hap1.p_ctg.fa \
    output hap1.meryl

meryl count k=16.42 /global/scratch/users/samridhichaturvedi/lonomia/hifiasm_data/lonomia_assembly1.bp.hap2.p_ctg.fa \
        output hap2.meryl
```

```bash
#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

# I ran busco from its own conda environment
module load anaconda3
conda activate merqury_env

#command to run busco
$MERQURY/merqury.sh hap1.meryl \
 /global/scratch/users/samridhichaturvedi/lonomia/hifiasm_data/lonomia_assembly1.bp.hap1.p_ctg.fa \
 hap1_merqury
 ```

# STEP 4: Post-assembly quality checks

## 1. BUSCO

```bash
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

# I ran busco from its own conda environment
module load anaconda3
conda activate busco_env

#command to run busco
## Haplotype 1 ###
#busco -i /global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap1.p_ctg/seqs/lonomia_assembly1.bp.hap1.p_ctg.purged.fa \
 #    -l insecta_odb10 -c 44 -m genome -o lonomia_hap1_purged -f

## Haplotype 2 ###
busco -i /global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap2.p_ctg/seqs/lonomia_assembly1.bp.hap2.p_ctg.purged.fa \
     -l insecta_odb10 -c 44 -m genome -o lonomia_hap2_purged -f
```

## 2. QUAST

https://github.com/ablab/quast

I cannot find the scripts I used for this but use the above link for running QUAST on both assemblies. 

## 3. Rerun Genometools seqstat

```bash
module load genometools
#hap 1
gt seqstat -contigs -genome 816376018 /global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap1.p_ctg/seqs/lonomia_assembly1.bp.hap1.p_ctg.purged.fa > purged_hap1_assembly.stats

#hap 2
gt seqstat -contigs -genome 816376018 /global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap2.p_ctg/seqs/lonomia_assembly1.bp.hap2.p_ctg.purged.fa > purged_hap2_assembly.stats
```


# STEP 4: Scaffolding

```bash
#!/bin/bash
#SBATCH --job-name=salsa
#SBATCH --account=fc_flyminer
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

. ~/.bashrc

# I ran busco from its own conda environment
module load anaconda3
conda activate salsa_env

#purged assembly
assembly=/global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap1.p_ctg/seqs/lonomia_assembly1.bp.hap1.p_ctg.purged.fa
#contig length file
length=/global/scratch/users/samridhichaturvedi/lonomia/purge_dups/lonomia_assembly1.bp.hap1.p_ctg/seqs/lonomia_assembly1.bp.hap1.p_ctg.purged.fa.fai
#bed file
bed=/global/scratch/users/samridhichaturvedi/lonomia/minimap/lonomia_sorted.bed

run_pipeline.py -a $assembly -l $length -b $bed -e DNASE -i 20 -p yes -o out
```

# STEP 5: Ran BlobToolKits on Salsa assemblies

```bash
module load anaconda3
conda create -n blobtools2 -c bioconda -c conda-forge python=3.6 docopt pyyaml ujson pysam tqdm nodejs seqtk

conda activate blobtools2

# Set variables (edit these!)
FASTA="/global/scratch/users/samridhichaturvedi/lonomia/salsa/out_hap1/scaffolds_FINAL.fasta"
BAM="/global/scratch/users/samridhichaturvedi/lonomia/minimap/salsa_assembly_mapping/salsa_assembly_hap1.bam"
BLASTN="/global/scratch/users/samridhichaturvedi/lonomia/blobtoolkit/blast/hap1/hap1_blast.out"
BUSCO_DIR="/global/scratch/users/samridhichaturvedi/lonomia/busco/purged/salsa/hap1_salsa/run_insecta_odb10/full_table.tsv"  # or full path to full_table.tsv
TAXDUMP="/global/scratch/users/samridhichaturvedi/lonomia/blobtoolkit/taxdump"
BLOBDIR="blobdir_hap1"

#create
blobtools create   --fasta "$FASTA"   --taxdump "$TAXDUMP"   "$BLOBDIR"

#add coverage
samtools index -c /global/scratch/users/samridhichaturvedi/lonomia/minimap/salsa_assembly_mapping/salsa_assembly_hap2.bam
blobtools add --cov "$BAM" "$BLOBDIR"

#add blastn hits
 blobtools add --hits "$BLASTN" --hits-cols 1=qseqid,2=qlen,3=sseqid,4=slen,5=length,6=pident,7=evalue,8=bitscore --taxrule bestsumorder --taxdump $TAXDUMP "$BLOBDIR"

#add busco
blobtools add --busco "$BUSCO_DIR" "$BLOBDIR"

#plots
blobtools view --plot --format svg --view snail blobdir_hap1
blobtools view --plot --format svg --view blob blobdir_hap1
blobtools view --plot --format svg --view cumulative blobdir_hap1
```