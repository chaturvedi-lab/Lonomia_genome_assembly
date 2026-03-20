# Genome Annotation Pipeline for *Lonomia casanarensis*
**PacBio HiFi primary haplotype (p_ctg), no RNA-seq**

This document describes a complete, reproducible pipeline to annotate the **primary haplotype assembly** of *Lonomia casanarensis* using PacBio HiFi data **without RNA-seq**, relying on repeat masking and protein homology–based gene prediction.

---

## Genome assembly
- **File**: `lonomia_assembly1.bp.hap1.p_ctg.fa`
- **Path**: `/project/samc/lonomia/lonomia_assembly1.bp.hap1.p_ctg.fa`
- **Assembly type**: HiFi, primary contigs (hifiasm `p_ctg`)
- **Taxon**: Lepidoptera

---

## Software environment

### Conda setup (recommended)
```bash
module load python
conda create -n lonomia_annot \
  -c conda-forge -c bioconda \
  seqkit cd-hit busco agat \
  python=3.10
conda activate lonomia_annot
```

### Install seqkit only
```bash
conda install -c bioconda seqkit
```

---

## Directory structure
```text
/project/samc/lonomia/
├── 00_inputs/
├── 00_qc/
├── 01_repeats/
├── 02_proteins/
├── 03_models/
├── 04_function/
├── 05_eval/
└── logs/
```

---

## Step 1: Assembly QC
```bash
seqkit stats \
  /project/samc/lonomia/lonomia_assembly1.bp.hap1.p_ctg.fa \
  > 00_qc/primary.seqkit.stats.txt
```

---

## Step 2: Repeat annotation and soft-masking

### Build repeat library
```bash
#before submitting the job create a output directory
mkdir 01_repeats

#!/bin/bash
#SBATCH -J repmod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=120G
#SBATCH -t 168:00:00
#SBATCH -p single
#SBATCH -A loni_philenor05
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=schaturvedi@tulane.edu
#SBATCH -o logs/repmod_%j.out
#SBATCH -e logs/repmod_%j.err

set -euo pipefail

module load python
source activate repeats

GENOME=./lonomia_assembly1.bp.hap1.p_ctg.fa
DBDIR=./01_repeats
DBNAME=${DBDIR}/Lcas_primary_db

mkdir -p logs 01_repeats

# Build DB (only needs to be done once)
if [ ! -f "${DBNAME}.nhr" ] && [ ! -f "${DBNAME}.nsq" ]; then
  BuildDatabase -name "$DBNAME" "$GENOME"
fi

# Run RepeatModeler (RepeatModeler2 uses -pa, not -threads)
RepeatModeler \
  -database "$DBNAME" \
  -threads ${SLURM_CPUS_PER_TASK} \
  -LTRStruct \
  > 01_repeats/repeatmodeler.log 2>&1
```

The output folder that has the files is RM_1647426.MonMar21551332026


### Soft-mask genome
```bash
#!/bin/bash
#SBATCH -J repmask
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=120G
#SBATCH -t 168:00:00
#SBATCH -p single
#SBATCH -A loni_philenor05
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=schaturvedi@tulane.edu
#SBATCH -o logs/repmask_%j.out
#SBATCH -e logs/repmask_%j.err

set -euo pipefail

module purge
module load python || true
source activate repeats

# ---- inputs ----
GENOME="lonomia_assembly1.bp.hap1.p_ctg.fa"

# RepeatModeler output library (adjust path)
# Typical name: <DBNAME>-families.fa
LIB="01_repeats/Lcas_primary_db-families.fa"

# ---- output ----
#mkdir 02_repeatmasker_custom
OUTDIR="02_repeatmasker_custom"
#mkdir -p "$OUTDIR" logs

# ---- temp ----
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
mkdir -p "$TMPDIR"

# ---- threads ----
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# ---- run ----
RepeatMasker \
  -pa "$THREADS" \
  -lib "$LIB" \
  -dir "$OUTDIR" \
  -gff -xsmall -a \
  "$GENOME"
```

---

## Step 3: Protein evidence

1. Install NCBI datasets CLI
```bash
conda create -n ncbi_datasets -c conda-forge ncbi-datasets-cli unzip
conda activate ncbi_datasets
datasets --version
```
2. Download protein FASTA for a few Lepidoptera references

Pick current, non-suppressed accessions. Here are safe ones:

Bombyx mori: GCF_000151625.1

Danaus plexippus (RefSeq): GCF_009731565.1

Manduca sexta (newer GenBank assembly): GCA_014839805.1

You can add more later, but these are enough to get BRAKER moving.

```bash
mkdir -p /project/samc/lonomia/02_proteins/ncbi
cd /project/samc/lonomia/02_proteins/ncbi

# download protein FASTA + annotation package (datasets)
datasets download genome accession GCF_000151625.1 --include protein --filename Bombyx_mori.zip
datasets download genome accession GCF_009731565.1 --include protein --filename Danaus_plexippus.zip
datasets download genome accession GCA_014839805.1 --include protein --filename Manduca_sexta.zip

# unzip
unzip -o Bombyx_mori.zip -d Bombyx_mori
unzip -o Danaus_plexippus.zip -d Danaus_plexippus
unzip -o Manduca_sexta.zip -d Manduca_sexta
```
I also downloaded the proteins and transcripts from hylesia_metabus and automeris_io datasets.

```bash
#automeris
wget https://zenodo.org/records/10697861/files/Automeris_io_braker2_augustus_pep.aa

#hlyesia
wget https://bipaa.genouest.org/sp/hylesia_metabus/download/hylesia_metabus/assembly_1.0/annotation_Helixer/proteins_fa/hmetabus_ass1.0_annotHelixer_proteins.fa
```

3. Extract protein fastas and combine

```bash
cd /project/samc/lonomia

# gather all protein FASTAs from the unzipped datasets
find 02_proteins/ -name "*protein.faa" -type f > 02_proteins/protein_files.txt

# concatenate into one evidence file
cat $(cat 02_proteins/protein_files.txt) > 02_proteins/lepidoptera_proteins.raw.faa

#combine to master file
cat Automeris_io_braker2_augustus_pep.aa hmetabus_ass1.0_annotHelixer_proteins.fa lepidoptera_proteins.raw.faa > saturniid_proteins_all.raw.faa

#sanity check
grep -c "^>" 02_proteins/saturniid_proteins_all.raw.faa
head -n 2 02_proteins/saturniid_proteins_all.raw.faa
```

```bash
cd-hit \
  -i 02_proteins/saturniid_proteins_all.raw.faa \
  -o 02_proteins/saturniid_proteins_all.cdhit95.faa \
  -c 0.95 -n 5 -T 32 -M 0
```

---

## Step 4: Gene prediction (BRAKER3)
```bash
#1. Pull the container. singularity pull docker://teambraker/braker3:latest

##2. Copy the augustus config to a writable directory and set variables. 

WORKDIR=/global/scratch/users/samridhichaturvedi/lonomia/annotation/03_braker/braker_primary
AUG_CFG="$WORKDIR/augustus_config"
mkdir -p "$AUG_CFG"

#3. Run singularity like this .. 
singularity exec --bind "$WORKDIR:$WORKDIR" braker3_latest.sif cp -r /opt/Augustus/config/. "$AUG_CFG/"

#4. Run braker by pointing to the config
singularity exec --home "$WORKDIR" --bind "$WORKDIR:$WORKDIR" braker3_latest.sif braker.pl --genome=/global/scratch/users/samridhichaturvedi/lonomia/annotation/02_repeatmasker_custom/lonomia_assembly1.bp.hap1.p_ctg.fa.masked \
   --prot_seq=/global/scratch/users/samridhichaturvedi/lonomia/annotation/02_proteins/saturniid_proteins_all.cdhit95.faa \
   --softmasking \
   --workingdir="$WORKDIR" \
   --AUGUSTUS_CONFIG_PATH "$AUG_CFG" \
   --threads=32  \
   --gff3
```
#5. Summarize Braker3 gff using AGAT

```bash
conda activate braker3-env
conda install -c bioconda agat

agat_sp_statistics.pl \
  --gff 03_braker/braker_primary/braker.gff3 \
  -o 06_stats/braker_agat_stats.txt
---

## Step 5: BUSCO
```bash
conda activate lonomia_annot

busco \
-i 03_braker/braker_primary/braker.aa \
-l lepidoptera_odb10 \
-m proteins \
-o 04_busco \
-c 32

#pulling unique AAs and redoing BUSCO
cd-hit \
  -i braker.longest.aa \
  -o braker.longest.cdhit99.aa \
  -c 0.99 -n 5

#rerun BUSCO on the output file
```

---

## Step 6: Functional annotation

### InterProScan
```bash
./interproscan/interproscan-5.77-108.0/interproscan.sh \
 -i 03_braker/braker_primary/braker.clean.aa \
 -dp -iprlookup \
 -f tsv \
 -cpu 32 \
 -b 05_interproscan/lonomia_interpro
```

Summarizing interproscan output for manuscript

```bash
#1. count number of protein hits
cut -f1 05_interproscan/lonomia_interpro.tsv | sort | uniq | wc -l #32013

#compare to braker3
grep -c "^>" 03_braker/braker_primary/braker.clean.aa #37088

#percent annotated
annot=$(cut -f1 05_interproscan/lonomia_interpro.tsv | sort | uniq | wc -l)
total=$(grep -c "^>" 03_braker/braker_primary/braker.clean.aa)
echo "scale=2; 100*$annot/$total" | bc #86.31

## 2. Count unique InterPro accessions
cut -f12 05_interproscan/lonomia_interpro.tsv | grep -v "^-$" | sort | uniq | wc -l #12620

## 3. Count unique Pfam domains
awk -F'\t' '$4=="Pfam"{print $5}' 05_interproscan/lonomia_interpro.tsv | sort | uniq | wc -l #5921

## 4. Count proteins IPR terms
awk -F'\t' '$12 != "-" {print $1}' 05_interproscan/lonomia_interpro.tsv | sort | uniq | wc -l

## 5. Top 20 most frequest IPR terms
awk -F'\t' '$12 != "-" {print $12"\t"$13}' 05_interproscan/lonomia_interpro.tsv \
  | sort | uniq -c | sort -nr | head -20

## 6. Top 20 most frequent Pfam domains 
awk -F'\t' '$4=="Pfam"{print $5"\t"$6}' 05_interproscan/lonomia_interpro.tsv \
>   | sort | uniq -c | sort -nr | head -20

```


### eggNOG-mapper
```bash
conda activate braker3-env
download_eggnog_data.py   -y   -f   --data_dir 05_eggnog
```

Run eggnog:

```bash
export EGGNOG_DATA_DIR=/global/scratch/users/samridhichaturvedi/lonomia/annotation/05_eggnog

emapper.py \
  -i 03_braker/braker_primary/braker.clean.aa \
  -m diamond \
  --data_dir $EGGNOG_DATA_DIR \
  --cpu 32 \
  -o 05_eggnog/lonomia_eggnog
```

Summarizing eggnog mapper output

```bash
ann=05_eggnog/lonomia_eggnog.emapper.annotations
faa=03_braker/braker_primary/braker.clean.aa

total=$(grep -c "^>" $faa)
annot=$(awk -F'\t' '!/^#/ && $1 != "" {a[$1]=1} END{print length(a)}' $ann)
go=$(awk -F'\t' '!/^#/ && $10 != "-" && $10 != "" {a[$1]=1} END{print length(a)}' $ann)
kegg_ko=$(awk -F'\t' '!/^#/ && $12 != "-" && $12 != "" {a[$1]=1} END{print length(a)}' $ann)
kegg_path=$(awk -F'\t' '!/^#/ && $13 != "-" && $13 != "" {a[$1]=1} END{print length(a)}' $ann)
desc=$(awk -F'\t' '!/^#/ && $8 != "-" && $8 != "" {a[$1]=1} END{print length(a)}' $ann)
pfam_prot=$(awk -F'\t' '!/^#/ && $NF != "-" && $NF != "" {a[$1]=1} END{print length(a)}' $ann)
uniq_pfam=$(awk -F'\t' '!/^#/ && $NF != "-" && $NF != "" {n=split($NF,a,","); for(i=1;i<=n;i++) if(a[i] != "" && a[i] != "-") print a[i]}' $ann | sort | uniq | wc -l)

echo "Total proteins: $total" #37088
echo "Proteins with eggNOG annotations: $annot" #30815
echo "Percent annotated: $(echo "scale=2; 100*$annot/$total" | bc)%" #83.08%
echo "Proteins with functional descriptions: $desc" #26997
echo "Proteins with GO terms: $go" #18965
echo "Proteins with KEGG KO terms: $kegg_ko" #18015
echo "Proteins with KEGG pathways: $kegg_path" #11279
echo "Proteins with Pfam annotations: $pfam_prot" #26364
echo "Unique Pfam domains: $uniq_pfam" #4930
```
Get top hits
```bash
##functional hits
awk -F'\t' '
!/^#/ && $8 != "-" && $8 != "" {print $8}
' 05_eggnog/lonomia_eggnog.emapper.annotations \
| sort | uniq -c | sort -nr | head -20

#kegg pathways
awk -F'\t' '
!/^#/ && $13 != "-" && $13 != "" {
  n=split($13,a,",")
  for(i=1;i<=n;i++) print a[i]
}' 05_eggnog/lonomia_eggnog.emapper.annotations \
| sort | uniq -c | sort -nr | head -20
---

## Final merge to create a final file

```bash
python merge_annotation_tables.py \
  --gtf 03_braker/braker_primary/braker.gtf \
  --interpro 05_interproscan/lonomia_interpro.tsv \
  --eggnog 05_eggnog/lonomia_eggnog.emapper.annotations \
  --proteins 03_braker/braker_primary/braker.clean.aa \
  --out 06_final_annotation/lonomia_annotation_merged.tsv
```


## maker annotation summary for manuscript
```bash
python make_annotation_table.py \
  --gtf 03_braker/braker_primary/braker.gtf \
  --fasta 03_braker/braker_primary/braker.clean.aa \
  --busco BUSCO_braker.nr.aa/short_summary*.txt \
  --interpro 05_interproscan/lonomia_interpro.tsv \
  --eggnog 05_eggnog/lonomia_eggnog.emapper.annotations \
  --out 06_stats/table1_annotation_summary.tsv
  ```