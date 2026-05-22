#!/usr/bin/env bash

# ============================================================
# DATABASE DOWNLOAD SCRIPT — Viral Metagenomics Pipeline
# ============================================================
#
# This script downloads and prepares all reference databases
# required by the pipeline:
#
#   1. Kraken2 PlusPF database  → taxonomic classification
#   2. Host genome (pig)        → host read depletion
#   3. PhiX174 genome           → sequencing control depletion
#   4. Viral RefSeq genomes     → Layer 1 (Bowtie2 alignment)
#   5. Viral protein database   → Layer 2 (DIAMOND blastx)
#
# USAGE:
#   1. Set DB_DIR to your preferred database directory
#   2. Set tool environment paths (QC_ENV, KRAKEN_ENV)
#   3. Run: bash 01_download_databases.sh
#
# NOTE ON URLs:
#   NCBI FTP and AWS S3 URLs may change over time as new genome
#   assemblies and database versions are released. If a download
#   fails, check the following sources for updated URLs:
#     Kraken2 databases : https://benlangmead.github.io/aws-indexes/k2
#     NCBI genomes      : https://www.ncbi.nlm.nih.gov/datasets/genome/
#     NCBI RefSeq viral : https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
#
# DISK SPACE REQUIREMENTS (approximate):
#   Kraken2 PlusPF  : ~90 GB (extracted)
#   Pig genome      : ~3 GB
#   Human genome    : ~3 GB
#   Mosquito genome : ~2 GB
#   PhiX174 genome  : <1 MB
#   Viral RefSeq    : ~5 GB
#   DIAMOND protein : ~3 GB
#   TOTAL           : ~106 GB
#
# NOTE:
#   The Kraken2 database is large. If your system RAM is smaller
#   than the database size, use --memory-mapping when running
#   kraken2. Placing the database on an SSD is strongly recommended
#   for reasonable classification speeds.
# ============================================================

# === CONFIGURE THESE BEFORE RUNNING ===
DB_DIR="/path/to/databases"           # destination for all databases
QC_ENV="/path/to/qc_env/bin"         # bowtie2, samtools, diamond
KRAKEN_ENV="/path/to/kraken2_env/bin" # kraken2
THREADS=32

# ============================================================
# 1. Kraken2 PlusPF database
# ============================================================
# Contains: bacteria, archaea, viruses, plasmids, fungi,
#           protozoa, and human genome (RefSeq)
#
# Available pre-built databases:
#   Standard (~8 GB):  https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
#   PlusPF   (~69 GB): https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz
#   PlusPFP  (~87 GB): https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240904.tar.gz
#
# More options: https://benlangmead.github.io/aws-indexes/k2

KRAKEN_DB_DIR="$DB_DIR/kraken2_db"
mkdir -p "$KRAKEN_DB_DIR"

echo "=== Downloading Kraken2 PlusPF database ==="
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz \
    -O "$KRAKEN_DB_DIR/k2_pluspf_20240904.tar.gz"

echo "=== Extracting Kraken2 database ==="
tar -xzf "$KRAKEN_DB_DIR/k2_pluspf_20240904.tar.gz" -C "$KRAKEN_DB_DIR"

# Remove compressed file after extraction to save disk space
rm "$KRAKEN_DB_DIR/k2_pluspf_20240904.tar.gz"

# ============================================================
# 2. Host genome — Sus scrofa (pig)
# ============================================================
# Source: NCBI RefSeq (Sscrofa11.1)
# Used for: Bowtie2 host depletion (step 5B)

PIG_DIR="$DB_DIR/pig_genome"
mkdir -p "$PIG_DIR"

echo "=== Downloading pig genome (Sscrofa11.1) ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz" \
    -O "$PIG_DIR/sscrofa11.1.fna.gz"

echo "=== Building Bowtie2 index for pig genome ==="
"$QC_ENV/bowtie2-build" \
    "$PIG_DIR/sscrofa11.1.fna.gz" \
    "$PIG_DIR/sscrofa_index" \
    --threads "$THREADS"

# ============================================================
# 3. PhiX174 control genome
# ============================================================
# Source: NCBI RefSeq (NC_001422.1)
# Used for: removing Illumina sequencing control reads

PHIX_DIR="$DB_DIR/phix_genome"
mkdir -p "$PHIX_DIR"

echo "=== Downloading PhiX174 genome ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
    -O "$PHIX_DIR/phix174.fna.gz"

echo "=== Building Bowtie2 index for PhiX174 ==="
"$QC_ENV/bowtie2-build" \
    "$PHIX_DIR/phix174.fna.gz" \
    "$PHIX_DIR/phix_index" \
    --threads "$THREADS"

# ============================================================
# 4. Host genome — Homo sapiens (human) [optional]
# ============================================================
# Source: NCBI RefSeq (GRCh38)
# Used for: human read depletion (uncomment in pipeline if needed)

HUMAN_DIR="$DB_DIR/human_genome"
mkdir -p "$HUMAN_DIR"

echo "=== Downloading human genome (GRCh38) ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
    -O "$HUMAN_DIR/hg38.fna.gz"

echo "=== Building Bowtie2 index for human genome ==="
"$QC_ENV/bowtie2-build" \
    "$HUMAN_DIR/hg38.fna.gz" \
    "$HUMAN_DIR/hg38_index" \
    --threads "$THREADS"

# ============================================================
# 5. Host genome — Aedes aegypti (mosquito)
# ============================================================
# Source: NCBI RefSeq (AaegL5.0)
# Used for: mosquito host read depletion in the mosquito virome pipeline

MOSQUITO_DIR="$DB_DIR/mosquito_genome"
mkdir -p "$MOSQUITO_DIR"

echo "=== Downloading Aedes aegypti genome (AaegL5.0) ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz" \
    -O "$MOSQUITO_DIR/aedes_aegypti.fna.gz"

echo "=== Building Bowtie2 index for Aedes aegypti ==="
"$QC_ENV/bowtie2-build" \
    "$MOSQUITO_DIR/aedes_aegypti.fna.gz" \
    "$MOSQUITO_DIR/aedes_index" \
    --threads "$THREADS"

# ============================================================
# 6. Viral RefSeq genomes — Layer 1 Bowtie2 database
# ============================================================
# Source: NCBI RefSeq viral genomes
# Used for: nucleotide-level viral classification (Bowtie2)

VIRAL_DIR="$DB_DIR/viral_db"
mkdir -p "$VIRAL_DIR"

echo "=== Downloading viral RefSeq genomes ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz" \
    -O "$VIRAL_DIR/ref_vir_all.fna.gz"

# Download additional viral genome files if available
# Check https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ for all parts

echo "=== Building Bowtie2 index for viral RefSeq ==="
"$QC_ENV/bowtie2-build" \
    "$VIRAL_DIR/ref_vir_all.fna.gz" \
    "$VIRAL_DIR/ref_vir_index" \
    --threads "$THREADS"

# ============================================================
# 5. Viral protein database — Layer 2 DIAMOND database
# ============================================================
# Source: NCBI RefSeq viral proteins
# Used for: protein-level viral classification (DIAMOND blastx)

echo "=== Downloading viral RefSeq proteins ==="
wget -c "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz" \
    -O "$VIRAL_DIR/viral_proteins.faa.gz"

echo "=== Building DIAMOND database ==="
"$QC_ENV/diamond" makedb \
    --in "$VIRAL_DIR/viral_proteins.faa.gz" \
    --db "$VIRAL_DIR/viral_proteins.dmnd" \
    --threads "$THREADS"

# ============================================================
echo ""
echo "=== ALL DATABASES READY ==="
echo ""
echo "  Kraken2 DB      : $KRAKEN_DB_DIR"
echo "  Pig index       : $PIG_DIR/sscrofa_index"
echo "  Human index     : $HUMAN_DIR/hg38_index"
echo "  Mosquito index  : $MOSQUITO_DIR/aedes_index"
echo "  PhiX index      : $PHIX_DIR/phix_index"
echo "  Viral index     : $VIRAL_DIR/ref_vir_index"
echo "  DIAMOND DB      : $VIRAL_DIR/viral_proteins.dmnd"
