#!/usr/bin/env bash

# ============================================================
# INSTALLATION SCRIPT — Viral Metagenomics Pipeline
# ============================================================
#
# This script creates all conda/mamba environments required
# to run the viral metagenomics pipeline.
#
# REQUIREMENTS:
#   - conda, mamba, or micromamba installed and in PATH
#   - Internet connection to download packages
#
# USAGE:
#   Replace CONDA_CMD below with your package manager:
#     conda   → use if you have Anaconda/Miniconda
#     mamba   → recommended for faster installs (conda-forge)
#     micromamba → lightweight alternative to mamba
#
#   Then run:
#     bash 00_install_environments.sh
#
# NOTE:
#   Each environment is isolated — tools are called using
#   their full path (e.g. ~/micromamba/envs/qc_env/bin/fastqc)
#   so environments do not need to be activated during the pipeline.
# ============================================================

CONDA_CMD="mamba"   # change to: conda or micromamba if needed

# ============================================================
# ENV 1: qc_env
# Tools: QC, trimming, alignment, assembly, and format utilities
# ============================================================
$CONDA_CMD create -n qc_env -c bioconda -c conda-forge -y \
    fastqc \
    fastp \
    multiqc \
    bowtie2 \
    samtools \
    seqkit \
    diamond \
    spades \
    pigz \
    htslib

# ============================================================
# ENV 2: kraken2_env
# Tools: taxonomic classification and interactive Krona plots
# ============================================================
$CONDA_CMD create -n kraken2_env -c bioconda -c conda-forge -y \
    kraken2 \
    krakentools \
    krona \
    biopython

# After installing krona, update its taxonomy database:
#   $(conda info --base)/envs/kraken2_env/opt/krona/updateTaxonomy.sh
# This step requires internet access and may take several minutes.

# ============================================================
# ENV 3: quast_env
# Tools: assembly quality assessment
# ============================================================
$CONDA_CMD create -n quast_env -c bioconda -c conda-forge -y \
    quast

# ============================================================
# ENV 4: bcftools_env
# Tools: variant calling and consensus genome generation
# ============================================================
$CONDA_CMD create -n bcftools_env -c bioconda -c conda-forge -y \
    bcftools

# ============================================================
# ENV 5: sratools
# Tools: downloading reads from NCBI SRA
# ============================================================
$CONDA_CMD create -n sratools -c bioconda -c conda-forge -y \
    sra-tools

# ============================================================
echo ""
echo "=== ALL ENVIRONMENTS INSTALLED ==="
echo ""
echo "  qc_env       : fastqc, fastp, multiqc, bowtie2, samtools, seqkit, diamond, spades, pigz"
echo "  kraken2_env  : kraken2, krakentools, krona, biopython"
echo "  quast_env    : quast"
echo "  bcftools_env : bcftools"
echo "  sratools     : prefetch, fasterq-dump"
echo ""
echo "IMPORTANT: Update Krona taxonomy before running the pipeline."
echo "  See the comment in ENV 2 section above."
