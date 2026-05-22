# Entrenamiento_Abbott_Prociencia

Este entrenamiento en el equipo de Abbott Core Diagnostics fue financiado por PROCIENCIA bajo el contrato PE501098965-2025 - "Fortalecimiento de capacidades bioinformáticas en análisis metagenómico para el diagnóstico molecular de patógenos en contextos clínicos y de salud pública"

---

## Viral Metagenomics Pipeline

A bioinformatics pipeline for viral discovery in metagenomic sequencing data. It applies a three-layer classification strategy to characterize known viruses, divergent viruses, and novel uncharacterized sequences (dark matter) from paired-end Illumina reads.

### Pipeline overview

| Step | Description |
|------|-------------|
| 1 | Download raw reads from NCBI SRA |
| 2 | FastQC quality assessment on raw reads |
| 3 | Adapter trimming and quality filtering (fastp) |
| 4 | FastQC on trimmed reads + MultiQC report |
| 5A | Kraken2 taxonomic classification → Krona plot |
| 5B | Host genome + PhiX depletion (Bowtie2) |
| 5C | Microbial depletion with Kraken2 (keep viruses + unclassified) |
| 6 | Three-layer viral classification (Bowtie2 + DIAMOND) |
| 6B | Comprehensive Krona plot combining all classification layers |
| 7 | Viral genome assembly (SPAdes --metaviral) |
| 8 | Assembly QC: QUAST, coverage depth, and consensus genome |

### Three-layer viral classification

- **Layer 1** — Bowtie2 alignment against viral RefSeq genomes: detects known viruses at the nucleotide level
- **Layer 2** — DIAMOND blastx against viral proteins: detects divergent viruses detectable only at the protein level
- **Layer 3** — Reads unclassified by both layers: dark matter, potential novel viruses

### Expected outputs

- `quality_control/multiqc/multiqc_report.html` — QC summary for all samples
- `krona/*.kreport` — Kraken2 classification reports
- `krona/*_krona.html` — Interactive Krona plots (pre-depletion overview)
- `krona/*_comprehensive_krona.html` — Comprehensive Krona combining all layers
- `viral_detection/abundance/*_summary.txt` — Read counts per layer per sample
- `viral_detection/abundance/*_all_viruses.txt` — All detected viruses (Layer 1)
- `viral_detection/layer2_divergent/*_diamond_summary.tsv` — Divergent virus table (Layer 2)
- `assembly/` — SPAdes assemblies per virus and dark matter
- `assembly_qc/*/coverage.txt` — Genome coverage per assembled virus
- `assembly_qc/*/consensus.fasta` — Consensus genome sequences
- `diversity_results/layer1_alpha_diversity.csv` — Alpha diversity metrics for Layer 1
- `diversity_results/layer2_alpha_diversity.csv` — Alpha diversity metrics for Layer 2
- `diversity_results/virome_composition_3layers.pdf` — Proportional bar chart of all layers
- `diversity_results/beta_bray_curtis.csv` — Bray-Curtis dissimilarity matrix (≥2 samples)
- `diversity_results/beta_pcoa_bray_curtis.pdf` — PCoA ordination plot (≥2 samples)

### Scripts

| File | Description |
|------|-------------|
| `00_install_environments.sh` | Creates all required conda/mamba environments |
| `01_download_databases.sh` | Downloads and indexes all reference databases |
| `02_virome_pipeline.sh` | Generic pipeline — configure host and paths before running |
| `03_example_pig_virome.sh` | Example run with pig serum virome (SRR38194065) |
| `04_diversity_analysis.R` | Alpha and beta diversity analysis for multi-sample viromes |

### Requirements

- conda, mamba, or micromamba
- ~106 GB disk space for databases
- Recommended: SSD for Kraken2 database (performance-critical)
- RAM: at least 32 GB; 64+ GB recommended for SPAdes assembly

