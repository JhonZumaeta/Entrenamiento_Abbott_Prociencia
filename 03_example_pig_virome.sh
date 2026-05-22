#!/usr/bin/env bash

# ============================================================
# VIRAL METAGENOMICS PIPELINE — Pig serum virome
# ============================================================
#
# Steps:
#   1.  Download reads from SRA
#   2.  FastQC on raw reads
#   3.  Adapter trimming with fastp
#   4.  FastQC on trimmed reads + MultiQC report
#   5A. Kraken2 classification (fastp reads) → Krona plot
#   5B. Host depletion: Pig + PhiX (Bowtie2)
#   5C. Kraken2 depletion (remove bacteria/archaea/fungi/plants)
#   6.  Three-layer viral classification (Bowtie2 + DIAMOND)
#   7.  Genome assembly (SPAdes --metaviral)
#   8.  Assembly QC (QUAST + coverage + consensus)
#
# Tool environments:
#   qc_env       : fastqc, fastp, multiqc, bowtie2, samtools,
#                  seqkit, diamond, spades.py, pigz
#   kraken2_env  : kraken2, krakentools, krona
#   quast_env    : quast
#   bcftools_env : bcftools
#   sratools     : prefetch, fasterq-dump
# ============================================================

# === BASE PATHS ===
BASE_DIR="/home/cayetano/Disco_4TB/training/Internship_Abbott/pig_pipeline_result"
ACCESSION_FILE="$BASE_DIR/pig_accesions.txt"

# === DATA DIRECTORIES ===
FASTQ_DIR="$BASE_DIR/fastq"
QC_DIR="$BASE_DIR/quality_control"
HOST_DIR="$BASE_DIR/host_depletion"
VIRAL_DIR="$BASE_DIR/viral_detection"
ASSEMBLY_DIR="$BASE_DIR/assembly"
QC_ASSEMBLY_DIR="$BASE_DIR/assembly_qc"
KRONA_DIR="$BASE_DIR/krona"

# === DATABASE DIRECTORIES ===
DB_DIR="/home/cayetano/Disco_4TB/training/Internship_Abbott/databases"
PIG_GENOME_DIR="$DB_DIR/pig_genome"
HUMAN_GENOME_DIR="$DB_DIR/human_genome"
PHIX_GENOME_DIR="$DB_DIR/phix_genome"
VIRAL_DB_DIR="$DB_DIR/viral_db"
KRAKEN_DB="/home/cayetano/Downloads/JHON/Internship_Abbott/kraken2_db"

# === REFERENCE GENOME INDEXES ===
PIG_INDEX="$PIG_GENOME_DIR/sscrofa_index"
HUMAN_INDEX="$HUMAN_GENOME_DIR/hg38_index"
PHIX_INDEX="$PHIX_GENOME_DIR/phix_index"

# === VIRAL DATABASES ===
VIRAL_INDEX="$VIRAL_DB_DIR/ref_vir_index"
DIAMOND_DB="$VIRAL_DB_DIR/viral_proteins.dmnd"

# === TOOL PATHS ===
QC_ENV="/home/cayetano/micromamba/envs/qc_env/bin"
SRA_ENV="/home/cayetano/micromamba/envs/sratools/bin"
KRAKEN_ENV="/home/cayetano/micromamba/envs/kraken2_env/bin"
QUAST="/home/cayetano/micromamba/envs/quast_env/bin/quast"
BCFTOOLS="/home/cayetano/micromamba/envs/bcftools_env/bin/bcftools"

# === PARAMETERS ===
THREADS=32
THREADS_HIGH=36
MIN_READS=500

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================
mkdir -p \
    "$FASTQ_DIR" \
    "$QC_DIR/fastqc_raw" \
    "$QC_DIR/fastp" \
    "$QC_DIR/fastqc_trimmed" \
    "$QC_DIR/multiqc" \
    "$HOST_DIR" \
    "$VIRAL_DIR/layer1_known" \
    "$VIRAL_DIR/layer2_divergent" \
    "$VIRAL_DIR/layer3_unclassified" \
    "$VIRAL_DIR/abundance" \
    "$ASSEMBLY_DIR" \
    "$QC_ASSEMBLY_DIR" \
    "$KRONA_DIR"

# ============================================================
# STEP 1: Download reads from SRA
# ============================================================
echo "=== STEP 1: Downloading reads from SRA ==="

SRA_TMP="$BASE_DIR/sra_tmp"
mkdir -p "$SRA_TMP"

"$SRA_ENV/prefetch" \
    --option-file "$ACCESSION_FILE" \
    --max-size 200G \
    --progress \
    --verify yes \
    -O "$SRA_TMP"

for SRA_FILE in "$SRA_TMP"/*/*.sra; do
    SAMPLE=$(basename "$(dirname "$SRA_FILE")")
    echo "  Converting: $SAMPLE"
    "$SRA_ENV/fasterq-dump" "$SRA_FILE" \
        --split-files \
        --outdir "$FASTQ_DIR" \
        --threads "$THREADS" \
        --progress
    rm -rf "$SRA_TMP/$SAMPLE"
done

rmdir "$SRA_TMP" 2>/dev/null || true

echo "  Compressing FASTQs..."
"$QC_ENV/pigz" -p "$THREADS" "$FASTQ_DIR"/*.fastq
echo "  Done: raw FASTQs in $FASTQ_DIR"

# ============================================================
# STEP 2: FastQC on raw reads
# ============================================================
echo "=== STEP 2: FastQC (raw reads) ==="

"$QC_ENV/fastqc" "$FASTQ_DIR"/*.fastq.gz \
    --outdir "$QC_DIR/fastqc_raw" \
    --threads "$THREADS"

# ============================================================
# STEP 3: Adapter trimming with fastp
# ============================================================
echo "=== STEP 3: Adapter trimming (fastp) ==="

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"

    echo "  Processing: $SAMPLE"
    "$QC_ENV/fastp" \
        -i "$R1" -I "$R2" \
        -o "$QC_DIR/fastp/${SAMPLE}_1_clean.fastq.gz" \
        -O "$QC_DIR/fastp/${SAMPLE}_2_clean.fastq.gz" \
        --html "$QC_DIR/fastp/${SAMPLE}_fastp.html" \
        --json "$QC_DIR/fastp/${SAMPLE}_fastp.json" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --low_complexity_filter \
        --dedup
done

# ============================================================
# STEP 4: FastQC on trimmed reads + MultiQC report
# ============================================================
echo "=== STEP 4: FastQC (trimmed reads) + MultiQC ==="

"$QC_ENV/fastqc" "$QC_DIR/fastp"/*_clean.fastq.gz \
    --outdir "$QC_DIR/fastqc_trimmed" \
    --threads "$THREADS"

"$QC_ENV/multiqc" \
    "$QC_DIR/fastqc_raw" \
    "$QC_DIR/fastp" \
    "$QC_DIR/fastqc_trimmed" \
    --outdir "$QC_DIR/multiqc" \
    --filename "multiqc_report" \
    --title "Virome Pipeline QC" \
    --force

echo "  Report: $QC_DIR/multiqc/multiqc_report.html"

# ============================================================
# STEP 5A: Kraken2 classification (fastp reads) → Krona plot
# ============================================================
echo "=== STEP 5A: Kraken2 classification + Krona ==="

for R1 in "$QC_DIR/fastp"/*_1_clean.fastq.gz; do
    SAMPLE=$(basename "$R1" _1_clean.fastq.gz)
    R2="$QC_DIR/fastp/${SAMPLE}_2_clean.fastq.gz"

    echo "  Classifying: $SAMPLE"
    "$KRAKEN_ENV/kraken2" \
        --db "$KRAKEN_DB" \
        --paired "$R1" "$R2" \
        --output "$KRONA_DIR/${SAMPLE}.kraken2" \
        --report "$KRONA_DIR/${SAMPLE}.kreport" \
        --threads "$THREADS" \
        --gzip-compressed \
        --memory-mapping

    awk '$4 !~ /^S/ || $2 >= 50' \
        "$KRONA_DIR/${SAMPLE}.kreport" \
        > "$KRONA_DIR/${SAMPLE}.kreport_filtered"

    "$KRAKEN_ENV/ktImportTaxonomy" \
        -t 5 -m 3 \
        "$KRONA_DIR/${SAMPLE}.kreport_filtered" \
        -o "$KRONA_DIR/${SAMPLE}_krona.html"

    echo "  Krona: $KRONA_DIR/${SAMPLE}_krona.html"
done

# ============================================================
# STEP 5B: Host depletion (Pig + PhiX)
# ============================================================
echo "=== STEP 5B: Host depletion ==="

_deplete_host() {
    local IN_R1="$1" IN_R2="$2" INDEX="$3"
    local OUT_BAM="$4" OUT_R1="$5" OUT_R2="$6" LOG="$7"

    "$QC_ENV/bowtie2" \
        -x "$INDEX" -1 "$IN_R1" -2 "$IN_R2" \
        --very-sensitive --threads "$THREADS" \
        2> "$LOG" \
    | "$QC_ENV/samtools" view -bS -@ "$THREADS_HIGH" \
    | "$QC_ENV/samtools" sort -o "$OUT_BAM" -@ "$THREADS_HIGH"

    "$QC_ENV/samtools" index "$OUT_BAM"

    "$QC_ENV/samtools" view -b -f 12 -F 256 "$OUT_BAM" \
    | "$QC_ENV/samtools" fastq \
        -1 "$OUT_R1" -2 "$OUT_R2" \
        -0 /dev/null -s /dev/null \
        -@ "$THREADS_HIGH"
}

for R1 in "$QC_DIR/fastp"/*_1_clean.fastq.gz; do
    SAMPLE=$(basename "$R1" _1_clean.fastq.gz)
    R2="$QC_DIR/fastp/${SAMPLE}_2_clean.fastq.gz"
    TMP="$HOST_DIR/${SAMPLE}_tmp"
    mkdir -p "$TMP"

    echo "  Processing: $SAMPLE"

    # Remove pig reads
    _deplete_host "$R1" "$R2" \
        "$PIG_INDEX" \
        "$TMP/${SAMPLE}_pig.bam" \
        "$TMP/${SAMPLE}_1_nopig.fastq.gz" \
        "$TMP/${SAMPLE}_2_nopig.fastq.gz" \
        "$TMP/${SAMPLE}_pig.log"

    # Remove human reads (uncomment if needed)
    # _deplete_host \
    #     "$TMP/${SAMPLE}_1_nopig.fastq.gz" \
    #     "$TMP/${SAMPLE}_2_nopig.fastq.gz" \
    #     "$HUMAN_INDEX" \
    #     "$TMP/${SAMPLE}_human.bam" \
    #     "$TMP/${SAMPLE}_1_nohuman.fastq.gz" \
    #     "$TMP/${SAMPLE}_2_nohuman.fastq.gz" \
    #     "$TMP/${SAMPLE}_human.log"

    # Remove PhiX reads
    _deplete_host \
        "$TMP/${SAMPLE}_1_nopig.fastq.gz" \
        "$TMP/${SAMPLE}_2_nopig.fastq.gz" \
        "$PHIX_INDEX" \
        "$TMP/${SAMPLE}_phix.bam" \
        "$TMP/${SAMPLE}_1_nophix.fastq.gz" \
        "$TMP/${SAMPLE}_2_nophix.fastq.gz" \
        "$TMP/${SAMPLE}_phix.log"

    # ============================================================
    # STEP 5C: Kraken2 depletion (remove non-viral reads)
    # ============================================================
    echo "  Kraken2 depletion: $SAMPLE"

    "$KRAKEN_ENV/kraken2" \
        --db "$KRAKEN_DB" \
        --paired \
        "$TMP/${SAMPLE}_1_nophix.fastq.gz" \
        "$TMP/${SAMPLE}_2_nophix.fastq.gz" \
        --output "$HOST_DIR/${SAMPLE}_kraken2.out" \
        --report "$TMP/${SAMPLE}_kraken2.report" \
        --threads "$THREADS" \
        --gzip-compressed \
        --memory-mapping

    # Extract virus (taxid 10239) + unclassified (taxid 0) reads
    "$KRAKEN_ENV/python" "$KRAKEN_ENV/extract_kraken_reads.py" \
        -k "$HOST_DIR/${SAMPLE}_kraken2.out" \
        --report "$TMP/${SAMPLE}_kraken2.report" \
        -s1 "$TMP/${SAMPLE}_1_nophix.fastq.gz" \
        -s2 "$TMP/${SAMPLE}_2_nophix.fastq.gz" \
        -o "$HOST_DIR/${SAMPLE}_1_microbial.fastq.gz" \
        -o2 "$HOST_DIR/${SAMPLE}_2_microbial.fastq.gz" \
        -t 10239 0 \
        --include-children \
        --fastq-output

    # Save data for comprehensive Krona
    "$QC_ENV/samtools" view -c -F 4 "$TMP/${SAMPLE}_pig.bam" > "$HOST_DIR/${SAMPLE}_pig_count.txt"
    cp "$TMP/${SAMPLE}_kraken2.report" "$HOST_DIR/${SAMPLE}_kraken2.report"

    echo "  Done: $SAMPLE → microbial reads ready"
done

# ============================================================
# STEP 6: Three-layer viral classification
# ============================================================
echo "=== STEP 6: Viral classification (3 layers) ==="

for R1 in "$HOST_DIR"/*_1_microbial.fastq.gz; do
    SAMPLE=$(basename "$R1" _1_microbial.fastq.gz)
    R2="$HOST_DIR/${SAMPLE}_2_microbial.fastq.gz"

    echo "  Processing: $SAMPLE"

    # --- Layer 1: Bowtie2 vs viral reference genomes ---
    if [[ ! -f "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam" ]]; then
        "$QC_ENV/bowtie2" \
            -x "$VIRAL_INDEX" -1 "$R1" -2 "$R2" \
            --very-sensitive --threads "$THREADS" \
            2> "$VIRAL_DIR/layer1_known/${SAMPLE}_bowtie2.log" \
        | "$QC_ENV/samtools" view -bS -@ "$THREADS_HIGH" \
        | "$QC_ENV/samtools" sort \
            -o "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam" \
            -@ "$THREADS_HIGH"

        "$QC_ENV/samtools" index "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam"

        # Extract reads that aligned (known viruses)
        "$QC_ENV/samtools" view -b -F 12 \
            "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam" \
        | "$QC_ENV/samtools" fastq \
            -1 "$VIRAL_DIR/layer1_known/${SAMPLE}_1_known.fastq.gz" \
            -2 "$VIRAL_DIR/layer1_known/${SAMPLE}_2_known.fastq.gz" \
            -0 /dev/null -s /dev/null -@ "$THREADS_HIGH"

        # Extract reads that did not align (input for Layer 2)
        "$QC_ENV/samtools" view -b -f 12 -F 256 \
            "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam" \
        | "$QC_ENV/samtools" fastq \
            -1 "$VIRAL_DIR/layer2_divergent/${SAMPLE}_1_unaligned.fastq.gz" \
            -2 "$VIRAL_DIR/layer2_divergent/${SAMPLE}_2_unaligned.fastq.gz" \
            -0 /dev/null -s /dev/null -@ "$THREADS_HIGH"
    fi

    # Abundance table for Layer 1
    "$QC_ENV/samtools" idxstats \
        "$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam" \
        > "$VIRAL_DIR/abundance/${SAMPLE}_abundance.txt"

    awk '$3 > 0 && $1 != "*" {print $1, $2, $3}' \
        "$VIRAL_DIR/abundance/${SAMPLE}_abundance.txt" | sort -k3 -rn \
        > "$VIRAL_DIR/abundance/${SAMPLE}_all_viruses.txt"

    awk -v min="$MIN_READS" '$3 >= min && $1 != "*" {print $1, $3}' \
        "$VIRAL_DIR/abundance/${SAMPLE}_abundance.txt" | sort -k2 -rn \
        > "$VIRAL_DIR/abundance/${SAMPLE}_top_viruses.txt"

    echo "  Viruses detected: $(wc -l < "$VIRAL_DIR/abundance/${SAMPLE}_all_viruses.txt")"
    echo "  Viruses for assembly (>=$MIN_READS reads): $(wc -l < "$VIRAL_DIR/abundance/${SAMPLE}_top_viruses.txt")"

    # --- Layer 2: DIAMOND blastx vs viral proteins ---
    zcat "$VIRAL_DIR/layer2_divergent/${SAMPLE}_1_unaligned.fastq.gz" \
         "$VIRAL_DIR/layer2_divergent/${SAMPLE}_2_unaligned.fastq.gz" \
    | "$QC_ENV/diamond" blastx \
        --db "$DIAMOND_DB" \
        --query - \
        --out "$VIRAL_DIR/layer2_divergent/${SAMPLE}_diamond.tsv" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
        --threads "$THREADS" \
        --sensitive \
        --evalue 1e-5 \
        --max-target-seqs 1

    # DIAMOND summary table
    awk -F'\t' 'BEGIN{OFS="\t"; print "virus_name","pident","evalue","read_count"}
    {
        v=$7; if(!(v in n)){n[v]=0; p[v]=0; e[v]=$5}
        n[v]++; p[v]+=$3; if($5<e[v]) e[v]=$5
    }
    END{for(v in n) printf "%s\t%.2f\t%s\t%d\n", v, p[v]/n[v], e[v], n[v]}' \
        "$VIRAL_DIR/layer2_divergent/${SAMPLE}_diamond.tsv" \
    | sort -t$'\t' -k4 -rn \
    > "$VIRAL_DIR/layer2_divergent/${SAMPLE}_diamond_summary.tsv"

    # --- Layer 3: Extract reads not classified by DIAMOND ---
    cut -f1 "$VIRAL_DIR/layer2_divergent/${SAMPLE}_diamond.tsv" | sort -u \
        > "$VIRAL_DIR/layer2_divergent/${SAMPLE}_classified_ids.txt"

    TOTAL_L2=$(zcat "$VIRAL_DIR/layer2_divergent/${SAMPLE}_1_unaligned.fastq.gz" \
        | awk 'NR%4==1' | wc -l)
    CLASSIFIED=$(wc -l < "$VIRAL_DIR/layer2_divergent/${SAMPLE}_classified_ids.txt")
    UNCLASSIFIED=$((TOTAL_L2 - CLASSIFIED))
    LAYER1=$(awk '{sum+=$3} END {print sum}' "$VIRAL_DIR/abundance/${SAMPLE}_all_viruses.txt")

    printf "sample\tlayer1_known\tlayer2_divergent\tlayer3_unclassified\n" \
        > "$VIRAL_DIR/abundance/${SAMPLE}_summary.txt"
    printf "%s\t%s\t%s\t%s\n" "$SAMPLE" "$LAYER1" "$CLASSIFIED" "$UNCLASSIFIED" \
        >> "$VIRAL_DIR/abundance/${SAMPLE}_summary.txt"

    echo "  Layer 1 reads : $LAYER1"
    echo "  Layer 2 reads : $CLASSIFIED"
    echo "  Layer 3 reads : $UNCLASSIFIED"

    # Extract unclassified read IDs
    zcat "$VIRAL_DIR/layer2_divergent/${SAMPLE}_1_unaligned.fastq.gz" \
        | awk 'NR%4==1 {print substr($1,2)}' | sort \
        > "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_all_ids.txt"

    comm -23 \
        "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_all_ids.txt" \
        "$VIRAL_DIR/layer2_divergent/${SAMPLE}_classified_ids.txt" \
        > "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_unclassified_ids.txt"

    "$QC_ENV/seqkit" grep \
        -f "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_unclassified_ids.txt" \
        "$VIRAL_DIR/layer2_divergent/${SAMPLE}_1_unaligned.fastq.gz" \
        -o "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_1_unclassified.fastq.gz"

    "$QC_ENV/seqkit" grep \
        -f "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_unclassified_ids.txt" \
        "$VIRAL_DIR/layer2_divergent/${SAMPLE}_2_unaligned.fastq.gz" \
        -o "$VIRAL_DIR/layer3_unclassified/${SAMPLE}_2_unclassified.fastq.gz"

    echo "  Done: $SAMPLE"
done

# ============================================================
# STEP 6B: Comprehensive Krona plot
# ============================================================
echo "=== STEP 6B: Comprehensive Krona ==="

for R1 in "$HOST_DIR"/*_1_microbial.fastq.gz; do
    SAMPLE=$(basename "$R1" _1_microbial.fastq.gz)
    KREPORT="$HOST_DIR/${SAMPLE}_kraken2.report"
    VIRAL_BAM="$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam"
    LAYER3_R1="$VIRAL_DIR/layer3_unclassified/${SAMPLE}_1_unclassified.fastq.gz"

    echo "  Building Krona: $SAMPLE"

    PIG_READS=$(cat "$HOST_DIR/${SAMPLE}_pig_count.txt")
    BACTERIA=$(awk '$5 == 2    {print $2; exit}' "$KREPORT")
    ARCHAEA=$(awk  '$5 == 2157 {print $2; exit}' "$KREPORT")
    FUNGI=$(awk    '$5 == 4751 {print $2; exit}' "$KREPORT")
    VIRUS_READS=$("$QC_ENV/samtools" view -c -F 12 "$VIRAL_BAM")
    DARK_MATTER=$(zcat "$LAYER3_R1" | awk 'NR%4==1' | wc -l)

    KRONA_TXT="$KRONA_DIR/${SAMPLE}_comprehensive.txt"
    {
        [[ "$PIG_READS"   -ge 50 ]] && echo -e "$PIG_READS\tEukaryota\tAnimalia\tSus scrofa"
        [[ "$BACTERIA"    -ge 50 ]] && echo -e "$BACTERIA\tBacteria"
        [[ "$ARCHAEA"     -ge 50 ]] && echo -e "$ARCHAEA\tArchaea"
        [[ "$FUNGI"       -ge 50 ]] && echo -e "$FUNGI\tEukaryota\tFungi"
        [[ "$VIRUS_READS" -ge 50 ]] && echo -e "$VIRUS_READS\tViruses"
        [[ "$DARK_MATTER" -ge 50 ]] && echo -e "$DARK_MATTER\tUnclassified"
    } > "$KRONA_TXT"

    "$KRAKEN_ENV/ktImportText" \
        "$KRONA_TXT" \
        -o "$KRONA_DIR/${SAMPLE}_comprehensive_krona.html"

    echo "  Krona: $KRONA_DIR/${SAMPLE}_comprehensive_krona.html"
done

# ============================================================
# STEP 7: Genome assembly (known viruses + Layer 3)
# ============================================================
echo "=== STEP 7: Genome assembly ==="

# Known viruses (Layer 1)
for TOP_VIRUSES in "$VIRAL_DIR/abundance"/*_top_viruses.txt; do
    SAMPLE=$(basename "$TOP_VIRUSES" _top_viruses.txt)
    BAM="$VIRAL_DIR/layer1_known/${SAMPLE}_viral.bam"

    if [[ ! -s "$TOP_VIRUSES" ]]; then
        echo "  No viruses to assemble: $SAMPLE"
        continue
    fi

    echo "  Sample: $SAMPLE"

    while read -r VIRUS READS; do
        VIRUS_SAFE=$(echo "$VIRUS" | tr '/' '_' | tr '.' '_')
        OUT_DIR="$ASSEMBLY_DIR/${SAMPLE}_${VIRUS_SAFE}"

        if [[ -d "$OUT_DIR/spades" ]]; then
            echo "  Assembly exists, skipping: $VIRUS"
            continue
        fi

        mkdir -p "$OUT_DIR/reads"

        echo "  Extracting reads: $VIRUS ($READS reads)"
        "$QC_ENV/samtools" view -b "$BAM" "$VIRUS" \
        | "$QC_ENV/samtools" fastq \
            -1 "$OUT_DIR/reads/R1.fastq.gz" \
            -2 "$OUT_DIR/reads/R2.fastq.gz" \
            -0 /dev/null -s /dev/null \
            -@ "$THREADS_HIGH"

        echo "  Assembling: $VIRUS"
        "$QC_ENV/spades.py" \
            --metaviral \
            -1 "$OUT_DIR/reads/R1.fastq.gz" \
            -2 "$OUT_DIR/reads/R2.fastq.gz" \
            -o "$OUT_DIR/spades" \
            --threads "$THREADS" \
            --memory 64

    done < "$TOP_VIRUSES"
done

# Layer 3 unclassified reads
for R1 in "$VIRAL_DIR/layer3_unclassified"/*_1_unclassified.fastq.gz; do
    SAMPLE=$(basename "$R1" _1_unclassified.fastq.gz)
    R2="$VIRAL_DIR/layer3_unclassified/${SAMPLE}_2_unclassified.fastq.gz"
    OUT_DIR="$ASSEMBLY_DIR/${SAMPLE}_unclassified"

    if [[ -d "$OUT_DIR/spades" ]]; then
        echo "  Assembly exists, skipping: ${SAMPLE}_unclassified"
        continue
    fi

    echo "  Assembling unclassified: $SAMPLE"
    mkdir -p "$OUT_DIR"

    "$QC_ENV/spades.py" \
        --metaviral \
        -1 "$R1" -2 "$R2" \
        -o "$OUT_DIR/spades" \
        --threads "$THREADS" \
        --memory 64
done

# ============================================================
# STEP 8: Assembly QC (QUAST + coverage + consensus genome)
# ============================================================
echo "=== STEP 8: Assembly QC ==="

if [[ ! -f "$VIRAL_DB_DIR/ref_vir_all.fna.bgz" ]]; then
    echo "  Converting viral DB to bgzip for random access..."
    zcat "$VIRAL_DB_DIR/ref_vir_all.fna.gz" \
        | "$QC_ENV/bgzip" -c > "$VIRAL_DB_DIR/ref_vir_all.fna.bgz"
    "$QC_ENV/samtools" faidx "$VIRAL_DB_DIR/ref_vir_all.fna.bgz"
fi

for SPADES_DIR in "$ASSEMBLY_DIR"/*/spades; do
    VIRUS_SAMPLE=$(basename "$(dirname "$SPADES_DIR")")
    SCAFFOLDS="$SPADES_DIR/scaffolds.fasta"
    OUT_DIR="$QC_ASSEMBLY_DIR/$VIRUS_SAMPLE"
    READS_R1="$(dirname "$SPADES_DIR")/reads/R1.fastq.gz"
    READS_R2="$(dirname "$SPADES_DIR")/reads/R2.fastq.gz"

    VIRUS_ACC_RAW=$(echo "$VIRUS_SAMPLE" | grep -oP 'NC_\d+_\d+' | head -1)

    # Skip unclassified assemblies — no reference genome available
    if [[ -z "$VIRUS_ACC_RAW" ]]; then
        echo "  Skipping (no reference): $VIRUS_SAMPLE"
        continue
    fi

    VIRUS_ACC="${VIRUS_ACC_RAW%_*}.${VIRUS_ACC_RAW##*_}"
    SAMPLE_ID=$(echo "$VIRUS_SAMPLE" | grep -oP 'SRR\d+|ERR\d+')

    mkdir -p "$OUT_DIR"
    echo "  Processing: $VIRUS_SAMPLE"

    REF_FASTA="$OUT_DIR/reference.fasta"
    "$QC_ENV/samtools" faidx "$VIRAL_DB_DIR/ref_vir_all.fna.bgz" \
        "$VIRUS_ACC" > "$REF_FASTA"

    "$QUAST" "$SCAFFOLDS" \
        -r "$REF_FASTA" \
        -o "$OUT_DIR/quast" \
        --threads "$THREADS"

    "$QC_ENV/bowtie2-build" "$REF_FASTA" "$OUT_DIR/ref_index" \
        --threads "$THREADS"

    "$QC_ENV/bowtie2" \
        -x "$OUT_DIR/ref_index" \
        -1 "$READS_R1" -2 "$READS_R2" \
        --very-sensitive --threads "$THREADS" \
        2> "$OUT_DIR/alignment.log" \
    | "$QC_ENV/samtools" sort -o "$OUT_DIR/aligned.bam" -@ "$THREADS_HIGH"

    "$QC_ENV/samtools" index "$OUT_DIR/aligned.bam"

    "$QC_ENV/samtools" depth -a "$OUT_DIR/aligned.bam" > "$OUT_DIR/depth.txt"
    "$QC_ENV/samtools" coverage "$OUT_DIR/aligned.bam" > "$OUT_DIR/coverage.txt"

    echo "  Coverage:"
    cat "$OUT_DIR/coverage.txt"

    "$BCFTOOLS" mpileup -f "$REF_FASTA" "$OUT_DIR/aligned.bam" \
    | "$BCFTOOLS" call -mv -Oz -o "$OUT_DIR/variants.vcf.gz"

    "$BCFTOOLS" index "$OUT_DIR/variants.vcf.gz"

    "$BCFTOOLS" consensus \
        -f "$REF_FASTA" "$OUT_DIR/variants.vcf.gz" \
        > "$OUT_DIR/consensus.fasta"

    sed -i "s/^>.*/>$SAMPLE_ID\_$VIRUS_ACC/" "$OUT_DIR/consensus.fasta"

    echo "  Done: $VIRUS_SAMPLE → $OUT_DIR/consensus.fasta"
done

# ============================================================
echo ""
echo "=== PIPELINE COMPLETE ==="
echo "  Raw FASTQs          : $FASTQ_DIR"
echo "  QC report           : $QC_DIR/multiqc/multiqc_report.html"
echo "  Krona plots         : $KRONA_DIR/"
echo "  Comprehensive Krona : $KRONA_DIR/*_comprehensive_krona.html"
echo "  DIAMOND tables      : $VIRAL_DIR/layer2_divergent/*_diamond_summary.tsv"
echo "  Host-depleted reads : $HOST_DIR"
echo "  Viral detection     : $VIRAL_DIR/abundance/"
echo "  Genome assemblies   : $ASSEMBLY_DIR"
echo "  Assembly QC         : $QC_ASSEMBLY_DIR"
