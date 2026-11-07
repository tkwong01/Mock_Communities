#!/bin/bash
#SBATCH --job-name=UnitigCaller_RefFree
#SBATCH --output=/cluster/tufts/hussainlab/tkwong01/mock_reads/logs/UnitigCaller_%A_%a.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# Exit immediately if any command fails
set -e
set -o pipefail

# --- Directories ---
WORK_DIR="/cluster/tufts/hussainlab/tkwong01/mock_reads"
cd "$WORK_DIR"

THREADS=4

# Loop through each community size setup (2 to 6 species)
for NSPS in 2 3 4 5 6; do
    echo "=================================================="
    echo "Processing ${NSPS}-Species Communities (Reference-Free)"
    echo "=================================================="

    COMM_PREFIX="community_${NSPS}sp"
    UNITIG_DIR="$WORK_DIR/${NSPS}sp_unitigs"
    MAPPING_DIR="$WORK_DIR/${NSPS}sp_mapping"
    
    mkdir -p "$UNITIG_DIR" "$MAPPING_DIR"

    # --- Step 1: Create Manifest of Time-Point FASTQ Files ---
    READS_LIST="$UNITIG_DIR/reads_list_${NSPS}sp.txt"
    > "$READS_LIST"
    
    for tp in {1..30}; do
        FASTQ="${COMM_PREFIX}_T${tp}_reads.fastq"
        if [ -f "$FASTQ" ]; then
            echo "$(realpath "$FASTQ")" >> "$READS_LIST"
        fi
    done

    # --- Step 2: Run unitig-caller De Novo Solely from Reads ---
    echo "--- Running reference-free unitig-caller for ${NSPS}sp ---"

    unitig-caller --call \
        --reads "$READS_LIST" \
        --kmer 31 \
        --rtab \
        --pyseer \
        --out "$UNITIG_DIR/uc_${NSPS}sp" \
        --threads "$THREADS"

    echo "Finished unitig-caller for ${NSPS}sp. Outputs generated: uc_${NSPS}sp.rtab and uc_${NSPS}sp.pyseer"

    # --- Step 3 (Optional): Build Bowtie2 Index using unitig-caller sequences ---
    # unitig-caller outputs the called unitig sequences in the .pyseer file (Column 1).
    # You can extract them into a FASTA file if you want to map back reads for depth counts:
    UNITIG_FASTA="$UNITIG_DIR/unitigs_${NSPS}sp.fasta"
    
    if [ -f "$UNITIG_DIR/uc_${NSPS}sp.pyseer" ]; then
        awk 'NR>1 {print ">unitig_" NR-1 "\n" $1}' "$UNITIG_DIR/uc_${NSPS}sp.pyseer" > "$UNITIG_FASTA"
        
        echo "Building Bowtie2 index from de novo unitigs..."
        bowtie2-build --threads "$THREADS" "$UNITIG_FASTA" "$UNITIG_DIR/unitig_idx" > /dev/null 2>&1

        # --- Step 4: Map Time-Point FASTQs to Get Depth/Counts Matrix ---
        SUMMARY_MATRIX="$WORK_DIR/timeseries_counts_${NSPS}sp.tsv"
        grep "^>" "$UNITIG_FASTA" | sed 's/>//' > "$WORK_DIR/unitig_ids_${NSPS}.txt"
        
        HEADER="Unitig_ID"
        for tp in {1..30}; do
            HEADER="${HEADER}\tT${tp}"
        done
        echo -e "$HEADER" > "$SUMMARY_MATRIX"

        cp "$WORK_DIR/unitig_ids_${NSPS}.txt" "$WORK_DIR/col_unitigs_${NSPS}.txt"

        for tp in {1..30}; do
            FASTQ="${COMM_PREFIX}_T${tp}_reads.fastq"
            SAM="$MAPPING_DIR/T${tp}.sam"
            SORTED_BAM="$MAPPING_DIR/T${tp}.sorted.bam"
            COUNTS_FILE="$MAPPING_DIR/T${tp}_counts.txt"

            if [ ! -f "$FASTQ" ]; then
                awk '{print "0"}' "$WORK_DIR/unitig_ids_${NSPS}.txt" > "$COUNTS_FILE"
            else
                bowtie2 -x "$UNITIG_DIR/unitig_idx" -U "$FASTQ" -p "$THREADS" -S "$SAM" --very-fast > /dev/null 2>&1
                samtools view -bS "$SAM" | samtools sort -o "$SORTED_BAM"
                samtools index "$SORTED_BAM"
                samtools idxstats "$SORTED_BAM" | awk '{print $3}' > "$COUNTS_FILE"
                rm -f "$SAM"
            fi

            paste "$WORK_DIR/col_unitigs_${NSPS}.txt" "$COUNTS_FILE" > "$WORK_DIR/temp_merged_${NSPS}.txt"
            mv "$WORK_DIR/temp_merged_${NSPS}.txt" "$WORK_DIR/col_unitigs_${NSPS}.txt"
        done

        cat "$WORK_DIR/col_unitigs_${NSPS}.txt" > "$SUMMARY_MATRIX"
        rm -f "$WORK_DIR/col_unitigs_${NSPS}.txt" "$WORK_DIR/unitig_ids_${NSPS}.txt"
        echo "Completed timeseries count matrix for ${NSPS}sp: $SUMMARY_MATRIX"
    fi
done

echo "All reference-free unitig analyses completed successfully!"
