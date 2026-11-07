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
COMM_DIR="$WORK_DIR/wgsim_communities" # Ensure it points to where your FASTQs live
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
        FASTQ="$COMM_DIR/${COMM_PREFIX}_T${tp}_reads.fastq"
        if [ -f "$FASTQ" ]; then
            echo "$FASTQ" >> "$READS_LIST"
        fi
    done

    if [ ! -s "$READS_LIST" ]; then
        echo "Error: No FASTQ files found for ${NSPS}sp in $COMM_DIR. Skipping."
        continue
    fi

    # --- Step 2: Run unitig-caller De Novo Solely from Reads ---
    echo "--- Running reference-free unitig-caller for ${NSPS}sp ---"

    unitig-caller --call \
        --reads "$READS_LIST" \
        --kmer 31 \
        --rtab \
        --pyseer \
        --out "$UNITIG_DIR/uc_${NSPS}sp" \
        --threads "$THREADS"

    echo "Finished unitig-caller for ${NSPS}sp."

    # --- Step 3: Build Bowtie2 Index using unitig-caller sequences ---
    UNITIG_FASTA="$UNITIG_DIR/unitigs_${NSPS}sp.fasta"
    
    if [ -f "$UNITIG_DIR/uc_${NSPS}sp.pyseer" ]; then
        awk 'NR>1 {print ">unitig_" NR-1 "\n" $1}' "$UNITIG_DIR/uc_${NSPS}sp.pyseer" > "$UNITIG_FASTA"
        
        echo "Building Bowtie2 index from de novo unitigs..."
        bowtie2-build --threads "$THREADS" "$UNITIG_FASTA" "$UNITIG_DIR/unitig_idx" > /dev/null 2>&1

        # --- Step 4: Map Time-Point FASTQs to Get Depth/Counts Matrix ---
        SUMMARY_MATRIX="$WORK_DIR/timeseries_counts_${NSPS}sp.tsv"
        NONZERO_MATRIX="$WORK_DIR/timeseries_counts_${NSPS}sp_nonzero.tsv"

        # 1. Extract clean unitig IDs and lengths
        grep "^>" "$UNITIG_FASTA" | sed 's/>//' > "$WORK_DIR/unitig_ids_${NSPS}.txt"
        awk '/^>/ {if (seq) print length(seq); printf "%s\n", substr($1,2); seq=""; next} {seq=seq $1} END {if (seq) print length(seq)}' "$UNITIG_FASTA" > "$WORK_DIR/unitig_lengths_${NSPS}.txt"

        # Initialize counts matrix with just the IDs and lengths
        cp "$WORK_DIR/unitig_ids_${NSPS}.txt" "$WORK_DIR/matrix_body_${NSPS}.txt"
        paste "$WORK_DIR/matrix_body_${NSPS}.txt" "$WORK_DIR/unitig_lengths_${NSPS}.txt" > "$WORK_DIR/matrix_merged_${NSPS}.txt"
        mv "$WORK_DIR/matrix_merged_${NSPS}.txt" "$WORK_DIR/matrix_body_${NSPS}.txt"

        # 2. Loop through each time point and append counts column
        for tp in {1..30}; do
            FASTQ="$COMM_DIR/${COMM_PREFIX}_T${tp}_reads.fastq"
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

            # Append this timepoint's counts as a new column
            paste "$WORK_DIR/matrix_body_${NSPS}.txt" "$COUNTS_FILE" > "$WORK_DIR/temp_merged_${NSPS}.txt"
            mv "$WORK_DIR/temp_merged_${NSPS}.txt" "$WORK_DIR/matrix_body_${NSPS}.txt"
        done

        # 3. Construct clean tab-separated header string
        HEADER="Unitig_ID\tLength"
        for tp in {1..30}; do
            HEADER="${HEADER}\tT${tp}"
        done

        # 4. Write full matrix with header
        echo -e "$HEADER" > "$SUMMARY_MATRIX"
        cat "$WORK_DIR/matrix_body_${NSPS}.txt" >> "$SUMMARY_MATRIX"

        # 5. Write non-zero filtered matrix (Columns 3 through 32 represent T1-T30)
        echo -e "$HEADER" > "$NONZERO_MATRIX"
        awk -F'\t' '
        {
            sum=0;
            for(i=3; i<=NF; i++) sum+=$i;
            if(sum > 0) print $0;
        }' "$WORK_DIR/matrix_body_${NSPS}.txt" >> "$NONZERO_MATRIX"

        # Cleanup temporary files
        rm -f "$WORK_DIR/matrix_body_${NSPS}.txt" "$WORK_DIR/unitig_ids_${NSPS}.txt" "$WORK_DIR/unitig_lengths_${NSPS}.txt"
        
        echo "Completed timeseries count matrices for ${NSPS}sp:"
        echo "  - Full: $SUMMARY_MATRIX"
        echo "  - Non-zero filtered: $NONZERO_MATRIX"
    fi
done

echo "All reference-free unitig analyses completed successfully!"
