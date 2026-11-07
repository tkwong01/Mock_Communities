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

        # --- Step 4: Map Time-Point FASTQs and Build Counts Matrix (Python-backed) ---
        SUMMARY_MATRIX="$WORK_DIR/timeseries_counts_${NSPS}sp.tsv"
        NONZERO_MATRIX="$WORK_DIR/timeseries_counts_${NSPS}sp_nonzero.tsv"
        MAPPING_SUMMARY_DIR="$WORK_DIR/mapping_temp_${NSPS}"
        mkdir -p "$MAPPING_SUMMARY_DIR"

        # 1. Run bowtie2 and samtools indexing for all 30 timepoints, saving count files
        for tp in {1..30}; do
            FASTQ="$COMM_DIR/${COMM_PREFIX}_T${tp}_reads.fastq"
            SAM="$MAPPING_DIR/T${tp}.sam"
            SORTED_BAM="$MAPPING_DIR/T${tp}.sorted.bam"
            COUNTS_FILE="$MAPPING_SUMMARY_DIR/T${tp}_counts.txt"

            if [ ! -f "$FASTQ" ]; then
                # If fastq missing, output zeros matching number of unitigs
                n_unitigs=$(grep -c "^>" "$UNITIG_FASTA")
                seq "$n_unitigs" | awk '{print "0"}' > "$COUNTS_FILE"
            else
                bowtie2 -x "$UNITIG_DIR"//unitig_idx -U "$FASTQ" -p "$THREADS" -S "$SAM" --very-fast > /dev/null 2>&1
                samtools view -bS "$SAM" | samtools sort -o "$SORTED_BAM"
                samtools index "$SORTED_BAM"
                samtools idxstats "$SORTED_BAM" | awk '{print $3}' > "$COUNTS_FILE"
                rm -f "$SAM"
            fi
        done

        # 2. Use Python to assemble, clean, calculate lengths, and filter zero-rows perfectly
        python3 -c '
import os

unitig_fasta = "'"$UNITIG_FASTA"'"
mapping_dir = "'"$MAPPING_SUMMARY_DIR"'"
out_full = "'"$SUMMARY_MATRIX"'"
out_nonzero = "'"$NONZERO_MATRIX"'"

# Parse unitig IDs and lengths from FASTA
unitigs = []
current_id = None
current_seq = []

with open(unitig_fasta) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                unitigs.append((current_id, len("".join(current_seq))))
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
if current_id:
    unitigs.append((current_id, len("".join(current_seq))))

# Load counts for T1 to T30
timepoints = [f"T{i}" for i in range(1, 31)]
counts_data = {}

for tp in timepoints:
    cfile = os.path.join(mapping_dir, f"{tp}_counts.txt")
    vals = []
    if os.path.exists(cfile):
        with open(cfile) as cf:
            for line in cf:
                vals.append(int(line.strip()))
    else:
        vals = [0] * len(unitigs)
    counts_data[tp] = vals

# Write matrices
header = ["Unitig_ID", "Length"] + timepoints

for target_path, filter_nonzero in [(out_full, False), (out_nonzero, True)]:
    with open(target_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for idx, (uid, ulen) in enumerate(unitigs):
            row_counts = [counts_data[tp][idx] for tp in timepoints]
            if filter_nonzero and sum(row_counts) == 0:
                continue
            row_str = [uid, str(ulen)] + [str(c) for c in row_counts]
            out.write("\t".join(row_str) + "\n")
print("Python matrix assembly and filtering complete.")
'

        # Cleanup temporary mapping counts directory
        rm -rf "$MAPPING_SUMMARY_DIR"
        
        echo "Completed timeseries count matrices for ${NSPS}sp:"
        echo "  - Full: $SUMMARY_MATRIX"
        echo "  - Non-zero filtered: $NONZERO_MATRIX"
    fi
done

echo "All reference-free unitig analyses completed successfully!"
