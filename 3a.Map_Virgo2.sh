#!/bin/bash
#SBATCH --job-name=CustomRef_Mapping
#SBATCH --output=/cluster/tufts/hussainlab/tkwong01/mock_reads/logs/CustomRef_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# Exit immediately if any command fails
set -e
set -o pipefail

# --- User Configuration ---
CUSTOM_FASTA="/cluster/tufts/hussainlab/DATABASES/VIRGO2/FastaFiles/VIRGO2.fa"

# Define index prefix right alongside the FASTA (or in a persistent shared folder if read-only)
# Using a persistent index directory in your work dir to avoid permission errors on the database folder:
INDEX_DIR="/cluster/tufts/hussainlab/tkwong01/mock_reads/shared_bowtie_index"
INDEX_PREFIX="$INDEX_DIR/VIRGO2_idx"

# --- Directories & Environment ---
source activate /cluster/tufts/hussainlab/tkwong01/envs/wgsimenv
WORK_DIR="/cluster/tufts/hussainlab/tkwong01/mock_reads"
COMM_DIR="$WORK_DIR/wgsim_communities"
cd "$WORK_DIR"

THREADS=8
mkdir -p "$INDEX_DIR"

# --- Step 1: Check Custom FASTA Existence ---
if [ ! -f "$CUSTOM_FASTA" ]; then
    echo "Error: Custom FASTA file not found at $CUSTOM_FASTA. Exiting."
    exit 1
fi

# --- Step 2: Build Bowtie2 Index ONCE (if it doesn't already exist) ---
if [ -f "${INDEX_PREFIX}.1.bt2" ] || [ -f "${INDEX_PREFIX}.1.bt2l" ]; then
    echo "Found existing Bowtie2 index at $INDEX_PREFIX. Skipping index building."
else
    echo "Index not found. Building Bowtie2 index for custom reference (this happens only once)..."
    bowtie2-build --threads "$THREADS" "$CUSTOM_FASTA" "$INDEX_PREFIX"
    echo "Index building complete."
fi

# Loop through each community size setup (2 to 6 species)
for NSPS in 2; do
    echo "=================================================="
    echo "Processing ${NSPS}-Species Communities (Custom Ref Mapping)"
    echo "=================================================="

    COMM_PREFIX="community_${NSPS}sp"
    CUSTOM_MAP_DIR="$WORK_DIR/${NSPS}sp_custom_mapping"
    mkdir -p "$CUSTOM_MAP_DIR"

    # --- Step 3: Map Time-Point FASTQs Across 30 Time Points ---
    CUSTOM_SUMMARY_MATRIX="$WORK_DIR/timeseries_counts_custom_${NSPS}sp.tsv"
    CUSTOM_NONZERO_MATRIX="$WORK_DIR/timeseries_counts_custom_${NSPS}sp_nonzero.tsv"
    CUSTOM_MAPPING_TEMP="$WORK_DIR/custom_mapping_temp_${NSPS}"
    mkdir -p "$CUSTOM_MAPPING_TEMP"

    for tp in {1..30}; do
        FASTQ="$COMM_DIR/${COMM_PREFIX}_T${tp}_reads.fastq"
        SAM="$CUSTOM_MAP_DIR/T${tp}.sam"
        SORTED_BAM="$CUSTOM_MAP_DIR/T${tp}.sorted.bam"
        COUNTS_FILE="$CUSTOM_MAPPING_TEMP/T${tp}_counts.txt"

        if [ ! -f "$FASTQ" ]; then
            n_refs=$(grep -c "^>" "$CUSTOM_FASTA")
            seq "$n_refs" | awk '{print "0"}' > "$COUNTS_FILE"
        else
            bowtie2 -x "$INDEX_PREFIX" \
                    -U "$FASTQ" \
                    -p "$THREADS" \
                    -S "$SAM" \
                    -N 0 \
                    -L 20 \
                    -D 20 \
                    -R 3 \
                    -k 10 \
                    --sam-no-qname-trunc \
                    --end-to-end \
                    --seed 343 \
                    --no-unal > /dev/null 2>&1

            samtools view -bS "$SAM" | samtools sort -o "$SORTED_BAM"
            samtools index "$SORTED_BAM"
            samtools idxstats "$SORTED_BAM" | awk '{print $3}' > "$COUNTS_FILE"
            rm -f "$SAM"
        fi
    done

    # --- Step 4: Python Assembly, Length Extraction, and Zero-Filtering ---
    python3 -c '
import os

custom_fasta = "'"$CUSTOM_FASTA"'"
mapping_dir = "'"$CUSTOM_MAPPING_TEMP"'"
out_full = "'"$CUSTOM_SUMMARY_MATRIX"'"
out_nonzero = "'"$CUSTOM_NONZERO_MATRIX"'"

# Parse reference IDs and sequence lengths from multi-FASTA
references = []
current_id = None
current_seq = []

with open(custom_fasta) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                references.append((current_id, len("".join(current_seq))))
            current_id = line[1:].split()[0]  # Take first word as reference ID
            current_seq = []
        else:
            current_seq.append(line)
if current_id:
    references.append((current_id, len("".join(current_seq))))

# Load depth/counts data for T1 to T30
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
        vals = [0] * len(references)
    counts_data[tp] = vals

# Write full and non-zero filtered matrices
header = ["Reference_ID", "Length"] + timepoints

for target_path, filter_nonzero in [(out_full, False), (out_nonzero, True)]:
    with open(target_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for idx, (rid, rlen) in enumerate(references):
            row_counts = [counts_data[tp][idx] for tp in timepoints]
            if filter_nonzero and sum(row_counts) == 0:
                continue
            row_str = [rid, str(rlen)] + [str(c) for c in row_counts]
            out.write("\t".join(row_str) + "\n")

print("Custom reference timeseries matrix generation complete.")
'

    # Cleanup temporary mapping counts directory
    rm -rf "$CUSTOM_MAPPING_TEMP"
    
    echo "Completed custom reference timeseries matrices for ${NSPS}sp:"
    echo "  - Full: $CUSTOM_SUMMARY_MATRIX"
    echo "  - Non-zero filtered: $CUSTOM_NONZERO_MATRIX"
done

echo "All custom multi-FASTA mapping workflows completed successfully!"
