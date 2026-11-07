#!/bin/bash
#SBATCH --job-name=ISS_Mock_Communities
#SBATCH --output=/cluster/tufts/hussainlab/tkwong01/mock_reads/logs/ISS_Community_%A_%a.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=2-6

# Exit immediately if any command fails
set -e
set -o pipefail

# --- Directories ---
FA_DIR="/cluster/tufts/hussainlab/tkwong01/15crispatus"
WORK_DIR="/cluster/tufts/hussainlab/tkwong01/mock_reads"
COMM_DIR="$WORK_DIR/iss_communities"
LOG_DIR="$WORK_DIR/logs"

mkdir -p "$COMM_DIR" "$LOG_DIR"
cd "$WORK_DIR"

# Number of species for this array task (SLURM_ARRAY_TASK_ID goes from 2 to 6)
NUM_SPECIES=$SLURM_ARRAY_TASK_ID

echo "=== Starting Task for Community Size: $NUM_SPECIES Species ==="

# --- Step 1: Switch to df_env2 and Generate Abundances ---
echo "--- Activating df_env2 for abundance generation ---"
export PATH="/cluster/tufts/hussainlab/tkwong01/envs/df_env2/bin:$PATH"

python3 1.Simulated_Communities.py \
    --fa_dir "$FA_DIR" \
    --output_dir "$COMM_DIR" \
    --num_species "$NUM_SPECIES" \
    --num_timepoints 30 \
    --std_dev 0.2 \
    --trend_type random \
    --seed 42

# --- Step 2: Switch to VIRGO2_env and Run InSilicoSeq ---
echo "--- Activating VIRGO2_env for InSilicoSeq read simulation ---"
export PATH="/cluster/tufts/hussainlab/tkwong01/envs/iss_env/bin:$PATH"

echo "--- Running InSilicoSeq for all generated time-points (${NUM_SPECIES} species) ---"

for abundance_file in "$COMM_DIR"/abundances_${NUM_SPECIES}sp_T*.txt; do
    # Extract time-point name (e.g., T1, T2, ..., T30) from filename
    tp_name=$(basename "$abundance_file" | sed -E 's/.*_(T[0-9]+)\.txt/\1/')
    OUT_PREFIX="$COMM_DIR/community_${NUM_SPECIES}sp_${tp_name}"
    expected_reads_file="${OUT_PREFIX}_reads.fastq"

    # Skip if reads already generated (makes the job resumable/idempotent)
    if [ -f "$expected_reads_file" ]; then
        echo "Reads for ${NUM_SPECIES}sp at $tp_name already exist. Skipping."
        continue
    fi

    # Calculate total target reads for ~5x coverage (~2.5 Mb per genome * NUM_SPECIES)
    TOTAL_BASES=$(( NUM_SPECIES * 2500000 ))
    TOTAL_READS=$(( (TOTAL_BASES * 5) / 150 ))

    echo "Simulating reads for ${NUM_SPECIES}sp at $tp_name (Target reads: $TOTAL_READS)..."
    iss generate \
        --genomes "$abundance_file" \
        --model novaseq \
        --n_reads "$TOTAL_READS" \
        --output "$OUT_PREFIX" \
        --seed 42
done

echo "Successfully completed mock community generation and read simulation for ${NUM_SPECIES}-species setup."
