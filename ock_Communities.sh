#!/bin/bash
#SBATCH --job-name=Mock_Communities
#SBATCH --output=/cluster/tufts/hussainlab/tkwong01/mock_reads/logs/Mock_Communities_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#
# Usage:
#   sbatch run_mock_communities.sh [NUM_DATASETS] ["SPECIES_LIST"] [BASE_SEED]
#
# Examples:
#   sbatch run_mock_communities.sh                    # 1 dataset, species 2 3 4 5, seed 42
#   sbatch run_mock_communities.sh 10                 # 10 datasets
#   sbatch run_mock_communities.sh 10 "2 3 5 10 15"   # custom community sizes
#   sbatch run_mock_communities.sh 10 "2 3 4 5" 100   # custom base seed
#
# Optional overrides via environment (e.g. sbatch --export=ALL,NUM_TIMEPOINTS=10 ...):
#   NUM_TIMEPOINTS (30)  STD_DEV (0.2)  TREND_TYPE (random)  COVERAGE (5)
#   GENOME_MB (2.5)      READ_LEN (150) KEEP_BAM (0)         KEEP_READS (1)
#
# Output layout:
#   $OUT_ROOT/dataset_01/
#       2sp/ abundances/  reads/  mapping/  community_2sp_summary.csv  community_2sp_graph.txt
#            timeseries_counts_2sp.tsv  timeseries_counts_2sp_nonzero.tsv
#       3sp/ ...
#   $OUT_ROOT/dataset_02/ ...
# Dataset i uses seed BASE_SEED + (i-1), so every dataset is different but reproducible.
# Re-submitting the same command resumes: finished reads/mappings are skipped.

set -euo pipefail

# ---------------- Arguments ----------------
NUM_DATASETS="${1:-1}"
SPECIES_LIST="${2:-2 3 4 5}"
BASE_SEED="${3:-42}"

: "${NUM_TIMEPOINTS:=30}"
: "${STD_DEV:=0.2}"
: "${TREND_TYPE:=random}"
: "${COVERAGE:=5}"
: "${GENOME_MB:=2.5}"
: "${READ_LEN:=150}"
: "${KEEP_BAM:=0}"
: "${KEEP_READS:=1}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---------------- Paths ----------------
FA_DIR="/cluster/tufts/hussainlab/tkwong01/15crispatus"
CUSTOM_FASTA="/cluster/tufts/hussainlab/DATABASES/VIRGO2/FastaFiles/VIRGO2.fa"
WORK_DIR="/cluster/tufts/hussainlab/tkwong01/mock_reads"
OUT_ROOT="$WORK_DIR/mock_datasets"
INDEX_DIR="$WORK_DIR/shared_bowtie_index"
INDEX_PREFIX="$INDEX_DIR/VIRGO2_idx"
GEN_SCRIPT="$WORK_DIR/1.Simulated_Communities.py"

# Python (pandas/numpy) from df_env2; wgsim/bowtie2/samtools from wgsimenv
PY="/cluster/tufts/hussainlab/tkwong01/envs/df_env2/bin/python3"
export PATH="/cluster/tufts/hussainlab/tkwong01/envs/wgsimenv/bin:$PATH"

mkdir -p "$WORK_DIR/logs" "$OUT_ROOT" "$INDEX_DIR"
cd "$WORK_DIR"

[ -f "$CUSTOM_FASTA" ] || { echo "Error: reference FASTA not found: $CUSTOM_FASTA"; exit 1; }
[ -f "$GEN_SCRIPT" ]   || { echo "Error: generator script not found: $GEN_SCRIPT"; exit 1; }

# ---------------- Bowtie2 index (built once) ----------------
if [ -f "${INDEX_PREFIX}.1.bt2" ] || [ -f "${INDEX_PREFIX}.1.bt2l" ]; then
    echo "Found existing Bowtie2 index at $INDEX_PREFIX."
else
    echo "Building Bowtie2 index (one time)..."
    bowtie2-build --threads "$THREADS" "$CUSTOM_FASTA" "$INDEX_PREFIX"
fi

echo "Datasets: $NUM_DATASETS | Species sizes: $SPECIES_LIST | Timepoints: $NUM_TIMEPOINTS | Base seed: $BASE_SEED"

# ---------------- Main loops ----------------
for (( d=1; d<=NUM_DATASETS; d++ )); do
    DS_NAME=$(printf "dataset_%02d" "$d")
    DS_SEED=$(( BASE_SEED + d - 1 ))

    for N in $SPECIES_LIST; do
        COMM_DIR="$OUT_ROOT/$DS_NAME/${N}sp"
        ABUN_DIR="$COMM_DIR/abundances"
        READS_DIR="$COMM_DIR/reads"
        MAP_DIR="$COMM_DIR/mapping"
        mkdir -p "$READS_DIR" "$MAP_DIR"

        echo "=================================================="
        echo " $DS_NAME | ${N} species | seed $DS_SEED"
        echo "=================================================="

        # ---- Step 1: abundances (skip if already generated) ----
        if [ -f "$ABUN_DIR/abundances_${N}sp_T${NUM_TIMEPOINTS}.txt" ]; then
            echo "Abundances already exist. Skipping generation."
        else
            "$PY" "$GEN_SCRIPT" \
                --fa_dir "$FA_DIR" \
                --output_dir "$COMM_DIR" \
                --num_species "$N" \
                --num_timepoints "$NUM_TIMEPOINTS" \
                --std_dev "$STD_DEV" \
                --trend_type "$TREND_TYPE" \
                --seed "$DS_SEED"
        fi

        TOTAL_READS=$(awk -v n="$N" -v mb="$GENOME_MB" -v c="$COVERAGE" -v l="$READ_LEN" \
                          'BEGIN{printf "%d", (n*mb*1e6*c)/l}')

        for (( tp=1; tp<=NUM_TIMEPOINTS; tp++ )); do
            ABUN_FILE="$ABUN_DIR/abundances_${N}sp_T${tp}.txt"
            FASTQ="$READS_DIR/community_${N}sp_T${tp}_reads.fastq"
            IDXSTATS="$MAP_DIR/T${tp}.idxstats.txt"
            BAM="$MAP_DIR/T${tp}.sorted.bam"

            # ---- Step 2: wgsim read simulation ----
            if [ -s "$IDXSTATS" ]; then
                : # mapping already finished, reads not needed
            elif [ ! -s "$FASTQ" ]; then
                TMP="$READS_DIR/tmp_T${tp}"
                rm -rf "$TMP"; mkdir -p "$TMP"
                while IFS=$'\t' read -r genome_path abundance; do
                    [ -z "$genome_path" ] && continue
                    genome_reads=$(awk -v t="$TOTAL_READS" -v a="$abundance" \
                                   'BEGIN{n=int(t*a); if (n<1 && a>0) n=1; print n}')
                    base=$(basename "$genome_path" .fa)
                    # varied wgsim seed per dataset/timepoint so reads differ across time-points
                    wgsim -1 "$READ_LEN" -N "$genome_reads" -S $(( DS_SEED*1000 + tp )) \
                          "$genome_path" "$TMP/${base}_R1.fq" /dev/null > /dev/null
                done < "$ABUN_FILE"
                cat "$TMP"/*_R1.fq > "$FASTQ"
                rm -rf "$TMP"
            fi

            # ---- Step 3: map to VIRGO2 and get per-reference counts ----
            if [ ! -s "$IDXSTATS" ]; then
                bowtie2 -x "$INDEX_PREFIX" -U "$FASTQ" -p "$THREADS" \
                        -N 0 -L 20 -D 20 -R 3 -k 10 \
                        --sam-no-qname-trunc --end-to-end --seed 343 --no-unal \
                        2> "$MAP_DIR/T${tp}.bowtie2.log" \
                  | samtools sort -@ "$THREADS" -o "$BAM" -
                samtools index "$BAM"
                samtools idxstats "$BAM" > "$IDXSTATS"
                [ "$KEEP_BAM" = "1" ] || rm -f "$BAM" "$BAM.bai"
                [ "$KEEP_READS" = "1" ] || rm -f "$FASTQ"
            fi
            echo "  $DS_NAME ${N}sp T${tp} done"
        done

        # ---- Step 4: assemble count matrices (lengths come from idxstats) ----
        "$PY" - "$MAP_DIR" "$NUM_TIMEPOINTS" \
                "$COMM_DIR/timeseries_counts_${N}sp.tsv" \
                "$COMM_DIR/timeseries_counts_${N}sp_nonzero.tsv" << 'PY'
import sys
import pandas as pd

map_dir, ntp, out_full, out_nz = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
ref, cols = None, {}
for t in range(1, ntp + 1):
    df = pd.read_csv(f"{map_dir}/T{t}.idxstats.txt", sep="\t", header=None,
                     names=["Reference_ID", "Length", "Mapped", "Unmapped"])
    df = df[df["Reference_ID"] != "*"].reset_index(drop=True)
    if ref is None:
        ref = df[["Reference_ID", "Length"]]
    cols[f"T{t}"] = df["Mapped"].to_numpy()
mat = pd.concat([ref, pd.DataFrame(cols)], axis=1)
mat.to_csv(out_full, sep="\t", index=False)
mat[mat.iloc[:, 2:].sum(axis=1) > 0].to_csv(out_nz, sep="\t", index=False)
print(f"Wrote {out_full} and {out_nz}")
PY
        echo "Completed $DS_NAME ${N}sp"
    done
done

echo "All $NUM_DATASETS mock datasets finished: $OUT_ROOT"
