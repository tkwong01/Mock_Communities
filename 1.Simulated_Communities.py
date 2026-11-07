# Generates longitudinal abundance files for ONE mock community (N species x T time-points).
# Called once per (dataset, species-count) by run_mock_communities.sh.
#
# Example:
# python 1.Simulated_Communities.py --fa_dir /path/to/15crispatus \
#     --output_dir ./mock_datasets/dataset_01/3sp --num_species 3 \
#     --num_timepoints 30 --std_dev 0.2 --trend_type random --seed 42
#
# Output layout inside --output_dir:
#   abundances/abundances_{N}sp_T{t}.txt   ISS/wgsim-style "path<TAB>abundance" files
#   community_{N}sp_summary.csv            ground-truth abundances (long format)
#   community_{N}sp_graph.txt              ASCII plot of the trajectories
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def generate_ascii_graph(df, num_species, num_timepoints, std_dev, trend_type, height=15, width=60):
    """ASCII line graph of strain abundances over time."""
    strains = df["Strain"].unique()
    symbols = ["*", "o", ".", "+", "x", "#", "@", "%", "&", "=", "-", "~", "^", "!", ":"]
    grid = [[" " for _ in range(width)] for _ in range(height)]

    for idx, strain in enumerate(strains):
        sym = symbols[idx % len(symbols)]
        for _, row in df[df["Strain"] == strain].sort_values("Time_Step").iterrows():
            col = int(round((int(row["Time_Step"]) / max(1, num_timepoints - 1)) * (width - 1)))
            row_idx = int(round((1.0 - min(max(float(row["Abundance"]), 0.0), 1.0)) * (height - 1)))
            grid[row_idx][col] = sym

    lines = [
        f"\n--- Abundance Trajectories ({num_species} Strains, {num_timepoints} TPs, Trend: {trend_type}, SD: {std_dev}) ---",
        f"Y-axis: Abundance [0.0 to 1.0] | X-axis: Time (T1 to T{num_timepoints})",
        "+" + "-" * width + "+",
    ]
    for r in range(height):
        y_val = 1.0 - (r / (height - 1))
        row_str = "".join(grid[r])
        if r in (0, height // 2, height - 1):
            lines.append(f"|{row_str}| {y_val:.2f}")
        else:
            lines.append(f"|{row_str}|")
    lines.append("+" + "-" * width + "+")
    lines.append("Legend (Strains & Symbols):")
    for idx, strain in enumerate(strains):
        lines.append(f"  [{symbols[idx % len(symbols)]}] {strain}")
    return "\n".join(lines)


def generate_time_series(fa_dir, output_dir, num_species=5, num_timepoints=30,
                         std_dev=0.2, trend_type="random", seed=42):
    rng = np.random.RandomState(seed)  # local RNG: no global state side effects
    fa_files = sorted(Path(fa_dir).glob("*.fa"))
    if len(fa_files) < num_species:
        raise ValueError(f"Requested {num_species} species, but only found {len(fa_files)} .fa files in {fa_dir}")

    out_path = Path(output_dir)
    abun_dir = out_path / "abundances"
    abun_dir.mkdir(parents=True, exist_ok=True)

    sel_idx = rng.choice(len(fa_files), size=num_species, replace=False)
    selected_fa = [fa_files[i] for i in sel_idx]
    strain_names = [f.stem.replace("_renamed", "") for f in selected_fa]
    print(f"Selected strains ({num_species}, seed={seed}): {strain_names}")

    base_log_abundances = rng.normal(0, 0.5, size=num_species)
    time_steps = np.arange(num_timepoints)
    trends = np.zeros((num_species, num_timepoints))
    for i in range(num_species):
        if trend_type == "growth":
            slope = 0.15 if i == 0 else -0.05
        elif trend_type == "decline":
            slope = -0.15 if i == 0 else 0.05
        elif trend_type == "cyclic":
            trends[i] = 0.3 * np.sin(2 * np.pi * time_steps / num_timepoints + (i * np.pi / 4))
            continue
        else:
            slope = rng.uniform(-0.1, 0.1)
        trends[i] = slope * time_steps

    records = []
    for t in time_steps:
        tp_name = f"T{t+1}"
        log_vals = base_log_abundances + trends[:, t] + rng.normal(0, std_dev, size=num_species)
        abundances = np.exp(log_vals) / np.exp(log_vals).sum()  # softmax -> sums to 1

        with open(abun_dir / f"abundances_{num_species}sp_{tp_name}.txt", "w") as f:
            for fa, abun in zip(selected_fa, abundances):
                f.write(f"{fa.resolve().as_posix()}\t{abun:.6f}\n")
                records.append({"Timepoint": tp_name, "Time_Step": int(t),
                                "Strain": fa.stem.replace("_renamed", ""), "Abundance": abun})

    df = pd.DataFrame(records)
    df.to_csv(out_path / f"community_{num_species}sp_summary.csv", index=False)

    graph = generate_ascii_graph(df, num_species, num_timepoints, std_dev, trend_type)
    print(graph)
    (out_path / f"community_{num_species}sp_graph.txt").write_text(graph + "\n")
    print(f"Wrote abundances, summary CSV and graph to {out_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate longitudinal mock community abundance files.")
    p.add_argument("--fa_dir", type=str, required=True, help="Folder containing .fa genome files")
    p.add_argument("--output_dir", type=str, required=True, help="Output directory for this community")
    p.add_argument("--num_species", type=int, default=4)
    p.add_argument("--num_timepoints", type=int, default=30)
    p.add_argument("--std_dev", type=float, default=0.2, help="Noise SD between time-points")
    p.add_argument("--trend_type", choices=["growth", "decline", "cyclic", "random"], default="random")
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()
    generate_time_series(a.fa_dir, a.output_dir, a.num_species, a.num_timepoints,
                         a.std_dev, a.trend_type, a.seed)
