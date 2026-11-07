
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 2     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 3     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 4     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 5     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 10     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
 # python 1.Simulated_Communities.py     --fa_dir "/cluster/tufts/hussainlab/tkwong01/15crispatus"     --output_dir "./iss_communities"     --num_species 15     --num_timepoints 30     --std_dev 0.2     --trend_type random     --seed 42
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def generate_ascii_graph(df, num_species, num_timepoints, std_dev, trend_type, height=15, width=60):
    """Generates an ASCII line graph for strain abundances over time, including std_dev and trend_type in the title."""
    strains = df["Strain"].unique()
    symbols = ["*", "o", ".", "+", "x", "#", "@", "%", "&", "=", "-", "~", "^", "!", ":"]
    
    # Initialize grid with spaces
    grid = [[" " for _ in range(width)] for _ in range(height)]
    
    # Map time steps to columns (0 to width-1)
    # Map abundances (0.0 to 1.0) to rows (height-1 to 0)
    for idx, strain in enumerate(strains):
        sym = symbols[idx % len(symbols)]
        strain_data = df[df["Strain"] == strain].sort_values("Time_Step")
        
        for _, row in strain_data.iterrows():
            t = int(row["Time_Step"])
            abun = float(row["Abundance"])
            
            # Compute grid coordinates
            col = int(round((t / max(1, num_timepoints - 1)) * (width - 1)))
            row_idx = int(round((1.0 - min(max(abun, 0.0), 1.0)) * (height - 1)))
            
            grid[row_idx][col] = sym

    # Build border and axes with updated title containing std and trend
    graph_lines = []
    graph_lines.append(f"\n--- Abundance Trajectories ({num_species} Strains, {num_timepoints} TPs, Trend: {trend_type}, SD: {std_dev}) ---")
    graph_lines.append(f"Y-axis: Abundance [0.0 to 1.0] | X-axis: Time (T1 to T{num_timepoints})")
    graph_lines.append("+" + "-" * width + "+")
    
    for r in range(height):
        # Calculate corresponding y-value label for reference edges
        y_val = 1.0 - (r / (height - 1))
        row_str = "".join(grid[r])
        if r == 0 or r == height // 2 or r == height - 1:
            graph_lines.append(f"|{row_str}| {y_val:.2f}")
        else:
            graph_lines.append(f"|{row_str}|")
            
    graph_lines.append("+" + "-" * width + "+")
    
    # Legend
    graph_lines.append("Legend (Strains & Symbols):")
    for idx, strain in enumerate(strains):
        sym = symbols[idx % len(symbols)]
        graph_lines.append(f"  [{sym}] {strain}")
        
    return "\n".join(graph_lines)


def generate_time_series(
    fa_dir: str,
    output_dir: str,
    num_species: int = 5,
    num_timepoints: int = 5,
    std_dev: float = 0.2,
    trend_type: str = "random",
    seed: int = 42,
):
    np.random.seed(seed)
    fa_path = Path(fa_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Gather available .fa files
    fa_files = sorted(list(fa_path.glob("*.fa")))
    if len(fa_files) < num_species:
        raise ValueError(
            f"Requested {num_species} species, but only found {len(fa_files)} .fa files in {fa_dir}"
        )

    # Select a deterministic or random subset of species
    selected_fa = np.random.choice(
        fa_files, size=num_species, replace=False
    )
    strain_names = [f.stem.replace("_renamed", "") for f in selected_fa]
    print(f"Selected strains ({num_species}): {strain_names}")

    # 2. Define baseline coefficients and time trends
    base_log_abundances = np.random.normal(0, 0.5, size=num_species)

    # Time trend slopes per strain
    time_steps = np.arange(num_timepoints)
    trends = np.zeros((num_species, num_timepoints))

    for i in range(num_species):
        if trend_type == "growth":
            slope = 0.15 if i == 0 else -0.05
        elif trend_type == "decline":
            slope = -0.15 if i == 0 else 0.05
        elif trend_type == "cyclic":
            trends[i] = 0.3 * np.sin(
                2 * np.pi * time_steps / num_timepoints + (i * np.pi / 4)
            )
            continue
        else:  # 'random' or default
            slope = np.random.uniform(-0.1, 0.1)

        trends[i] = slope * time_steps

    # 3. Simulate abundances across time-points
    records = []
    for t_idx, t in enumerate(time_steps):
        tp_name = f"T{t+1}"
        abundance_file = out_path / f"abundances_{num_species}sp_{tp_name}.txt"

        # Combine base + trend + Gaussian noise (std_dev)
        noise = np.random.normal(0, std_dev, size=num_species)
        log_vals = base_log_abundances + trends[:, t_idx] + noise

        # Convert log values to proportions using Softmax to sum to 1.0
        exp_vals = np.exp(log_vals)
        abundances = exp_vals / np.sum(exp_vals)

        # Write ISS-compatible abundance file for this time-point
        with open(abundance_file, "w") as f:
            for fa, abun in zip(selected_fa, abundances):
                clean_path = fa.resolve().as_posix().strip()
                f.write(f"{clean_path}\t{abun:.6f}\n")
                records.append(
                    {
                        "Timepoint": tp_name,
                        "Time_Step": t,
                        "Strain": fa.stem.replace("_renamed", ""),
                        "Abundance": abun,
                    }
                )

        print(
            f"Generated {tp_name} abundances -> saved to {abundance_file.name}"
        )

    # Save a summary CSV for your tracking/visualization
    df = pd.DataFrame(records)
    summary_csv = out_path / f"community_{num_species}sp_summary.csv"
    df.to_csv(summary_csv, index=False)
    print(f"Saved longitudinal summary table to {summary_csv}")

    # 4. Generate and save ASCII graph with parameters in title
    ascii_graph = generate_ascii_graph(df, num_species, num_timepoints, std_dev, trend_type)
    print(ascii_graph)
    
    graph_txt = out_path / f"community_{num_species}sp_graph.txt"
    with open(graph_txt, "w") as f:
        f.write(ascii_graph + "\n")
    print(f"Saved ASCII text graph to {graph_txt}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate longitudinal microbial community abundance files with ASCII graph."
    )
    parser.add_argument(
        "--fa_dir",
        type=str,
        default="/cluster/tufts/hussainlab/tkwong01/15crispatus",
        help="Path to folder containing .fa files",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./iss_communities",
        help="Output directory",
    )
    parser.add_argument(
        "--num_species",
        type=int,
        default=4,
        help="Number of random species to include",
    )
    parser.add_argument(
        "--num_timepoints",
        type=int,
        default=5,
        help="How many time-points to simulate over time",
    )
    parser.add_argument(
        "--std_dev",
        type=float,
        default=0.15,
        help="Standard deviation (noise level) between time-points",
    )
    parser.add_argument(
        "--trend_type",
        type=str,
        choices=["growth", "decline", "cyclic", "random"],
        default="growth",
        help="Time trend pattern across time-points",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility"
    )

    args = parser.parse_args()
    generate_time_series(
        fa_dir=args.fa_dir,
        output_dir=args.output_dir,
        num_species=args.num_species,
        num_timepoints=args.num_timepoints,
        std_dev=args.std_dev,
        trend_type=args.trend_type,
        seed=args.seed,
    )
