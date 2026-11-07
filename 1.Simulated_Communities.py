import argparse
from pathlib import Path
import numpy as np
import pandas as pd


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
    # Base coefficients for each strain
    base_log_abundances = np.random.normal(0, 0.5, size=num_species)

    # Time trend slopes per strain
    time_steps = np.arange(num_timepoints)
    trends = np.zeros((num_species, num_timepoints))

    for i in range(num_species):
        if trend_type == "growth":
            # Strain 0 grows over time, others decline or stay flat
            slope = 0.15 if i == 0 else -0.05
        elif trend_type == "decline":
            # Strain 0 declines over time
            slope = -0.15 if i == 0 else 0.05
        elif trend_type == "cyclic":
            # Oscillatory sinusoidal trend
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
                f.write(f"{fa.resolve()}\t{abun:.6f}\n")
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

    # Save a summary CSV for your own tracking/visualization
    df = pd.DataFrame(records)
    summary_csv = out_path / f"community_{num_species}sp_summary.csv"
    df.to_csv(summary_csv, index=False)
    print(f"Saved longitudinal summary table to {summary_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate longitudinal microbial community abundance files."
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
