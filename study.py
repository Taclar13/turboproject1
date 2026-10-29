"""Run only the agreed pressure-ratio sweep and its two optimizations."""
from dataclasses import asdict
from pathlib import Path
import csv
import json
import argparse
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycle import Inputs, InfeasibleCycle, solve_cycle


def sweep(p, lower=2.0, upper=40.0, count=401, gas=None):
    if not 1 < lower < upper or count < 3:
        raise ValueError("Require 1 < lower < upper and count >= 3")
    rows = []
    for ratio in np.linspace(lower, upper, count):
        try:
            rows.append({"valid": True, "reason": "", **solve_cycle(float(ratio), p, gas)})
        except InfeasibleCycle as error:
            rows.append({"valid": False, "reason": str(error), "pressure_ratio": float(ratio)})
    return rows


def optimize(p, objective, rows, gas=None):
    """Refine all sampled local maxima and compare with feasible endpoints.

    The coarse sweep locates peaks; bounded minimization refines -objective.
    The result is a maximum within the specified sampled/refined design domain,
    not a claim about all possible pressure ratios or real engine geometries.
    """
    if objective not in ("efficiency", "wnet"):
        raise ValueError("Only efficiency and wnet are optimization objectives")
    valid = [r for r in rows if r["valid"]]
    if not valid:
        raise ValueError("No feasible cycles in sweep")
    candidates = list(valid)
    # Refine boundaries next to infeasible grid points using feasibility bisection.
    for left, right in zip(rows, rows[1:]):
        if left["valid"] == right["valid"]:
            continue
        good = left if left["valid"] else right
        bad = right if left["valid"] else left
        a, b = good["pressure_ratio"], bad["pressure_ratio"]
        for _ in range(40):
            mid = (a + b) / 2
            try:
                point = solve_cycle(mid, p, gas)
                a, good = mid, point
            except InfeasibleCycle:
                b = mid
        candidates.append(good)
    for i in range(1, len(rows) - 1):
        left, mid, right = rows[i-1:i+2]
        if not all(r["valid"] for r in (left, mid, right)):
            continue
        if mid[objective] >= max(left[objective], right[objective]):
            result = minimize_scalar(
                lambda ratio: -solve_cycle(ratio, p, gas)[objective],
                bounds=(left["pressure_ratio"], right["pressure_ratio"]),
                method="bounded", options={"xatol": 1e-8})
            if not result.success:
                raise RuntimeError(result.message)
            candidates.append(solve_cycle(float(result.x), p, gas))
    best = dict(max(candidates, key=lambda row: row[objective]))
    best["at_search_boundary"] = any(
        abs(best["pressure_ratio"] - r["pressure_ratio"]) < 1e-6
        for r in (rows[0], rows[-1]))
    return best


def plot_results(rows, best_eta, best_work, output):
    valid = [r for r in rows if r["valid"]]
    ratio = np.array([r["pressure_ratio"] for r in valid])
    work = np.array([r["wnet"] / 1000 for r in valid])
    efficiency = np.array([100 * r["efficiency"] for r in valid])
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 140})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, y, ylabel, best in zip(
            axes, (efficiency, work), ("Thermal efficiency (%)", "Net specific work (kJ/kg air)"),
            (best_eta, best_work)):
        ax.plot(ratio, y, color="#175b83", linewidth=2)
        val = 100 * best["efficiency"] if best is best_eta else best["wnet"] / 1000
        ax.scatter([best["pressure_ratio"]], [val], color="#c54b28", zorder=3)
        ax.axvline(best["pressure_ratio"], color="#c54b28", linestyle="--", alpha=.5)
        ax.set(xlabel="Compressor total-pressure ratio", ylabel=ylabel,
               title=f"Maximum at pressure ratio {best['pressure_ratio']:.3f}")
        ax.grid(alpha=.2)
    fig.suptitle("Non-ideal air-standard Brayton cycle | fixed inputs")
    fig.savefig(output / "pressure_ratio_sweep.png", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    # bbox_inches protects axis labels and the colorbar when a viewer uses a
    # constrained preview pane (for example, VS Code's image preview).
    fig, ax = plt.subplots(figsize=(8.5, 6.25), constrained_layout=True)
    ax.plot(work, efficiency, color="#c7cdd2", linewidth=1, zorder=0)
    sc = ax.scatter(work, efficiency, c=ratio, cmap="viridis", s=16)
    fig.colorbar(sc, ax=ax, label="Compressor total-pressure ratio")
    for best, label, offset, marker in (
            (best_eta, "Maximum efficiency", (-155, -55), "*"),
            (best_work, "Maximum specific work", (-145, -65), "D")):
        x, y = best["wnet"] / 1000, 100 * best["efficiency"]
        ax.scatter(x, y, s=110, marker=marker, color="#c54b28", edgecolor="white", zorder=4)
        ax.annotate(f"{label}\nPR = {best['pressure_ratio']:.3f}", (x, y),
                    xytext=offset, textcoords="offset points", fontsize=10,
                    arrowprops={"arrowstyle": "-", "color": "#666666"})
    ax.set(xlabel="Net specific work (kJ/kg air)", ylabel="Thermal efficiency (%)",
           title="Efficiency versus specific work at each pressure ratio")
    ax.grid(alpha=.2)
    fig.savefig(output / "efficiency_work_tradeoff.png", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def run(config_path):
    config = json.loads(Path(config_path).read_text())
    p = Inputs(**config["inputs"])
    rows = sweep(p, **config["sweep"])
    best_eta = optimize(p, "efficiency", rows)
    best_work = optimize(p, "wnet", rows)
    output = Path(__file__).resolve().parent / "results"
    output.mkdir(exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (output / "sweep.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {"inputs": asdict(p), "sweep": config["sweep"],
               "valid_count": sum(r["valid"] for r in rows),
               "maximum_efficiency": best_eta, "maximum_specific_work": best_work}
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    plot_results(rows, best_eta, best_work, output)
    lines = ["# Computed baseline results", "", "Illustrative assumptions, not a calibrated engine.", "",
             "| Objective | Pressure ratio | Efficiency | Specific work (kJ/kg) | Shaft power (MW) |",
             "|---|---:|---:|---:|---:|"]
    for label, row in (("Maximum efficiency", best_eta), ("Maximum specific work", best_work)):
        lines.append(f"| {label} | {row['pressure_ratio']:.4f} | {row['efficiency']:.3%} | "
                     f"{row['wnet']/1000:.3f} | {row['power']/1e6:.4f} |")
    lines += ["", f"Feasible sweep points: {summary['valid_count']} / {len(rows)}.",
              "Optimization is restricted to the configured pressure-ratio interval.",
              "All working-fluid mass flows are equal; f is an equivalent fuel requirement."]
    (output / "RESULTS.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(Path(__file__).with_name("baseline.json")))
    run(parser.parse_args().config)
