"""Ambient-condition sensitivity study using the same Brayton-cycle solver.

This module deliberately reuses cycle.py and study.py.  It does not change
component efficiencies or introduce compressor maps: at every ambient condition
it evaluates a prescribed-pressure-ratio design study over the same PR domain.
"""
from dataclasses import asdict, replace
from pathlib import Path
import argparse
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cycle import Inputs
from study import sweep, optimize

R_DRY_AIR = 287.05  # J/(kg K), ideal-gas density approximation for inlet air.


def density_ideal_air(temperature, pressure):
    """Return ideal-gas dry-air density [kg/m^3] at the stated total condition."""
    if not (math.isfinite(temperature) and math.isfinite(pressure)
            and temperature > 0 and pressure > 0):
        raise ValueError("Temperature and pressure must be finite and positive")
    return pressure / (R_DRY_AIR * temperature)


def volume_flow_for_reference(reference):
    """Define a notional inlet-flow capacity from the standard-day reference."""
    return reference["mass_flow_kg_s"] / density_ideal_air(
        reference["temperature_K"], reference["pressure_Pa"])


def inputs_at_ambient(base, temperature, pressure, volume_flow):
    """Preserve inlet volume flow and calculate the resulting air mass flow."""
    mass_flow = density_ideal_air(temperature, pressure) * volume_flow
    return replace(base, T1=float(temperature), Pamb=float(pressure), mass_flow=mass_flow)


def compact_point(row, prefix):
    """Select the metrics used in tables and plots from an optimized cycle point."""
    return {
        f"{prefix}_pressure_ratio": row["pressure_ratio"],
        f"{prefix}_efficiency": row["efficiency"],
        f"{prefix}_wnet": row["wnet"],
        f"{prefix}_power": row["power"],
        f"{prefix}_wc": row["wc"],
        f"{prefix}_q_fuel": row["q_fuel"],
        f"{prefix}_T4": row["T4"],
        f"{prefix}_fuel_flow": row["fuel_flow"],
    }


def evaluate_condition(base, temperature, pressure, volume_flow, pr_sweep):
    """Optimize one ambient point for efficiency and for net specific work."""
    p = inputs_at_ambient(base, temperature, pressure, volume_flow)
    rows = sweep(p, **pr_sweep)
    best_eta = optimize(p, "efficiency", rows)
    best_work = optimize(p, "wnet", rows)
    result = {
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "mass_flow_kg_s": p.mass_flow,
        "density_kg_m3": density_ideal_air(temperature, pressure),
        "valid_pressure_ratio_points": sum(row["valid"] for row in rows),
    }
    result.update(compact_point(best_eta, "eta_opt"))
    result.update(compact_point(best_work, "work_opt"))
    return result


def evaluate_grid(base, temperatures, pressures, volume_flow, pr_sweep):
    return [evaluate_condition(base, temperature, pressure, volume_flow, pr_sweep)
            for temperature in temperatures for pressure in pressures]


def write_csv(path, rows):
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_series(rows, variable, xlabel, output, name):
    x = np.array([row[variable] for row in rows])
    eta = np.array([100 * row["eta_opt_efficiency"] for row in rows])
    work = np.array([row["work_opt_wnet"] / 1000 for row in rows])
    power_eta = np.array([row["eta_opt_power"] / 1e6 for row in rows])
    power_work = np.array([row["work_opt_power"] / 1e6 for row in rows])
    pr_eta = np.array([row["eta_opt_pressure_ratio"] for row in rows])
    pr_work = np.array([row["work_opt_pressure_ratio"] for row in rows])
    wc = np.array([row["work_opt_wc"] / 1000 for row in rows])
    qfuel = np.array([row["work_opt_q_fuel"] / 1000 for row in rows])
    t4 = np.array([row["work_opt_T4"] for row in rows])
    mdot = np.array([row["mass_flow_kg_s"] for row in rows])

    fig, axes = plt.subplots(3, 2, figsize=(11, 11), constrained_layout=True)
    panels = (
        (axes[0, 0], x, eta, "Peak thermal efficiency (%)", "#175b83"),
        (axes[0, 1], x, work, "Peak net specific work (kJ/kg air)", "#175b83"),
        (axes[1, 0], x, power_eta, "Shaft power at efficiency optimum (MW)", "#175b83"),
        (axes[1, 1], x, power_work, "Shaft power at work optimum (MW)", "#c54b28"),
        (axes[2, 0], x, t4, "Turbine-exit temperature at work optimum (K)", "#175b83"),
        (axes[2, 1], x, mdot, "Air mass flow at fixed inlet volume flow (kg/s)", "#175b83"),
    )
    for ax, xdata, ydata, ylabel, color in panels:
        ax.plot(xdata, ydata, color=color, linewidth=2)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.grid(alpha=.2)

    ax = axes[0, 0].twinx()
    ax.plot(x, pr_eta, color="#c54b28", linestyle="--", label="PR at efficiency optimum")
    ax.set_ylabel("Compressor pressure ratio")
    ax.legend(loc="lower left", fontsize=8)
    ax = axes[0, 1].twinx()
    ax.plot(x, pr_work, color="#c54b28", linestyle="--", label="PR at work optimum")
    ax.set_ylabel("Compressor pressure ratio")
    ax.legend(loc="lower right", fontsize=8)
    ax = axes[2, 0].twinx()
    ax.plot(x, wc, color="#c54b28", linestyle="--", label="Compressor work")
    ax.plot(x, qfuel, color="#665191", linestyle=":", label="Fuel heat input")
    ax.set_ylabel("Specific energy (kJ/kg air)")
    ax.legend(loc="best", fontsize=8)

    fig.suptitle(f"Non-ideal Brayton cycle | {name}", fontsize=14)
    fig.savefig(output / f"{name.lower().replace(' ', '_')}.png", bbox_inches="tight", pad_inches=.15)
    plt.close(fig)


def plot_pressure_series(rows, output):
    """Plot only the pressure-dependent quantities at fixed ambient temperature.

    Specific cycle quantities are analytically invariant to ambient pressure in
    this prescribed-ratio model.  Plotting optimizer differences at 1e-10--1e-6
    scale would manufacture visual noise, so they are reported as constants.
    """
    pressure_kpa = np.array([row["pressure_Pa"] / 1000 for row in rows])
    density = np.array([row["density_kg_m3"] for row in rows])
    mdot = np.array([row["mass_flow_kg_s"] for row in rows])
    power_eta = np.array([row["eta_opt_power"] / 1e6 for row in rows])
    power_work = np.array([row["work_opt_power"] / 1e6 for row in rows])
    change = np.array([row["work_opt_power_change_pct_from_reference"] for row in rows])
    first = rows[0]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    axes[0, 0].plot(pressure_kpa, density, color="#175b83", linewidth=2)
    axes[0, 0].set(xlabel="Ambient total pressure (kPa)",
                    ylabel="Inlet density (kg/m³)",
                    title="Ideal-gas inlet density")
    axes[0, 1].plot(pressure_kpa, mdot, color="#175b83", linewidth=2)
    axes[0, 1].set(xlabel="Ambient total pressure (kPa)",
                    ylabel="Air mass flow (kg/s)",
                    title="Fixed inlet volume-flow capacity")
    axes[1, 0].plot(pressure_kpa, power_eta, color="#175b83", linewidth=2,
                    label="At efficiency optimum")
    axes[1, 0].plot(pressure_kpa, power_work, color="#c54b28", linewidth=2,
                    label="At specific-work optimum")
    axes[1, 0].set(xlabel="Ambient total pressure (kPa)", ylabel="Shaft power (MW)",
                    title="Installed shaft power")
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].plot(pressure_kpa, change, color="#c54b28", linewidth=2)
    axes[1, 1].axhline(0, color="#666666", linewidth=1)
    axes[1, 1].set(xlabel="Ambient total pressure (kPa)",
                    ylabel="Change from 101.325 kPa (%)",
                    title="Power change at specific-work optimum")
    for ax in axes.flat:
        ax.grid(alpha=.2)

    invariant = (
        "Specific cycle performance is invariant\n"
        f"at fixed Tₐₘᵦ = {first['temperature_K']:.2f} K\n"
        f"ηₜₕ,max = {100 * first['eta_opt_efficiency']:.3f}%, PR = {first['eta_opt_pressure_ratio']:.4f}\n"
        f"wₙₑₜ,max = {first['work_opt_wnet'] / 1000:.3f} kJ/kg, PR = {first['work_opt_pressure_ratio']:.4f}\n"
        "Pₐₘᵦ → ρ → ṁ → shaft power"
    )
    axes[1, 1].text(.03, .97, invariant, transform=axes[1, 1].transAxes,
                    ha="left", va="top", fontsize=7.8,
                    bbox={"boxstyle": "round,pad=.4", "facecolor": "#f2f2f2", "edgecolor": "#aaaaaa"})
    fig.suptitle("Ambient-pressure sensitivity | fixed inlet volume-flow capacity", fontsize=14)
    fig.savefig(output / "ambient_pressure_sweep.png", bbox_inches="tight", pad_inches=.15)
    plt.close(fig)


def plot_contours(rows, temperatures, pressures, output):
    shape = (len(temperatures), len(pressures))
    def grid(key, scale=1):
        return np.array([row[key] * scale for row in rows]).reshape(shape)

    fields = (
        (grid("eta_opt_efficiency", 100), "Peak thermal efficiency (%)"),
        (grid("eta_opt_pressure_ratio"), "PR at efficiency optimum"),
        (grid("work_opt_wnet", 1e-3), "Peak net specific work (kJ/kg air)"),
        (grid("work_opt_pressure_ratio"), "PR at work optimum"),
        (grid("work_opt_power", 1e-6), "Shaft power at work optimum (MW)"),
        (grid("mass_flow_kg_s"), "Air mass flow (kg/s)"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(11, 12), constrained_layout=True)
    for ax, (values, label) in zip(axes.flat, fields):
        contour = ax.contourf(pressures / 1000, temperatures, values, levels=16, cmap="viridis")
        fig.colorbar(contour, ax=ax, label=label)
        ax.set(xlabel="Ambient total pressure (kPa)", ylabel="Ambient total temperature (K)")
    fig.suptitle("Environmental sensitivity | fixed inlet volume-flow capacity", fontsize=14)
    fig.savefig(output / "environmental_contours.png", bbox_inches="tight", pad_inches=.15)
    plt.close(fig)


def run(config_path):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    base_config = json.loads((config_path.parent / config["base_config"]).read_text())
    base = Inputs(**base_config["inputs"])
    reference = config["reference_condition"]
    pr_sweep = config["pressure_ratio_sweep"]
    volume_flow = volume_flow_for_reference(reference)
    temperatures = np.linspace(config["temperature_sweep"]["lower"],
                               config["temperature_sweep"]["upper"],
                               config["temperature_sweep"]["count"])
    pressures = np.linspace(config["pressure_sweep"]["lower"],
                            config["pressure_sweep"]["upper"],
                            config["pressure_sweep"]["count"])
    reference_row = evaluate_condition(
        base, reference["temperature_K"], reference["pressure_Pa"], volume_flow, pr_sweep)
    temperature_rows = evaluate_grid(base, temperatures, [reference["pressure_Pa"]], volume_flow, pr_sweep)
    pressure_rows = evaluate_grid(base, [reference["temperature_K"]], pressures, volume_flow, pr_sweep)
    grid_rows = evaluate_grid(base, temperatures, pressures, volume_flow, pr_sweep)
    for row in temperature_rows + pressure_rows + grid_rows:
        row["work_opt_power_change_pct_from_reference"] = 100 * (
            row["work_opt_power"] / reference_row["work_opt_power"] - 1)

    output = Path(__file__).resolve().parent / "results" / "environmental"
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "temperature_sweep.csv", temperature_rows)
    write_csv(output / "pressure_sweep.csv", pressure_rows)
    write_csv(output / "environmental_grid.csv", grid_rows)
    plot_series(temperature_rows, "temperature_K", "Ambient total temperature (K)", output,
                "Ambient temperature sweep")
    plot_pressure_series(pressure_rows, output)
    plot_contours(grid_rows, temperatures, pressures, output)

    summary = {
        "base_cycle_inputs": asdict(base),
        "reference_condition": reference,
        "inlet_volume_flow_m3_s": volume_flow,
        "reference_result": reference_row,
        "temperature_sweep": config["temperature_sweep"],
        "pressure_sweep": config["pressure_sweep"],
        "pressure_ratio_sweep": pr_sweep,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    lines = ["# Environmental-sensitivity results", "",
             "The cycle retains fixed component efficiencies and a prescribed compressor pressure ratio.",
             "A constant inlet volume-flow capacity is defined by 10 kg/s at the standard-day reference.", "",
             f"Reference inlet volume flow: {volume_flow:.4f} m³/s", "",
             "| Reference objective | PR | Efficiency | Specific work (kJ/kg air) | Shaft power (MW) |",
             "|---|---:|---:|---:|---:|"]
    for label, prefix in (("Maximum efficiency", "eta_opt"), ("Maximum specific work", "work_opt")):
        lines.append(f"| {label} | {reference_row[f'{prefix}_pressure_ratio']:.4f} | "
                     f"{reference_row[f'{prefix}_efficiency']:.3%} | "
                     f"{reference_row[f'{prefix}_wnet']/1000:.3f} | "
                     f"{reference_row[f'{prefix}_power']/1e6:.4f} |")
    lines += ["", "This is an environmental sensitivity study, not map-based off-design engine modeling."]
    (output / "RESULTS.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(Path(__file__).with_name("environmental.json")))
    run(parser.parse_args().config)
