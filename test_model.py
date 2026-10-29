"""Verification against analytical thermodynamics and independent numerics.

Run: python -m unittest -v
These tests verify implementation; they do not validate against a real engine.
"""
from dataclasses import replace
import math
import unittest
import numpy as np
from scipy.integrate import quad
from properties import Air, ConstantAir, T_BREAK
from cycle import Inputs, solve_cycle, InfeasibleCycle
from study import sweep, optimize
from ambient_study import (density_ideal_air, volume_flow_for_reference,
                           inputs_at_ambient, evaluate_condition)


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.p = Inputs()

    def test_ideal_analytical_states_work_efficiency(self):
        gas = ConstantAir()
        p = replace(self.p, eta_c=1, eta_t=1, eta_b=1, eta_m=1,
                    pi_in=1, pi_b=1, pi_ex=1)
        for ratio in (2, 5, 10, 20, 40):
            with self.subTest(ratio=ratio):
                row = solve_cycle(ratio, p, gas)
                a = (1.4 - 1) / 1.4
                T2, T4 = p.T1 * ratio**a, p.TIT / ratio**a
                w = gas.cp_value * (p.TIT - T4 - T2 + p.T1)
                self.assertAlmostEqual(row["T2"], T2, places=6)
                self.assertAlmostEqual(row["T4"], T4, places=6)
                self.assertAlmostEqual(row["wnet"], w, places=4)
                self.assertAlmostEqual(row["efficiency"], 1 - ratio**(-a), places=9)

    def test_ideal_optimum_and_efficiency_endpoint(self):
        gas = ConstantAir()
        p = replace(self.p, eta_c=1, eta_t=1, eta_b=1, eta_m=1,
                    pi_in=1, pi_b=1, pi_ex=1)
        rows = sweep(p, gas=gas)
        optimum = optimize(p, "wnet", rows, gas)
        expected = (p.TIT / p.T1)**(1.4 / (2 * (1.4 - 1)))
        self.assertAlmostEqual(optimum["pressure_ratio"], expected, places=4)
        best_eta = optimize(p, "efficiency", rows, gas)
        self.assertEqual(best_eta["pressure_ratio"], 40)
        self.assertTrue(best_eta["at_search_boundary"])

    def test_lossy_constant_property_analytical_case(self):
        p, gas, ratio = self.p, ConstantAir(), 12.0
        row = solve_cycle(ratio, p, gas)
        a = gas.R / gas.cp_value
        T2 = p.T1 * (1 + (ratio**a - 1) / p.eta_c)
        expansion = p.pi_in * p.pi_b * p.pi_ex * ratio
        T4 = p.TIT * (1 - p.eta_t * (1 - expansion**(-a)))
        expected = p.eta_m * gas.cp_value * (p.TIT - T4) - gas.cp_value * (T2 - p.T1)
        self.assertAlmostEqual(row["wnet"], expected, places=4)

    def test_property_integrals_against_quadrature(self):
        gas = Air()
        for T in (500, 999, 1000, 1001, 1400, 2500, 3500):
            points = [T_BREAK] if T > T_BREAK else None
            h = quad(gas.cp, gas.Tmin, T, points=points, epsabs=1e-6)[0]
            s = quad(lambda x: gas.cp(x) / x, gas.Tmin, T, points=points, epsabs=1e-8)[0]
            self.assertAlmostEqual(h, gas.h(T), places=5)
            self.assertAlmostEqual(s, gas.s0(T), places=7)
        self.assertLess(abs(gas.h(1000 + 1e-7) - gas.h(1000 - 1e-7)), 1e-3)
        self.assertTrue(1000 < gas.cp(300) < 1020)
        self.assertTrue(1.39 < gas.gamma(300) < 1.41)

    def test_independent_root_solvers(self):
        for ratio in (2, 12, 25, 40):
            brent = solve_cycle(ratio, self.p)
            bisect = solve_cycle(ratio, self.p, method="bisection")
            for name in ("T2s", "T2", "T4s", "T4", "efficiency", "wnet"):
                self.assertAlmostEqual(brent[name], bisect[name], delta=1e-4)

    def test_energy_entropy_and_ordering_over_sweep(self):
        for ratio in np.linspace(2, 40, 41):
            row = solve_cycle(float(ratio), self.p)
            self.assertLess(row["T1"], row["T2s"])
            self.assertLess(row["T2s"], row["T2"])
            self.assertLess(row["T2"], row["T3"])
            self.assertLess(row["T4s"], row["T4"])
            self.assertLess(row["T4"], row["T3"])
            self.assertGreater(row["ds_c"], 0)
            self.assertGreater(row["ds_t"], 0)
            for name in ("residual_energy_c", "residual_energy_t", "residual_energy_b", "residual_cycle_energy"):
                self.assertLess(abs(row[name]) / row["q_fuel"], 1e-9)
            for name in ("residual_entropy_c", "residual_entropy_t"):
                self.assertLess(abs(row[name]), 1e-7)

    def test_mass_flow_scaling_and_no_fuel_mass_multiplier(self):
        base = solve_cycle(12, self.p)
        double = solve_cycle(12, replace(self.p, mass_flow=20))
        changed_lhv = solve_cycle(12, replace(self.p, LHV=25e6))
        self.assertAlmostEqual(double["power"], 2 * base["power"])
        self.assertEqual(double["efficiency"], base["efficiency"])
        self.assertEqual(changed_lhv["wnet"], base["wnet"])
        self.assertAlmostEqual(changed_lhv["f"], 2 * base["f"])

    def test_invalid_inputs_and_infeasible_points(self):
        with self.assertRaises(ValueError):
            solve_cycle(10, replace(self.p, eta_c=1.1))
        with self.assertRaises(ValueError):
            solve_cycle(10, replace(self.p, T1=199.99))  # below selected data range
        for ratio in (1.01, 1000):
            with self.assertRaises(InfeasibleCycle):
                solve_cycle(ratio, self.p)
        rows = sweep(self.p, lower=1.001, upper=200, count=50)
        self.assertTrue(any(not r["valid"] for r in rows))
        best = optimize(self.p, "wnet", rows)
        self.assertGreater(best["wnet"], 0)

    def test_optima_against_dense_grid(self):
        rows = sweep(self.p)
        dense = [solve_cycle(float(r), self.p) for r in np.linspace(2, 40, 2001)]
        for objective in ("efficiency", "wnet"):
            best = optimize(self.p, objective, rows)
            grid = max(dense, key=lambda row: row[objective])
            self.assertGreaterEqual(best[objective] + 1e-9, grid[objective])
            self.assertLess(abs(best["pressure_ratio"] - grid["pressure_ratio"]), .02)

    def test_constant_volume_flow_environmental_logic(self):
        reference = {"temperature_K": 288.15, "pressure_Pa": 101325.0,
                     "mass_flow_kg_s": 10.0}
        volume_flow = volume_flow_for_reference(reference)
        self.assertAlmostEqual(
            density_ideal_air(288.15, 101325.0) * volume_flow, 10.0, places=10)
        hot = inputs_at_ambient(self.p, 323.15, 101325.0, volume_flow)
        low_pressure = inputs_at_ambient(self.p, 288.15, 70000.0, volume_flow)
        self.assertLess(hot.mass_flow, reference["mass_flow_kg_s"])
        self.assertLess(low_pressure.mass_flow, reference["mass_flow_kg_s"])
        short_sweep = {"lower": 2.0, "upper": 12.0, "count": 31}
        standard = evaluate_condition(self.p, 288.15, 101325.0, volume_flow, short_sweep)
        hotter = evaluate_condition(self.p, 323.15, 101325.0, volume_flow, short_sweep)
        low_pressure = evaluate_condition(self.p, 288.15, 70000.0, volume_flow, short_sweep)
        high_pressure = evaluate_condition(self.p, 288.15, 105000.0, volume_flow, short_sweep)
        self.assertLess(hotter["work_opt_power"], standard["work_opt_power"])
        self.assertEqual(standard["valid_pressure_ratio_points"], 31)
        for name in ("eta_opt_efficiency", "work_opt_wnet"):
            self.assertAlmostEqual(low_pressure[name], high_pressure[name], places=8)
        for name in ("eta_opt_pressure_ratio", "work_opt_pressure_ratio"):
            self.assertAlmostEqual(low_pressure[name], high_pressure[name], places=5)
        self.assertAlmostEqual(high_pressure["work_opt_power"] / low_pressure["work_opt_power"],
                               105000 / 70000, places=10)


if __name__ == "__main__":
    unittest.main()
