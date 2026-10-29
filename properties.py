"""Thermally perfect air: 79% N2 / 21% O2 by mole; SI mass-specific units.

NASA7 cp coefficients from Cantera data/gri30.yaml (see data_sources.md).
Only cp coefficients are needed: h and s0 are exact integrals anchored at
200 K and integrated piecewise across 1000 K. No numerical quadrature.
"""
import math

T_MIN, T_BREAK, T_MAX = 200.0, 1000.0, 3500.0
R_UNIVERSAL = 8314.46261815324  # J/(kmol K)
M_AIR = 0.79 * 28.014 + 0.21 * 31.998  # kg/kmol
R_AIR = R_UNIVERSAL / M_AIR

# NASA7 fits with a shared valid range of 200--3500 K for N2 and O2.
# This lower bound permits the environmental study's 273.15 K cases without
# extrapolating the former 300 K GRI-Mech nitrogen fit.
N2_LOW = (3.53100528, -1.23660987e-4, -5.02999433e-7, 2.43530612e-9, -1.40881235e-12)
N2_HIGH = (2.95257626, 1.39690040e-3, -4.92631603e-7, 7.86010195e-11, -4.60755204e-15)
O2_LOW = (3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12)
O2_HIGH = (3.66096065, 6.56365811e-4, -1.41149485e-7, 2.05797658e-11, -1.29913248e-15)

# cp_mass = Ru / Mmix * sum(mole_fraction * cp_molar / Ru).
LOW = tuple(R_AIR * (0.79 * n + 0.21 * o) for n, o in zip(N2_LOW, O2_LOW))
HIGH = tuple(R_AIR * (0.79 * n + 0.21 * o) for n, o in zip(N2_HIGH, O2_HIGH))


def primitive(T, coefficients, entropy=False):
    """Antiderivative of cp, or cp/T for the Gibbs entropy relation."""
    if entropy:
        return coefficients[0] * math.log(T) + sum(
            coefficients[i] * T**i / i for i in range(1, 5))
    return sum(a * T**(i + 1) / (i + 1) for i, a in enumerate(coefficients))


class Air:
    R = R_AIR
    Tmin, Tmax = T_MIN, T_MAX

    def check(self, T):
        if not math.isfinite(T) or not self.Tmin <= T <= self.Tmax:
            raise ValueError(f"Temperature {T} K outside [{self.Tmin}, {self.Tmax}] K")

    def cp(self, T):
        self.check(T)
        a = LOW if T <= T_BREAK else HIGH
        return sum(coefficient * T**i for i, coefficient in enumerate(a))

    def integral(self, T, entropy=False):
        self.check(T)
        stop = min(T, T_BREAK)
        value = primitive(stop, LOW, entropy) - primitive(T_MIN, LOW, entropy)
        if T > T_BREAK:
            value += primitive(T, HIGH, entropy) - primitive(T_BREAK, HIGH, entropy)
        return value

    def h(self, T):
        """J/kg; sensible enthalpy relative to the selected 200 K reference."""
        return self.integral(T)

    def s0(self, T):
        """J/(kg K); temperature part of entropy relative to 200 K."""
        return self.integral(T, entropy=True)

    def gamma(self, T):
        cp = self.cp(T)
        return cp / (cp - self.R)


class ConstantAir(Air):
    """Constant-property verification model, not the final analysis model."""
    def __init__(self, cp=1005.0, gamma=1.4):
        self.cp_value = cp
        self.R = cp * (1 - 1 / gamma)
        self.Tmin, self.Tmax = 100.0, 5000.0

    def cp(self, T):
        self.check(T)
        return self.cp_value

    def h(self, T):
        self.check(T)
        return self.cp_value * (T - T_MIN)

    def s0(self, T):
        self.check(T)
        return self.cp_value * math.log(T / T_MIN)
