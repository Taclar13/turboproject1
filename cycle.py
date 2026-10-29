"""Four-station non-ideal, variable-cp air-standard Brayton cycle.

All temperatures/pressures are TOTAL; work is positive input/output magnitude.
Units: K, Pa, J/kg, kg/s, W. Constant working-fluid mass through all components.
"""
from dataclasses import dataclass
import math
from scipy.optimize import brentq
from properties import Air


@dataclass(frozen=True)
class Inputs:
    T1: float = 300.0
    Pamb: float = 101325.0
    TIT: float = 1400.0
    eta_c: float = 0.85
    eta_t: float = 0.88
    eta_b: float = 0.99
    eta_m: float = 0.99  # turbine fluid-to-shaft efficiency, BEFORE compressor load
    pi_in: float = 0.99
    pi_b: float = 0.96
    pi_ex: float = 0.99
    LHV: float = 50e6
    mass_flow: float = 10.0


class InfeasibleCycle(ValueError):
    """Physical or property-domain failure; never silently return a valid point."""


def bisect_root(function, lower, upper, xtol=1e-8):
    """Readable independent root solver for checking SciPy's Brent method."""
    f_low, f_high = function(lower), function(upper)
    if f_low == 0:
        return lower
    if f_high == 0:
        return upper
    if f_low * f_high > 0:
        raise ValueError("Bisection requires a sign-changing bracket")
    for _ in range(200):
        middle = (lower + upper) / 2
        f_mid = function(middle)
        if f_mid == 0 or upper - lower <= xtol:
            return middle
        if f_low * f_mid < 0:
            upper = middle
        else:
            lower, f_low = middle, f_mid
    raise RuntimeError("Bisection failed to converge")


def temperature_root(function, lower, upper, method):
    if method not in ("brentq", "bisection"):
        raise ValueError("Unknown root method")
    if function(lower) * function(upper) > 0:
        raise InfeasibleCycle("Temperature root outside property-valid bracket")
    if method == "bisection":
        return bisect_root(function, lower, upper)
    return brentq(function, lower, upper, xtol=1e-8, rtol=1e-12)


def validate(p, gas):
    for name, value in vars(p).items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
        if name.startswith(("eta_", "pi_")) and value > 1:
            raise ValueError(f"{name} must be <= 1")
    gas.check(p.T1)
    gas.check(p.TIT)
    if p.TIT <= p.T1:
        raise ValueError("TIT must exceed inlet temperature")


def solve_cycle(pressure_ratio, p=None, gas=None, method="brentq"):
    p, gas = p or Inputs(), gas or Air()
    validate(p, gas)
    if not math.isfinite(pressure_ratio) or pressure_ratio <= 1:
        raise InfeasibleCycle("Compressor pressure ratio must exceed 1")

    # Boundary assumption: large ambient reservoir; exhaust KE is negligible.
    # pi_ex represents total-pressure loss to that reservoir, not a nozzle model.
    P1 = p.pi_in * p.Pamb
    P2 = pressure_ratio * P1
    P3 = p.pi_b * P2
    P4 = p.Pamb / p.pi_ex
    if P3 <= P4:
        raise InfeasibleCycle("No turbine expansion available after pressure losses")

    h1, h3 = gas.h(p.T1), gas.h(p.TIT)
    # Gibbs + ideal gas + isentropic reference: ds = ds0 - R ln(Pout/Pin).
    entropy_c = lambda T: gas.s0(T) - gas.s0(p.T1) - gas.R * math.log(pressure_ratio)
    T2s = temperature_root(entropy_c, p.T1, gas.Tmax, method)
    # Definition of compressor isentropic efficiency; First Law.
    wc = (gas.h(T2s) - h1) / p.eta_c
    h2 = h1 + wc
    if h2 >= h3:
        raise InfeasibleCycle("Compressor exit reaches/exceeds prescribed TIT")
    energy_c = lambda T: gas.h(T) - h2
    T2 = T2s if p.eta_c == 1 else temperature_root(energy_c, T2s, p.TIT, method)

    # Air-standard combustor: equivalent fuel energy, NO added working-fluid mass.
    q_added = h3 - h2
    q_fuel = q_added / p.eta_b
    fuel_air_ratio = q_fuel / p.LHV

    entropy_t = lambda T: gas.s0(T) - gas.s0(p.TIT) - gas.R * math.log(P4 / P3)
    T4s = temperature_root(entropy_t, gas.Tmin, p.TIT, method)
    # Definition of turbine isentropic efficiency; First Law.
    wt = p.eta_t * (h3 - gas.h(T4s))
    h4 = h3 - wt
    energy_t = lambda T: gas.h(T) - h4
    T4 = T4s if p.eta_t == 1 else temperature_root(energy_t, T4s, p.TIT, method)

    wnet = p.eta_m * wt - wc
    if wnet <= 0:
        raise InfeasibleCycle("Non-positive net shaft work")
    eta = wnet / q_fuel
    if not 0 < eta < 1:
        raise InfeasibleCycle("Thermal efficiency outside physical bounds")

    # Independent recomputation from solved temperatures, rather than targets.
    hc, ht = gas.h(T2) - h1, h3 - gas.h(T4)
    q_exhaust = gas.h(T4) - h1  # relative to inlet reservoir enthalpy
    mechanical_loss = (1 - p.eta_m) * ht
    burner_loss = (1 - p.eta_b) * q_fuel
    balance = q_fuel - wnet - q_exhaust - mechanical_loss - burner_loss
    return dict(
        pressure_ratio=pressure_ratio, T1=p.T1, T2s=T2s, T2=T2,
        T3=p.TIT, T4s=T4s, T4=T4, P1=P1, P2=P2, P3=P3, P4=P4,
        wc=wc, wt=wt, wnet=wnet, efficiency=eta, q_added=q_added,
        q_fuel=q_fuel, f=fuel_air_ratio, power=p.mass_flow * wnet,
        fuel_flow=p.mass_flow * fuel_air_ratio,
        ds_c=gas.s0(T2) - gas.s0(p.T1) - gas.R * math.log(P2 / P1),
        ds_t=gas.s0(T4) - gas.s0(p.TIT) - gas.R * math.log(P4 / P3),
        residual_entropy_c=entropy_c(T2s), residual_entropy_t=entropy_t(T4s),
        residual_energy_c=hc - wc, residual_energy_t=ht - wt,
        residual_energy_b=gas.h(p.TIT) - gas.h(T2) - p.eta_b * q_fuel,
        residual_cycle_energy=balance,
    )
