# Sources

## Cycle framework

Klaus Brun and Rainer Kurz, *Introduction to Gas Turbine Theory*, Solar Turbines,
third edition (2019), Chapter 1. Checked against the user-provided PDF:

- pp. 30-35: specific-work/efficiency distinction and process diagrams.
- pp. 36-39: ideal component/cycle framework.
- pp. 40-41: non-ideal efficiencies, pressure losses, variable gas properties,
  static-to-total corrections, and numerical solution discussion.

The code uses our four-station numbering and one equivalent turbine. It does
not reproduce the later split gas-producer/power-turbine architecture. The Gibbs
integral formulation is our derivation from thermodynamics, not a claim that the
book supplies these exact code equations. Solar's PDF is not redistributed.

Official book landing page: https://catsolar.my.salesforce-sites.com/gatedcontent

## Property coefficients

NASA7 coefficients for N2 and O2 are embedded in `properties.py`. They follow
the NASA Glenn thermodynamic-polynomial database reported in:

B. J. McBride, M. J. Zehe, and S. Gordon, *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*, NASA/TP-2002-211556
(2002).

Only the first five NASA7 coefficients are needed here:

```text
cp_molar/Ru = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
```

The selected N2 and O2 fits use 200--1000 K and 1000--3500 K ranges. The
mixture therefore uses the common 200--3500 K range, split at 1000 K. The
200 K lower bound supports the ambient-temperature study down to 273.15 K
without extrapolation.

The absolute enthalpy/entropy integration constants (NASA7 a6/a7) are unused.
Instead, integrate cp piecewise from a common 200 K reference. This prevents
enthalpy/entropy jumps from rounded fit constants; cp may have a tiny fit mismatch.

Composition: 0.79 N2 / 0.21 O2 by mole, no argon or water.
Molar masses: 28.014 / 31.998 kg/kmol.
Ru = 8314.46261815324 J/(kmol K).
cp_mass = Ru/Mmix * sum(y_i * cp_molar_i/Ru).

NASA7 documentation:
https://cantera.org/stable/reference/thermo/species-thermo.html

These are NASA7 coefficients. No NASA9 coefficients are used. The Cantera
license remains in the repository because the original project version used
Cantera-distributed coefficient data; the present coefficient provenance is the
NASA Glenn report above.
