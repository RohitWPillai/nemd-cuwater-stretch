#!/usr/bin/env python3
"""Build-side density calibration check.

Reads the rho(z) accumulation written by build.in (cuw_build_density.profile),
applies the gate's measurement contract (raw 0.5 A bins -> 1.5 A centred moving
average; central plateau = z_mid +/- 7.5 A), and reports the plateau density
against the shipping gate 0.997 +/- 0.015 g/cm3. Instructor-side; run after
`lmp_serial -in build.in` (or a candidate override) from this folder.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lammps_io import read_profile

M_H2O = 29.915          # g/cm3 per O-atom/A^3: 18.0153 g/mol / 6.02214e23 * 1e24
RHO_GATE = 0.997        # g/cm3, central-plateau target (300 K)
RHO_TOL = 0.015         # g/cm3, gate half-width
PLATEAU_HALF = 7.5      # A, plateau half-width about z_mid

def main():
    prof = "cuw_build_density.profile"
    if not os.path.exists(prof):
        sys.exit(f"{prof} not found - run `lmp_serial -in build.in` (with nprof > 0) first.")
    z, n, nO = read_profile(prof)          # bin centre (A), O count, O number density (A^-3)
    rho_raw = nO * M_H2O

    # 1.5 A centred moving average of the 0.5 A bins (3-bin boxcar, interior bins only)
    rho_s = (rho_raw[:-2] + rho_raw[1:-1] + rho_raw[2:]) / 3.0
    z_s = z[1:-1]

    # channel midpoint from the density-weighted centroid of the water slab
    # (walls are pinned and symmetric, so this equals the geometric gap midpoint)
    z_mid = (rho_raw * z).sum() / rho_raw.sum()
    plat = (z_s >= z_mid - PLATEAU_HALF) & (z_s <= z_mid + PLATEAU_HALF)
    rho_plat = rho_s[plat].mean()
    n_tot = n.sum()                        # time-averaged O count = molecule count

    dev = rho_plat - RHO_GATE
    ok = abs(dev) <= RHO_TOL
    n_next = round(n_tot * RHO_GATE / rho_plat)

    print("\nBuild calibration: central plateau density")
    print(f"    molecules in profile   N = {n_tot:.1f}")
    print(f"    channel midpoint       z_mid = {z_mid:.2f} A")
    print(f"    plateau rho (smoothed) = {rho_plat:.4f} g/cm3 "
          f"({plat.sum()} bins, z_mid +/- {PLATEAU_HALF} A)")
    print(f"    gate {RHO_GATE} +/- {RHO_TOL} g/cm3: "
          f"{'PASS' if ok else 'FAIL'} (deviation {dev:+.4f})")
    if not ok:
        print(f"    proportional nwater for {RHO_GATE}: {n_next}")

if __name__ == "__main__":
    main()
