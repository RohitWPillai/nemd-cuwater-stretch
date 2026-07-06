#!/usr/bin/env python3
"""Shipping gates for data.cuw_channel (SPEC section 5.5).

Reads the outputs of validate.in (a 20 ps NVE+SHAKE window at 1 fs, no
thermostats) and scores the four gates:
  1. total-energy drift <= 2 meV/atom over the window (linear fit x window);
  2. no vacuum pocket on the smoothed rho(z): no bin < 0.5 x plateau in the
     central 60 % of the water gap, and no contiguous span >= 2 A of bins
     < 0.5 g/cm3 more than 6 A from each wall face;
  3. walls crystalline: bath-layer Cu MSD < 0.5 A2 over the window;
  4. central plateau rho = 0.997 +/- 0.015 g/cm3 on the smoothed profile
     (z_mid +/- 7.5 A, section 5.4 contract).
Instructor-side; run after `lmp_serial -in validate.in` from this folder.
Exits 1 if any gate fails.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lammps_io import read_profile, read_params, read_timeseries

import numpy as np

M_H2O = 29.915          # g/cm3 per O-atom/A^3: 18.0153 g/mol / 6.02214e23 * 1e24
DRIFT_GATE = 2.0        # meV/atom over the window
MSD_GATE = 0.5          # A^2 per atom
POCKET_FRAC = 0.5       # vacuum threshold as a fraction of the plateau (central test)
POCKET_ABS = 0.5        # g/cm3, absolute vacuum threshold (span test)
POCKET_SPAN = 2.0       # A, minimum contiguous low span that counts as a pocket
WALL_MARGIN = 6.0       # A, layering zone excluded from the span test
PLATEAU_HALF = 7.5      # A, plateau half-width about z_mid (section 5.4 contract)
PLAT_TARGET = 0.997     # g/cm3, central-plateau density target
PLAT_TOL = 0.015        # g/cm3, plateau tolerance

def main():
    for f in ("cuw_validate_params.txt", "cuw_validate_energy.dat",
              "cuw_validate_msd.dat", "cuw_validate_density.profile"):
        if not os.path.exists(f):
            sys.exit(f"{f} not found - run `lmp_serial -in validate.in` first.")
    p = read_params("cuw_validate_params.txt")
    natoms, dt = p["natoms"], p["dt"]
    zfb, zft = p["zface_bot"], p["zface_top"]

    # --- gate 1: energy drift (fit slope x window; robust to fluctuation) ---
    step, etot, _ = read_timeseries("cuw_validate_energy.dat")
    t = step * dt                                        # ps
    slope = np.polyfit(t, etot, 1)[0]                    # eV/ps
    window = t[-1] - t[0]
    drift = slope * window * 1000.0 / natoms             # meV/atom over the window
    ptp = (etot.max() - etot.min()) * 1000.0 / natoms    # peak-to-peak, for the record
    ok1 = abs(drift) <= DRIFT_GATE

    # --- gate 2: vacuum pockets on the smoothed profile ---
    z, n, nO = read_profile("cuw_validate_density.profile")
    rho_raw = nO * M_H2O
    rho_s = (rho_raw[:-2] + rho_raw[1:-1] + rho_raw[2:]) / 3.0   # 1.5 A boxcar
    z_s = z[1:-1]
    binw = z[1] - z[0]
    z_mid = (rho_raw * z).sum() / rho_raw.sum()
    rho_plat = rho_s[np.abs(z_s - z_mid) <= PLATEAU_HALF].mean()

    gap = zft - zfb
    central = np.abs(z_s - 0.5 * (zfb + zft)) <= 0.3 * gap
    rho_min_c = rho_s[central].min()
    ok2a = rho_min_c >= POCKET_FRAC * rho_plat

    inner = (z_s >= zfb + WALL_MARGIN) & (z_s <= zft - WALL_MARGIN)
    low = rho_s[inner] < POCKET_ABS
    run, longest = 0, 0
    for flag in low:
        run = run + 1 if flag else 0
        longest = max(longest, run)
    span = longest * binw
    ok2b = span < POCKET_SPAN

    # --- gate 3: wall crystallinity (bath-layer Cu; the anchors are pinned) ---
    _, msd_bot, msd_top = read_timeseries("cuw_validate_msd.dat")
    msd_max = max(msd_bot.max(), msd_top.max())
    ok3 = msd_max < MSD_GATE

    # --- gate 4: central plateau density (section 5.4 measurement contract) ---
    ok4 = abs(rho_plat - PLAT_TARGET) <= PLAT_TOL

    print("\nShipping gates for data.cuw_channel (20 ps NVE+SHAKE at 1 fs)")
    print(f"    1. energy drift        = {drift:+.4f} meV/atom over {window:.1f} ps "
          f"(peak-to-peak {ptp:.4f}): {'PASS' if ok1 else 'FAIL'} (gate |drift| <= {DRIFT_GATE})")
    print(f"    2a. central 60% of gap: min rho = {rho_min_c:.4f} g/cm3 vs "
          f"0.5 x plateau {rho_plat:.4f} = {POCKET_FRAC * rho_plat:.4f}: "
          f"{'PASS' if ok2a else 'FAIL'}")
    print(f"    2b. longest span < {POCKET_ABS} g/cm3 beyond {WALL_MARGIN} A from the faces "
          f"(z {zfb + WALL_MARGIN:.2f}..{zft - WALL_MARGIN:.2f} A) = {span:.1f} A: "
          f"{'PASS' if ok2b else 'FAIL'} (gate < {POCKET_SPAN} A)")
    print(f"    3. bath-layer Cu MSD   = {msd_max:.4f} A2 max "
          f"(bot {msd_bot.max():.4f} / top {msd_top.max():.4f}): "
          f"{'PASS' if ok3 else 'FAIL'} (gate < {MSD_GATE})")
    print(f"    4. central plateau rho = {rho_plat:.4f} g/cm3: "
          f"{'PASS' if ok4 else 'FAIL'} (gate {PLAT_TARGET} +/- {PLAT_TOL})")
    all_ok = ok1 and ok2a and ok2b and ok3 and ok4
    print(f"    shipping verdict: {'PASS - ship data.cuw_channel' if all_ok else 'FAIL - do not ship'}")
    if not all_ok:
        sys.exit(1)

if __name__ == "__main__":
    main()
