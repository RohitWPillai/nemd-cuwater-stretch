#!/usr/bin/env python3
"""Day 1 stretch - interfacial density / layering in the Cu/water channel.

Reads the density profile written by density.in, converts the O-atom counts to
g/cm3 of water, reports the central plateau and the near-wall layering (both
read on a 1.5 A centred moving average of the raw 0.5 A bins), and saves a
plot of rho(z) to cuw_density.png. Run after `lmp_serial -in density.in`.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from lammps_io import read_profile, read_params, read_timeseries

# Figure overlap gate: figcheck.py is in the instructor tree (references/tools,
# two directories above the package), not in stretch/ - a shipped copy of
# stretch/ runs without it and saves the figure unchecked.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                                "references", "tools"))
try:
    import figcheck
except ImportError:
    figcheck = None

M_H2O = 18.01528 / 6.02214076e23   # g per water molecule
A3_TO_CM3 = 1.0e-24                # cm3 per A3


def smooth_1p5(rho):
    """Centred 3-bin moving average (1.5 A at the contract 0.5 A bin width)."""
    sm = rho.copy()
    sm[1:-1] = (rho[:-2] + rho[1:-1] + rho[2:]) / 3.0
    return sm


def peaks_in(z, sm, zlo, zhi):
    """Layering peaks in [zlo, zhi]: local maxima of the smoothed profile that
    clear 1.15x BOTH adjacent minima (each minimum searched to the previous/
    next local maximum, or to the window edge). Returns [(z, rho), ...]."""
    w = (z >= zlo) & (z <= zhi)
    if not w.any():
        return []
    s, e = w.argmax(), len(z) - 1 - w[::-1].argmax()
    idx = [i for i in range(max(s, 1), min(e, len(sm) - 2) + 1)
           if sm[i] > sm[i - 1] and sm[i] >= sm[i + 1]]
    out = []
    for k, i in enumerate(idx):
        lo = idx[k - 1] if k > 0 else s
        hi = idx[k + 1] if k + 1 < len(idx) else e
        if sm[i] >= 1.15 * sm[lo:i + 1].min() and sm[i] >= 1.15 * sm[i:hi + 1].min():
            out.append((z[i], sm[i]))
    return out


def plot(z, raw, sm, rho_bulk, zfb, zft, binw, out="cuw_density.png"):
    """Save rho(z) with z vertical, so it reads like the channel side-on."""
    try:
        import matplotlib
        if not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")     # headless: save the PNG, no window
        import matplotlib.pyplot as plt
    except ImportError:
        print("    (matplotlib not found - skipping the plot)")
        return
    m = (z > zfb - 1.5) & (z < zft + 1.5)     # the channel, walls cropped
    plt.plot(raw[m], z[m], color="#bbbbbb", lw=0.8,
             label=r"raw $%.1f\ \mathrm{\AA}$ bins" % binw)
    plt.plot(sm[m], z[m], lw=1.4,
             label=r"$%.1f\ \mathrm{\AA}$ smoothed" % (3.0 * binw))
    # mark the two values the sheet reads off the profile
    plt.axvline(rho_bulk, color="#999999", ls="--", lw=0.9,
                label=r"$\rho_{\mathrm{bulk}} \approx %.3f$ g/cm$^3$" % rho_bulk)
    ipk = sm.argmax()
    plt.plot(sm[ipk], z[ipk], "o", color="#B23A2E", ms=6,
             label=r"$\rho_{\mathrm{peak}} \approx %.2f$ g/cm$^3$" % sm[ipk])
    plt.axhline(zfb, color="#999999", ls=":", lw=0.9, label="Cu wall faces")
    plt.axhline(zft, color="#999999", ls=":", lw=0.9)
    plt.legend(fontsize=8, frameon=False, loc="center right")
    plt.xlabel(r"water density $\rho(z)$ (g/cm$^3$)")
    plt.ylabel(r"$z$ ($\mathrm{\AA}$)")
    plt.title("density profile across the Cu/water channel")
    plt.tight_layout()
    fig = plt.gcf()
    if figcheck is not None:
        figcheck.savefig(fig, out, strict=True, dpi=150)
    else:
        fig.savefig(out, dpi=150)
    print(f"    plot -> {out}")
    if os.environ.get("DISPLAY"):     # ssh -X: also pop the figure up on screen
        try:
            plt.show()
        except KeyboardInterrupt:      # Ctrl-C with the window up: fine, the file is already saved
            pass


def main():
    for f in ("cuw_density.profile", "cuw_params.txt"):
        if not os.path.exists(f):
            sys.exit(f"{f} not found. If you submitted with ../submit.sh, wait for the "
                     "job (squeue --me) to finish; otherwise run `lmp_serial -in density.in` first.")
    z, n, nden = read_profile("cuw_density.profile")
    par = read_params("cuw_params.txt")
    zfb, zft = par["zface_bot"], par["zface_top"]

    rho = nden * M_H2O / A3_TO_CM3          # O atoms per A3 -> g/cm3 of water
    sm = smooth_1p5(rho)

    # central plateau: the 15 A band around the channel midplane
    zmid = 0.5 * (zfb + zft)
    plat = (z >= zmid - 7.5) & (z <= zmid + 7.5)
    rho_bulk = sm[plat].mean()

    # layering: peaks within 10 A of each wall face; the first peak is the
    # one nearest the face
    pk_bot = peaks_in(z, sm, zfb, zfb + 10.0)
    pk_top = peaks_in(z, sm, zft - 10.0, zft)
    first_bot = pk_bot[0][1] if pk_bot else float("nan")
    first_top = pk_top[-1][1] if pk_top else float("nan")
    c_bot, c_top = first_bot / rho_bulk, first_top / rho_bulk

    watT = read_timeseries("cuw_density_T.dat", col=1)
    t_ps = par["nprod"] * par["dt"]

    print("\nStretch sheet 1: density")
    print(f"    central plateau      rho_bulk ~ {rho_bulk:.3f} g/cm3")
    print(f"    first-layer peak     rho_peak ~ {first_bot:.3f} g/cm3 (bottom) / "
          f"{first_top:.3f} g/cm3 (top)")
    print(f"    contrast vs plateau  {c_bot:.2f}x / {c_top:.2f}x   "
          f"(peaks within 10 A of each face: {len(pk_bot)} / {len(pk_top)})")
    print(f"    fluid temperature    {watT.mean():.1f} K over the {t_ps:.0f} ps production "
          f"window  (wall baths at {par['Tbot']:.0f} K)")
    if min(c_bot, c_top) > 1.2:
        print("    -> the water is NOT uniform: it stacks into layers against the copper wall.")
    else:
        print("    -> no layering left: the water meets the wall essentially uniform.")
    if par["Tbot"] > 300.0:
        print(f"    -> wall baths at {par['Tbot']:.0f} K, not 300 K: compare the contrast "
              "values with your 300 K run - hotter water smears the layering.")
    if par["binw"] != 0.5:
        print(f"    -> binw = {par['binw']:.1f} A, not 0.5 A: the smoothing window is "
              f"{3.0 * par['binw']:.1f} A, and the plateau/peak/contrast values above assume "
              "0.5 A raw bins with 1.5 A smoothing (the default binw).")
    plot(z, rho, sm, rho_bulk, zfb, zft, par["binw"])


if __name__ == "__main__":
    main()
