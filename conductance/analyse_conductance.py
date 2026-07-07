#!/usr/bin/env python3
"""Day 1 stretch - interfacial thermal conductance (Kapitza) at the Cu/water wall.

Reads the cumulative bath energy tallies and the per-bin kinetic ingredients
written by conductance.in and reports, per wall: the rate energy enters or
leaves through that bath and the interfacial temperature jump dT (the wall face
against the fluid conduction line, with the adhered first water layer excluded
from the fit). Real water reaches its steady thermal profile slowly, so a run
at the student defaults shows only the TRANSIENT and CANNOT measure G: the
conductance G = J / dT is quoted from the shipped reference alone, not from
this run. The summary flags the run as not steady whenever the two tallies
disagree by more than 30 %. The shipped reference/ folder holds a
200 ps run at the same gradient; the summary quotes its converged jumps and G,
and the figure overlays its T(z) and its steady tallies against this run's.
Saves the two-panel figure to cuw_conductance.png.
Run after `lmp_serial -in conductance.in`.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from lammps_io import read_profile_multi, read_params, read_timeseries

# Figure overlap gate: figcheck.py is in the instructor tree (references/tools,
# two directories above the package), not in stretch/ - a shipped copy of
# stretch/ runs without it and saves the figure unchecked.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                                "references", "tools"))
try:
    import figcheck
except ImportError:
    figcheck = None

KB = 8.617333e-5        # Boltzmann constant (eV/K)
MVV2E = 1.0364269e-4    # (g/mol)(A/ps)^2 -> eV (LAMMPS metal-units mvv2e)
EVPS_A2_TO_MW_M2 = 1.602177e7   # eV/ps per A^2 -> MW/m^2
ADHERED = 4.0           # A of water excluded from the fit at each face: the first
                        # adhered layer (density sheet: first peak ~2.4 A off the
                        # face, first minimum ~3.9 A) is bound to the wall and does
                        # not sit on the bulk conduction line


def corrected_T_profile(zk, ck, vals, dof):
    """Peculiar T(z) from the per-bin kinetic ingredients (cuw_ke*.profile):
    each bin's mass-weighted streaming velocity is subtracted and the kinetic
    energy divided over `dof` DOF per atom - 2 for SHAKE-rigid water (6 per
    molecule, not 9), 3 for Cu. Returns (z, count, T) over the occupied bins
    (window-mean count > 0.5: a bin visited by a fraction of an atom carries
    no measurable temperature)."""
    msq, mvx, mvy, mvz, mm = (vals[:, i] for i in range(5))
    occ = ck > 0.5
    pec = msq[occ] - (mvx[occ] ** 2 + mvy[occ] ** 2 + mvz[occ] ** 2) / mm[occ]
    return zk[occ], ck[occ], MVV2E * pec / (dof * KB)


def fit_extrap(z, T, z0, name):
    """OLS line through (z, T) -> (T at z0, SE of that value, slope, R^2).
    The SE is the prediction-mean error sqrt(sigma^2*(1/n + (z0-zbar)^2/Sxx)),
    so extrapolating further past the data costs more error. The same fit
    scores this run and the shipped reference."""
    if len(z) < 3:
        raise SystemExit(f"conductance fit: fewer than 3 occupied bins in the {name} window - "
                         "check the run completed and the binning/geometry.")
    n = len(z)
    sxx = np.sum((z - z.mean()) ** 2)
    if sxx <= 0:                              # all bins at one z: no gradient (se would /0)
        raise SystemExit(f"conductance fit: all {name} z coordinates are identical - "
                         "check the binning/geometry.")
    s, c = np.polyfit(z, T, 1)
    resid = T - (s * z + c)
    sig2 = np.sum(resid ** 2) / (n - 2)
    se = float(np.sqrt(sig2 * (1.0 / n + (z0 - z.mean()) ** 2 / sxx)))
    ss = np.sum((T - T.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 1e-6 else float("nan")  # flat/all-equal -> unresolved
    return s * z0 + c, se, s, r2


def rate_fit(t, F, nb=10):
    """Window-fitted exchange rate: the OLS slope of a cumulative bath tally
    against time (eV/ps), with the SE taken over nb equal sub-window increment
    rates (a cumulative curve integrates its own noise, so the naive fit SE
    would understate the error)."""
    s = np.polyfit(t, F, 1)[0]
    nper = len(t) // nb
    if nper < 2:                                # too few rows to split into nb blocks:
        return float(s), float("nan")          # a real SE is unavailable (do not fabricate 0)
    blocks = [(F[(i + 1) * nper - 1] - F[i * nper]) / (t[(i + 1) * nper - 1] - t[i * nper])
              for i in range(nb)]
    return float(s), float(np.std(blocks, ddof=1) / np.sqrt(nb))


def G_of(P, seP, j, sej, area):
    """Conductance at one interface: G = J / dT with J = |P| / area, errors in
    quadrature. Returns (None, None) unless the jump is resolved (positive and
    at least 3x its fit error) - dividing J by an unresolved dT gives noise,
    not a conductance."""
    if j <= 0 or j < 3.0 * sej or abs(P) < 1e-12:   # unresolved jump OR no bath power -> no G
        return None, None
    G = abs(P) / area * EVPS_A2_TO_MW_M2 / j
    return G, G * float(np.hypot(seP / abs(P), sej / j))


def score_run(rdir, name):
    """One run directory -> everything the summary quotes: per-wall energy
    rates and their imbalance, the corrected T(z) on both species, the fluid
    conduction fit, the interfacial jumps and (where resolved) G. `name` tags
    the run in fit error messages."""
    par = read_params(os.path.join(rdir, "cuw_params.txt"))
    zfb, zft = par["zface_bot"], par["zface_top"]
    zdb, zdt = par["zsplit_bot"], par["zsplit_top"]
    bothot = par["Tbot"] > par["Ttop"]
    sgn = 1.0 if bothot else -1.0

    # bath tallies -> energy in at the hot wall / out at the cold wall. The
    # Langevin tally counts energy EXTRACTED from the atoms as positive, so
    # the feeding hot bath runs negative and its sign is flipped for display.
    step, fb, ft = read_timeseries(os.path.join(rdir, "cuw_heat.dat"))
    t = step * par["dt"]
    Pb, seb = rate_fit(t, fb)
    Pt, sett = rate_fit(t, ft)
    denom = 0.5 * (abs(Pb) + abs(Pt))          # floor: both bath powers ~0 (unresolved window)
    imb = abs(Pb + Pt) / denom if denom > 1e-12 else float("nan")

    # corrected thermometers: water on 6 DOF per molecule, Cu on the full 3
    # DOF per atom; the pinned anchor planes (outside zsplit) are dropped
    zk, ck, vals = read_profile_multi(os.path.join(rdir, "cuw_ke.profile"), 5)
    zw, cw, Tw = corrected_T_profile(zk, ck, vals, 2.0)
    zk2, ck2, vals2 = read_profile_multi(os.path.join(rdir, "cuw_ke_wall.profile"), 5)
    zc, cc, Tc = corrected_T_profile(zk2, ck2, vals2, 3.0)

    # fluid conduction line: one fit across the interior, the adhered layer
    # excluded at each face, extrapolated back to both faces
    m = (zw >= zfb + ADHERED) & (zw <= zft - ADHERED)
    Twb, seWb, slope, r2 = fit_extrap(zw[m], Tw[m], zfb, name + " fluid-interior")
    Twt, seWt, _, _ = fit_extrap(zw[m], Tw[m], zft, name + " fluid-interior")

    # wall faces: per-bath fits over the thermostatted planes only
    mb = (zc > zdb + 0.25) & (zc < zfb + 0.75)
    mt = (zc > zft - 0.75) & (zc < zdt - 0.25)
    Tcb, seCb, _, _ = fit_extrap(zc[mb], Tc[mb], zfb, name + " bottom-bath")
    Tct, seCt, _, _ = fit_extrap(zc[mt], Tc[mt], zft, name + " top-bath")

    # jumps, signed so that heat flowing wall->water->wall reads positive
    jb, sejb = sgn * (Tcb - Twb), float(np.hypot(seCb, seWb))
    jt, sejt = sgn * (Twt - Tct), float(np.hypot(seCt, seWt))
    Gb, seGb = G_of(Pb, seb, jb, sejb, par["area"])
    Gt, seGt = G_of(Pt, sett, jt, sejt, par["area"])

    # the water bins handed to the plot are the fit's: inside the adhered
    # layer the density is structured, so a bin slices molecules
    # systematically (H shells without their O) and the per-atom thermometer
    # only reads a temperature where that slicing averages out - the
    # part-molecule bins sit 30-80 K off the local conduction line
    mid = np.abs(zw - 0.5 * (zfb + zft)) <= 7.5
    return {"par": par, "t": t, "t_ps": par["nprod"] * par["dt"],
            "hot": "bottom" if bothot else "top",
            "cold": "top" if bothot else "bottom",
            "e_in": -(fb if bothot else ft), "e_out": ft if bothot else fb,
            "rate_in": -(Pb if bothot else Pt), "se_in": seb if bothot else sett,
            "rate_out": Pt if bothot else Pb, "se_out": sett if bothot else seb,
            "imb": imb, "zw": zw[m], "Tw": Tw[m],
            "zbath": zc[mb | mt], "Tbath": Tc[mb | mt],
            "slope": slope, "r2": r2, "Twb": Twb, "Twt": Twt, "Tcb": Tcb, "Tct": Tct,
            "jb": jb, "sejb": sejb, "jt": jt, "sejt": sejt,
            "Gb": Gb, "seGb": seGb, "Gt": Gt, "seGt": seGt,
            "t_mid": float(np.sum(cw[mid] * Tw[mid]) / np.sum(cw[mid]))}


def load_reference():
    """Read reference/follow/ (a shipped 200 ps run at the default gradient)
    and score it exactly like this run. Returns the score_run dict, or None
    when the folder or its files are absent."""
    rdir = os.path.join(_HERE, "reference", "follow")
    need = ["cuw_heat.dat", "cuw_ke.profile", "cuw_ke_wall.profile", "cuw_params.txt"]
    if not all(os.path.exists(os.path.join(rdir, f)) for f in need):
        return None
    return score_run(rdir, "reference")


def plot(run, ref=None, out="cuw_conductance.png"):
    """Two panels. Left: T(z) across the channel (z vertical, side-on) with
    the interfacial jump marked at each wall - the wall-face temperature
    (circle) against the fluid conduction line extrapolated to the same face
    (square). Right: the cumulative bath tallies as energy in and energy out
    against time. ref = a load_reference() dict -> its FITTED steady answer
    underlies both panels in green: the converged conduction line with its
    wall-face circles on the left, the steady in/out rates on the right (the
    raw reference bins carry the same part-molecule and random-walk noise as
    this run's, only smaller - the fits are what a longer run buys)."""
    try:
        import matplotlib
        if not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")     # headless: save the PNG, no window
        import matplotlib.pyplot as plt
    except ImportError:
        print("    (matplotlib not found - skipping the plot)")
        return
    BLUE, RED, GREEN, GREY = "#1F3A5F", "#B23A2E", "#4A7C59", "#888888"
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.4, 4.0))
    par = run["par"]
    zfb, zft = par["zface_bot"], par["zface_top"]

    # left: T(z), z vertical. The gap between each red circle (wall face) and
    # blue square (fluid line at the same face) is the Kapitza jump dT. The
    # water scatter shows the bins the fit uses (see the score_run note on
    # the part-molecule bins inside the adhered layer).
    if ref is not None:
        rfb, rft = ref["par"]["zface_bot"], ref["par"]["zface_top"]
        zr = np.linspace(rfb, rft, 50)
        axL.plot(ref["slope"] * (zr - rfb) + ref["Twb"], zr, color=GREEN,
                 lw=1.3, ls="--",
                 label=f"shipped reference ({ref['t_ps']:.0f} ps, fitted)")
        axL.plot([ref["Tcb"], ref["Tct"]], [rfb, rft], "o", color=GREEN,
                 ms=6, ls="none")
    axL.scatter(run["Tw"], run["zw"], s=12, color=GREY, alpha=0.65,
                label=r"$T(z)$ bins (water + baths)")
    axL.scatter(run["Tbath"], run["zbath"], s=12, color=GREY, alpha=0.65)
    zl = np.linspace(zfb, zft, 50)
    tl = run["slope"] * (zl - zfb) + run["Twb"]
    axL.plot(tl, zl, color=BLUE, lw=1.4, label="fluid conduction fit")
    axL.axhline(zfb, color="#aaaaaa", ls=":", lw=0.8)
    axL.axhline(zft, color="#aaaaaa", ls=":", lw=0.8)
    for Tcf, Twf, zf in ((run["Tcb"], run["Twb"], zfb), (run["Tct"], run["Twt"], zft)):
        axL.plot([Tcf], [zf], "o", color=RED, ms=7, zorder=5)
        axL.plot([Twf], [zf], "s", color=BLUE, ms=6, zorder=5)
        axL.annotate("", xy=(Tcf, zf), xytext=(Twf, zf),
                     arrowprops=dict(arrowstyle="<->", color=RED, lw=1.0))
    xmin = min(run["Tw"].min(), run["Tbath"].min())
    xmax = max(run["Tw"].max(), run["Tbath"].max())
    axL.text(xmin, zfb + 0.7, r"$\Delta T_{\rm bot} = %+.1f$ K" % run["jb"],
             color=RED, fontsize=8, ha="left", va="bottom")
    axL.text(xmax, zft - 0.7, r"$\Delta T_{\rm top} = %+.1f$ K" % run["jt"],
             color=RED, fontsize=8, ha="right", va="top")
    axL.plot([], [], "o", color=RED, label="wall face")
    axL.plot([], [], "s", color=BLUE, label="fluid at face")
    axL.set_xlabel(r"$T(z)$ (K)")
    axL.set_ylabel(r"$z$ ($\mathrm{\AA}$)")
    axL.set_title("temperature jump at each wall")
    axL.legend(fontsize=7, frameon=False, loc="center left")

    # right: cumulative energy in at the hot wall vs out at the cold wall.
    # In a steady run the two curves climb together; while the run is still
    # transient they disagree - that gap is the printed imbalance.
    axR.plot(run["t"], run["e_in"], color=RED, lw=1.4,
             label=f"in at the hot ({run['hot']}) wall")
    axR.plot(run["t"], run["e_out"], color=BLUE, lw=1.4,
             label=f"out at the cold ({run['cold']}) wall")
    if ref is not None:
        tl2 = np.array([0.0, run["t"][-1]])
        axR.plot(tl2, ref["rate_in"] * tl2, color=GREEN, lw=1.2,
                 label="reference steady rate (in / out)")
        axR.plot(tl2, ref["rate_out"] * tl2, color=GREEN, lw=1.2, ls="--")
    axR.set_xlabel(r"$t$ (ps)")
    axR.set_ylabel(r"exchanged energy $E(t)$ (eV)")
    axR.set_title("energy in vs energy out  (imbalance %.0f%%)" % (run["imb"] * 100))
    axR.legend(fontsize=7, frameon=True, facecolor="white", framealpha=0.9,
               edgecolor="none", loc="upper left")

    fig.tight_layout()
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
    for f in ("cuw_heat.dat", "cuw_ke.profile", "cuw_ke_wall.profile", "cuw_params.txt"):
        if not os.path.exists(f):
            sys.exit(f"{f} not found. If you submitted with ../submit.sh, wait for the "
                     "job (squeue --me) to finish; otherwise run `lmp_serial -in "
                     "conductance.in` first.")
    par = read_params("cuw_params.txt")
    this_ps = par.get("nprod", 0) * par.get("dt", 0)   # this run's length, .get so a corrupt params can't KeyError

    # Reference FIRST: loaded before ANY this-run value (a params key OR a file read), so a
    # corrupt/incompatible this-run params or a truncated file still shows the converged answer.
    # A corrupt shipped reference degrades to no reference, not an abort.
    try:
        ref = load_reference()
    except (SystemExit, OSError, KeyError, ValueError):
        ref = None

    print("\nStretch sheet 4: interfacial conductance")

    # ---- this run's own short, transient measurement; any failure (no gradient,
    #      a degenerate fit window) degrades to the reference-only summary ----
    run = None
    try:
        if abs(par["Tbot"] - par["Ttop"]) < 1e-9:
            raise SystemExit("this run set Tbot == Ttop, so there is no thermal gradient and no "
                             "conductance to\n      measure; the shipped reference is shown below. "
                             "Rerun conductance.in with its\n      default gradient (Tbot 330, "
                             "Ttop 270).")
        run = score_run(".", "this run's")
        print(f"    wall baths           Tbot = {par['Tbot']:.0f} K / Ttop = {par['Ttop']:.0f} K   "
              f"(the {run['hot']} wall is hot)")
        print(f"    energy in / out      in {run['rate_in']:+.3f} +/- {run['se_in']:.3f} eV/ps at "
              f"the {run['hot']} wall,")
        print(f"                         out {run['rate_out']:+.3f} +/- {run['se_out']:.3f} eV/ps at "
              f"the {run['cold']} wall")
        print(f"    imbalance            {run['imb'] * 100:.1f} %   (energy in equals energy out "
              "in a steady run)")
        if run["imb"] > 0.30:
            print("      -> the two tallies have not converged - this run is shorter than the")
            print("         steady-state time. The jumps below are a snapshot of the transient,")
            print("         not the steady answer.")
        print(f"    conduction fit       dT/dz = {run['slope']:+.3f} K/A across the fluid interior "
              f"(R^2 = {run['r2']:.3f};")
        print(f"                         the adhered first {ADHERED:.0f} A of water at each face "
              "excluded)")
        print(f"    temperature jumps    dT_bot = {run['jb']:+.1f} +/- {run['sejb']:.1f} K, "
              f"dT_top = {run['jt']:+.1f} +/- {run['sejt']:.1f} K")
        # This run does NOT measure G. A student-length run is transient (SPEC 6.4:
        # "a student-length run CANNOT measure G"), so G = J / dT is quoted only from
        # the shipped reference below, not from this run - even when both jumps
        # happen to clear 3x their fit error, that ratio is a transient snapshot, not
        # the steady conductance. The message says which case this run fell in.
        if sum(g is not None for g in (run["Gb"], run["Gt"])) < 2:
            print("      (a jump that does not clear 3x its fit error is not resolved: G = J / dT")
            print("       there would divide by noise, so it is not quoted for this window)")
        else:
            print("      (both jumps clear 3x their fit error here, but a student-length run is")
            print("       still transient, so this run's G = J / dT is a snapshot of the transient,")
            print("       not the measured conductance - G is quoted from the shipped reference below)")
        print(f"    water temperature    {run['t_mid']:.1f} K at mid-channel over the "
              f"{run['t_ps']:.0f} ps production window")
    except (SystemExit, ZeroDivisionError, FloatingPointError, ValueError, KeyError, OSError) as e:
        # any numeric degeneracy OR an incompatible/truncated this-run file degrades to the
        # reference, never a fatal abort (true programmer errors like AttributeError still surface)
        print(f"    this run did not produce a usable conductance measurement ({e});")
        print("      the shipped reference below is the value to read.")
        run = None

    # shipped reference tier - ALWAYS shown, whatever this run did above
    if ref is not None:
        rpar = ref["par"]
        print(f"    shipped reference ({ref['t_ps']:.0f} ps of averaging after "
              f"{rpar['nequil'] * rpar['dt']:.0f} ps of settling; this run: "
              f"{this_ps:.0f} ps):")
        print(f"      energy in / out = {ref['rate_in']:+.3f} / {ref['rate_out']:+.3f} eV/ps  "
              f"(imbalance {ref['imb'] * 100:.1f} %)")
        print(f"      dT_bot = {ref['jb']:+.1f} +/- {ref['sejb']:.1f} K, "
              f"dT_top = {ref['jt']:+.1f} +/- {ref['sejt']:.1f} K")
        if ref["Gb"] is not None and ref["Gt"] is not None:
            print(f"      G_bot = {ref['Gb']:.0f} +/- {ref['seGb']:.0f} MW/m2K, "
                  f"G_top = {ref['Gt']:.0f} +/- {ref['seGt']:.0f} MW/m2K   "
                  "(G = J / dT per wall)")
        if run is not None and not (np.isclose(par.get("Tbot"), rpar["Tbot"]) and
                                    np.isclose(par.get("Ttop"), rpar["Ttop"])):
            print(f"      (no shipped reference at Tbot = {par.get('Tbot'):g} K, "
                  f"Ttop = {par.get('Ttop'):g} K - the figure shows this run alone)")
            ref = None
    if run is not None:
        plot(run, ref=ref)
    else:
        print("    (no figure for this run - it produced no usable temperature profile; the "
              "shipped reference values above are the ones to read)")


if __name__ == "__main__":
    main()
