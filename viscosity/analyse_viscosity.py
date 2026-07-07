#!/usr/bin/env python3
"""Day 1 stretch - shear viscosity of the confined water.

Reads the Couette velocity profile, the water-only shear-stress trace and the
per-bin kinetic ingredients written by viscosity.in and reports
eta = |p_xz| / (dvx/dz) in mPa*s: the water shear stress over the central
shear rate. The temperature panel is the peculiar T(z) - each bin's streaming
velocity subtracted and the kinetic energy divided over the 6 DOF a
SHAKE-rigid molecule keeps - so it reads thermal energy, not flow energy.
The shipped reference/ folder holds a 200 ps run at the same drive; the
summary quotes its converged eta, and the figure overlays its profile and
its mean stress against this run's running mean.
Saves the three-panel figure to cuw_viscosity.png.
Run after `lmp_serial -in viscosity.in`.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from lammps_io import read_profile, read_profile_multi, read_params, read_timeseries

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
BAR_PER_PS_TO_MPAS = 1.0e-4   # eta [mPa*s] = p_xz [bar] / slope [1/ps] * 1e-4


def corrected_T_profile(zk, ck, vals):
    """Peculiar water T(z) from the per-bin kinetic ingredients in
    cuw_ke.profile: each bin's mass-weighted streaming velocity is subtracted
    and the kinetic energy is divided over 2 DOF per atom (6 per SHAKE-rigid
    molecule, not 9). Returns (z, count, T) over the occupied bins (window-mean
    count > 0.5, the same threshold the vx bins use: a bin visited by a
    fraction of an atom carries no measurable temperature)."""
    msq, mvx, mvy, mvz, mm = (vals[:, i] for i in range(5))
    occ = ck > 0.5
    pec = msq[occ] - (mvx[occ] ** 2 + mvy[occ] ** 2 + mvz[occ] ** 2) / mm[occ]
    return zk[occ], ck[occ], MVV2E * pec / (2.0 * KB)


def t_mid_of(zo, co, tbin, zmid):
    """Count-weighted mean of the peculiar T(z) over |z - zmid| <= 7.5 A
    (the mid-channel window every stretch case shares)."""
    mid = np.abs(zo - zmid) <= 7.5
    return float(np.sum(co[mid] * tbin[mid]) / np.sum(co[mid]))


def central_fit(z, vx, name="cuw_vx.profile"):
    """Unweighted OLS line through the central half of the occupied span
    (away from the structured near-wall layers) -> slope, intercept, R^2,
    SE(slope). The same fit scores this run and the shipped reference."""
    if len(z) < 3:                            # empty/near-empty profile: guard min()/max() below
        raise SystemExit(f"viscosity fit: fewer than 3 occupied bins in {name} - check the run "
                         "completed and wrote a filled channel.")
    zc, W = 0.5 * (z.min() + z.max()), z.max() - z.min()
    cen = (z > zc - 0.25 * W) & (z < zc + 0.25 * W)
    if np.count_nonzero(cen) < 3:
        raise SystemExit(f"viscosity fit: fewer than 3 occupied bins in the central half of "
                         f"{name} - check the run completed and the binning/geometry.")
    zf, vf = z[cen], vx[cen]
    s, c = np.polyfit(zf, vf, 1)
    resid = vf - (s * zf + c)
    ss_tot = np.sum((vf - vf.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else float("nan")
    se_s = np.sqrt(np.sum(resid ** 2) / (len(zf) - 2) / np.sum((zf - zf.mean()) ** 2))
    return s, c, r2, se_s


def block_se(x, nb=10):
    """Standard error of the mean of a correlated series: split into nb equal
    blocks and take the SE of the block means (samples 0.1 ps apart are not
    independent, so the naive SE would understate the error)."""
    if len(x) < 2:
        return float("nan")                    # no SE from a single sample
    nper = len(x) // nb
    if nper < 1:
        return float(np.std(x, ddof=1) / np.sqrt(len(x)))
    b = x[:nb * nper].reshape(nb, nper).mean(axis=1)
    return float(b.std(ddof=1) / np.sqrt(nb))


def eta_of(par, z, vx, pxz_series, name="cuw_vx.profile"):
    """One run -> everything the summary quotes: the central fit, the mean
    stress with its block SE, and eta with the two errors in quadrature."""
    s, c, r2, se_s = central_fit(z, vx, name=name)
    pxz = float(pxz_series.mean())
    se_p = block_se(pxz_series)
    # floor both divisors: a near-zero shear rate (a flat/reversed short run) or a
    # near-zero mean stress (the zero-drive cross-check) is a legitimate outcome and
    # must not raise ZeroDivisionError before the reference is shown.
    eta = abs(pxz) / max(abs(s), 1e-6) * BAR_PER_PS_TO_MPAS
    se_eta = eta * float(np.hypot(se_p / max(abs(pxz), 1e-12), se_s / max(abs(s), 1e-6)))
    return {"s": s, "c": c, "r2": r2, "se_s": se_s,
            "pxz": pxz, "se_p": se_p, "eta": eta, "se_eta": se_eta}


def load_reference():
    """Read reference/follow/ (a shipped 200 ps run at the default drive) and
    score it exactly like this run. Returns a dict with the fit, eta and the
    corrected mid-channel T, or None when the folder or its files are absent."""
    rdir = os.path.join(_HERE, "reference", "follow")
    need = ["cuw_vx.profile", "cuw_stress.dat", "cuw_ke.profile", "cuw_params.txt"]
    if not all(os.path.exists(os.path.join(rdir, f)) for f in need):
        return None
    rpar = read_params(os.path.join(rdir, "cuw_params.txt"))
    z, n, vx = read_profile(os.path.join(rdir, "cuw_vx.profile"))
    m = n > 0.5
    z, vx = z[m], vx[m]
    p = read_timeseries(os.path.join(rdir, "cuw_stress.dat"), col=1)
    r = eta_of(rpar, z, vx, p, name="reference/follow/cuw_vx.profile")
    zk, ck, vals = read_profile_multi(os.path.join(rdir, "cuw_ke.profile"), 5)
    zo, co, tb = corrected_T_profile(zk, ck, vals)
    r["t_mid"] = t_mid_of(zo, co, tb, 0.5 * (rpar["zface_bot"] + rpar["zface_top"]))
    r.update({"par": rpar, "z": z, "vx": vx})
    return r


def plot(z, vx, fit, zfb, zft, Tbot, zT, T, t_ps, pxz_series, ref=None,
         out="cuw_viscosity.png"):
    """Three panels, side by side (the day1 layout): the velocity profile with
    the central fit gives the shear rate; the peculiar T(z) shows whether
    viscous heating has humped the channel; the p_xz trace with its running
    mean shows how far the stress average has settled. ref = a
    load_reference() dict -> its converged profile underlies the measured
    points and its mean stress is the dashed line in the right panel."""
    try:
        import matplotlib
        if not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")     # headless: save the PNG, no window
        import matplotlib.pyplot as plt
    except ImportError:
        print("    (matplotlib not found - skipping the plot)")
        return
    BLUE, RED, GREEN, GREY = "#1F3A5F", "#B23A2E", "#4A7C59", "#888888"
    fig, (axL, axM, axR) = plt.subplots(1, 3, figsize=(10.5, 3.6))

    # left: vx(z), z vertical (side-on like the channel), central fit in red
    if ref is not None:
        t_ref = ref["par"]["nprod"] * ref["par"]["dt"]
        axL.plot(ref["vx"], ref["z"], color=GREEN, lw=1.3,
                 label=f"shipped reference ({t_ref:.0f} ps)")
    axL.scatter(vx, z, s=14, color=BLUE, alpha=0.75, label=r"measured $v_x(z)$")
    if fit is not None:                        # the central fit only if this run fitted
        zl = np.linspace(zfb, zft, 50)
        axL.plot(fit["s"] * zl + fit["c"], zl, color=RED, lw=1.4, label="central fit")
    axL.set_xlabel(r"$v_x(z)$ ($\mathrm{\AA}$/ps)")
    axL.set_ylabel(r"$z$ ($\mathrm{\AA}$)")
    axL.set_title(r"velocity profile $v_x(z)$")
    axL.legend(fontsize=7, frameon=False, loc="upper left")

    # centre: peculiar T(z). Flat at a gentle shear; shear hard and viscous
    # heating humps it (hot centre) - the channel is no longer isothermal.
    axM.scatter(zT, T, s=14, color=RED, alpha=0.85)
    axM.axvline(zfb, color="#aaaaaa", ls=":", lw=0.8)
    axM.axvline(zft, color="#aaaaaa", ls=":", lw=0.8)
    axM.axhline(Tbot, color=GREY, ls="--", lw=0.8)
    axM.set_xlabel(r"$z$ ($\mathrm{\AA}$)")
    axM.set_ylabel(r"fluid temperature $T(z)$ (K)")
    axM.set_title(r"peculiar temperature $T(z)$")

    # right: instantaneous p_xz (noisy) vs its running mean. The running mean
    # should flatten once eta has converged; if it is still drifting, it has not.
    t = np.linspace(t_ps / len(pxz_series), t_ps, len(pxz_series))
    run = np.cumsum(pxz_series) / np.arange(1, len(pxz_series) + 1)
    axR.plot(t, pxz_series, color=GREY, lw=0.7, alpha=0.7, label="instantaneous")
    axR.plot(t, run, color=RED, lw=1.6, label=r"running mean $%+.0f$ bar" % run[-1])
    if ref is not None:
        axR.axhline(ref["pxz"], color=GREEN, ls="--", lw=1.2,
                    label=r"reference mean $%+.0f$ bar" % ref["pxz"])
    axR.set_xlabel(r"$t$ (ps)")
    axR.set_ylabel(r"$p_{xz}$ (bar)")
    eta_txt = (r"$\eta \approx %.2f$ mPa$\cdot$s" % fit["eta"]) if fit is not None else "reference shown"
    axR.set_title(r"shear stress $p_{xz}$  (%s)" % eta_txt)
    axR.legend(fontsize=7, frameon=True, facecolor="white", framealpha=0.9,
               edgecolor="none", loc="upper right")

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
    for f in ("cuw_vx.profile", "cuw_stress.dat", "cuw_ke.profile", "cuw_params.txt"):
        if not os.path.exists(f):
            sys.exit(f"{f} not found. If you submitted with ../submit.sh, wait for the "
                     "job (squeue --me) to finish; otherwise run `lmp_serial -in "
                     "viscosity.in` first.")
    par = read_params("cuw_params.txt")

    # Reference FIRST: loaded before ANY this-run value (a params key OR a file read), so a
    # corrupt/incompatible this-run params or a truncated/incomplete output file still shows
    # the converged answer. A corrupt shipped reference degrades to no overlay, not an abort.
    try:
        ref = load_reference()
    except (SystemExit, OSError, KeyError, ValueError):
        ref = None
    t_ps = par.get("nprod", 0) * par.get("dt", 0)

    print("\nStretch sheet 3: viscosity")

    # ---- this run's own short, noisy measurement; any failure - a degenerate fit,
    #      no drive, an incompatible params, OR an incomplete output file - degrades
    #      to the reference-only summary, never an abort. ----
    z = vx = pxz_series = zo = tbin = fit = None
    vwall = zfb = zft = None
    ns_slope = 0.0
    t_heated = False
    try:
        vwall, zfb, zft = par["vwall"], par["zface_bot"], par["zface_top"]
        h = zft - zfb
        ns_slope = 2 * vwall / h             # slope if the water stuck to the walls
        z, n, vx = read_profile("cuw_vx.profile")
        m = n > 0.5                          # keep only bins that hold atoms
        z, vx = z[m], vx[m]
        pxz_series = read_timeseries("cuw_stress.dat", col=1)
        # peculiar T(z) (corrected thermometer) + its mid-channel mean
        zk, ck, vals = read_profile_multi("cuw_ke.profile", 5)
        zo, co, tbin = corrected_T_profile(zk, ck, vals)
        t_mid = t_mid_of(zo, co, tbin, 0.5 * (zfb + zft))
        t_heated = t_mid - par["Tbot"] > 10.0
        if vwall == 0.0:
            raise SystemExit("this run set vwall = 0, so there is no shear rate to divide the "
                             "stress by;\n      the shipped reference is shown below. Rerun "
                             "viscosity.in with its default vwall (0.5 A/ps).")
        fit = eta_of(par, z, vx, pxz_series)
        v_face_bot = fit["s"] * zfb + fit["c"]   # fluid velocity extrapolated to each wall face
        v_face_top = fit["s"] * zft + fit["c"]
        print(f"    central shear rate   dvx/dz = {fit['s']:+.5f} (A/ps)/A   "
              f"(R^2 = {fit['r2']:.3f})")

        # Sanity of the short run - all NON-FATAL (the reference below is always shown).
        tol = 1.15
        flat = fit["s"] < 0.25 * ns_slope
        # a gross overshoot (not from viscous heating) means the fit is inconsistent with the
        # wall speed - eta is unusable, so treat it like flat: suppress the number, defer to ref.
        severe = (not flat) and (not t_heated) and (abs(fit["s"]) > 2.0 * ns_slope
                 or max(abs(v_face_bot), abs(v_face_top)) > 2.0 * vwall)
        if flat:
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            print("    NOTE: the central slope is far below the no-slip line (flat or reversed);")
            print("      on a short run this is usually noise. eta = |p_xz|/|dvx/dz| is unreliable")
            print("      here - read it as indicative only and rely on the shipped reference below.")
        elif severe:
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            print("    NOTE: the fitted slope/wall-face speed is grossly inconsistent with the")
            print("      imposed wall speed - this short-run eta is unusable. Rely on the reference.")
        elif (not t_heated) and (abs(fit["s"]) > tol * ns_slope
                                 or max(abs(v_face_bot), abs(v_face_top)) > tol * vwall):
            # noisy-but-driving (skip once viscous heating sets in - its own note fires below)
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            print("    NOTE: the short run's shear rate is noisy (the wetting copper pins the first")
            print("      water layer, so a straight-line fit over-reads the near-wall gradient). The")
            print("      walls are driving correctly; this just makes the single-run eta below")
            print("      noisier - rely on the shipped reference for the converged value.")

        print(f"    mean shear stress    p_xz = {fit['pxz']:+.1f} +/- {fit['se_p']:.1f} bar "
              f"over the {t_ps:.0f} ps production window")
        if flat or severe:
            # eta = |p_xz|/|s| is meaningless when the shear is ~0/reversed or the fit is
            # inconsistent with the wall; the floored divisor would print an artefact - suppress it.
            print("    viscosity            not meaningfully determined on this run (the fit is")
            print("      flat/reversed or inconsistent with the wall - see the NOTE); use the")
            print("      shipped reference below.")
        else:
            print(f"    viscosity            eta = |p_xz| / (dvx/dz) = "
                  f"{fit['eta']:.2f} +/- {fit['se_eta']:.2f} mPa*s")
            print("    (single-run values move with the seed and the run length - "
                  + ("compare the" if ref is not None else "watch the running"))
            print("     shipped reference below and the running mean in the figure.)"
                  if ref is not None else "     mean in the figure.)")
        print(f"    water temperature    {t_mid:.1f} K mid-channel over the {t_ps:.0f} ps "
              f"production window  (wall baths at {par['Tbot']:.0f} K)")
        if t_heated:
            print("      -> more than 10 K above the bath: viscous heating has set in - the")
            print("         channel is no longer isothermal, so a single eta no longer")
            print("         describes it (read T(z) in the figure).")
    except (SystemExit, ZeroDivisionError, FloatingPointError, ValueError, KeyError, OSError) as e:
        print(f"    {e}")

    # shipped reference tier - ALWAYS shown, whatever this run did above
    overlay = ref
    if ref is not None:
        t_ref = ref["par"]["nprod"] * ref["par"]["dt"]
        print(f"    shipped reference ({t_ref:.0f} ps of averaging; this run: {t_ps:.0f} ps):")
        print(f"      eta = {ref['eta']:.2f} +/- {ref['se_eta']:.2f} mPa*s at "
              f"vwall {ref['par']['vwall']:g} A/ps  (R^2 = {ref['r2']:.3f}, "
              f"mid-channel {ref['t_mid']:.1f} K)")
        if vwall is not None and not (np.isclose(par.get("eps_sl"), ref["par"]["eps_sl"]) and
                                      np.isclose(vwall, ref["par"]["vwall"])):
            print(f"      (no shipped reference at eps_sl = {par.get('eps_sl'):g} eV, "
                  f"vwall = {vwall:g} A/ps - the figure shows this run alone)")
            overlay = None
    if z is not None and zo is not None and zfb is not None:
        plot(z, vx, fit, zfb, zft, par.get("Tbot", 300.0), zo, tbin, t_ps, pxz_series, ref=overlay)
    else:
        print("    (no figure - this run's output was incomplete; the reference above is the value)")


if __name__ == "__main__":
    main()
