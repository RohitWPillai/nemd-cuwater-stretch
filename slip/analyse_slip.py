#!/usr/bin/env python3
"""Day 1 stretch - slip length at the Cu/water interface.

Reads the Couette velocity profile written by slip.in and reports the slip
length b: how far past the wall face the central linear profile must
extrapolate to reach the wall speed. b=0 -> the water sticks (no-slip);
b>0 -> it slips. The two walls are symmetric (equal T and eps_sl), so a
single b describes both: the profile runs from -vwall to +vwall over
(h + 2b), so b = vwall/|s| - h/2, with h the face-to-face channel width.
The shipped reference/ folder holds 100 ps runs at both eps_sl values;
the summary quotes their converged b, and the figure overlays the
reference profile that matches this run's eps_sl and vwall.
Saves vx(z) with the extrapolated fit to cuw_slip.png.
Run after `lmp_serial -in slip.in`.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from lammps_io import read_profile, read_profile_multi, read_params

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


def corrected_T_mid(zk, ck, vals, zmid):
    """Mid-channel water temperature from the per-bin kinetic ingredients in
    cuw_ke.profile: each bin's mass-weighted streaming velocity is subtracted
    and the kinetic energy is divided over 2 DOF per atom (6 per SHAKE-rigid
    molecule, not 9), then the bins with |z - zmid| <= 7.5 A are count-weighted
    into one number."""
    msq, mvx, mvy, mvz, mm = (vals[:, i] for i in range(5))
    occ = ck > 0.5                            # a bin visited by < 1 atom carries no temperature
    pec = msq[occ] - (mvx[occ] ** 2 + mvy[occ] ** 2 + mvz[occ] ** 2) / mm[occ]
    tbin = MVV2E * pec / (2.0 * KB)
    zo, co = zk[occ], ck[occ]
    mid = np.abs(zo - zmid) <= 7.5
    return float(np.sum(co[mid] * tbin[mid]) / np.sum(co[mid]))


def central_fit(z, vx, name="cuw_vx.profile"):
    """Unweighted OLS line through the central half of the occupied span
    (away from the structured near-wall layers) -> slope, intercept, R^2,
    SE(slope). The same fit scores this run and the shipped reference."""
    if len(z) < 3:                            # empty/near-empty profile: guard min()/max() below
        raise SystemExit(f"slip fit: fewer than 3 occupied bins in {name} - check the run "
                         "completed and wrote a filled channel.")
    zc, W = 0.5 * (z.min() + z.max()), z.max() - z.min()
    cen = (z > zc - 0.25 * W) & (z < zc + 0.25 * W)
    if np.count_nonzero(cen) < 3:
        raise SystemExit(f"slip fit: fewer than 3 occupied bins in the central half of "
                         f"{name} - check the run completed and the binning/geometry.")
    zf, vf = z[cen], vx[cen]
    s, c = np.polyfit(zf, vf, 1)
    resid = vf - (s * zf + c)
    ss_tot = np.sum((vf - vf.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else float("nan")
    se_s = np.sqrt(np.sum(resid ** 2) / (len(zf) - 2) / np.sum((zf - zf.mean()) ** 2))
    return s, c, r2, se_s


def load_reference(tier):
    """Read reference/<tier>/ (a shipped 100 ps run; tier = "follow" or
    "push") and fit it exactly like this run. Returns the fit, b and SE(b)
    in a dict, or None when the folder or its files are absent."""
    rdir = os.path.join(_HERE, "reference", tier)
    vxf, parf = os.path.join(rdir, "cuw_vx.profile"), os.path.join(rdir, "cuw_params.txt")
    if not (os.path.exists(vxf) and os.path.exists(parf)):
        return None
    rpar = read_params(parf)
    z, n, vx = read_profile(vxf)
    m = n > 0.5
    z, vx = z[m], vx[m]
    s, c, r2, se_s = central_fit(z, vx, name=f"reference/{tier}/cuw_vx.profile")
    h = rpar["zface_top"] - rpar["zface_bot"]
    b = rpar["vwall"] / abs(s) - h / 2
    se_b = rpar["vwall"] * se_s / s ** 2
    return {"par": rpar, "z": z, "vx": vx, "s": s, "c": c, "r2": r2, "b": b, "se_b": se_b}


def plot(z, vx, s, c, zfb, zft, vwall, b, ref=None, this_ps=None, out="cuw_slip.png"):
    """Save the velocity profile vx(z) (z vertical, side-on like the channel)
    with the central Couette fit extrapolated to the wall faces, so the slip
    length b reads off as the offset between the fit and the wall speed.
    ref = a load_reference() dict -> its converged profile is drawn under
    the measured points (the shipped 100 ps run vs this one)."""
    try:
        import matplotlib
        if not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")     # headless: save the PNG, no window
        import matplotlib.pyplot as plt
    except ImportError:
        print("    (matplotlib not found - skipping the plot)")
        return
    BLUE, RED, GREEN = "#1F3A5F", "#B23A2E", "#4A7C59"
    fig, ax = plt.subplots(figsize=(4.0, 4.2))
    bb = b if (b is not None and np.isfinite(b)) else 0.0    # figure span when the fit is missing
    zl = np.linspace(min(zfb - 1.0, zfb - bb - 0.5), max(zft + 1.0, zft + bb + 0.5), 80)
    ax.plot(2 * vwall / (zft - zfb) * (zl - 0.5 * (zfb + zft)), zl,
            color="#999999", ls="--", lw=1.0, label=r"no-slip reference ($b=0$)")
    if ref is not None:
        t_ref = ref["par"]["nprod"] * ref["par"]["dt"]
        ax.plot(ref["vx"], ref["z"], color=GREEN, lw=1.3,
                label=f"shipped reference ({t_ref:.0f} ps)")
    meas_lbl = (r"measured $v_x(z)$ (%.0f ps)" % this_ps) if this_ps else r"measured $v_x(z)$"
    ax.scatter(vx, z, s=14, color=BLUE, alpha=0.75, label=meas_lbl)
    if s is not None:                          # the central fit + b arrow only if this run fitted
        ax.plot(s * zl + c, zl, color=RED, lw=1.4, label="central fit, extrapolated")
    ax.axhline(zfb, color="#888888", ls=":", lw=0.9)
    ax.axhline(zft, color="#888888", ls=":", lw=0.9)
    ax.axvline(-vwall, color="#bbbbbb", ls="--", lw=0.7)
    ax.axvline(vwall, color="#bbbbbb", ls="--", lw=0.7)
    # mark b: the distance past the top wall face where the fit reaches +vwall
    if s is not None and b is not None and np.isfinite(b) and abs(s) > 1e-12 and abs(b) < (zft - zfb):
        z_hit = (vwall - c) / s
        ax.annotate("", xy=(vwall, zft), xytext=(vwall, z_hit),
                    arrowprops=dict(arrowstyle="<->", color=RED, lw=1.0))
        x_line = s * 0.5 * (zft + z_hit) + c      # fit line at the label height
        ax.text(min(vwall, x_line) - 0.18 * vwall, 0.5 * (zft + z_hit),
                r"$b \approx %.1f\ \mathrm{\AA}$" % b,
                color=RED, ha="right", va="center", fontsize=8)
    ax.set_xlabel(r"$v_x(z)$ ($\mathrm{\AA}$/ps)")
    ax.set_ylabel(r"$z$ ($\mathrm{\AA}$)")
    title_b = (r"$b = %+.2f\ \mathrm{\AA}$" % b) if (b is not None and np.isfinite(b)) else "reference shown"
    ax.set_title(r"Couette velocity profile  (%s)" % title_b)
    # lower right: the only corner the fit line, the b arrow (top, for b > 0)
    # and the measured points never reach
    ax.legend(fontsize=8, frameon=False, loc="lower right")
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
    for f in ("cuw_vx.profile", "cuw_ke.profile", "cuw_params.txt"):
        if not os.path.exists(f):
            sys.exit(f"{f} not found. If you submitted with ../submit.sh, wait for the "
                     "job (squeue --me) to finish; otherwise run `lmp_serial -in slip.in` first.")
    par = read_params("cuw_params.txt")

    # Reference FIRST: loaded before ANY this-run value (a params key OR a file read), so a
    # corrupt/incompatible this-run params or a truncated/incomplete output file still shows
    # the converged answer. A corrupt shipped reference degrades to no overlay, not an abort.
    try:
        ref = {t: load_reference(t) for t in ("follow", "push")}
    except (SystemExit, OSError, KeyError, ValueError):
        ref = {"follow": None, "push": None}
    tags = {"follow": "the wetting default", "push": "the Push value"}
    have = [t for t in ("follow", "push") if ref[t]]
    t_ps = par.get("nprod", 0) * par.get("dt", 0)

    print("\nStretch sheet 2: slip length")

    # ---- this run's own short, noisy measurement. Any failure - a degenerate fit,
    #      no drive, an incompatible params, OR an incomplete output file - degrades
    #      to the reference-only summary below, never an abort. ----
    z = vx = s = c = b = se_b = None
    vwall = zfb = zft = None
    try:
        vwall, zfb, zft = par["vwall"], par["zface_bot"], par["zface_top"]
        h = zft - zfb
        z, n, vx = read_profile("cuw_vx.profile")
        m = n > 0.5                          # keep only bins that hold atoms
        z, vx = z[m], vx[m]
        if vwall == 0.0:
            raise SystemExit("this run set vwall = 0, so there is no Couette slope to fit;\n"
                             "      the shipped reference is shown below. Rerun slip.in with its "
                             "default vwall (0.5 A/ps) for your own value.")
        s, c, r2, se_s = central_fit(z, vx)
        zk, ck, vals = read_profile_multi("cuw_ke.profile", 5)
        t_mid = corrected_T_mid(zk, ck, vals, 0.5 * (zfb + zft))
        ns_slope = 2 * vwall / h                 # slope if the water stuck to the walls
        v_face_bot = s * zfb + c                 # fluid velocity extrapolated to each wall face
        v_face_top = s * zft + c

        print(f"    central shear rate   dvx/dz = {s:+.5f} (A/ps)/A   (R^2 = {r2:.3f})")
        print(f"    no-slip reference    2*vwall/h = {ns_slope:+.5f} (A/ps)/A  (slope if the water")
        print("      stuck to the walls; slip makes the measured slope smaller)")
        print(f"    water temperature    {t_mid:.1f} K mid-channel over the {t_ps:.0f} ps production "
              f"window  (wall baths at {par['Tbot']:.0f} K)")
        if t_mid - par["Tbot"] > 10.0:
            print("      -> more than 10 K above the bath: the shear is viscously heating the water")

        # Sanity of the short run - all NON-FATAL (the reference below is always shown).
        # flat/reversed (s far below the no-slip line) and noisy-but-driving (slope or
        # wall-face overshoot) are both short-run fit noise for the pinned wetting layer,
        # not a setup error; a real problem only shows if it persists in the reference.
        tol = 1.15
        flat = s < 0.25 * ns_slope
        # a gross overshoot (slope or face speed more than ~2x expected) is not "near the
        # wall speed" - the fit is unusable, so treat it like flat: suppress b, defer to ref.
        severe = (not flat) and (abs(s) > 2.0 * ns_slope
                                 or max(abs(v_face_bot), abs(v_face_top)) > 2.0 * vwall)
        slope_hot = (not flat) and (not severe) and abs(s) > tol * ns_slope
        face_hot = (not flat) and (not severe) and max(abs(v_face_bot), abs(v_face_top)) > tol * vwall
        if flat:
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            print("    NOTE: the central slope is far below the no-slip line (flat or reversed). On")
            print("      a short run this is usually short-run/thermal noise, or a different")
            print("      machine's trajectory; it is a real setup problem only if it persists in")
            print("      the shipped reference overlay below. Read b as indicative only.")
        elif severe:
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            print("    NOTE: the fitted slope/wall-face speed is grossly inconsistent with the")
            print("      imposed wall speed - this short-run fit is unusable (check the run and the")
            print("      -var vwall value). Rely on the shipped reference below.")
        elif slope_hot or face_hot:
            print(f"    water velocity at the wall faces = {v_face_bot:+.3f} / {v_face_top:+.3f} "
                  f"A/ps  (wall speed -/+{vwall:g} A/ps)")
            what = []
            if slope_hot:
                what.append(f"the central slope runs {(abs(s) / ns_slope - 1) * 100:.0f}% above the "
                            "no-slip line")
            if face_hot:
                over = max(abs(v_face_bot), abs(v_face_top)) / vwall - 1.0
                what.append(f"the fitted wall-face speed overshoots the wall by {over * 100:.0f}%")
            print("    NOTE: " + "; and ".join(what) + ".")
            print("      The walls are driving correctly (the face speeds are near the wall speed),")
            print("      so this is short-run fit noise, not a setup error: the wetting copper pins")
            print("      the first water layer, and a straight-line fit over-reads the near-wall")
            print("      gradient on a short run. Read b below as indicative only and rely on the")
            print("      shipped reference for the converged value.")

        if flat or severe or abs(s) < 1e-6:
            print("    slip length          not meaningfully determined on this run (the fit is")
            print("      flat/reversed or inconsistent with the wall - see the NOTE); use the")
            print("      shipped reference below.")
        else:
            b = vwall / abs(s) - h / 2           # symmetric slip length (walls are equivalent)
            se_b = vwall * se_s / s ** 2         # SE(b) propagated through b = vwall/|s| - h/2
            print(f"    slip length          b = {b:+.2f} +/- {se_b:.2f} A  (symmetric walls -> "
                  "one b for both),")
            print(f"      measured to the innermost Cu plane on each side (z = {zfb:.2f} / "
                  f"{zft:.2f} A).")
            if b > 5.0:
                print("    -> b > 5 A: the water slides along the copper - expected when eps_sl is")
                print("       lowered (the Push run), not for the wetting default.")
            elif b < -5.0:
                print("    -> b < -5 A: the effective no-slip plane is well inside the wall face -")
                print("       with the wetting default this means the fit picked up the pinned")
                print("       first layer; check the profile figure before reading the sign.")
            else:
                print("    -> |b| <= 5 A: near no-slip - the wetting copper holds the water to")
                print("       within a few A of the stick condition.")
        if par["eps_sl"] != 0.0256:
            print(f"    -> eps_sl = {par['eps_sl']:g} eV, not the wetting 0.0256 eV: compare b with "
                  "your default run -")
            print("       a weaker O-Cu attraction lets the water slide further.")
    except (SystemExit, ZeroDivisionError, FloatingPointError, ValueError, KeyError, OSError) as e:
        print(f"    {e}")

    # shipped reference tier (the two-tier pattern: this short run vs the
    # converged answer) - ALWAYS shown, whatever this run did above
    if have:
        t_ref = ref[have[0]]["par"]["nprod"] * ref[have[0]]["par"]["dt"]
        print(f"    shipped reference ({t_ref:.0f} ps of averaging; this run: {t_ps:.0f} ps):")
        for t in have:
            r = ref[t]
            print(f"      eps_sl {r['par']['eps_sl']:g} eV ({tags[t]}): "
                  f"b = {r['b']:+.2f} +/- {r['se_b']:.2f} A  (R^2 = {r['r2']:.3f})")
        if len(have) == 2:
            sep = ref["push"]["b"] - ref["follow"]["b"]
            se_sep = float(np.hypot(ref["follow"]["se_b"], ref["push"]["se_b"]))
            print(f"      -> weakening the O-Cu attraction moves b by {sep:+.1f} A "
                  f"({sep / se_sep:.0f}x its standard error): wettability controls slip.")
    overlay = None
    for t in have:
        r = ref[t]
        if (vwall is not None and np.isclose(par.get("eps_sl"), r["par"]["eps_sl"]) and
                np.isclose(vwall, r["par"]["vwall"])):
            overlay = r
            break
    if have and overlay is None and vwall is not None:
        print(f"      (no shipped reference at eps_sl = {par.get('eps_sl'):g} eV, "
              f"vwall = {vwall:g} A/ps - the figure shows this run alone)")
    if z is not None and zfb is not None:
        plot(z, vx, s, c, zfb, zft, vwall, b, ref=overlay, this_ps=t_ps)
    else:
        print("    (no figure - this run's output was incomplete; the reference above is the value)")


if __name__ == "__main__":
    main()
