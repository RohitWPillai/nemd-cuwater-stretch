"""Shared output-parsing helpers for the Cu/water channel analysers.

These readers are the shared output parsers for the analysers, kept in
one place so the four
`analyse_<case>.py` scripts parse LAMMPS output identically. The physics
(conductance G, the Couette/slip fit, the viscosity) stays in each case folder.

Import it with the same reach-up convention you already use to run a case
(`cd density && python analyse_density.py`, which `include ../shared_setup.lmp`
mirrors on the LAMMPS side). At the top of each analyser:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from lammps_io import read_profile, read_params, read_timeseries
"""
import getpass
import os
import sys
import tempfile

import numpy as np

# Point matplotlib's cache at a writable dir BEFORE any analyser imports matplotlib,
# so a non-writable $HOME on HPC does not warn and rebuild the font cache (~10 s) each run.
# Per-user path: on a shared login node a fixed /tmp name is owned by whoever ran first,
# and everyone else gets a not-writable warning on every run.
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(tempfile.gettempdir(), f"cuw_mplcache_{getpass.getuser()}"))


def _last_chunk_block(lines, ncol=4):
    """Parse pre-filtered (comment/blank-stripped) `fix ave/chunk` lines and return
    the LAST complete block as an (n_chunks, ncol) array, or None if no complete block
    is present (a trailing block still being written is tolerated and ignored).
    Each block = a header `timestep n_chunks [total]`, then n_chunks rows
    `chunk Coord1 Ncount value [value ...]`. Used by read_profile and
    read_profile_multi, which fail loud on None."""
    last, i = None, 0
    while i < len(lines):
        try:
            n = int(lines[i].split()[1])               # header: timestep n_chunks [total]
        except (ValueError, IndexError):
            break                                      # garbled/partial header -> stop
        if n <= 0 or i + 1 + n > len(lines):
            break                                      # trailing block not fully written
        rows, ok = [], True
        for k in range(n):
            c = lines[i + 1 + k].split()
            if len(c) < ncol:                          # row truncated at a column boundary
                ok = False
                break
            try:
                rows.append([float(x) for x in c[:ncol]])
            except ValueError:                         # a value truncated mid-number
                ok = False
                break
        if not ok:
            break                                      # incomplete trailing block -> keep last complete
        last = np.array(rows)
        i += 1 + n
    return last


def read_profile(path):
    """LAMMPS `fix ave/chunk` file -> (z, Ncount, value) of the LAST complete block.

    Handles single- or multi-block files. Each block is a 2+-token header
    (timestep n_chunks [total_count]; only n_chunks is used) followed by n_chunks
    rows: chunk Coord1 Ncount value. Returns the spatial coordinate, the per-bin
    count, and the binned quantity (density / v_x / temperature, by case)."""
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#") and l.strip()]
    a = _last_chunk_block(lines)
    if a is None:
        sys.exit(f"{path}: no complete data block - the run may still be writing, was killed "
                 f"mid-write, or is too short. Wait for the job to finish (squeue --me) or "
                 f"rerun after ./clean.sh.")
    return a[:, 1], a[:, 2], a[:, 3]                                    # Coord1, Ncount, value


def read_profile_multi(path, nval):
    """LAMMPS `fix ave/chunk` file carrying nval binned quantities -> (z, Ncount,
    values) of the LAST complete block, with values shaped (n_chunks, nval).

    Same block structure as read_profile (header + n_chunks rows of
    `chunk Coord1 Ncount value ...`), for fixes that average several quantities
    at once (e.g. the per-bin kinetic ingredients of the water thermometer:
    v_msq v_mvx v_mvy v_mvz v_mm)."""
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#") and l.strip()]
    a = _last_chunk_block(lines, 3 + nval)
    if a is None:
        sys.exit(f"{path}: no complete data block - the run may still be writing, was killed "
                 f"mid-write, or is too short. Wait for the job to finish (squeue --me) or "
                 f"rerun after ./clean.sh.")
    return a[:, 1], a[:, 2], a[:, 3:]


def read_params(path):
    """Read a `key value` text file (cuw_params.txt) -> {key: float}.

    Lines whose second token is non-numeric are skipped, so header/comment
    lines do no harm."""
    p = {}
    with open(path) as f:
        for l in f:
            t = l.split()
            if len(t) >= 2:
                try:
                    p[t[0]] = float(t[1])
                except ValueError:
                    pass
    return p


def read_timeseries(path, col=None):
    """LAMMPS `fix ave/time` log.

    col=None -> (col0, col1, col2) arrays (timestep + first two quantities),
                exiting if there are < 2 data rows (a slope needs >= 2 points).
    col=int  -> just that column as a 1-D array (col 1 = the first quantity).

    Short lines, a final line left un-terminated (a run still writing / killed
    mid-write), and tokens cut mid-number are skipped rather than raising."""
    need = 3 if col is None else col + 1                # columns this branch reads
    with open(path) as f:
        text = f.read()
    file_lines = text.split("\n")
    if text and not text.endswith("\n"):               # final line written without a newline =
        file_lines = file_lines[:-1]                   # a run still writing / killed mid-write -> drop it
    rows = []
    for l in file_lines:
        if l.startswith("#") or not l.strip():
            continue
        t = l.split()
        if len(t) < need:                              # short/truncated line
            continue
        try:
            rows.append([float(x) for x in t[:need]])
        except ValueError:                             # a value cut mid-number (e.g. '2.7e') -> skip
            continue
    a = np.array(rows)
    if a.ndim != 2 or len(a) < 1:
        sys.exit(f"{path}: no complete data rows - the run may still be writing, was killed "
                 f"mid-write, or is too short. Wait for the job to finish (squeue --me) or "
                 f"rerun after ./clean.sh.")
    if col is not None:
        return a[:, col]
    if len(a) < 2:
        sys.exit(f"{path}: too few data rows for a slope - run longer.")
    return a[:, 0], a[:, 1], a[:, 2]
