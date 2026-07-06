# build/ — construction of data.cuw_channel

`build.in` constructs the Cu(100)/water channel and runs the staged equilibration that
writes `data.cuw_channel` into this folder. The validated copy ships one directory up as
`stretch/data.cuw_channel`; students read that file and do not run `build.in`.

Run from this folder:

```
lmp_serial -in build.in
```

Probe run (minimise + 1000 steps at 2 fs, ~25 s, temporary data file):

```
lmp_serial -in build.in -var nrelax 1000 -var nequil 0 -var nprof 0 -var dataout probe_data.tmp
```

## Geometry

- Walls: 6x6 lateral FCC(100) unit cells at a = 3.615 A (Lx = Ly = 21.69 A); 6 (100)
  planes per wall (72 atoms each, 864 Cu total). The outer 2 planes of each wall are the
  rigid anchor, the inner 4 the Langevin bath. Bottom wall is atom type 3, top wall type 4.
- Water gap, wall surface plane to wall surface plane: 19 plane spacings = 34.34 A
  (the ~35 A target snapped to the Cu lattice).
- Water fill: TIP4P/2005 molecules (`h2o.mol`) inserted with random orientations on an
  FCC seeding lattice (a = 4.338 A = Lx/5, 700 sites), then trimmed to exactly `nwater`
  molecules by an even spread over molecule IDs. The deck errors out if the trimmed count
  differs from `nwater` or the wall count differs from 864.

## Stages

1. Minimise (SHAKE off; the harmonic bond/angle constants hold the water geometry).
2. `nrelax` steps at 2 fs: Langevin baths on the wall bath layers plus a plain Langevin
   on the water (relaxation aid, unfixed afterwards).
3. `nequil` steps at 2 fs: bath-layer Langevin only, water NVE + SHAKE.
4. `nprof` steps at 2 fs: rho(z) accumulation into `cuw_build_density.profile`
   (O-atom number density n_O, raw 0.5 A bins; rho[g/cm3] = 29.915 * n_O[A^-3]),
   then `write_data ${dataout}`.

## Knobs (override with -var NAME value)

| knob | default | meaning |
|---|---|---|
| nwater | 520 | water molecule count (density calibration sets this) |
| eps_sl | 0.0256 | O-Cu LJ epsilon (eV) |
| Twall | 300.0 | bath temperature (K) |
| nrelax | 25000 | stage-2 steps (2 fs): bath + water Langevin |
| nequil | 75000 | stage-3 steps (2 fs): bath Langevin only |
| nprof | 10000 | stage-4 steps (2 fs): rho(z) accumulation; multiple of 10 |
| seed | 90187 | RNG seed |
| dataout | data.cuw_channel | write_data target |

## Notes

- Atom IDs in the written data file have gaps: the trim deletes whole molecules, and an
  ID reset is not applied because the tip4p pair/kspace styles locate each O's two H
  atoms at O-ID+1 and O-ID+2 — reordering breaks that. `velocity ... create` on this
  system therefore uses `loop geom` (the default `loop all` needs consecutive IDs).
- The EAM file sets the Cu mass (63.55) at `pair_coeff` time, overriding the `mass`
  commands for types 3 and 4.
- Stage 4 resets the timestep to 0 before defining its `ave/chunk` fix: the fix writes
  only at timestep multiples of its window (`nprof`), counted from step 0, so without
  the reset any `nrelax`/`nequil` override whose sum is not a multiple of `nprof`
  produces an empty profile file.

## Density calibration series

Gate on the shipped configuration: central-plateau rho = 0.997 +/- 0.015 g/cm3 at 300 K
(1.5 A moving average of the raw profile, plateau = z_mid +/- 7.5 A). Candidate runs are
recorded here.

| nwater | equil length (ps) | plateau rho (g/cm3) | in gate? |
|---|---|---|---|
| 540 | 40 (candidate) | 1.0356 | no (+0.039) |
| 520 | 40 (candidate) | 1.0036 | yes (+0.007) |
| 520 | 200 (final, shipped) | 1.0047 | yes (+0.008) |

Candidate runs: 25 ps stage-2 + 15 ps stage-3 at 2 fs, then the 15 ps rho(z) accumulation
(`-var nrelax 12500 -var nequil 7500 -var nprof 7500`). The proportional update
nwater -> round(nwater * 0.997 / rho_plateau) converged in one step. The shipped
configuration is the final 200 ps equilibration at nwater 520 (last row; default knobs);
`python3 calibrate.py` reports the gate from `cuw_build_density.profile`.

## Shipping validation

`validate.in` re-reads the written data file and runs a 20 ps NVE + SHAKE window at 1 fs
with no thermostats; `validate.py` scores the shipping gates on its outputs and exits 1
on any failure:

```
lmp_serial -in validate.in
python3 validate.py
```

| gate | measured (shipped file) | limit |
|---|---|---|
| total-energy drift (fit x window) | -0.003 meV/atom over 20 ps | <= 2 meV/atom |
| min smoothed rho, central 60 % of gap | 0.856 g/cm3 | >= 0.5 x plateau (0.502) |
| low-density span beyond 6 A from the faces | 0.0 A | < 2 A of rho < 0.5 g/cm3 |
| bath-layer Cu MSD | 0.059 A2 | < 0.5 A2 |

Window-mean water temperature under NVE: 299.9 K. The MSD gate is scored on the bath
layers only: the anchor planes are pinned and would dilute an all-Cu average. Wall-face
positions for the vacuum test come from the innermost wall atom on each side
(`bound()` at the start of the window; written to `cuw_validate_params.txt`).
The copy shipped as `../data.cuw_channel` is the validated file, byte-identical.
