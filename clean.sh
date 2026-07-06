#!/bin/bash
# clean.sh -- remove all regeneratable Cu/water-channel output files so you can
# start a fresh run.
#
# Removes, from this top-level folder and the four measurement folders:
#   cuw_*             simulation outputs + plots (profiles, .dat, params,
#                     trajectories, cuw_*.png)
#   log.lammps        LAMMPS run logs
#   slurm-*.out       Cirrus/SLURM job logs
#   __pycache__/      python bytecode caches
#
# It NEVER touches source (*.in / *.py / *.lmp / *.sh), the shipped inputs
# (data.cuw_channel, Cu.lammps.eam), or the reference/ folders.
#
#   ./clean.sh        remove the regeneratable files
#   ./clean.sh -n     dry run: list what WOULD be removed (deletes nothing)
#
set -u
cd "$(dirname "$0")" || exit 1

dry=0
[ "${1:-}" = "-n" ] && dry=1

for d in . density slip viscosity conductance; do
  [ -d "$d" ] || continue
  if [ "$dry" = 1 ]; then
    find "$d" -maxdepth 1 -type f \( -name 'cuw_*' -o -name 'log.lammps' -o -name 'slurm-*.out' -o -name '.DS_Store' \) -print
    find "$d" -maxdepth 1 -type d -name '__pycache__' -print
  else
    find "$d" -maxdepth 1 -type f \( -name 'cuw_*' -o -name 'log.lammps' -o -name 'slurm-*.out' -o -name '.DS_Store' \) -delete
    find "$d" -maxdepth 1 -type d -name '__pycache__' -exec rm -rf {} +
  fi
done

if [ "$dry" = 1 ]; then
  echo "(dry run -- nothing deleted)"
else
  echo "cleaned: generated Cu/water-channel files removed (source files kept)."
fi
