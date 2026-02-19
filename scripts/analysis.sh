#!/bin/bash
set -Eeuo pipefail

declare -A LiPF6inACN=(
  [cutoffin]="5.3" ## inner cutoff for domain, A
  [cutoffout]="9.1" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="0" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)
declare -A LiPF6inH2O=(
  [cutoffin]="5.3" ## inner cutoff for domain, A
  [cutoffout]="9.1" ## outer cutoff for domain, A
  [boxX]="28.57" ## box size A
  [boxY]="28.57" ## box size A
  [boxZ]="28.57" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="0" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="3 1 7"
)

measures=(preprocess processTraj rate conductivity ionionprob mechanism)
echo "Choose a measurement : "
select measure in "${measures[@]}"; do
  [[ -n "${measure:-}" ]] && break
  echo "Invalid choice. Try again."
done


targets=(LiPF6inACN LiPF6inH2O)

load_params() {
  local tgt=$1
  declare -n P="${tgt}"
  if [[ ${#P[@]} -eq 0 ]]; then
    echo "ERROR: unknown target ' ${tgt}. Known: ${targets[*]}" >&2
    return 1
  fi

  for k in "${!P[@]}"; do
    printf -v "$k" "%s" "${P[$k]}"
  done
}

echo "Choose a target :"
select target in "${targets[@]}"; do
  [[ -n "${targets:-}" ]] && break
  echo "Invalid choice. Try again."
done

echo "Selected target trajectory : $target"
load_params "$target"

density="05"

fields=( "00" "09" "19")

echo "Choose a field to analyse :"
select field in "${fields[@]}"; do
  [[ -n "${fields:-}" ]] && break
  echo "Invalid choice. Try again."
done

if [[ $measure == "preprocess" ]]; then
  echo "LAMMPS dcd to binary trajectory in : ./data/traj/$target/"
  mkdir -p ../data/traj/$target
  ../bin/preprocess $target $density $field 

elif [[ $measure == "processTraj" ]]; then
  input="../data/traj/$target/"
  echo "Extract cluster-based nearest-counterion distance trajectories from : $input"

  mkdir -p ../data/cnnDist
  mkdir -p ../data/cnnAngle
  mkdir -p ../data/cnnId
  mkdir -p ../data/cnnCid
  ../bin/processTraj $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime

elif [[ $measure == "rate" ]]; then
  input="../data/cnnDist/$target"
  echo "Compute TPT variables and rate from : $target"
  mkdir -p ../results/rate
  ../bin/rate $target $density $field $cutoffin $cutoffout $timestep

elif [[ $measure == "conductivity" ]]; then
  input="../data/traj/$target"
  echo "Measure conductivity from trajectory file : $input"
  mkdir -p ../results/conductivity
  ../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime

elif [[ $measure == "ionionprob" ]]; then
  input="../data/cnnDist/$target"
  echo "Compute CBNC distance probability distributione : $input"
  mkdir -p ../results/ionionprob
  ../bin/ionionprob $target $density $field 

elif [[ $measure == "mechanism" ]]; then
  echo "Analysis on reactive trajectories using all cnn data"
  mkdir -p ../results/mechanism
  ../bin/mechanism $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime
fi
