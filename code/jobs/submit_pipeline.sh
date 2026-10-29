#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$root"
target=vwh-all
force=
dry_run=false
while (($#)); do
  case "$1" in
    --target|--force)
      (($# >= 2)) || { echo "Missing value for $1" >&2; exit 1; }
      if [[ "$1" == --target ]]; then target=$2; else force=$2; fi
      shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done
case "$target" in
  vwh-all|conversion|panels|panel-reports|evals|eval-reports|decomposition|decomp-reports) ;;
  private-figure|prepare-figure|figure-1a) ;;
  *) echo "Unknown target: $target" >&2; exit 1 ;;
esac
export DATA_PATH="${DATA_PATH:-$root/data}"
export OUTPUT_PATH="${OUTPUT_PATH:-$root/output}"
export REPLICATION_LIBRARY="${REPLICATION_LIBRARY:-$root/.replication/library}"
# Planning does not install packages, execute science, or contact Slurm.
plan=$(make --no-print-directory -n -s PLAN=1 FORCE="$force" "$target" | sed -n 's/^: //p')
if $dry_run; then printf '%s\n' "$plan"; exit 0; fi
source code/jobs/init_server.sh "${APM_MAX_THREADS:-1}"
"${R:-Rscript}" --vanilla setup.R
[[ -n "$plan" ]] || { echo 'All requested stages are complete.'; exit 0; }
logdir="$OUTPUT_PATH/.replication/logs"
mkdir -p "$logdir"
record="$OUTPUT_PATH/.replication/jobs-$(date -u +%Y%m%dT%H%M%SZ)-$$.tsv"
declare -A ids=()
while read -r stage dependencies; do
  [[ -n "$stage" ]] || continue
  # Public rendering has no private dependency. A queued private preparation
  # supplies that dependency only in this explicit private submission.
  if [[ "$stage" == figure-1a && -n "${ids[prepare-figure]:-}" ]]; then dependencies=prepare-figure; fi
  dep_ids=()
  for dep in $dependencies; do
    [[ -z "${ids[$dep]:-}" ]] || dep_ids+=("${ids[$dep]}")
  done
  args=(--parsable "--chdir=$root" "--output=$logdir/$stage-%j.out" "--error=$logdir/$stage-%j.err")
  if ((${#dep_ids[@]})); then
    joined=$(IFS=:; echo "${dep_ids[*]}")
    args+=("--dependency=afterok:$joined")
  fi
  mode=resume
  if [[ ",$force," == *,all,* || ",$force," == *,"$stage",* ]]; then mode=fresh; fi
  reply=$(sbatch "${args[@]}" "code/slurm/${stage//-/_}.slurm" "$root" "$mode")
  [[ "$reply" =~ ^[0-9]+(\;[^[:space:]]+)?$ ]] || { echo "Invalid sbatch response: $reply" >&2; exit 1; }
  ids[$stage]=${reply%%;*}
  printf '%s\t%s\n' "$stage" "${ids[$stage]}" | tee -a "$record"
done <<< "$plan"
echo "Submitted job list: $record"
