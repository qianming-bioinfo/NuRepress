#!/usr/bin/env bash
set -Eeuo pipefail

print_help() {
  cat <<'EOF'
Usage:
  bash call_nucleosomes.sh -i INPUTS -o OUT_DIR [options] -- [DANPOS options...]

Required:
  -i, --input        Input data in alias:path format, comma-separated.
                     Example:
                     "sample1:/path/sample1.bam,sample2:/path/sample2.bam"

  -o, --outdir       Output directory

Optional:
  -b, --background   Background/input data in alias:path format, comma-separated.
                     Example:
                     "sample1:/path/sample1_input.bam,sample2:/path/sample2_input.bam"

                     Behavior:
                     1) If one background entry is provided, it will be used as the
                        single background for all experimental datasets.
                     2) If multiple background entries are provided, aliases must match
                        the input aliases used in the run, and the wrapper will convert
                        them to DANPOS-native exp:bg mapping.

  --compare-pairs    Treatment comparison pairs in alias:alias format, comma-separated.
                     Examples:
                     "sample1:sample2"
                     "koA:wt,koB:wt"

                     Notes:
                     1) The order is preserved exactly as given to DANPOS.
                     2) If omitted, all inputs are passed to DANPOS as independent
                        datasets separated by commas.
                     3) If provided, the wrapper converts them to DANPOS-native
                        path:path pairs.

  -x, --command      DANPOS subcommand
                     Default: dpos

  -p, --danpos       DANPOS executable name or full path
                     Default: danpos.py
                     If a plain command name is given, it will be resolved from PATH.

  -t, --threads      Threads for samtools index
                     Default: value of SLURM_CPUS_PER_TASK or 1

  -c, --config       Config file with one key=value per line
                     Supported keys:
                       input=
                       background=
                       compare_pairs=
                       outdir=
                       command=
                       danpos=
                       threads=
                       skip_index=
                       log_dir=

  --log-dir          Log directory
                     Default: <OUT_DIR>/logs

  --skip-index       Skip BAM quickcheck and BAM index generation

  -h, --help         Show this help message

Pass-through:
  All arguments after -- are passed directly to DANPOS unchanged.

Examples:
  # 1) Run multiple datasets independently
  bash call_nucleosomes.sh \
    -i "sample1:/path/sample1.bam,sample2:/path/sample2.bam" \
    -o /path/to/output_dir \
    -- \
    -m 1 --mifrsz 80 --mafrsz 800 --extend 70

  # 2) Compare one treatment pair
  bash call_nucleosomes.sh \
    -i "sample1:/path/sample1.bam,sample2:/path/sample2.bam" \
    --compare-pairs "sample1:sample2" \
    -o /path/to/output_dir \
    -- \
    -m 1 --mifrsz 80 --mafrsz 800 --extend 70

  # 3) Compare one treatment pair, with matched input/background for each dataset
  bash call_nucleosomes.sh \
    -i "sample1:/path/sample1.bam,sample2:/path/sample2.bam" \
    -b "sample1:/path/sample1_input.bam,sample2:/path/sample2_input.bam" \
    --compare-pairs "sample1:sample2" \
    -o /path/to/output_dir \
    -- \
    -m 1 --mifrsz 80 --mafrsz 800 --extend 70
EOF
}

COMMAND="dpos"
DANPOS_CMD="danpos.py"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

INPUT_SPEC=""
BG_SPEC=""
COMPARE_PAIRS=""
OUT_DIR=""
CONFIG_FILE=""
LOG_DIR=""
SKIP_INDEX=0

DANPOS_ARGS=()
HAS_DANPOS_ARGS=0

RAW_ARGS=("$@")

declare -A INPUT_PATH_BY_ALIAS=()
declare -A BG_PATH_BY_ALIAS=()
declare -A USED_INPUT_ALIAS_SEEN=()

INPUT_ALIAS_ORDER=()
BG_ALIAS_ORDER=()
USED_INPUT_ALIAS_ORDER=()
COMPARE_LEFT_ALIASES=()
COMPARE_RIGHT_ALIASES=()

INPUT_COUNT=0
BG_COUNT=0
USED_INPUT_COUNT=0
COMPARE_PAIR_COUNT=0

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "$s"
}

join_by() {
  local delim="$1"
  shift || true
  local first=1
  local x
  for x in "$@"; do
    if [[ $first -eq 1 ]]; then
      printf '%s' "$x"
      first=0
    else
      printf '%s%s' "$delim" "$x"
    fi
  done
}

quote_cmd() {
  local out=""
  local arg
  for arg in "$@"; do
    out+=" $(printf '%q' "$arg")"
  done
  printf '%s\n' "${out# }"
}

resolve_danpos_cmd() {
  local cmd="$1"

  if [[ "$cmd" == */* ]]; then
    [[ -f "$cmd" ]] || die "DANPOS executable not found: $cmd"
    [[ -x "$cmd" ]] || die "DANPOS executable is not executable: $cmd"
    printf '%s\n' "$cmd"
  else
    command -v "$cmd" >/dev/null 2>&1 || die "DANPOS command not found in PATH: $cmd"
    command -v "$cmd"
  fi
}

load_config() {
  local cfg="$1"
  [[ -f "$cfg" ]] || die "Config file not found: $cfg"

  local line key value
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="$(trim "$line")"
    [[ -z "$line" ]] && continue
    [[ "${line:0:1}" == "#" ]] && continue
    [[ "$line" == *=* ]] || die "Invalid config line: $line"

    key="$(trim "${line%%=*}")"
    value="$(trim "${line#*=}")"

    case "$key" in
      input) INPUT_SPEC="$value" ;;
      background) BG_SPEC="$value" ;;
      compare_pairs) COMPARE_PAIRS="$value" ;;
      outdir) OUT_DIR="$value" ;;
      command) COMMAND="$value" ;;
      danpos) DANPOS_CMD="$value" ;;
      threads) THREADS="$value" ;;
      log_dir) LOG_DIR="$value" ;;
      skip_index)
        case "$value" in
          1|true|TRUE|yes|YES) SKIP_INDEX=1 ;;
          0|false|FALSE|no|NO) SKIP_INDEX=0 ;;
          *) die "Invalid skip_index value in config: $value" ;;
        esac
        ;;
      *) die "Unsupported config key: $key" ;;
    esac
  done < "$cfg"
}

parse_alias_path_list() {
  local spec="$1"
  local kind="$2"

  local -a items=()
  local item alias path

  IFS=',' read -r -a items <<< "$spec"

  for item in "${items[@]}"; do
    item="$(trim "$item")"
    [[ -z "$item" ]] && continue
    [[ "$item" == *:* ]] || die "Each ${kind} entry must be alias:path, got: $item"

    alias="$(trim "${item%%:*}")"
    path="$(trim "${item#*:}")"

    [[ -n "$alias" ]] || die "Empty alias in ${kind} entry: $item"
    [[ -n "$path" ]] || die "Empty path in ${kind} entry: $item"

    if [[ "$kind" == "input" ]]; then
      [[ -z "${INPUT_PATH_BY_ALIAS[$alias]+_}" ]] || die "Duplicate input alias: $alias"
      INPUT_PATH_BY_ALIAS["$alias"]="$path"
      INPUT_ALIAS_ORDER+=("$alias")
      INPUT_COUNT=$((INPUT_COUNT + 1))
    elif [[ "$kind" == "background" ]]; then
      [[ -z "${BG_PATH_BY_ALIAS[$alias]+_}" ]] || die "Duplicate background alias: $alias"
      BG_PATH_BY_ALIAS["$alias"]="$path"
      BG_ALIAS_ORDER+=("$alias")
      BG_COUNT=$((BG_COUNT + 1))
    else
      die "Internal error: unsupported parse kind: $kind"
    fi
  done
}

add_used_input_alias() {
  local alias="$1"
  if [[ -z "${USED_INPUT_ALIAS_SEEN[$alias]+_}" ]]; then
    USED_INPUT_ALIAS_SEEN["$alias"]=1
    USED_INPUT_ALIAS_ORDER+=("$alias")
    USED_INPUT_COUNT=$((USED_INPUT_COUNT + 1))
  fi
}

parse_compare_pairs_and_build_used_aliases() {
  local spec="$1"

  local -a pairs=()
  local pair left right
  local alias

  COMPARE_LEFT_ALIASES=()
  COMPARE_RIGHT_ALIASES=()
  USED_INPUT_ALIAS_ORDER=()
  USED_INPUT_ALIAS_SEEN=()
  USED_INPUT_COUNT=0
  COMPARE_PAIR_COUNT=0

  if [[ -z "$spec" ]]; then
    if (( INPUT_COUNT > 0 )); then
      for alias in "${INPUT_ALIAS_ORDER[@]}"; do
        add_used_input_alias "$alias"
      done
    fi
    return 0
  fi

  IFS=',' read -r -a pairs <<< "$spec"

  for pair in "${pairs[@]}"; do
    pair="$(trim "$pair")"
    [[ -z "$pair" ]] && continue
    [[ "$pair" == *:* ]] || die "Each compare pair must be alias:alias, got: $pair"

    left="$(trim "${pair%%:*}")"
    right="$(trim "${pair#*:}")"

    [[ -n "$left" ]] || die "Empty left alias in compare pair: $pair"
    [[ -n "$right" ]] || die "Empty right alias in compare pair: $pair"

    [[ -n "${INPUT_PATH_BY_ALIAS[$left]+_}" ]] || die "Unknown input alias in compare pair: $left"
    [[ -n "${INPUT_PATH_BY_ALIAS[$right]+_}" ]] || die "Unknown input alias in compare pair: $right"

    COMPARE_LEFT_ALIASES+=("$left")
    COMPARE_RIGHT_ALIASES+=("$right")
    COMPARE_PAIR_COUNT=$((COMPARE_PAIR_COUNT + 1))

    add_used_input_alias "$left"
    add_used_input_alias "$right"
  done

  (( COMPARE_PAIR_COUNT > 0 )) || die "--compare-pairs was provided but no valid pair was parsed"
}

ensure_bai() {
  local bam="$1"

  [[ -f "$bam" ]] || die "BAM not found: $bam"

  if ! samtools quickcheck -v "$bam" >/dev/null 2>&1; then
    samtools quickcheck -v "$bam" || true
    die "BAM failed samtools quickcheck: $bam"
  fi

  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    samtools index -@ "$THREADS" "$bam"
  fi

  [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] || die "Failed to create BAM index for: $bam"
}

ensure_all_bai() {
  local -A seen=()
  local alias bam

  if (( INPUT_COUNT > 0 )); then
    for alias in "${INPUT_ALIAS_ORDER[@]}"; do
      bam="${INPUT_PATH_BY_ALIAS[$alias]}"
      if [[ -z "${seen[$bam]+_}" ]]; then
        ensure_bai "$bam"
        seen["$bam"]=1
      fi
    done
  fi

  if (( BG_COUNT > 0 )); then
    for alias in "${BG_ALIAS_ORDER[@]}"; do
      bam="${BG_PATH_BY_ALIAS[$alias]}"
      if [[ -z "${seen[$bam]+_}" ]]; then
        ensure_bai "$bam"
        seen["$bam"]=1
      fi
    done
  fi
}

build_danpos_input_spec() {
  local -a items=()
  local alias idx left right

  if (( COMPARE_PAIR_COUNT == 0 )); then
    for alias in "${INPUT_ALIAS_ORDER[@]}"; do
      items+=("${INPUT_PATH_BY_ALIAS[$alias]}")
    done
  else
    for ((idx=0; idx<COMPARE_PAIR_COUNT; idx++)); do
      left="${COMPARE_LEFT_ALIASES[$idx]}"
      right="${COMPARE_RIGHT_ALIASES[$idx]}"
      items+=("${INPUT_PATH_BY_ALIAS[$left]}:${INPUT_PATH_BY_ALIAS[$right]}")
    done
  fi

  join_by "," "${items[@]}"
}

build_danpos_bg_spec() {
  local -a items=()
  local alias

  if (( BG_COUNT == 0 )); then
    printf '%s\n' ""
    return 0
  fi

  if (( BG_COUNT == 1 )); then
    printf '%s\n' "${BG_PATH_BY_ALIAS[${BG_ALIAS_ORDER[0]}]}"
    return 0
  fi

  for alias in "${USED_INPUT_ALIAS_ORDER[@]}"; do
    [[ -n "${BG_PATH_BY_ALIAS[$alias]+_}" ]] || die "Missing background for input alias used in run: $alias"
    items+=("${INPUT_PATH_BY_ALIAS[$alias]}:${BG_PATH_BY_ALIAS[$alias]}")
  done

  join_by "," "${items[@]}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      INPUT_SPEC="$2"
      shift 2
      ;;
    -b|--background)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      BG_SPEC="$2"
      shift 2
      ;;
    --compare-pairs|--pairs)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      COMPARE_PAIRS="$2"
      shift 2
      ;;
    -o|--outdir)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      OUT_DIR="$2"
      shift 2
      ;;
    -x|--command)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      COMMAND="$2"
      shift 2
      ;;
    -p|--danpos)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      DANPOS_CMD="$2"
      shift 2
      ;;
    -t|--threads)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      THREADS="$2"
      shift 2
      ;;
    -c|--config)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      CONFIG_FILE="$2"
      shift 2
      ;;
    --log-dir)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      LOG_DIR="$2"
      shift 2
      ;;
    --skip-index)
      SKIP_INDEX=1
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    --)
      shift
      if [[ $# -gt 0 ]]; then
        DANPOS_ARGS=("$@")
        HAS_DANPOS_ARGS=1
      fi
      break
      ;;
    *)
      die "Unknown wrapper option: $1"
      ;;
  esac
done

if [[ -n "$CONFIG_FILE" ]]; then
  load_config "$CONFIG_FILE"
fi

[[ -n "$INPUT_SPEC" ]] || die "-i/--input is required"
[[ -n "$OUT_DIR" ]] || die "-o/--outdir is required"

parse_alias_path_list "$INPUT_SPEC" "input"

if [[ -n "$BG_SPEC" ]]; then
  parse_alias_path_list "$BG_SPEC" "background"
fi

parse_compare_pairs_and_build_used_aliases "$COMPARE_PAIRS"

DANPOS_EXE="$(resolve_danpos_cmd "$DANPOS_CMD")"
DANPOS_INPUT_SPEC="$(build_danpos_input_spec)"
DANPOS_BG_SPEC="$(build_danpos_bg_spec)"

mkdir -p "$OUT_DIR"

if [[ -z "$LOG_DIR" ]]; then
  LOG_DIR="${OUT_DIR}/logs"
fi
mkdir -p "$LOG_DIR"

RUN_TS="$(date '+%Y%m%d_%H%M%S')"
WRAPPER_STDOUT_LOG="${LOG_DIR}/wrapper.${RUN_TS}.out.log"
WRAPPER_STDERR_LOG="${LOG_DIR}/wrapper.${RUN_TS}.err.log"
RUN_METADATA="${LOG_DIR}/run_metadata.${RUN_TS}.txt"

exec > >(tee -a "$WRAPPER_STDOUT_LOG") 2> >(tee -a "$WRAPPER_STDERR_LOG" >&2)

if [[ "$HAS_DANPOS_ARGS" -eq 1 ]]; then
  DANPOS_EXTRA_ARGS_STR="$(quote_cmd "${DANPOS_ARGS[@]}")"
else
  DANPOS_EXTRA_ARGS_STR=""
fi

echo "run_timestamp=${RUN_TS}"
echo "hostname=$(hostname || true)"
echo "workdir=$(pwd)"
echo "outdir=${OUT_DIR}"
echo "log_dir=${LOG_DIR}"
echo "command=${COMMAND}"
echo "danpos_cmd=${DANPOS_CMD}"
echo "danpos_exe=${DANPOS_EXE}"
echo "threads=${THREADS}"
echo "skip_index=${SKIP_INDEX}"
echo "input_spec=${INPUT_SPEC}"
echo "background_spec=${BG_SPEC:-None}"
echo "compare_pairs=${COMPARE_PAIRS:-None}"
echo "danpos_input_spec=${DANPOS_INPUT_SPEC}"
echo "danpos_bg_spec=${DANPOS_BG_SPEC:-None}"
echo "wrapper_args=$(quote_cmd "${RAW_ARGS[@]}")"
echo "danpos_extra_args=${DANPOS_EXTRA_ARGS_STR}"

{
  echo "run_timestamp=${RUN_TS}"
  echo "hostname=$(hostname || true)"
  echo "workdir=$(pwd)"
  echo "outdir=${OUT_DIR}"
  echo "log_dir=${LOG_DIR}"
  echo "command=${COMMAND}"
  echo "danpos_cmd=${DANPOS_CMD}"
  echo "danpos_exe=${DANPOS_EXE}"
  echo "threads=${THREADS}"
  echo "skip_index=${SKIP_INDEX}"
  echo "input_spec=${INPUT_SPEC}"
  echo "background_spec=${BG_SPEC:-None}"
  echo "compare_pairs=${COMPARE_PAIRS:-None}"
  echo "danpos_input_spec=${DANPOS_INPUT_SPEC}"
  echo "danpos_bg_spec=${DANPOS_BG_SPEC:-None}"
  echo "wrapper_args=$(quote_cmd "${RAW_ARGS[@]}")"
  echo "danpos_extra_args=${DANPOS_EXTRA_ARGS_STR}"
  echo "danpos_path=$(command -v "$DANPOS_CMD" 2>/dev/null || true)"
  echo "samtools_path=$(command -v samtools || true)"
} > "$RUN_METADATA"

if [[ "$SKIP_INDEX" -eq 0 ]]; then
  ensure_all_bai
fi

CMD=(
  "$DANPOS_EXE" "$COMMAND"
  "$DANPOS_INPUT_SPEC"
)

if [[ -n "$DANPOS_BG_SPEC" ]]; then
  CMD+=(-b "$DANPOS_BG_SPEC")
fi

CMD+=(-o "$OUT_DIR")

if [[ "$HAS_DANPOS_ARGS" -eq 1 ]]; then
  CMD+=("${DANPOS_ARGS[@]}")
fi

echo "full_command=$(quote_cmd "${CMD[@]}")" | tee -a "$RUN_METADATA"
"${CMD[@]}"