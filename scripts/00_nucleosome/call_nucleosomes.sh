#!/usr/bin/env bash
set -Eeuo pipefail

print_help() {
  cat <<'EOF'
Usage:
  bash call_nucleosomes.sh -i CHIP_PAIRS -b BG_PAIRS -o OUT_DIR [options] -- [DANPOS options...]

Required:
  -i, --input        DANPOS sample input string
                     Example:
                     "sample1_histone_mark1:/path/to/sample1_histone_mark1.bam,sample2_histone_mark1:/path/to/sample2_histone_mark1.bam"

  -b, --background   DANPOS background input string
                     Example:
                     "sample1_histone_mark1:/path/to/sample1_input.bam,sample2_histone_mark1:/path/to/sample2_input.bam"

  -o, --outdir       Output directory

Optional:
  -x, --command      DANPOS subcommand
                     Default: dpos

  -p, --danpos       Path to danpos.py
                     Default: danpos.py

  -t, --threads      Threads for samtools index
                     Default: value of SLURM_CPUS_PER_TASK or 1

  -c, --config       Config file with one key=value per line
                     Supported keys:
                       input=
                       background=
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
  All arguments after -- are passed directly to danpos.py unchanged.

Examples:
  bash call_nucleosomes.sh \
    -i "sample1_histone_mark1:/path/to/sample1_histone_mark1.bam,sample2_histone_mark1:/path/to/sample2_histone_mark1.bam" \
    -b "sample1_histone_mark1:/path/to/sample1_input.bam,sample2_histone_mark1:/path/to/sample2_input.bam" \
    -o /path/to/output_dir \
    -x dpos \
    -p /path/to/danpos.py \
    -t 32 \
    -- \
    -m 1 --mifrsz 80 --mafrsz 250 --extend 70

  bash call_nucleosomes.sh \
    -c /path/to/run.conf \
    -- \
    -m 1 --mifrsz 80 --mafrsz 250 --extend 70
EOF
}

COMMAND="dpos"
DANPOS_PY="danpos.py"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
INPUT_PAIRS=""
BG_PAIRS=""
OUT_DIR=""
CONFIG_FILE=""
LOG_DIR=""
SKIP_INDEX=0
DANPOS_ARGS=()
RAW_ARGS=("$@")

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
      input) INPUT_PAIRS="$value" ;;
      background) BG_PAIRS="$value" ;;
      outdir) OUT_DIR="$value" ;;
      command) COMMAND="$value" ;;
      danpos) DANPOS_PY="$value" ;;
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

extract_bams_from_pairs() {
  local pairs="$1"
  awk -F',' '
    {
      for (i = 1; i <= NF; i++) {
        x = $i
        sub(/^[^:]*:/, "", x)
        print x
      }
    }
  ' <<< "$pairs"
}

quote_cmd() {
  local out=""
  local arg
  for arg in "$@"; do
    out+=" $(printf '%q' "$arg")"
  done
  printf '%s\n' "${out# }"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      INPUT_PAIRS="$2"
      shift 2
      ;;
    -b|--background)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      BG_PAIRS="$2"
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
      DANPOS_PY="$2"
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
      DANPOS_ARGS+=("$@")
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

[[ -n "$INPUT_PAIRS" ]] || die "-i/--input is required"
[[ -n "$BG_PAIRS" ]] || die "-b/--background is required"
[[ -n "$OUT_DIR" ]] || die "-o/--outdir is required"

if [[ "$DANPOS_PY" != "danpos.py" ]]; then
  [[ -f "$DANPOS_PY" ]] || die "danpos.py not found: $DANPOS_PY"
fi

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

echo "run_timestamp=${RUN_TS}"
echo "hostname=$(hostname || true)"
echo "workdir=$(pwd)"
echo "outdir=${OUT_DIR}"
echo "log_dir=${LOG_DIR}"
echo "command=${COMMAND}"
echo "danpos_py=${DANPOS_PY}"
echo "threads=${THREADS}"
echo "skip_index=${SKIP_INDEX}"
echo "input_pairs=${INPUT_PAIRS}"
echo "background_pairs=${BG_PAIRS}"
echo "wrapper_args=$(quote_cmd "${RAW_ARGS[@]}")"
echo "danpos_extra_args=$(quote_cmd "${DANPOS_ARGS[@]}")"

{
  echo "run_timestamp=${RUN_TS}"
  echo "hostname=$(hostname || true)"
  echo "workdir=$(pwd)"
  echo "outdir=${OUT_DIR}"
  echo "log_dir=${LOG_DIR}"
  echo "command=${COMMAND}"
  echo "danpos_py=${DANPOS_PY}"
  echo "threads=${THREADS}"
  echo "skip_index=${SKIP_INDEX}"
  echo "input_pairs=${INPUT_PAIRS}"
  echo "background_pairs=${BG_PAIRS}"
  echo "wrapper_args=$(quote_cmd "${RAW_ARGS[@]}")"
  echo "danpos_extra_args=$(quote_cmd "${DANPOS_ARGS[@]}")"
  echo "python_path=$(command -v python || true)"
  echo "samtools_path=$(command -v samtools || true)"
} > "$RUN_METADATA"

if [[ "$SKIP_INDEX" -eq 0 ]]; then
  while IFS= read -r bam; do
    [[ -n "$bam" ]] && ensure_bai "$bam"
  done < <(
    {
      extract_bams_from_pairs "$INPUT_PAIRS"
      extract_bams_from_pairs "$BG_PAIRS"
    } | sort -u
  )
fi

CMD=(
  python "$DANPOS_PY" "$COMMAND"
  "$INPUT_PAIRS"
  -b "$BG_PAIRS"
  -o "$OUT_DIR"
  "${DANPOS_ARGS[@]}"
)

echo "full_command=$(quote_cmd "${CMD[@]}")" | tee -a "$RUN_METADATA"
"${CMD[@]}"
