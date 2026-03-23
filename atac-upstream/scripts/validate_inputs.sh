#!/usr/bin/env bash
set -euo pipefail

TYPE=""
INPUT_FILE=""
DRY_RUN=0
ERRORS=0

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok() { printf "${GREEN}[OK]${NC} %s\n" "$1"; }
warn() { printf "${YELLOW}[WARN]${NC} %s\n" "$1"; }
err() { printf "${RED}[FAIL]${NC} %s\n" "$1"; }

usage() {
  cat <<'EOF'
Usage: validate_inputs.sh --type <nfcore|encode> --input <file> [--dry-run]

Validate ATAC-seq pipeline input files.

Options:
  --type     Validation mode: nfcore or encode
  --input    Path to input file
  --dry-run  Print validation plan without execution
  -h, --help Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)
      TYPE="$2"
      shift 2
      ;;
    --input)
      INPUT_FILE="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      err "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$TYPE" || -z "$INPUT_FILE" ]]; then
  err "Both --type and --input are required."
  usage
  exit 1
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  warn "Dry-run mode enabled; no validation executed."
  printf "[DRY-RUN] Would validate '%s' as type '%s'\n" "$INPUT_FILE" "$TYPE"
  exit 0
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  err "Input file not found: $INPUT_FILE"
  exit 1
fi

validate_nfcore() {
  local header
  header=$(awk 'NF && $0 !~ /^#/ { print; exit }' "$INPUT_FILE")
  if [[ -z "$header" ]]; then
    err "No CSV header found in $INPUT_FILE"
    ERRORS=$((ERRORS + 1))
    return
  fi

  IFS=',' read -r c1 c2 c3 c4 _extra <<<"$header"
  if [[ "$c1" != "sample" || "$c2" != "fastq_1" || "$c3" != "fastq_2" || "$c4" != "replicate" ]]; then
    err "Invalid CSV header. Expected: sample,fastq_1,fastq_2,replicate"
    ERRORS=$((ERRORS + 1))
    return
  fi
  ok "CSV header is valid"

  local line_no=0
  while IFS= read -r line; do
    line_no=$((line_no + 1))

    if [[ -z "$line" || "$line" == \#* ]]; then
      continue
    fi

    if [[ "$line" == "$header" ]]; then
      continue
    fi

    IFS=',' read -r sample fastq1 fastq2 replicate extra <<<"$line"
    if [[ -z "$sample" || -z "$fastq1" || -z "$fastq2" || -z "$replicate" ]]; then
      err "Line $line_no has empty required fields"
      ERRORS=$((ERRORS + 1))
      continue
    fi

    if [[ -n "${extra:-}" ]]; then
      err "Line $line_no has extra columns"
      ERRORS=$((ERRORS + 1))
    fi

    if [[ ! -f "$fastq1" ]]; then
      err "Line $line_no missing fastq_1 file: $fastq1"
      ERRORS=$((ERRORS + 1))
    fi

    if [[ ! -f "$fastq2" ]]; then
      err "Line $line_no missing fastq_2 file: $fastq2"
      ERRORS=$((ERRORS + 1))
    fi
  done <"$INPUT_FILE"
}

validate_encode() {
  if ! command -v jq >/dev/null 2>&1; then
    err "jq is required to validate ENCODE JSON inputs"
    ERRORS=$((ERRORS + 1))
    return
  fi

  if ! jq empty "$INPUT_FILE" >/dev/null 2>&1; then
    err "Invalid JSON syntax: $INPUT_FILE"
    ERRORS=$((ERRORS + 1))
    return
  fi
  ok "JSON syntax is valid"

  local missing=0
  for key in '.atac.title' '.atac.genome_tsv' '.atac.paired_end'; do
    if ! jq -e "$key" "$INPUT_FILE" >/dev/null 2>&1; then
      err "Missing required key: $key"
      ERRORS=$((ERRORS + 1))
      missing=1
    fi
  done

  if [[ "$missing" -eq 0 ]]; then
    ok "Required ENCODE keys are present"
  fi
}

case "$TYPE" in
  nfcore)
    validate_nfcore
    ;;
  encode)
    validate_encode
    ;;
  *)
    err "Invalid --type '$TYPE'. Use 'nfcore' or 'encode'."
    exit 1
    ;;
esac

if [[ "$ERRORS" -eq 0 ]]; then
  ok "Validation passed for $INPUT_FILE"
  exit 0
fi

err "Validation failed with $ERRORS error(s)."
exit 1
