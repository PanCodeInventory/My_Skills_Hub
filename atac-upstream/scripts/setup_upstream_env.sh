#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="upstream"
DRY_RUN=0
AUTO_YES=0

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok() { printf "${GREEN}[OK]${NC} %s\n" "$1"; }
warn() { printf "${YELLOW}[WARN]${NC} %s\n" "$1"; }
fail() { printf "${RED}[FAIL]${NC} %s\n" "$1"; }

usage() {
  cat <<'EOF'
Usage: setup_upstream_env.sh [--dry-run] [-y]

Recreate the "upstream" mamba environment for ATAC-seq upstream analysis.

Options:
  --dry-run  Print commands without executing
  -y         Auto-confirm environment recreation
  -h, --help Show this help message
EOF
}

run_cmd() {
  local cmd="$1"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf "[DRY-RUN] %s\n" "$cmd"
  else
    eval "$cmd"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -y)
      AUTO_YES=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      fail "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

if ! command -v mamba >/dev/null 2>&1; then
  fail "mamba is not available in PATH."
  exit 1
fi

if [[ "$AUTO_YES" -ne 1 ]]; then
  printf "This will remove and recreate environment '%s'. Continue? [y/N]: " "$ENV_NAME"
  read -r reply
  if [[ ! "$reply" =~ ^[Yy]$ ]]; then
    warn "Aborted by user."
    exit 0
  fi
fi

run_cmd "mamba env remove -n ${ENV_NAME} -y || true"
run_cmd "mamba create -n ${ENV_NAME} -y python=3.11 nextflow nf-core openjdk=17 pip git jq yq samtools bedtools macs2 multiqc fastqc"
run_cmd "mamba run -n ${ENV_NAME} pip install --upgrade pip"
run_cmd "mamba run -n ${ENV_NAME} pip install caper croo looper pypiper eido polars-lts-cpu pyarrow"

if [[ "$DRY_RUN" -eq 1 ]]; then
  ok "Dry run complete."
  exit 0
fi

verify_cmd() {
  local label="$1"
  local cmd="$2"
  if eval "$cmd" >/dev/null 2>&1; then
    ok "$label"
  else
    fail "$label"
    return 1
  fi
}

VERIFY_FAILED=0

verify_cmd "python installed" "mamba run -n ${ENV_NAME} python --version" || VERIFY_FAILED=1
verify_cmd "nextflow installed" "mamba run -n ${ENV_NAME} nextflow -version" || VERIFY_FAILED=1
verify_cmd "nf-core installed" "mamba run -n ${ENV_NAME} nf-core --version" || VERIFY_FAILED=1
verify_cmd "java installed" "mamba run -n ${ENV_NAME} java -version" || VERIFY_FAILED=1
verify_cmd "git installed" "mamba run -n ${ENV_NAME} git --version" || VERIFY_FAILED=1
verify_cmd "jq installed" "mamba run -n ${ENV_NAME} jq --version" || VERIFY_FAILED=1
verify_cmd "yq installed" "mamba run -n ${ENV_NAME} yq --version" || VERIFY_FAILED=1
verify_cmd "samtools installed" "mamba run -n ${ENV_NAME} samtools --version" || VERIFY_FAILED=1
verify_cmd "bedtools installed" "mamba run -n ${ENV_NAME} bedtools --version" || VERIFY_FAILED=1
verify_cmd "macs2 installed" "mamba run -n ${ENV_NAME} macs2 --version" || VERIFY_FAILED=1
verify_cmd "multiqc installed" "mamba run -n ${ENV_NAME} multiqc --version" || VERIFY_FAILED=1
verify_cmd "fastqc installed" "mamba run -n ${ENV_NAME} fastqc --version" || VERIFY_FAILED=1
verify_cmd "caper installed" "mamba run -n ${ENV_NAME} caper --version" || VERIFY_FAILED=1
verify_cmd "croo installed" "mamba run -n ${ENV_NAME} croo --version" || VERIFY_FAILED=1
verify_cmd "looper installed" "mamba run -n ${ENV_NAME} looper --version" || VERIFY_FAILED=1
verify_cmd "pypiper installed" "mamba run -n ${ENV_NAME} python -c 'import pypiper'" || VERIFY_FAILED=1
verify_cmd "eido installed" "mamba run -n ${ENV_NAME} python -c 'import eido'" || VERIFY_FAILED=1
verify_cmd "polars-lts-cpu installed" "mamba run -n ${ENV_NAME} python -c 'import polars; print(polars.__version__)'" || VERIFY_FAILED=1
verify_cmd "pyarrow installed" "mamba run -n ${ENV_NAME} python -c 'import pyarrow; print(pyarrow.__version__)'" || VERIFY_FAILED=1

if [[ "$VERIFY_FAILED" -eq 0 ]]; then
  ok "Environment '${ENV_NAME}' is ready."
  exit 0
fi

fail "One or more verification checks failed."
exit 1
