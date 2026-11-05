#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/backup_before_revision.sh [--include-caches]

Creates a verified, timestamped pre-revision archive under backups/.
By default it includes source code, configuration, documentation, tests, results,
models and processed data, but excludes virtual environments and large re-downloadable
raw/interim caches. Use --include-caches to include data/raw and data/interim too.
EOF
}

INCLUDE_CACHES=0
for arg in "$@"; do
  case "$arg" in
    --include-caches) INCLUDE_CACHES=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$ROOT"
if [[ ! -f pyproject.toml || ! -d src/cml_repairnet ]]; then
  echo "ERROR: run this script from the CML-RepairNet project root." >&2
  echo "Current directory: $ROOT" >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP_DIR="$ROOT/backups"
ARCHIVE="$BACKUP_DIR/CML_RepairNet_pre_revision_${STAMP}.tar.gz"
META_DIR="$BACKUP_DIR/pre_revision_${STAMP}_metadata"
mkdir -p "$BACKUP_DIR" "$META_DIR"

printf 'project_root=%s\ncreated_at=%s\ninclude_caches=%s\n' \
  "$ROOT" "$(date --iso-8601=seconds)" "$INCLUDE_CACHES" > "$META_DIR/backup_info.txt"
python --version > "$META_DIR/python_version.txt" 2>&1 || true
python -m pip freeze > "$META_DIR/pip_freeze.txt" 2>&1 || true
find . -maxdepth 4 -type f \
  ! -path './.git/*' ! -path './.venv/*' ! -path './venv/*' \
  ! -path './backups/*' -printf '%P\n' | sort > "$META_DIR/file_inventory_before_revision.txt"

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git status --short > "$META_DIR/git_status_short.txt" || true
  git status > "$META_DIR/git_status.txt" || true
  git rev-parse HEAD > "$META_DIR/git_head.txt" || true
  git diff --binary > "$META_DIR/uncommitted_changes.patch" || true
  git log -n 20 --oneline --decorate > "$META_DIR/git_log_last20.txt" || true
  git bundle create "$META_DIR/repository.bundle" --all >/dev/null 2>&1 || true
else
  echo "not_a_git_worktree" > "$META_DIR/git_status.txt"
fi

INCLUDES=(
  pyproject.toml README.md LICENSE Makefile environment.yml CITATION.cff
  src scripts tests configs docs results models data/processed
  manuscript supplementary figures tables
)
if [[ "$INCLUDE_CACHES" -eq 1 ]]; then
  INCLUDES+=(data/raw data/interim data/external)
fi

EXISTING=()
for path in "${INCLUDES[@]}"; do
  [[ -e "$path" ]] && EXISTING+=("$path")
done
if [[ ${#EXISTING[@]} -eq 0 ]]; then
  echo "ERROR: no project paths were found for backup." >&2
  exit 1
fi

# Metadata is copied into a temporary path under the project for one atomic archive,
# then removed. backups/ itself is explicitly excluded to avoid recursive archives.
TMP_META=".pre_revision_metadata_${STAMP}"
rm -rf "$TMP_META"
cp -a "$META_DIR" "$TMP_META"
trap 'rm -rf "$TMP_META"' EXIT

tar \
  --exclude='./.git' --exclude='./.venv' --exclude='./venv' \
  --exclude='./backups' --exclude='*/__pycache__' --exclude='*.pyc' \
  -czf "$ARCHIVE" "${EXISTING[@]}" "$TMP_META"

tar -tzf "$ARCHIVE" >/dev/null
sha256sum "$ARCHIVE" > "${ARCHIVE}.sha256"
stat -c '%n %s bytes' "$ARCHIVE" > "$META_DIR/archive_size.txt"
cp "${ARCHIVE}.sha256" "$META_DIR/"

printf '\nPASS: pre-revision backup created and verified.\n'
printf 'Archive: %s\n' "$ARCHIVE"
printf 'SHA-256: %s\n' "$(cut -d' ' -f1 "${ARCHIVE}.sha256")"
printf 'Metadata: %s\n' "$META_DIR"
printf '\nDo not delete this backup until the revised pipeline and manuscript are frozen.\n'
