#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL_FILE="${SCRIPT_DIR}/publically_available_data.txt"
OUT_DIR="${SCRIPT_DIR}/data"

mkdir -p "${OUT_DIR}"

while IFS= read -r url || [[ -n "${url}" ]]; do
  [[ -z "${url}" || "${url}" =~ ^# ]] && continue
  wget --continue --no-verbose --show-progress \
    --directory-prefix="${OUT_DIR}" \
    "${url}"
done < "${URL_FILE}"
