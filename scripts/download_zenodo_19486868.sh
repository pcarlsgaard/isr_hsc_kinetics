#!/usr/bin/env bash
# Download Zenodo record 19486868 (Milsom lab clonal hematopoiesis kinetics; DKFZ Academic
# Non-Commercial Research License - do not redistribute). Verifies Zenodo-published MD5s.
# Usage: scripts/download_zenodo_19486868.sh [manuscript|scrnaseq|all]
set -euo pipefail
what="${1:-manuscript}"
dest="$(cd "$(dirname "$0")/.." && pwd)/data/raw/zenodo_19486868"
base="https://zenodo.org/api/records/19486868/files"
mkdir -p "$dest"
fetch() { # name md5
    local out="$dest/$1"
    if [[ ! -f "$out" ]] || ! echo "$2  $out" | md5sum -c --status; then
        curl -fSL --retry 5 -C - "$base/$1/content" -o "$out"
    fi
    echo "$2  $out" | md5sum -c
    sha256sum "$out"
}
[[ "$what" == manuscript || "$what" == all ]] && fetch manuscript.zip 6ff4eab9aa4cb179ec048388d669663a
[[ "$what" == scrnaseq   || "$what" == all ]] && fetch scRNAseqData.zip 1b23a898cb7da5a2771b5cfac5696194
exit 0
