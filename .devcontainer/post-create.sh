#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
git lfs install --local
bash scripts/check-environment.sh
