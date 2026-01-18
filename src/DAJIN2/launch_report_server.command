#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "Python is not available in PATH."
    exit 1
fi

PORT=8000
open "http://127.0.0.1:${PORT}/report.html"
"${PYTHON}" -m http.server "${PORT}"
