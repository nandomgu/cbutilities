#!/usr/bin/env bash
# Connect to the Host defined in ssh-remote/config (default: myserver).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${ROOT}/config"
HOST="${1:-myserver}"

if [[ ! -f "${CONFIG}" ]]; then
  echo "Missing ${CONFIG}" >&2
  echo "Copy the example and edit it:" >&2
  echo "  cp ${ROOT}/config.example ${CONFIG}" >&2
  exit 1
fi

exec ssh -F "${CONFIG}" "${HOST}"
