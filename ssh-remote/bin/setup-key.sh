#!/usr/bin/env bash
# Generate a project SSH key (ed25519) under ~/.ssh if one does not exist.
set -euo pipefail

KEY="${HOME}/.ssh/id_ed25519"
COMMENT="${1:-$(whoami)@$(hostname)-cbutilities-ssh-remote}"

mkdir -p "${HOME}/.ssh"
chmod 700 "${HOME}/.ssh"

if [[ -f "${KEY}" ]]; then
  echo "Key already exists: ${KEY}"
  echo "Public key:"
  cat "${KEY}.pub"
  exit 0
fi

ssh-keygen -t ed25519 -f "${KEY}" -C "${COMMENT}" -N ""
chmod 600 "${KEY}"
chmod 644 "${KEY}.pub"

echo "Created ${KEY}"
echo "Add this public key to the server (~/.ssh/authorized_keys):"
echo
cat "${KEY}.pub"
