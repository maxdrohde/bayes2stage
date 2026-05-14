#!/bin/bash

set -e

cd "$(dirname "$0")"

HOSTS_FILE="hosts.txt"
REMOTE_SIM_DIR="~/r_package_development/bayes2stage/inst/SIMULATION_CODE"

# Extract remote hosts from hosts.txt (skip localhost ":", comments, empty lines)
REMOTES=$(grep -v '^#' "$HOSTS_FILE" | grep -v '^[[:space:]]*$' | grep -v '/:$' | grep -v '^:$' | sed 's|^[0-9]*/||')

echo "========================================"
echo "Cleaning results and logs"
echo "========================================"

# Clean local
echo ""
echo "→ local"
echo "  Deleting results/, logs/, plots/, processed_data/..."
rm -rf results/ logs/ plots/ processed_data/
echo "  Done."

# Clean remotes
for host in $REMOTES; do
    echo ""
    echo "→ $host"
    echo "  Deleting results/ and logs/..."
    ssh -o ConnectTimeout=10 "$host" "rm -rf $REMOTE_SIM_DIR/results $REMOTE_SIM_DIR/logs"
    echo "  Done."
done

echo ""
echo "========================================"
echo "All machines cleaned!"
echo "========================================"
