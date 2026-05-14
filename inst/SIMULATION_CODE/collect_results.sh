#!/bin/bash

set -e

cd "$(dirname "$0")"

HOSTS_FILE="hosts.txt"
REMOTE_SIM_DIR="~/r_package_development/bayes2stage/inst/SIMULATION_CODE"

# Extract remote hosts from hosts.txt (skip localhost ":", comments, empty lines)
REMOTES=$(grep -v '^#' "$HOSTS_FILE" | grep -v '^[[:space:]]*$' | grep -v '/:$' | grep -v '^:$' | sed 's|^[0-9]*/||')

if [ -z "$REMOTES" ]; then
    echo "No remote hosts found in $HOSTS_FILE"
    echo "Local results already in ./results/"
    exit 0
fi

echo "========================================"
echo "Collecting results from remote machines"
echo "========================================"

for host in $REMOTES; do
    (
        echo "→ $host"
        if ! rsync -a --whole-file -e "ssh -o Compression=no" "$host:$REMOTE_SIM_DIR/results/" "./results/" 2>&1; then
            echo "  WARNING: rsync failed for results from $host"
        fi
        if ! rsync -a --whole-file -e "ssh -o Compression=no" "$host:$REMOTE_SIM_DIR/logs/" "./logs/" 2>&1; then
            echo "  WARNING: rsync failed for logs from $host"
        fi
        echo "  $host done"
    ) &
done
wait

echo ""
echo "========================================"
echo "Collection complete!"
echo ""
echo "Results: ./results/"
echo "Logs:    ./logs/"
echo ""
echo "Next: Rscript process_results.R"
echo "========================================"
