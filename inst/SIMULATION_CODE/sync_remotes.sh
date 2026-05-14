#!/bin/bash

set -e
set -o pipefail

cd "$(dirname "$0")"

HOSTS_FILE="hosts.txt"
PROJECT_DIR="/Users/max/r_package_development/bayes2stage"
REMOTE_PROJECT_DIR="~/r_package_development/bayes2stage"

# Extract remote hosts from hosts.txt (skip localhost ":", comments, empty lines)
REMOTES=$(grep -v '^#' "$HOSTS_FILE" | grep -v '^[[:space:]]*$' | grep -v '/:$' | grep -v '^:$' | sed 's|^[0-9]*/||')

echo "========================================"
echo "Checking R versions on all machines"
echo "========================================"

# Get local R version
LOCAL_R_VERSION=$(Rscript --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
echo "→ local: R $LOCAL_R_VERSION"

# Check remote R versions
VERSION_MISMATCH=0
for host in $REMOTES; do
    REMOTE_R_VERSION=$(ssh -o ConnectTimeout=10 "$host" "eval \$(/usr/libexec/path_helper -s); Rscript --version 2>&1" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
    if [ -z "$REMOTE_R_VERSION" ]; then
        echo "→ $host: R NOT INSTALLED"
        VERSION_MISMATCH=1
    elif [ "$REMOTE_R_VERSION" != "$LOCAL_R_VERSION" ]; then
        echo "→ $host: R $REMOTE_R_VERSION (MISMATCH)"
        VERSION_MISMATCH=1
    else
        echo "→ $host: R $REMOTE_R_VERSION"
    fi
done

if [ "$VERSION_MISMATCH" -eq 1 ]; then
    echo ""
    echo "ERROR: R version mismatch detected. All machines must have R $LOCAL_R_VERSION"
    exit 1
fi

echo ""
echo "========================================"
echo "Syncing packages on all machines (renv)"
echo "========================================"

#-------------------------------------------------------------------------------
# Local: restore renv and install bayes2stage
#-------------------------------------------------------------------------------
echo ""
echo "→ local"
echo "  Restoring renv packages..."
if ! (cd "$PROJECT_DIR" && Rscript -e "renv::restore(prompt = FALSE)"); then
    echo "  ERROR: renv restore failed locally"
    exit 1
fi
echo "  Installing bayes2stage..."
if ! (cd "$PROJECT_DIR" && Rscript -e "renv::install('.')"); then
    echo "  ERROR: bayes2stage install failed locally"
    exit 1
fi
echo "  Done."

#-------------------------------------------------------------------------------
# Remote installations
#-------------------------------------------------------------------------------
if [ -z "$REMOTES" ]; then
    echo ""
    echo "No remote hosts found in $HOSTS_FILE"
    echo "========================================"
    echo "Local sync complete!"
    echo "========================================"
    exit 0
fi

for host in $REMOTES; do
    echo ""
    echo "→ $host"

    # Ensure parent directory exists on remote
    ssh -o ConnectTimeout=10 "$host" "mkdir -p $REMOTE_PROJECT_DIR"

    # Sync project (exclude results/logs to avoid overwriting remote work)
    echo "  Syncing project files..."
    rsync -az --delete \
        --exclude 'results/' \
        --exclude 'logs/' \
        --exclude '.git/' \
        "$PROJECT_DIR/" "$host:$REMOTE_PROJECT_DIR/"

    # Restore renv and install bayes2stage
    echo "  Restoring renv packages..."
    if ! ssh -o ConnectTimeout=10 "$host" "eval \$(/usr/libexec/path_helper -s); cd $REMOTE_PROJECT_DIR && Rscript -e \"renv::restore(prompt = FALSE)\""; then
        echo "  ERROR: renv restore failed on $host"
        exit 1
    fi

    echo "  Installing bayes2stage..."
    if ! ssh -o ConnectTimeout=10 "$host" "eval \$(/usr/libexec/path_helper -s); cd $REMOTE_PROJECT_DIR && Rscript -e \"renv::install('.')\""; then
        echo "  ERROR: bayes2stage install failed on $host"
        exit 1
    fi

    echo "  Done."
done

echo ""
echo "========================================"
echo "All machines synced and ready!"
echo "========================================"
