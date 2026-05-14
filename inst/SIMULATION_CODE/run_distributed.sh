#!/bin/bash

set -e
set -o pipefail

#-------------------------------------------------------------------------------
# CONFIGURATION
#-------------------------------------------------------------------------------
HOSTS_FILE="hosts.txt"

# Number of iterations per simulation setting
J_MAX=10000

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------
# Compute derived values
#-------------------------------------------------------------------------------
N_SETTINGS=$(Rscript -e "source('simulation_config.R'); cat(nrow(get_simulation_grid()))")

mkdir -p "results"
mkdir -p "logs"

echo "========================================"
echo "Simulation (DISTRIBUTED)"
echo "========================================"
echo "Hosts file:       $HOSTS_FILE"
echo "Settings:         1 to $N_SETTINGS"
echo "Iterations:       $J_MAX per setting"
echo "Total tasks:      $(( N_SETTINGS * J_MAX )) jobs"
echo "========================================"
echo ""

#-------------------------------------------------------------------------------
# Sync results bidirectionally (prevents duplicates on restart)
#-------------------------------------------------------------------------------
REMOTE_SIM_DIR="~/r_package_development/bayes2stage/inst/SIMULATION_CODE"
REMOTES=$(grep -v '^#' "$HOSTS_FILE" | grep -v '^[[:space:]]*$' | grep -v '/:$' | grep -v '^:$' | sed 's|^[0-9]*/||')

if [ -n "$REMOTES" ]; then
    echo "Syncing results..."
    for host in $REMOTES; do
        # Pull: remote → local
        rsync -a --ignore-existing "$host:$REMOTE_SIM_DIR/results/" "./results/" 2>/dev/null || true
        # Push: local → remote
        rsync -a --ignore-existing "./results/" "$host:$REMOTE_SIM_DIR/results/" 2>/dev/null || true
    done
    echo "Done."
    echo ""
fi

#-------------------------------------------------------------------------------
# Run distributed simulations
#-------------------------------------------------------------------------------
TOTAL_JOBS=$(( N_SETTINGS * J_MAX ))

parallel \
    --sshloginfile "$HOSTS_FILE" \
    --eta \
    --shuf \
    --retries 3 \
    --joblog "logs/parallel.log" \
    --halt soon,fail=5% \
    "eval \$(/usr/libexec/path_helper -s); cd ~/r_package_development/bayes2stage/inst/SIMULATION_CODE && mkdir -p logs && Rscript run_sim.R {1} {2} > logs/output_{1}_{2}.log 2>&1" \
    ::: $(seq 1 "$N_SETTINGS") \
    ::: $(seq 1 "$J_MAX") \
    2>&1 | perl -e '
        $| = 1;
        $total = '"$TOTAL_JOBS"';
        local $/ = "\r";
        while (<STDIN>) {
            # Parse: ETA: Ns Left: N AVG: N.Ns  host1:run/done/pct/avg  host2:run/done/pct/avg
            next unless /ETA:\s*(\d+)s/;
            my $eta = $1;

            # Sum running and done across all hosts (format: hostname:run/done/pct%/avg)
            my ($running, $done) = (0, 0);
            while (/([\w.-]*):(\d+)\/(\d+)\/\d+%/g) {
                $running += $2;
                $done += $3;
            }

            # Derive left from same source as done/running
            my $left = $total - $done - $running;

            my $pct = $total > 0 ? int($done * 100 / $total) : 0;
            my $eta_fmt = sprintf("%d:%02d:%02d", $eta/3600, ($eta%3600)/60, $eta%60);
            printf "\r  %d/%d (%2d%%) | Running: %2d | Left: %d | ETA: %s     ",
                   $done, $total, $pct, $running, $left, $eta_fmt;
        }
        print "\n";
    '

echo ""
echo "========================================"
echo "Simulation complete!"
echo ""
echo "Next steps:"
echo "  1. Collect results:  ./collect_results.sh"
echo "  2. Process results:  Rscript process_results.R"
echo "========================================"
