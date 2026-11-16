#!/usr/bin/env bash
set -euo pipefail

PACKAGE_ROOT="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
IMAGE="stressme-final:v1"
PARENT_IMAGE="queensysbio/stressme:v1.1"

docker version >/dev/null
docker image inspect "$PARENT_IMAGE" >/dev/null 2>&1 || docker pull "$PARENT_IMAGE"
docker image inspect "$IMAGE" >/dev/null 2>&1 || docker build -t "$IMAGE" "$PACKAGE_ROOT"

mkdir -p "$PACKAGE_ROOT/rerun_outputs"

run_condition() {
    local name="$1"
    local config="$2"
    local output="$PACKAGE_ROOT/rerun_outputs/$name"
    local container_name
    container_name="final_me_$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')"
    mkdir -p "$output"
    echo "Running $name..."
    docker run --rm \
        --name "$container_name" \
        -v "$PACKAGE_ROOT/inputs:/mfa:ro" \
        -v "$PACKAGE_ROOT/source_code:/scripts:ro" \
        -v "$PACKAGE_ROOT/run_configs:/run_configs:ro" \
        -v "$output:/mount_point" \
        "$IMAGE" \
        bash -lc "python /scripts/launch_from_json.py /run_configs/$config > /mount_point/run.log 2>&1"
    echo "Completed $name: $output"
}

run_condition "FM01_NH4Cl" "FM01_NH4Cl_args.json"
run_condition "FM01_Formamide" "FM01_Formamide_args.json"
echo "All finalized simulations completed."

