param(
    [string]$PackageRoot = $PSScriptRoot
)

$ErrorActionPreference = "Stop"
$image = "stressme-final:v1"
$parentImage = "queensysbio/stressme:v1.1"
$inputs = Join-Path $PackageRoot "inputs"
$source = Join-Path $PackageRoot "source_code"
$configs = Join-Path $PackageRoot "run_configs"
$outputs = Join-Path $PackageRoot "rerun_outputs"

docker version | Out-Null
if ($LASTEXITCODE -ne 0) { throw "Docker is not available or the daemon is not running." }

if (-not (docker image inspect $parentImage 2>$null)) {
    docker pull $parentImage
    if ($LASTEXITCODE -ne 0) { throw "Failed to pull $parentImage" }
}
if (-not (docker image inspect $image 2>$null)) {
    docker build -t $image $PackageRoot
    if ($LASTEXITCODE -ne 0) { throw "Failed to build $image" }
}

New-Item -ItemType Directory -Force -Path $outputs | Out-Null
$runs = @(
    @{ Name = "FM01_NH4Cl"; Config = "FM01_NH4Cl_args.json" },
    @{ Name = "FM01_Formamide"; Config = "FM01_Formamide_args.json" }
)

foreach ($run in $runs) {
    $out = Join-Path $outputs $run.Name
    New-Item -ItemType Directory -Force -Path $out | Out-Null
    Write-Host "Running $($run.Name)..."
    docker run --rm `
        --name ("final_me_" + $run.Name.ToLower()) `
        -v "${inputs}:/mfa:ro" `
        -v "${source}:/scripts:ro" `
        -v "${configs}:/run_configs:ro" `
        -v "${out}:/mount_point" `
        $image `
        bash -lc "python /scripts/launch_from_json.py /run_configs/$($run.Config) > /mount_point/run.log 2>&1"
    if ($LASTEXITCODE -ne 0) {
        throw "$($run.Name) failed. Inspect $out\run.log"
    }
    Write-Host "Completed $($run.Name): $out"
}

Write-Host "All finalized simulations completed."

