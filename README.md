# Finalized StressME two-condition package

This package reproduces the finalized ME-model comparison between **FM01 +
NH4Cl** and **FM01 + Formamide**.

- **FM01 + NH4Cl:** the reverse direction of `LDH_D` and the forward direction
  of `LDH_D2` are blocked to remove the diagnosed internal lactate-
  dehydrogenase cycle.
- **FM01 + Formamide:** the finalized MFA0831-updated state is retained
  unchanged. Ammonium uptake and secretion are closed.

Both models use the updated MFA workbook, experimental exchange bounds,
experimental fixed growth rates, the betaxanthin pathway, and the overflow-
share constraints recorded in their `config.json` files.

## 1. System requirements

### Required

- A 64-bit computer with at least 4 CPU cores.
- At least 16 GB RAM recommended. Run the two conditions sequentially on
  systems with limited memory.
- At least 12 GB free disk space for Docker and simulation outputs.
- Docker Engine 24 or newer, or a recent Docker Desktop release.
- Internet access for the initial Docker image pull and one-time image build.

### Operating-system setup

- **Windows 10/11:** install Docker Desktop and enable the WSL2 backend. Start
  Docker Desktop before running `run_models.ps1`.
- **Linux:** install Docker Engine and ensure the current user can run
  `docker`, or prefix the documented commands with `sudo`.
- **macOS:** Docker Desktop can be used. On Apple Silicon, Docker may emulate
  the x86-64 image and run substantially slower. Add
  `--platform linux/amd64` to the Docker build/run commands if required.

No host Python installation is required to solve the ME models. The legacy
COBRAme/ecolime/StressME stack and qMINOS solver are supplied by the pinned
Docker image.

## 2. Software versions

The supplied `Dockerfile` builds `stressme-final:v1` from:

- `queensysbio/stressme:v1.1`
- Python and scientific packages supplied by that image
- `openpyxl==3.0.10`, added for reading the updated MFA workbook under the
  image's Python version

The parent image contains COBRA 0.5.11, COBRAme, ecolime, FoldME, AcidifyME,
OxidizeME, solvemepy/qMINOS, and the prebuilt StressME model.

Optional local figure analysis requires Python 3.10 or newer with:

```text
pandas
numpy
matplotlib
scipy
scikit-learn
sympy
openpyxl
```

These optional host packages are not needed for model solving because all
final result tables and figures are already included.

## 3. Package contents

```text
FINAL_ME_FM01_NH4CL_FORMAMIDE/
|-- Dockerfile
|-- run_models.ps1                 Windows launcher
|-- run_models.sh                  Linux/macOS/Git-Bash launcher
|-- inputs/                        MFA workbook and exchange-bound JSON files
|-- run_configs/                   Exact command arguments for both conditions
|-- source_code/                   Model runner and analysis scripts
|-- model_snapshots/               Exact configured models and qMINOS solutions
|-- results/
|   |-- complete_ME_results.xlsx   Consolidated workbook
|   `-- csv/                       Raw, combined and figure-traceability CSVs
|-- figures/                       Final carbon and cofactor PNGs
`-- FILE_MANIFEST_SHA256.csv       File sizes and SHA-256 checksums
```

The `model_snapshots/` directories contain the exact configured model pickle,
solution pickle, qMINOS basis, arguments, configuration and original run log
for each finalized condition.

## 4. Quick start

Open a terminal in this package directory.

Confirm that Docker is running:

```bash
docker version
```

### Windows PowerShell

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\run_models.ps1
```

### Linux, macOS or Git Bash

```bash
chmod +x ./run_models.sh
./run_models.sh
```

The launcher performs these steps:

1. verifies that Docker is available;
2. pulls `queensysbio/stressme:v1.1` when necessary;
3. builds the pinned `stressme-final:v1` derivative image;
4. creates `rerun_outputs/FM01_NH4Cl/` and
   `rerun_outputs/FM01_Formamide/`;
5. runs the two models sequentially;
6. reports each container exit code and leaves the complete `run.log` in its
   output directory.

Existing rerun output directories are not silently deleted. Move or rename
them before rerunning if a clean directory is required.

## 5. Exact finalized model designs

### FM01 + NH4Cl

- Temperature: 37 degrees C
- External pH parameter: 7.0
- ROS multiplier: 1
- Nitrogen source: ammonium
- Fixed growth rate: 0.4678 h^-1
- MFA sheet: `2. fmdA_M9R`
- MFA scale: 7.713728121938614
- MFA tolerance: +/-10%
- Experimental exchanges:
  `exchange_A_fmdA-am_biobounds.json`
- Betaxanthin pathway added; betaxanthin MFA reaction bounds skipped
- Overflow percentage constraints enabled
- Reverse `LDH_D` blocked
- Forward `LDH_D2` blocked

### FM01 + Formamide

- Temperature: 37 degrees C
- External pH parameter: 7.0
- ROS multiplier: 1
- Nitrogen source: formamide
- Direct FmdA formamide-to-formate pathway added
- Formamide uptake limit: 20 mmol gDW^-1 h^-1
- Ammonium exchange closed in both directions
- Fixed growth rate: 0.1295 h^-1
- MFA sheet: `4. fmdA_M9 formamide`
- MFA scale: 24.082458337347077
- MFA tolerance: +/-10%
- Experimental exchanges:
  `exchange_A_fmdA-fm_biobounds.json`
- Betaxanthin pathway added; betaxanthin MFA reaction bounds skipped
- Overflow percentage constraints enabled
- No LDH directional block is applied

The complete machine-readable argument lists are in `run_configs/`. The
complete applied settings and resulting phenotype are in each finalized
`config.json`.

## 6. Expected output files

Each rerun condition should produce:

- `config.json`: complete applied settings, constraint details and phenotype
- `configured_model.pickle`: exact configured ME model
- `solution.pickle`: qMINOS solution vector
- `qminos_basis.npy`: reusable solver basis
- `fluxes.csv`: all reported reaction fluxes with equations
- `proteome.csv`: predicted gene/protein mass fractions
- `phenotype.csv`: phenotype summary
- `mfa_comparison.csv`: mapped MFA values and solved model fluxes
- `run.log`: complete solver log

The original finalized solutions provide these reference values:

| Condition | Growth (h^-1) | Glucose uptake | Acetate secretion |
|---|---:|---:|---:|
| FM01 + NH4Cl | 0.4678 | 12.4367 | 1.1051 |
| FM01 + Formamide | 0.1295 | 4.0940 | 3.4766 |

Flux units are mmol gDW^-1 h^-1 unless stated otherwise. Very small numerical
differences can occur across Docker/CPU environments, but the run must remain
qMINOS feasible with the fixed growth rate and constraints satisfied.

## 7. Result and figure traceability

`results/complete_ME_results.xlsx` contains:

- full fluxes, proteome predictions and phenotypes for both conditions;
- flattened model configuration;
- MFA mapping and model/MFA comparisons;
- carbon-category totals and complete reaction-level carbon breakdown;
- complete reaction-level cofactor production and consumption;
- exact ATP, NADH, NADPH and ubiquinol records displayed in the figures.

The same information is available as machine-readable CSV files:

- `results/csv/simulation_outputs/`: original per-condition tables
- `results/csv/combined_results/`: stacked two-condition tables
- `results/csv/figure_traceability/carbon_balance/`: carbon figure inputs
- `results/csv/figure_traceability/cofactor_figures/`: complete and plotted
  cofactor data, including reaction names, equations, flux contributions and
  associated proteome mass fractions

For the carbon balance, use `carbon_balance_reaction_breakdown.csv` to resolve
every plotted category into its component exchange reactions. For cofactor
figures, the `*_figure_data.csv` files identify the plotted production and
consumption records, while `ME_cofactor_reaction_balance_long.csv` retains all
contributions, including values below the plotting cutoff.

## 8. Recomputing cofactor balances from new runs

The extraction script runs inside the StressME image because it needs to load
the configured ME-model pickles. Arrange or map the two rerun directories as
`sim_01` (NH4Cl) and `sim_02` (formamide), then run:

```bash
docker run --rm \
  -v "ABSOLUTE_PATH_TO_RESULTS:/results" \
  -v "ABSOLUTE_PATH_TO_PACKAGE/source_code:/scripts:ro" \
  stressme-final:v1 \
  python /scripts/extract_me_cofactor_balance.py
```

The script writes the complete steady-state cofactor reaction ledger and
turnover summary under the mounted results directory. Production and
consumption totals should match at steady state up to numerical tolerance.

## 9. Troubleshooting

### Docker daemon is unavailable

If `docker version` reports only client information, start Docker Desktop or
the Linux Docker service and retry.

### Container exits with code 137

The system ran out of memory. Run one condition at a time, close other
memory-intensive applications, and increase Docker's memory allocation.

### Permission denied on Linux

Run the launcher with `sudo`, or configure the current user for the Docker
group. Ensure the package and output directories are writable.

### Apple Silicon image error

Add `--platform linux/amd64` to `docker build` and `docker run`. Emulation may
be slow.

### Solver reports infeasible

Verify that the three files in `inputs/` are unchanged, that the correct
condition-specific JSON is used, and that the argument files in `run_configs/`
have not been edited. Compare the new `config.json` with the supplied finalized
snapshot.

## 10. Integrity and provenance

`FILE_MANIFEST_SHA256.csv` provides the size and SHA-256 digest of every
delivered file. Recompute hashes with `Get-FileHash -Algorithm SHA256` on
PowerShell or `sha256sum` on Linux/macOS. Simulation and analysis decisions are
documented in `METHODS_AND_PROVENANCE.md` and `PUBLICATION_METHODS.md`.

