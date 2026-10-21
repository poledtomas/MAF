# Analysis configuration

Configuration files can use YAML, JSON or TOML. YAML requires the `yaml`
extra; TOML requires Python 3.11 or later. JSON needs no additional parser.
Relative paths are resolved against the configuration file.

```yaml
input: /path/to/simulations
output: "results_{run_pattern}/"
run_pattern: "run_*"
workers: 1
group_by: []
sampling:
  mode: event
observables:
  dndy:
    enabled: true
    pdg: [211, -211, 2212]
    y_bins: {min: -2, max: 2, count: 40}
```

`run_pattern` is a glob filter for run directory names. It must select the
intended energy and parameters. `group_by: []` combines all selected files;
keep different energies and parameter sets separate. Event numbers may
repeat across files. `workers` is the number of inputs processed
concurrently. Cached contributions are stored in `output/cache/`.
`{run_pattern}` in the output path is replaced with a sanitized filter name
without glob characters.

## Particle selection and binning

| Short configuration key | Meaning |
|---|---|
| `pdg: [211, -211, 2212]` | a separate result for each species |
| `pdg: [all_charged]` | one combined result for all species with nonzero charge |
| omitted `pdg` | all species, including neutral particles |
| `charged: true` / `false` | charged / neutral particles only |
| `pt_range: [a, b]` | pT acceptance in GeV |
| `y_range: [a, b]` | rapidity y acceptance |
| `eta_range: [a, b]` | pseudorapidity η acceptance |
| `boundary: open` | strict cuts a < x < b; the default is a ≤ x < b |
| `y_bins`, `eta_bins`, `pt_bins`, `mt_bins` | bin edges or `{min, max, count}` |
| `enabled: false` | skip this block |
| `kind` | calculation module when the block name identifies a variant |
| `name` | custom base name for the output |

Cuts are combined by intersection. y and η are different quantities for
massive particles. Histogram bins always use [lower, upper], including
when `boundary: open` is set. The final upper edge is excluded. Bins must
not extend beyond an explicit acceptance on the same axis. Frame settings
are `frame`, `input_frame` and, when needed, `cm_rapidity`; see
[physics.md](physics.md).

Numeric PDG lists in the short format expand into separate results, such
as `dndy_211.dat` and `dndy_-211.dat`. The charged alias produces one file,
`dndy_all_charged.dat`. `pdg: [all]` is not supported; omit PDG to select all
species. In the explicit list format, `selection.pdg` defines a combined
selection. Species expansion applies only to the short mapping format.

## Available observables

| Module / kind | Result |
|---|---|
| [`pt_spectra` / `pt`](pt-spectra.md) | dN/dpT per event, or an invariant spectrum with `invariant: true` |
| [`mt_spectra` / `mt`](pt-spectra.md) | dN/dmT; requires one reference `mass` in GeV |
| [`mt_minus_m0`, `invariant_mt`](pt-spectra.md) | mT−m₀ or invariant mT spectrum |
| [`dndy`](rapidity.md) | dN/dy per event |
| [`net_proton`](rapidity.md) | d(Np−Npbar)/dy per event; omit PDG |
| [`dndeta`](rapidity.md) | dNch/dη per event; charged particles only |
| [`yields` / `yield`](yields-means-ratios.md) | mean multiplicity per event |
| [`charged_multiplicity`](yields-means-ratios.md) | mean charged multiplicity |
| [`mean_pt`, `mean_mt_minus_m0`](yields-means-ratios.md) | particle-weighted means |
| [`ratio`](yields-means-ratios.md) | ratio of ensemble yields, with explicit `denominator_pdg` |
| [`v1`, `v2`, `vn`](flow.md) | flow with explicit `flow.method`; `vn` requires `harmonic` |

Scalar observables do not accept bins. Invariant spectra require an
explicit rapidity window. Normalization, units and statistical conventions
are documented in [physics.md](physics.md); RP/EP settings are in
[flow.md](flow.md).

## Plot appearance

```yaml
plot:
  rc_params:
    figure.figsize: [7.0, 5.5]
    figure.dpi: 150
    font.size: 12
    axes.labelsize: 13
    legend.fontsize: 9
    axes.grid: false
    lines.markersize: 7
    errorbar.capsize: 3
```

These keys are Matplotlib rcParams. Styling changes do not affect the
numerical analysis or cache validity. Axis quantities and units are read
from DAT headers. To create PDFs:

```bash
maf plot --all --config analysis.yaml --format pdf
```

## Multiple variants of one calculation

The block name identifies the result. If it differs from a built-in module
name, specify `kind`, for example `v1_proton_y` with `kind: v1`. This allows
different cuts, axes or methods in one analysis. YAML keys must be unique.
Numeric PDG IDs add a species suffix; `all_charged` adds `_all_charged`.

A block with `enabled: false` is skipped. If `enabled` is omitted, the block
is active; at least one observable must be enabled. Within `flow.reference`,
use the full selection keys `pt`, `y`, `eta`, rather than the short keys
`pt_range`, `y_range`, `eta_range`. A reference is one combined selection,
even when it contains multiple numeric PDG IDs.

## Command-line controls

`maf --config analysis.yaml` and `maf run --config analysis.yaml` run the
same analysis. The installed `maf` command, `python -m maf` and
`python main.py` accept the same arguments. Use `maf --help` for the complete
command list and `maf <command> --help` for command-specific options.
`maf --version` prints the package version without loading a configuration.

| Option | Effect |
|---|---|
| `--dry-run` | validate the configuration and list selected runs without writing output |
| `--workers N` | override the configured number of concurrent input files |
| `--no-resume` | recompute contributions instead of using valid cache entries |
| `--fail-fast` | stop after the first failed input |
| `--verbose` | enable debug logging; accepted before or after a command |

Normal analysis progress reports count completed files, freshly analyzed
files, cache hits and failures, with a percentage and elapsed seconds.
A completed file can be successful or failed; inspect the failure count and
`failures.json` before using a partial result. A successful run returns exit
code 0, a command failure returns 1, and an analysis with saved partial
results returns 2. Argument errors use argparse's exit code 2.

Progress messages use stderr. JSON from `--dry-run` and `discover` stays on
stdout so it can be redirected to a file. A dry run validates settings and
discovers inputs; it does not scan every event or validate the entire OSCAR
contents.
