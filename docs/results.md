# Results, DAT tables and accumulators

The configuration's `output` setting selects the output directory.
For a single physical ensemble, results are stored directly in its
`results/` subdirectory:

```text
output/
├── index.json
├── failures.json
├── results/
│   ├── summary.json
│   ├── accumulators.npz
│   └── <observable>.dat
├── plots/
└── cache/
```

When `group_by` creates multiple groups, each has its own
`results/group-<id>/`. `index.json` lists current accumulators and
summaries; older files from other configurations are not automatically
selected for plotting. `failures.json` lists run-processing failures
and is an empty list after a successful analysis. Partial results are
marked in diagnostics and messages. The `plot` command creates `plots/`.

## DAT tables

Each result has its own `.dat` file. Numeric columns are:

| Column | Meaning |
|---|---|
| 1 | bin center or scalar-value index |
| 2 | lower bin edge |
| 3 | upper bin edge |
| 4 | normalized value |
| 5 | statistical error |

Header lines start with `#` and record the observable name, axis
quantities and units, analysis signature, complete JSON calculation
definition and column descriptions. For example, flow quantities are
`v_1{RP}`, `v_2{EP}`, `v_2{2}` or `v_2{4}`, with unit `1`. Scalar results
have one row with index 0 and `nan` bin boundaries. Undefined values
and errors are also written as `nan`.

Dimensionless quantities retain unit `1` in DAT/JSON metadata. Plot labels
omit `[1]` while keeping physical units such as `[GeV]` and GeV⁻¹.
Canonical quantities are typeset with Matplotlib MathText: subscripts such
as v₁, v₂, Nch and pT, Greek symbols, fractions and unit exponents use
mathematical notation. This does not require an external LaTeX installation.
Custom prose labels are preserved. DAT/JSON quantities and units stay in
their original machine-readable form.

Plots read values and axes from DAT when the table is present. Signatures
and edges are checked against the summary. Missing required headers,
signature mismatches or inconsistent edges are errors. Covariance is
stored in JSON/NPZ rather than extra DAT columns.

## summary.json

The summary records the configuration, software versions, group metadata,
event and independent-unit counts, sources and diagnostics. The optional
result-level `run_pattern` records the simulation filter used for the
legend. It is also retained in NPZ metadata and DAT headers as a JSON string;
older files without it remain readable. Each
observable contains:

| Field | Contents |
|---|---|
| `definition` | kind, particles, cuts, bins, units and method settings |
| `value`, `stat_error` | central values and statistical errors |
| `covariance` | covariance of output bins |
| `normalization` | applied normalization factors |
| `axes`, `units`, `bins` | quantities, units and edges |
| `auxiliary` | for example, EP resolution or reference c₂/c₄ |

Infinite and undefined values are JSON `null`. Source records retain
input identity and fingerprint, metadata, OSCAR headers and parser
diagnostics, sampling and selected-event counts.

`unknown_impact` counts events without a known impact parameter, for
example a value of −1. They are valid inputs when no impact-parameter
cut is requested. A geometry cut requires a known value.

## accumulators.npz

NPZ is a checkpoint of raw statistics for independent units, loaded
without pickle (`allow_pickle=False`). The schema version is 1.

| Field | Contents |
|---|---|
| `metadata` | JSON with schema, version, configuration, definitions, sources and unit identities |
| `a` | raw numerators, shape (units, features) |
| `b` | raw denominators with the same shape |
| `w` | sum of event weights within each unit |
| `events` | number of complete selected events within each unit |

Features concatenate observables in configuration order. Ordinary
histograms use one feature per bin. EP additionally stores the subevent
correlation; cumulants store reference and differential correlations.
Metadata `feature_layout` names each feature. The raw feature count may
therefore differ from the exported bin count.

In hydro mode, raw contributions are divided by the sum of weights within
the hydro unit during estimation, after any cross-file merge. Otherwise,
unequal afterburner oversampling would change hydro-event ensemble weights.

## Cache and merging

`cache/` contains atomic input checkpoints. Cache keys include the file
fingerprint, configuration, software version and metadata. Changes to
numerical definitions or inputs therefore invalidate affected cached
contributions. Changes to `plot.rc_params` do not affect cache validity.
`--no-resume` forces recomputation.

The `merge` command combines compatible NPZ files with the same definition
and sampling. It combines raw contributions rather than averaging final
spectra, errors or resolutions. Duplicate source identities are rejected.
Physically identical samples copied under different paths cannot be
reliably detected, so provide disjoint inputs.

Supported older results from 0.1.0 through 0.4.1 remain readable according
to the explicit version list in the loader. New analyses include the
current version in cache signatures. Older v1/v2 definitions infer the
harmonic from `kind`; EP definitions without `empty_bins` retain
`exclude`. Loading does not recompute their numerical results.

See the [main README](../README.md) for plotting commands.
