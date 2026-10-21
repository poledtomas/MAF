# Calculation modules and program structure

Users select calculations in the configuration. Each observable has a
small calculation function; file reading, selections, statistics and
output are shared. Installation and common commands are covered in the
[main README](../README.md).

## Program flow

1. `main.py` loads the job; `config.py` expands short configurations into explicit definitions.
2. `analysis.py` discovers runs and streams complete OSCAR events.
3. `observables/definition.py` selects particles and calls enabled modules.
4. `statistics/aggregation.py` retains raw contributions from independent units.
5. Aggregation is followed by value, error and covariance estimation and DAT/JSON/NPZ export.
6. `plotting/` later plots saved results without reading events again.

Each input file is read once during analysis. Individual modules receive
data from the same event. Cached contributions are stored per file;
final results are grouped by physical metadata. These two levels serve
different purposes.

## Module overview

| Files in `src/maf/observables/` | Calculation |
|---|---|
| `pt_spectra.py`, `mt_spectra.py` | pT and mT histograms and normalization |
| `dndy.py`, `dndeta.py`, `net_proton.py` | rapidity and pseudorapidity distributions |
| `yields.py`, `mean_pt.py`, `mean_mt.py`, `ratios.py` | scalar observables |
| `flow_reaction_plane.py`, `v1.py`, `v2.py` | flow relative to a known RP |
| `flow_event_plane.py`, `flow_resolution.py` | reconstructed EP and resolution |
| `flow_cumulants.py` | two- and four-particle Q correlations |
| `flow_settings.py` | configuration and shared flow-method adapters |
| `definition.py` | validation, selections, dispatch and normalization |

## Using a function directly

```python
from maf.observables import dndy

counts, event_count = dndy.calculate(
    rapidity=[-0.4, 0.2, 0.3],
    bins=[-1, 0, 1],
)
# counts = [1, 2], event_count = 1
```

`calculate(...)` receives data from **one event** and returns a numerator
and denominator. A yield returns the particle count and 1; mean-pT
returns the pT sum and particle count; particle-weighted RP returns cosine
sums and counts in each bin. Ensemble averaging is performed later using
the combined raw contributions.

EP also retains the subevent correlation; cumulants retain reference and
differential correlations. The number of stored features may differ from
the output bin count. NPZ `feature_layout` records their order.

## Adding an observable

Add a small file with `calculate` and, when needed, `normalization`.
In `definition.py`, add the import, the name to `_BUILTINS`, parameter
validation and a call from `calculate_event`. A nonlinear estimator also
needs an appropriate `estimate` implementation and accumulator layout.

Short names and aliases are expanded in `config.py`. `io/dat.py` defines
axis quantities and units for DAT and plots, keeping numerical export
and figure labels consistent.

Ordinary modules use plain functions and NumPy without inheritance or
decorators. `register_observable` is also available for custom Python API
extensions. Document physics definitions and analytical examples in the
relevant `docs/` page; numerical checks should verify the actual estimator.
