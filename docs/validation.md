# MAF 0.5.0 validation

The current suite contains **141 tests; all passed**. Validation used
Python 3.9.12, NumPy 2.0.2, SciPy 1.13.1,
Matplotlib 3.9.4 and PyYAML 6.0.3.
Production simulations were not recomputed.

From the project root, in an environment with the required dependencies:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Tests use synthetic events and do not require access to production data.
Tests requiring optional dependencies are skipped when those dependencies
are unavailable.

## RP and EP versus y and η

Seven tests in `tests/test_flow_rapidity.py` check:

- general harmonic order, a known nonzero RP and the v1/v2 aliases;
- all 14 outputs of the four supplied configurations for n=1,2,3;
- direct scalar sums of cosines, reference vectors and subevent correlations;
- charged POI, neutral reference particles and self-subtraction only for overlap;
- the difference between rapidity and pseudorapidity for massive particles;
- zero versus excluded empty bins and ineligible subevents;
- settings validation, YAML/JSON agreement, NPZ/DAT and loading 0.4.1 definitions.

The test reference is an independent scalar calculation. RP/EP definitions
and physics assumptions are documented in [flow.md](flow.md).

A separate installed-command check with `maf` 0.5.0 covered input discovery,
analysis and plotting: four synthetic OSCAR files, 20 events, two workers,
one combined ensemble, 14 DAT tables and 14 PNG files. DAT quantities
and units for y/η and v₂{RP}/v₂{EP} were also verified. The check's inputs
and outputs were stored in a temporary directory outside the project.

## Q cumulants

Expressions are compared against explicit enumeration of distinct ordered
pairs and quadruplets: all POI/RFP selections up to four particles and
random selections up to eight particles, including disjoint and partially
overlapping selections, for n=1,2,3. Additional checks cover low
multiplicities, large combinatorial weights, c₂/c₄ signs and differential
estimators. A proton/pion-reference integration test checks tuple weighting,
joint jackknife and covariance, merging and NPZ. Definitions and source
equations are in [q-cumulants.md](q-cumulants.md).

## Input, particles and statistics

Other tests check OSCAR headers, units, incomplete events, independent
y/η/pT cuts, known and unknown geometry, and physics normalization.
`all_charged` is checked against actual electric charge, including
antiparticles, neutral particles and multiply charged nuclei. Numeric
PDG lists still expand into separate outputs in the short configuration.

Aggregation checks cover unequal multiplicities, empty events, sampling
hierarchies and weights, including covariance and joint EP-resolution
uncertainty. The suite also covers cache, resuming, saving, merging, DAT
and plots. Older NPZ versions, including 0.4.0/0.4.1, remain readable.
New analyses include 0.5.0 in cache signatures to distinguish them from
older calculations.

## Simulation labels

Plot legends use the saved `run_pattern`. Legacy outputs can use the
configuration pattern or their result directory label. Checks cover saved
labels, explicit Python API label overrides, cached analysis and persistence
through NPZ, JSON and DAT. The run filter is result metadata; it does not
change the numerical cache identity.

## CLI and packaging

The CLI checks cover complete root and plot help, version reporting without
a configuration, explicit `run` and implicit `--config` execution, verbosity
before and after a command, dry-run JSON and the absence of dry-run output
files. Serial and parallel runs verify fresh/cache/failure counts, forced
recomputation, unchanged numerical summaries and checkpoint reuse. Partial
results retain exit code 2; fail-fast reports the stopped file count.

A wheel was built with isolated build dependencies and loaded outside the
source tree. Its version, author, repository/documentation/issue URLs and
bundled MIT license were checked, together with version agreement in
`CITATION.cff`. Wheel help/version and all three source launchers (`maf`,
`python -m maf`, `python main.py`) were checked.

The package and CLI version are defined in `src/maf/_version.py`;
`pyproject.toml` reads that value when building. Update the citation version
in `CITATION.cff` when preparing a new release. Build tools require
setuptools 77 or later for SPDX license metadata; runtime still supports
Python 3.9 or later.
