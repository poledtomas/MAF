# Shared physics and statistics conventions

This page describes conventions shared by all observables. Individual
definitions and configuration examples are listed in the
[observable documentation index](README.md).

## Input and units

MAF reads final particle states in OSCAR2013/OSCAR2013Extended format.
Columns are mapped by header names. Four-momenta and masses use GeV with
c=1; header units are checked when present. The parser validates finite
values, nonnegative mass and energy, integer PDG IDs, charges and particle
IDs, and E≈√(m²+p²). Default tolerances of rtol=5×10⁻⁴ and atol=10⁻⁶ GeV
accommodate rounded ASCII input. Extra columns remain available, but
primary-particle classifications or universal decay corrections are
not inferred from them.

An accepted event must be complete, confirmed by an end marker, with its
final output block. A complete empty event is valid. Incomplete or invalid
records follow the configured reader policy and appear in diagnostics.
Initial states must not be mistaken for final hadronic samples even when
their header format is identical.

## Kinematics and cuts

```text
pT = sqrt(px²+py²)
y = 0.5 [ln(E+pz)−ln(E−pz)]
η = asinh(pz/pT)
mT = sqrt(m₀²+pT²).
```

y and η are not interchangeable for massive particles. At pT=0, a particle
moving along the beam has infinite η, while a particle at rest has
undefined η. Such particles remain available to observables that do not
require finite η. Histograms and cuts on a coordinate require that
coordinate to be finite. Azimuth is undefined at pT=0, so these particles
do not contribute to flow.

`y_range`, `eta_range` and `pt_range` acceptances are combined by
intersection. Cuts and histogram bins use [lower, upper) by default.
`boundary: open` changes acceptance cuts to strict inequalities without
changing histogram bin conventions. The final histogram upper edge is
excluded. Explicit acceptances must not truncate the declared bin width.

## Reference frames

The default `frame: input` uses the input momenta. Converting rapidity
between the laboratory and center-of-mass frames requires a known
`input_frame` and an explicit center-of-mass rapidity, `cm_rapidity`:

```text
y_CM = y_lab − Y_CM
y_lab = y_CM + Y_CM.
```

Y_CM is not inferred from an energy label. pT is unchanged by a
longitudinal boost. η cannot be shifted by subtracting Y_CM; η cuts and
distributions in a different frame are rejected because they would
require a full momentum boost.

## Ensemble and normalization

Selected files from one physical configuration form a common ensemble.
Raw event contributions are summed before final ratios and
normalizations are formed. An equally weighted mean of normalized
file spectra is generally a different estimator. A complete event
without selected particles contributes zero yield and remains in
the yield normalization.

Count distributions are divided by actual bin widths. Invariant spectra
use exact annular areas and an explicit Δy. Mean values and particle
ratios are ratios of ensemble sums. RP/EP weighting and empty-bin
policies are part of the configuration. Specific formulas are given
in the individual observable guides.

## Independent statistical units

Each event contributes a numerator A and denominator B for each quantity.
Positive event weights may be supplied through a `weight VALUE` field
in the end marker or through `sampling.events["event:ensemble"].weight`
in a manifest. Interaction weights from collision histories are not
event weights and are not used.

| `sampling.mode` | Meaning |
|---|---|
| `event` | individual events are explicitly independent |
| `hydro` | afterburner samples of one hydro event form one independent unit |
| `run` | an entire run is an independent block; the default mode records this assumption |

For `hydro`, supply `hydro_ids` or an `events` mapping with a hydro ID for
each event/ensemble. Contributions are first averaged with event weights
within each hydro event. Hydro events then receive equal ensemble weight
regardless of their oversample count. Separate hydro prior weights are
not implemented.

Hydro IDs are namespaced to the source file by default. For one hydro
production split across files, use a common `sampling.namespace` and
reuse IDs only for genuinely identical initial conditions. This
relationship is not inferred from run names or numbers. In `run` mode,
contributions are summed within a run, and central values retain event
weighting even for runs of unequal length.

## Statistical uncertainties

Uncertainties use a paired delete-one-independent-unit jackknife. It
recomputes the complete estimator, including denominators, EP resolution
or reference cumulants, preserving covariance between bins and related
contributions. For an unweighted mean, it gives the sample standard
deviation divided by √N.

At least two independent units and a valid estimator in every relevant
jackknife replicate are required. Otherwise, affected errors and
covariance elements are undefined: `null` in JSON and `nan` in DAT.
Errors from small or strongly unequal blocks can be unstable.

Within-observable bin covariance is stored in JSON/NPZ. Cross-observable
covariance is not exported by default, but paired raw contributions
allow it to be computed. A Poisson approximation is not substituted.
These errors quantify finite sampling under the declared hierarchy;
model systematics and centrality-selection uncertainty are not included.

## Metadata and experimental comparison

Energy, collision system and centrality must have explicit meanings from
metadata, a manifest or configured naming rules. A centrality label does
not itself select events. Use an appropriate preselected production or
a supported event cut; percentile-centrality calibration is not automatic.

Comparisons require compatible particles, acceptance, frame,
normalization, units and final-particle conventions. Feed-down and
detector corrections are not applied automatically. Supported count-bin
rebinning must preserve existing edges; flow is not rebinned through
that interface.
