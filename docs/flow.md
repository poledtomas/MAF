# Flow vₙ: reaction plane and event plane

MAF computes flow versus rapidity y, pseudorapidity η or transverse
momentum pT.
## Harmonic order, particles and axis

| Setting | Meaning |
|---|---|
| `kind: v1` | first harmonic, n=1 |
| `kind: v2` | second harmonic, n=2 |
| `kind: vn`, `harmonic: 3` | any positive integer harmonic |
| `axis: y`, `y_bins` | rapidity dependence |
| `axis: eta`, `eta_bins` | pseudorapidity dependence |
| `axis: pt`, `pt_bins` | pT dependence |
| `pdg: [2212]` | protons |
| `pdg: [211, -211, 2212]` | separate results for each species |
| `pdg: [all_charged]` | one combined result for all charged species |

The block name identifies the output; `kind` and, where applicable,
`harmonic` determine n. If `harmonic` is specified for `v1` or `v2`, it
must match that kind. Particles with pT=0 have undefined azimuth and do
not enter flow calculations. Rapidity and pseudorapidity are separate
axes; their cuts are combined only when both are specified.

## Proton v₁{RP}(y)

Place the following block under `observables` in the configuration:

```yaml
observables:
  v1_proton_y:
    enabled: true
    kind: v1
    pdg: [2212]
    axis: y
    y_bins: {min: -1, max: 1, count: 10}
    pt_range: [0.2, 1.6]
    eta_range: [-1, 1]
    boundary: open
    flow:
      method: reaction_plane
      angle: 0.0
      weighting: event
```

This computes cos[n(φ−ΨRP)] relative to a **known** reaction plane.
`angle` is in radians. Zero applies to a production aligned with the known
plane ΨRP=0; the impact-parameter magnitude does not determine its azimuth.
For different known planes in individual events, replace `angle` with
`metadata_key: reaction_plane` and supply each event's angle through
manifest metadata. A missing angle is an error.

`weighting: event` first averages particles within each occupied event
bin, then averages those event means. An empty bin in an event is omitted.
`weighting: particle` divides the combined cosine sum by the combined
particle count. The estimators differ for unequal multiplicities.
The RP default is `particle`; this example explicitly selects `event`.
RP requires neither plane reconstruction nor a resolution correction.

The example produces `v1_proton_y_2212.dat`, with quantity `v_1{RP}`.
For v₂, change `kind` to `v2` and rename the block. For higher n, use
`kind: vn` and `harmonic`. Multiple blocks can run together.

## Identified-hadron v₂{EP}(y)

Selected particles of interest are called POI. The independent
`flow.reference` selection defines the particles used to reconstruct
the plane.

```yaml
observables:
  v2_ep_y_identified:
    enabled: true
    kind: vn
    harmonic: 2
    pdg: [211, -211, 321, -321, 2212, -2212]
    charged: true
    axis: y
    y_range: [-1, 1]
    pt_range: [0.2, 1.6]
    boundary: open
    y_bins: {min: -1, max: 1, count: 10}
    flow:
      method: event_plane
      weighting: event
      empty_bins: exclude
      subevents: alternating
      reference:
        pdg: [all_charged]
        eta: [-1, 1]
        pt: [0.2, 1.6]
        boundary: open
```

POI and reference selections may differ. In this example, the POI have
no η cut, while the reference does. Omitting `pdg` from the reference
includes neutral particles. A numeric PDG list within the reference
forms one combined selection and does not expand into separate planes.

### Reconstruction and autocorrelations

Each event provides the pT-weighted vector Qₙ=Σj pT,j exp(i n φj),
with plane angle Ψₙ=arg(Qₙ)/n. For each POI, its own vector is subtracted
if it belongs to the reference. A POI outside the reference leaves Q
unchanged. Its observed contribution is cos[n(φPOI−Ψₙ)] relative to
this corrected plane.

The reference is split by even and odd indices in the original particle
list (`subevents: alternating`). The correlation between planes A and B
gives

```text
C = ⟨cos[n(ΨA−ΨB)]⟩,
Rsub = √C,
R(χsub) = Rsub,
Rfull = R(√2 χsub),
vₙ{EP} = vₙ,obs / Rfull.
```

The resolution function is

```text
R(χ) = √π/(2√2) χ exp(−χ²/4) [I₀(χ²/4) + I₁(χ²/4)].
```

The inversion uses SciPy from the `flow` installation extra. Resolution
is obtained from the ensemble correlation, rather than averaging
resolutions from individual files. `summary.json` stores the correlation
and resolution under `auxiliary`. This convention is described by
[Poskanzer and Voloshin, arXiv:nucl-ex/9805001](https://arxiv.org/abs/nucl-ex/9805001).

### Weighting and empty bins

| Setting | Ensemble mean of observed flow |
|---|---|
| `weighting: event`, `empty_bins: exclude` | eligible events with a valid contribution to the bin |
| `weighting: event`, `empty_bins: zero` | all eligible events; an empty bin contributes zero |
| `weighting: particle` | combined cosine sums divided by the number of valid POI |

The EP defaults are `weighting: event` and `empty_bins: exclude`.
`zero` is allowed only with event weighting. Including empty bins can
substantially change results for rare species; the policy is an explicit
part of the observable definition. The supplied rapidity configurations
use `zero`, while the example above demonstrates `exclude`.

An event with an empty or numerically zero Q in either subevent is
excluded from EP. A POI with zero corrected Q has no defined plane and
provides no valid contribution; the resulting empty bin then follows
the selected policy. C≤0 or missing resolution gives `null` in JSON and
`nan` in DAT. A zero vector is not assigned a plane angle of zero.

## Supplied y and η configurations

[flow_rapidity.yaml](../examples/flow_rapidity.yaml) and
[flow_rapidity.json](../examples/flow_rapidity.json) contain:

| Block | POI | Reference |
|---|---|---|
| `v2_rp_y_identified` | each PDG separately; −1<y<1, −1<η<1, 0.2<pT<1.6 | known RP |
| `v2_ep_y_identified` | each charged PDG separately; −1<y<1, 0.2<pT<1.6 | all species, −1<η<1, 0.2<pT<1.6 |
| `v2_ep_y` | all charged; −1<y<1, 0.2<pT<2 | all species, −1<η<1, 0.2<pT<2 |
| `v2_ep_eta` | all charged; −2<η<2, 0.2<pT<2 | all species, −2<η<2, 0.2<pT<2 |

pT is in GeV. The y axis has 10 bins; η has 25 bins. Species, cuts and
harmonic order can be changed. For the first harmonic, use a name such
as `v1_ep_y` and set `harmonic: 1`. For a pT dependence, use `axis: pt`
and `pt_bins`.

## Physics assumptions and statistics

Extrapolating subevent resolution to the full plane assumes statistically
comparable subevents. The even/odd split depends on particle order and
does not guarantee equal multiplicity. An η gap or spectator plane is
not introduced automatically. For n=1, the reference uses positive pT
weights without a sign(y) factor; RP and EP may therefore measure
different components of v₁.

`sampling` defines the independent statistical unit. The jackknife
recomputes observed flow and resolution jointly, preserving their
covariance. The centrality selection must match the analyzed production;
run names alone do not select events. Flow is not separately rebinned
through the experimental-comparison interface, which rebins count
distributions.

DAT distinguishes quantities such as `v_1{RP}` and `v_2{EP}`, with unit `1`.
Definitions retain the reference, n, weighting and empty-bin policy.
See [results.md](results.md) and the [main README](../README.md) for the
format and plotting. Modules: `flow_reaction_plane.py`,
`flow_event_plane.py` and `flow_resolution.py`; `v1.py` and `v2.py` provide
simple interfaces for a known RP.
