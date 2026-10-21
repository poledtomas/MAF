# Flow from two- and four-particle Q cumulants

The `cumulant_2` and `cumulant_4` methods return differential vₙ{2} and
vₙ{4} versus y, η or pT. They use direct Q-vector correlations following
[Bilandzic, Snellings and Voloshin, arXiv:1010.0233v2](https://arxiv.org/html/1010.0233v2),
Eqs. (16), (18), (24), (25), (28), (32), (34) and (35).
The implementation assumes unit particle weights and uniform azimuthal
acceptance.

## Configuring v₂{2}(pT)

```yaml
observables:
  v2_cumulant2:
    kind: v2
    enabled: true
    pdg: [all_charged]
    axis: pt
    eta_range: [-0.8, 0.8]
    pt_range: [0.2, 2.0]
    pt_bins: {min: 0.2, max: 2.0, count: 18}
    flow:
      method: cumulant_2
      weighting: tuple
      reference:
        pdg: [all_charged]
        eta: [-0.8, 0.8]
        pt: [0.2, 3.0]
```

Particles of interest (POI) are selected for each differential bin;
reference flow particles (RFP) define the reference correlations. The
reference must be explicit. Numeric PDG IDs within it form one combined
selection. POI and RFP can be identical, partially overlapping or disjoint.
Autocorrelations are subtracted using their actual intersection.

`weighting: tuple` means combinatorial event weighting. It does not assign
pT or detector-efficiency weights to individual particles. `event` and
`particle` weighting are not supported for cumulants.

## Configuring proton v₂{4}(y)

```yaml
observables:
  proton_v2_cumulant4_y:
    kind: v2
    enabled: true
    pdg: [2212]
    axis: y
    pt_range: [0.2, 2.0]
    y_bins: {min: -1, max: 1, count: 10}
    flow:
      method: cumulant_4
      weighting: tuple
      reference:
        pdg: [211, -211]
        eta: [-0.8, 0.8]
        pt: [0.2, 3.0]
```

This example uses disjoint proton POI and pion RFP. For identical
selections, match both species and cuts. Higher harmonics use `kind: vn`
and `harmonic: n`, where n is any positive integer.

## Correlations and weights

Each event contributes a raw numerator and the number of actual ordered
tuples with distinct particles. Let M be the number of RFP, m_p the number
of POI in the bin and m_q the size of their intersection. Denominators are:

```text
reference 2:       M(M−1)
reference 4:       M(M−1)(M−2)(M−3)
differential 2:    m_p M−m_q
differential 4:    (m_p M−3m_q)(M−1)(M−2).
```

A zero tuple weight contributes neither numerator nor denominator.
Reference pairs require M≥2 and quadruplets require M≥4. Differential
quadruplets can exist with three RFP and one distinct POI. Contributions
are summed over the ensemble before finalization; individual event vₙ
values or cumulants are not averaged.

## Final estimators

Double brackets denote correlations formed from combined weighted sums;
a prime denotes a differential correlation with a POI:

```text
c_n{2} = <<2>>
d_n{2}(bin) = <<2'>>
v_n{2}(bin) = d_n{2}(bin) / sqrt(c_n{2})

c_n{4} = <<4>> − 2 <<2>>²
d_n{4}(bin) = <<4'>> − 2 <<2'>> <<2>>
v_n{4}(bin) = −d_n{4}(bin) / (−c_n{4})^(3/4).
```

vₙ{2} requires cₙ{2}>0; vₙ{4} requires cₙ{4}<0. An invalid sign or missing
tuples gives an undefined result: `null` in JSON and `nan` in DAT.
Absolute values are not used to force a real result.

## Statistics and output

The jackknife over independent units defined by `sampling` recomputes
reference and differential correlations jointly, preserving their
covariance. The resulting vₙ is dimensionless and is not divided by
bin width.

DAT headers record quantities such as `v_2{2}` or `v_2{4}`, unit `1` and
the full selection definition. In `summary.json`, `auxiliary` contains
reference c₂ and, for the four-particle method, c₄. NPZ retains raw
contributions for merging.

Example outputs are `v2_cumulant2_all_charged.dat` and
`proton_v2_cumulant4_y_2212.dat`. See [results.md](results.md) and
[physics.md](physics.md) for the format and normalization conventions.
The calculation module is `observables/flow_cumulants.py`; tests verify
overlap and combinatorics by explicitly enumerating pairs and quadruplets.
