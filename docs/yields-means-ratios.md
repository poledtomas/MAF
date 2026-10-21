# Yields, mean values and particle ratios

These observables are scalars: each selection produces one value, so
histogram bins must be omitted. Numeric PDG lists expand by species,
except for ratios with explicit combined numerator and denominator
selections.

## Yield and charged multiplicity

```yaml
observables:
  yields:
    enabled: true
    pdg: [211, -211, 2212]
    y_range: [-0.5, 0.5]

  charged_multiplicity:
    enabled: true
    eta_range: [-1, 1]
```

`yields` is the short name for `kind: yield`. It returns ⟨N⟩, the mean
multiplicity per event within the configured acceptance.
`charged_multiplicity` always selects nonzero electric charge and returns
⟨Nch⟩. Both use unit `1`. A complete event without selected particles
contributes zero and remains in the event count.

For yields, add `density: y` or `density: eta` together with the matching
`y_range` or `eta_range` to divide by the window width. This gives the
average density over that window.

## Mean pT

```yaml
observables:
  mean_pt:
    enabled: true
    pdg: [211, -211, 2212]
    y_range: [-0.5, 0.5]
    pt_range: [0.2, 3.0]
```

The result is the particle-weighted ensemble mean:

```text
⟨pT⟩ = Σevents Σparticles pT / Σevents N.
```

Events with more selected particles contribute proportionally more than
to an equally weighted mean of event means. An event without selected
particles contributes zero to both numerator and denominator; the result
is undefined if the combined particle count is zero. The unit is GeV.

## Mean mT−m₀

```yaml
observables:
  proton_mean_mt:
    kind: mean_mt_minus_m0
    enabled: true
    pdg: [2212]
    mass: 0.938272
    y_range: [-0.5, 0.5]
```

With mT=√(m₀²+pT²), the result is the combined sum of (mT−m₀) divided by
the combined particle count. `mass` must be positive and correspond to
one selected species. The result is in GeV.

## Particle ratios

```yaml
observables:
  kaon_to_pion:
    kind: ratio
    enabled: true
    pdg: [321, -321]
    denominator_pdg: [211, -211]
    y_range: [-0.5, 0.5]
```

This example computes (K⁺+K⁻)/(π⁺+π⁻), the ratio of ensemble yields.
For this module, numerator PDG IDs are combined into one result, as are
denominator PDG IDs. The denominator uses the same remaining cuts as the
numerator. Individual event ratios are not averaged: those ratios would
be undefined for events with a zero denominator. A zero combined
denominator gives `null`/`nan`. The unit is `1`; the example produces
`kaon_to_pion.dat`.

## Uncertainties and export

The jackknife recomputes numerator and denominator jointly, preserving
their covariance. `sampling` defines the independent statistical unit;
see [physics.md](physics.md).

A scalar DAT table contains one row: index 0, `nan` bin boundaries,
the value and its statistical error. The header records the quantity
and unit.

Calculation modules: `observables/yields.py`, `mean_pt.py`, `mean_mt.py`
and `ratios.py`.
