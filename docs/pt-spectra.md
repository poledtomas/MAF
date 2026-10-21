# pT and mT spectra

These modules histogram particles in individual events and normalize their
combined contributions. Each species in a numeric `pdg` list receives a
separate spectrum in the short configuration format. Kinematic cuts are
independent settings. All momenta and masses use GeV with c=1.

## pT spectrum

```yaml
observables:
  pt_spectra:
    enabled: true
    pdg: [211, -211, 2212]
    y_range: [-0.5, 0.5]
    pt_bins: {min: 0.0, max: 3.0, count: 30}
    invariant: false
```

`invariant: false` corresponds to `kind: pt`. For a bin count Hᵢ, the result
is ⟨Hᵢ⟩/ΔpT: dN/dpT per event within the selected acceptance, in GeV⁻¹.
`y_range` selects particles but does not divide the result by the rapidity
window width. Add `density: y` to a non-invariant spectrum for d²N/(dpT dy).

## Invariant pT spectrum

```yaml
observables:
  proton_invariant_pt:
    kind: pt
    enabled: true
    pdg: [2212]
    y_range: [-0.1, 0.1]
    pt_range: [0.2, 2.0]
    pt_bins: {min: 0.2, max: 2.0, count: 36}
    boundary: open
    invariant: true
```

`invariant: true` corresponds to `kind: invariant_pt` and requires an
explicit `y_range`. The bin value is

```text
⟨Hᵢ⟩ / [π (pT_high² − pT_low²) Δy].
```

This is a finite-bin average of the invariant yield
(1/2πpT) d²N/(dpT dy) with the transverse annular measure. At the arithmetic
bin center, the denominator equals 2π pT_center ΔpT Δy. The unit is GeV⁻²;
the normalization is finite even for a bin starting at pT=0. Do not add a
separate `density` normalization to an invariant spectrum.

## mT and mT−m₀ spectra

Transverse mass is mT=√(m₀²+pT²). Set `mass` in GeV and select one particle
species, optionally including its antiparticle. `mass` is the explicit
reference mass used in this transformation and must be positive.

```yaml
observables:
  proton_mt:
    kind: mt
    enabled: true
    pdg: [2212]
    mass: 0.938272
    y_range: [-0.5, 0.5]
    mt_bins: [0.94, 1.0, 1.2, 1.5, 2.0]

  proton_mt_minus_m0:
    kind: mt_minus_m0
    enabled: true
    pdg: [2212]
    mass: 0.938272
    y_range: [-0.5, 0.5]
    bins: {min: 0.0, max: 2.0, count: 20}
```

`mt_spectra` is the short name for `kind: mt`. Ordinary spectra give
dN/dmT or dN/d(mT−m₀) per event, in GeV⁻¹. With `kind: mt`, use
`invariant: true`, or select `kind: invariant_mt` directly, to obtain
⟨Hᵢ⟩/[π(mT_high²−mT_low²)Δy] in GeV⁻². mT bins must not start below
the configured reference mass.

## Bins, uncertainties and files

Bins use [lower, upper); the final upper edge is excluded. Bins must not
extend beyond an explicit acceptance on the same axis. For mT, a pT cut
is checked against the transformed bin edges. Nonuniform bins use their
actual width or annular area.

A complete event with no selected particles contributes a zero histogram
and remains in the event count. Statistical errors use the joint event
contributions and the declared [sampling hierarchy](physics.md), preserving
fluctuations and covariance.

Example filenames are `pt_spectra_211.dat` and
`proton_invariant_pt_2212.dat`. See the [result format](results.md) for
export details. Calculation modules: `observables/pt_spectra.py` and
`observables/mt_spectra.py`.
