# Rapidity, pseudorapidity and net protons

Rapidity y=½ ln[(E+pz)/(E−pz)] and pseudorapidity η=asinh(pz/pT)
are different quantities. In general, y≠η for massive particles. Their
frames and behavior at zero pT are covered by the
[shared physics conventions](physics.md).

## dN/dy

```yaml
observables:
  dndy:
    enabled: true
    pdg: [211, -211, 321, -321, 2212, -2212]
    y_bins: {min: -3, max: 3, count: 60}
```

The module counts particles in each bin of each event. The result is
⟨Hᵢ⟩/Δy, or dN/dy per event. Both y and the density are dimensionless;
the DAT unit is `1`. Optional `pt_range` or `eta_range` cuts restrict the
particle population. No η cut is applied automatically.

Numeric PDG lists produce separate tables, such as `dndy_2212.dat`.
Use `pdg: [all_charged]` for combined charged dN/dy, or omit PDG to
include all species, including neutral particles.

## dNch/dη

```yaml
observables:
  dndeta:
    enabled: true
    eta_bins: {min: -5, max: 5, count: 50}
```

`dndeta` always selects particles with `charge != 0`. It returns ⟨Hᵢ⟩/Δη,
or charged dNch/dη per event, with unit `1`. An optional numeric PDG list
further restricts that charged selection. `charged: false` is invalid.
Infinite or undefined η values do not enter the histogram.

## Net-proton dN/dy

```yaml
observables:
  net_proton:
    enabled: true
    name: net_proton_dndy
    pt_range: [0.4, 2.0]
    y_bins: {min: -2, max: 2, count: 40}
```

The module selects protons (PDG 2212) and antiprotons (PDG −2212)
internally; omit `pdg`. Within each event, it subtracts the antiproton
histogram from the proton histogram, then normalizes by the combined
event count and rapidity bin width:

```text
d(Np−Npbar)/dy = ⟨Hp,i−Hpbar,i⟩ / Δy.
```

This preserves proton–antiproton covariance in the statistical error.
Net protons are distinct from net baryons: neutrons and other baryons
do not contribute. The example produces `net_proton_dndy.dat`, with unit `1`.

## Normalization and uncertainties

Complete empty events and events without the selected species contribute
zero and remain in the normalization. Contributions are first summed
over all selected files in one ensemble; normalized distributions from
individual files are not averaged with equal file weights. Statistical
errors and covariance follow the [sampling hierarchy](physics.md).

Bins use [lower, upper). Values are divided by the actual width of each
bin. `count` is the number of bins; there is one more edge than bins.
See [results.md](results.md) for DAT and plotting details.
Modules: `observables/dndy.py`, `dndeta.py` and `net_proton.py`.
