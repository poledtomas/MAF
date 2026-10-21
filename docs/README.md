# MAF documentation

Start with the [main README](../README.md) for installation, analysis and
plotting. The pages below describe individual quantities, their
definitions, configuration and outputs.

## Observables

| Document | Quantities |
|---|---|
| [pT and mT spectra](pt-spectra.md) | dN/dpT, dN/dmT, mT−m₀ and invariant spectra |
| [Rapidity and pseudorapidity](rapidity.md) | dN/dy, dNch/dη, d(Np−Npbar)/dy |
| [Yields, means and ratios](yields-means-ratios.md) | ⟨N⟩, ⟨Nch⟩, ⟨pT⟩, ⟨mT−m₀⟩ and particle ratios |
| [RP and EP flow](flow.md) | vₙ{RP}, vₙ{EP}; y, η, pT; identified and charged particles |
| [Q cumulants](q-cumulants.md) | differential vₙ{2} and vₙ{4} |

## Shared settings

| Document | Contents |
|---|---|
| [Configuration](configuration.md) | inputs, run filters, PDG IDs, cuts, bins and plot styling |
| [Physics and statistics](physics.md) | units, frames, normalization and sampling |
| [Result format](results.md) | DAT, JSON, NPZ, cache and merging |
| [Calculation modules](simple-modules.md) | program structure and adding an observable |
| [Validation](validation.md) | tests, analytical checks and synthetic inputs |

Example configurations are in [examples/](../examples/):
`observables.yaml` / `.json` demonstrate several calculations, while
`flow_rapidity.yaml` / `.json` provide four RP/EP combinations with y/η axes.
Adjust the input and run selection for your production.
