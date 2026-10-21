"""Create four tiny synthetic runs; no production or experimental data used."""
import argparse
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    centrality = {"type": "synthetic", "range": [0, 5], "unit": "percent"}
    manifest = {"defaults": {"sampling": {"mode": "event"}, "centrality": centrality}, "runs": {}}
    header = ("#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n"
              "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n# SMASH-synthetic\n")
    for parameter_index, strength in enumerate((0.2, 0.4)):
        for run_index in (1, 2):
            name = f"7.7_AuAu_demo{parameter_index}_run_{run_index}"
            path = root / "runs" / name
            path.mkdir(parents=True, exist_ok=True)
            blocks = [header]
            for event_id, count in enumerate((0, 2 + parameter_index, 4 + 2*parameter_index)):
                particles = []
                for ident in range(count):
                    particles.append(f"0 0 0 0 0.938 {math.hypot(0.938,0.4):.16g} 0.4 0 0 2212 {ident} 1\n")
                blocks.append(f"# event {event_id} ensemble 0 out {count}\n")
                blocks.extend(particles)
                blocks.append(f"# event {event_id} ensemble 0 end 0 impact 1.0 scattering_projectile_target yes\n")
            (path / "particle_lists.oscar").write_text("".join(blocks))
            manifest["runs"][name] = {"parameters": {"strength": strength}, "tags": ["synthetic-only"]}
    selection = {"pdg": [2212], "frame": "cm", "input_frame": "cm", "particle_definition": "inclusive_final"}
    definition = {"name": "proton_dndy", "kind": "dndy", "bins": [-1, 0, 1], "selection": selection, "units": "1"}
    config = {"group_by": ["energy", "collision_system", "parameters.strength"], "reader": {"errors": "strict"},
              "observables": [definition, {"name": "proton_mean_pt", "kind": "mean_pt", "selection": selection},
                              {"name": "proton_invariant_pt", "kind": "invariant_pt", "bins": [0, 0.5, 1, 2], "selection": dict(selection, y=[-0.5, 0.5])}]}
    (root / "analysis.json").write_text(json.dumps(config, indent=2) + "\n")
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (root / "analytic.csv").write_text("low,high,value,stat\n-1,0,0,0.1\n0,1,2,0.1\n")
    metadata = {"experiment": "Synthetic analytic validation only", "reference": "Analytic fixture; not measured experimental data",
                "collision_system": "Au+Au", "energy": {"value": 7.7, "unit": "GeV"}, "centrality": centrality, "definition": definition}
    (root / "analytic-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(root)


if __name__ == "__main__":
    main()
