from pathlib import Path
import numpy as np
from maf import Event, Particles

HEADER = "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n# SMASH-test\n"


def row(pt=1., pz=0., mass=1., pdg=2212, charge=1, ident=0):
    return f"0 0 0 0 {mass:.16g} {np.sqrt(mass**2+pt**2+pz**2):.16g} {pt:.16g} 0 {pz:.16g} {pdg} {ident} {charge}\n"


def block(event_id, rows, *, ensemble=0, end=True, weight=None):
    result=f"# event {event_id} ensemble {ensemble} out {len(rows)}\n"+"".join(rows)
    if end:
        result+=f"# event {event_id} ensemble {ensemble} end 0 impact 1.0"
        if weight is not None:
            result+=f" weight {weight}"
        result+="\n"
    return result


def write_run(root,name,events):
    run=Path(root)/name
    run.mkdir(parents=True,exist_ok=True)
    (run/"particle_lists.oscar").write_text(HEADER+"".join(events))
    return run


def event(count=1, *, pt=1., pz=0., mass=1., pdg=2212, charge=1):
    def a(value):
        return np.broadcast_to(np.asarray(value), (count,)).copy()
    return Event(Particles(a(np.sqrt(np.asarray(mass)**2+np.asarray(pt)**2+np.asarray(pz)**2)),a(pt),a(0.),a(pz),a(mass),a(pdg),a(charge)),0)
