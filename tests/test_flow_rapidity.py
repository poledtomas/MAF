"""Rapidity presets checked with scalar sums independent of the NumPy modules."""
import json
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest
import numpy as np
from maf import Event, Particles, Observable
from maf.config import expand_observables, load_job
from maf.io.dat import read_dat
from maf.observables import flow_event_plane, flow_reaction_plane
from maf.statistics.aggregation import Result

ROOT = Path(__file__).resolve().parents[1]


def sample(index=0):
    eta = np.array([-2.1, -1.8, -1.2, -.9, -.5, -.1, .1, .4, .9, 1.2, 1.8, 2.1])
    pt = np.array([.3, .4, .5, .6, .8, 1, 1.2, 1.5, 1.55, 1.9, 2.1, .1])
    phi = np.linspace(-.15, .35, len(pt)) + .07*index
    phi += .04*np.sin(np.arange(len(pt))*(index+1))
    mass = np.array([.938, .139, .938, .139, .494, .938, .938, .139, .494, .938, .139, .938])
    pdg = np.array([2212, 211, 2112, -211, 321, 2212, -2212, 111, -321, 2212, 211, 2212])
    charge = np.array([1, 1, 0, -1, 1, 1, -1, 0, -1, 1, 1, 1])
    pz = pt*np.sinh(eta)
    return Event(Particles(np.sqrt(mass**2+pt**2+pz**2), pt*np.cos(phi), pt*np.sin(phi),
                           pz, mass, pdg, charge), index)


def scalar_ep(event, obs):
    """Literal per-particle angles and sums; subtract only actual reference members."""
    p = event.particles
    n = obs.harmonic
    ref_eta = 2 if obs.axis == 'eta' else 1
    upper_pt = 1.6 if 'identified' in obs.name else 2
    phi = [math.atan2(y, x) for x, y in zip(p.px, p.py)]
    pt = [math.hypot(x, y) for x, y in zip(p.px, p.py)]
    eta = [math.asinh(z/t) for z, t in zip(p.pz, pt)]
    rapidity = [.5*math.log((e+z)/(e-z)) for e, z in zip(p.energy, p.pz)]
    ref = [i for i in range(len(p)) if abs(eta[i]) < ref_eta and .2 < pt[i] < upper_pt]
    def vector(indices):
        return sum(pt[i]*complex(math.cos(n*phi[i]), math.sin(n*phi[i])) for i in indices)
    a, b = vector([i for i in ref if i % 2 == 0]), vector([i for i in ref if i % 2])
    numerator, denominator = np.zeros(obs.size+1), np.zeros(obs.size+1)
    if abs(a) < 1e-12 or abs(b) < 1e-12:
        return numerator, denominator
    numerator[-1] = math.cos(math.atan2(a.imag, a.real)-math.atan2(b.imag, b.real))
    denominator[:] = 1  # ROOT rapidity EP: empty bins contribute zero.
    totals = [[] for _ in range(obs.size)]
    for i in range(len(p)):
        if not (.2 < pt[i] < upper_pt and p.charge[i] != 0):
            continue
        if obs.selection.get('pdg') and p.pdg[i] not in obs.selection['pdg']:
            continue
        x = eta[i] if obs.axis == 'eta' else rapidity[i]
        if not obs.bins[0] < x < obs.bins[-1]:
            continue
        plane = vector([j for j in ref if j != i])
        if abs(plane) < 1e-12:
            continue
        for k, (low, high) in enumerate(zip(obs.bins[:-1], obs.bins[1:])):
            if low <= x < high:
                totals[k].append(math.cos(n*phi[i]-math.atan2(plane.imag, plane.real)))
    for k, values in enumerate(totals):
        if values:
            numerator[k] = sum(values)/len(values)
    return numerator, denominator


class RapidityFlowTests(unittest.TestCase):
    def test_general_rp_harmonic_and_rotated_known_plane(self):
        px, py = np.cos([.1, .4, 1]), np.sin([.1, .4, 1])
        for n in (1, 2, 3, 5):
            obs = Observable('v', 'vn', harmonic=n, bins=[-1, 1],
                             flow={'method':'reaction_plane', 'angle':.3})
            a, b = flow_reaction_plane.calculate(px, py, [0, 0, 0], [-1, 1], .3, n)
            self.assertAlmostEqual(obs.estimate(a, b)[0], sum(math.cos(n*(x-.3)) for x in [.1, .4, 1])/3)
            for alias in ('v1', 'v2'):
                old = Observable('v', alias, bins=[-1, 1], flow={'method':'reaction_plane','angle':.3})
                self.assertEqual(old.harmonic, int(alias[1:]))

    def test_invalid_harmonics_and_empty_bin_settings(self):
        for n in (None, 0, -1, 1.5, True, '2'):
            with self.assertRaisesRegex(ValueError, 'harmonic'):
                Observable('v', 'vn', harmonic=n, bins=[-1, 1], flow={'method':'reaction_plane','angle':0})
        with self.assertRaisesRegex(ValueError, 'harmonic'):
            Observable('v', 'v1', harmonic=2, bins=[-1, 1], flow={'method':'reaction_plane','angle':0})
        for extra in ({'empty_bins':'wrong'}, {'empty_bins':'zero','weighting':'particle'}):
            with self.assertRaises(ValueError):
                Observable('v', 'vn', harmonic=2, bins=[-1, 1],
                           flow=dict(method='event_plane', reference={}, **extra))

    @unittest.skipUnless(importlib.util.find_spec('scipy'), 'Optional SciPy required')
    def test_reference_presets_and_scalar_ep_for_multiple_harmonics(self):
        document = json.loads((ROOT/'examples/flow_rapidity.json').read_text())
        definitions = expand_observables(document['observables'])
        self.assertEqual(len(definitions), 14)  # six RP, six EP species, two inclusive EP.
        for n in (1, 2, 3):
            for definition in definitions:
                definition = dict(definition, harmonic=n)
                obs = Observable(**definition)
                event = sample(n)
                a, b = obs.evaluate(event)
                if obs.flow['method'] == 'event_plane':
                    expected_a, expected_b = scalar_ep(event, obs)
                    np.testing.assert_allclose(a, expected_a, atol=1e-13)
                    np.testing.assert_array_equal(b, expected_b)
                    self.assertNotIn('charged', obs.flow['reference'])
                    self.assertNotIn('eta', obs.selection if obs.axis == 'y' else {})
                else:
                    p = event.particles
                    values = [math.cos(n*math.atan2(p.py[i],p.px[i])) for i in range(len(p))
                              if p.pdg[i] in obs.selection['pdg'] and abs(p.eta[i]) < 1
                              and .2 < p.pt[i] < 1.6 and -1 < p.rapidity[i] < 1]
                    total = a.sum()
                    # At most one test particle of each identified species occurs in each bin.
                    self.assertAlmostEqual(total, sum(values), places=12)

    def test_ep_empty_bins_zero_differs_from_occupied_only(self):
        args = (np.zeros(4), np.ones(4), np.full(4, .5), np.ones(4,bool), np.ones(4,bool), [-1,0,1], 2)
        a, b = flow_event_plane.calculate(*args, empty_bins='zero')
        np.testing.assert_array_equal(a, [0,1,1])
        np.testing.assert_array_equal(b, [1,1,1])
        a, b = flow_event_plane.calculate(*args)  # historical default retained
        np.testing.assert_array_equal(b, [0,1,1])
        args = (np.array([0.]), np.ones(1), np.zeros(1), np.ones(1,bool), np.ones(1,bool), [-1,0,1], 2)
        a, b = flow_event_plane.calculate(*args, empty_bins='zero')
        np.testing.assert_array_equal(b, [0,0,0])  # undefined subevent plane is not angle zero

    def test_y_and_eta_differ_for_massive_hadrons(self):
        event = sample()
        # Massive hadrons have different y and eta; the axes must not be interchanged.
        p = event.particles
        self.assertGreater(np.max(np.abs(p.eta-p.rapidity)), .3)
        for axis, expected in [('y',p.rapidity), ('eta',p.eta)]:
            obs = Observable('v', 'vn', harmonic=3, axis=axis, bins=[-2,-.5,0,.5,2],
                             flow={'method':'reaction_plane','angle':0})
            a, b = obs.evaluate(event)
            manual = flow_reaction_plane.calculate(p.px,p.py,expected,obs.bins,0,3)
            np.testing.assert_allclose(a, manual[0])
            np.testing.assert_allclose(b, manual[1])

    @unittest.skipUnless(importlib.util.find_spec('scipy'), 'Optional SciPy required')
    def test_roundtrip_dat_harmonic_method_and_old_checkpoint(self):
        obs = Observable('v3_ep_y', 'vn', harmonic=3, bins=[-1,0,1], selection={'charged':True},
                         flow={'method':'event_plane','reference':{'eta':[-1,1]},'empty_bins':'zero'})
        result = Result([obs], 'same', '{}')
        for i in range(5):
            event = sample(i)
            result.add(str(i),'event',event,[obs.evaluate(event)])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            result.save(path/'result.npz')
            loaded = Result.load(path/'result.npz')
            self.assertEqual(result.summary(), loaded.summary())
            loaded.export(path)
            _, data = read_dat(path/'v3_ep_y.dat')
            self.assertEqual(data['axes']['x'], {'quantity':'y','unit':'1'})
            self.assertEqual(data['axes']['y'], {'quantity':'v_3{EP}','unit':'1'})
            self.assertEqual(data['definition']['harmonic'],3)
            self.assertEqual(data['definition']['flow']['empty_bins'],'zero')
            # Real 0.4.1 definitions did not contain harmonic or empty_bins.
            oldobs = Observable('v', 'v2', bins=[-1,1], flow={'method':'event_plane','reference':{}})
            old = Result([oldobs], 'same', '{}')
            old.save(path/'old.npz')
            with np.load(path/'old.npz',allow_pickle=False) as archive:
                arrays = {k:archive[k].copy() for k in archive.files}
            metadata = json.loads(str(arrays['metadata']))
            metadata['version'] = '0.4.1'
            metadata['observables'][0].pop('harmonic')
            metadata['observables'][0]['flow'].pop('empty_bins')
            arrays['metadata'] = json.dumps(metadata)
            np.savez(path/'old.npz',**arrays)
            restored = Result.load(path/'old.npz').observables[0]
            self.assertEqual(restored.harmonic,2)
            self.assertEqual(restored.flow['empty_bins'],'exclude')

    @unittest.skipUnless(importlib.util.find_spec('yaml'), 'Optional YAML required')
    def test_yaml_and_json_examples_agree_and_main_configuration_is_valid(self):
        import yaml
        self.assertEqual(yaml.safe_load((ROOT/'examples/flow_rapidity.yaml').read_text()),
                         json.loads((ROOT/'examples/flow_rapidity.json').read_text()))
        load_job(ROOT/'analysis.yaml')


if __name__ == '__main__':
    unittest.main()
