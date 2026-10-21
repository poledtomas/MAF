"""Reference acceptances, tuple enumeration and nonlinear independent-unit errors."""
import copy
import itertools
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from maf import Event,Particles,Observable
from maf.config import expand_observables
from maf.core.selection import Selection
from maf.observables import flow_cumulants,flow_event_plane
from maf.statistics.aggregation import Result
from maf.statistics.errors import jackknife_derived


def event_with_angles(phi, pt=None, pdg=None):
    phi=np.asarray(phi,float)
    pt=np.ones(len(phi)) if pt is None else np.asarray(pt,float)
    return Event(Particles(np.hypot(1,pt),pt*np.cos(phi),pt*np.sin(phi),np.zeros(len(phi)),np.ones(len(phi)),
                           np.full(len(phi),2212) if pdg is None else np.asarray(pdg,int),np.ones(len(phi),int)),0)


def require_scipy(test):
    try:
        import scipy
    except ImportError:
        test.skipTest('Optional SciPy required for event-plane resolution')


class ReferenceConventions(unittest.TestCase):
    def test_exact_spectrum_and_charged_acceptance_preset(self):
        path=Path(__file__).resolve().parents[1]/'examples/observables.json'
        config=json.loads(path.read_text())
        definitions=expand_observables(config['observables'])
        spectra=[Observable(**d) for d in definitions if d['kind']=='invariant_pt']
        self.assertEqual([s.selection['pdg'] for s in spectra],[[211],[321],[2212]])
        for spectrum in spectra:
            self.assertEqual(spectrum._cut.y,[-.1,.1])
            np.testing.assert_allclose(spectrum.bins,np.linspace(.2,2,37))
            centers=(np.asarray(spectrum.bins[:-1])+spectrum.bins[1:])/2
            np.testing.assert_allclose(spectrum.normalization(),1/(2*np.pi*centers*.05*.2))
        charged=Observable(**next(d for d in definitions if d['kind']=='dndeta'))
        self.assertTrue(charged._cut.charged)
        np.testing.assert_allclose(charged.bins,np.linspace(-5,5,51))

    def test_open_acceptance_with_unchanged_default(self):
        event=event_with_angles([0,0,0],pt=[.2,.3,2])
        np.testing.assert_array_equal(Selection(pt=[.2,2],boundary='open').mask(event.particles),[False,True,False])
        np.testing.assert_array_equal(Selection(pt=[.2,2]).mask(event.particles),[True,True,False])

    def test_reaction_plane_event_vs_particle_weighting(self):
        events=[event_with_angles([0]),event_with_angles([np.pi]*3)]
        estimates=[]
        for weighting in ('event','particle'):
            obs=Observable('flow','v1',bins=[-1,1],flow={'method':'reaction_plane','angle':0,'weighting':weighting})
            result=Result([obs],'same','{}')
            for i,event in enumerate(events):
                result.add(str(i),'event',event,[obs.evaluate(event)])
            estimates.append(result.summary()['observables']['flow']['value'][0])
        np.testing.assert_allclose(estimates,[0,-.5],atol=1e-15)

    def test_invalid_uniform_binning_and_weighting(self):
        for bins in ({'min':0,'max':1,'count':0},{'min':0,'max':1,'count':2.5},{'min':0,'max':1}):
            with self.assertRaises(ValueError):
                expand_observables({'pt_spectra':{'pt_bins':bins}})
        with self.assertRaises(ValueError):
            Observable('v','v2',bins=[0,2],flow={'method':'reaction_plane','angle':0,'weighting':'unknown'})


class CumulantTupleTests(unittest.TestCase):
    def test_q_vectors_against_all_distinct_ordered_tuples(self):
        rng=np.random.default_rng(915)
        for multiplicity in range(7):
            phi=rng.uniform(-np.pi,np.pi,multiplicity)
            coordinate=np.arange(multiplicity)/max(1,multiplicity)
            reference=np.ones(multiplicity,bool)
            particles=np.arange(multiplicity)%2==0
            bins=[0,.5,1]
            for harmonic in (1,2):
                for four in (False,True):
                    a,b=flow_cumulants.calculate(phi,coordinate,reference,particles,bins,harmonic,four)
                    offset=2 if four else 1
                    for order,feature in ((2,0),(4,1)) if four else ((2,0),):
                        tuples=list(itertools.permutations(range(multiplicity),order))
                        def cosine(indices):
                            half=order//2
                            angle=sum(phi[i] for i in indices[:half])-sum(phi[i] for i in indices[half:])
                            return np.cos(harmonic*angle)
                        self.assertAlmostEqual(b[feature],len(tuples))
                        self.assertAlmostEqual(a[feature],sum(map(cosine,tuples)),places=9)
                        for index,(low,high) in enumerate(zip(bins[:-1],bins[1:])):
                            chosen=[t for t in tuples if particles[t[0]] and low<=coordinate[t[0]]<high]
                            position=offset+index+(2 if order==4 else 0)
                            self.assertAlmostEqual(b[position],len(chosen))
                            self.assertAlmostEqual(a[position],sum(map(cosine,chosen)),places=9)

    def test_four_particle_weights_do_not_overflow_32_bit(self):
        count=1000
        a,b=flow_cumulants.calculate(np.zeros(count),np.zeros(count),np.ones(count,bool),np.ones(count,bool),[-1,1],2,True)
        self.assertEqual(b[1],count*(count-1)*(count-2)*(count-3))
        self.assertGreater(b[1],2**31)
        np.testing.assert_allclose(flow_cumulants.finalize(a,b,True),[1])

    def test_fully_aligned_cumulants_equal_one(self):
        for four in (False,True):
            a,b=flow_cumulants.calculate(np.zeros(6),np.zeros(6),np.ones(6,bool),np.ones(6,bool),[-1,1],2,four)
            np.testing.assert_allclose(flow_cumulants.finalize(a,b,four),[1])

    def test_invalid_cumulant_signs_and_missing_tuples(self):
        self.assertTrue(np.isnan(flow_cumulants.finalize([-.5,.2],[1,1])[0]))
        self.assertTrue(np.isnan(flow_cumulants.finalize([.1,.2,.1,.1],[1,1,1,1],four=True)[0]))
        for count in range(4):
            a,b=flow_cumulants.calculate(np.zeros(count),np.zeros(count),np.ones(count,bool),np.ones(count,bool),[-1,1],2,True)
            self.assertEqual(b[1],0)
            self.assertTrue(np.isfinite(a).all())
            self.assertTrue(np.isfinite(b).all())

    def test_reference_overlap_and_bin_selection(self):
        a,b=flow_cumulants.calculate(np.zeros(3),np.zeros(3),np.array([True,True,False]),np.ones(3,bool),[-1,1],2)
        np.testing.assert_allclose(a,b)
        self.assertEqual(b[1],4)
        obs=Observable('v','v2',axis='pt',bins=[.2,2],flow={'method':'cumulant_2','reference':{'pt':[.2,2]}})
        # A particle outside the plotted bins is not a particle of interest.
        a,b=obs.evaluate(event_with_angles([0,0,0],pt=[.1,.5,1]))
        np.testing.assert_allclose(obs.estimate(a,b),[1])


class EventPlaneTests(unittest.TestCase):
    def test_autocorrelation_removed_before_azimuthal_mean(self):
        phi=np.array([0,np.pi/4,0,np.pi/4])
        a,b=flow_event_plane.calculate(phi,np.ones(4),np.zeros(4),np.ones(4,bool),np.ones(4,bool),[-1,1],2)
        self.assertAlmostEqual(a[0]/b[0],1/np.sqrt(5))
        self.assertNotAlmostEqual(a[0]/b[0],1/np.sqrt(2))

    def test_single_particle_and_empty_subevents_excluded(self):
        a,b=flow_event_plane.calculate(np.zeros(1),np.ones(1),np.zeros(1),np.ones(1,bool),np.ones(1,bool),[-1,1],2)
        np.testing.assert_array_equal(a,[0,0])
        np.testing.assert_array_equal(b,[0,0])

    def test_resolution_counts_do_not_depend_on_first_bin_occupancy(self):
        a,b=flow_event_plane.calculate(np.zeros(4),np.ones(4),np.full(4,1.5),np.ones(4,bool),np.ones(4,bool),[0,1,2],2)
        np.testing.assert_array_equal(b,[0,1,1])
        np.testing.assert_allclose(a,[0,1,1])

    def test_resolution_inversion_and_undefined_values(self):
        require_scipy(self)
        from maf.observables.flow_resolution import resolution,full_resolution
        chi=np.array([.1,1,3,10])
        np.testing.assert_allclose(full_resolution(resolution(chi)**2),resolution(np.sqrt(2)*chi),rtol=1e-10)
        self.assertTrue(np.isnan(full_resolution(np.array([-1,0,np.nan]))).all())
        self.assertEqual(full_resolution(1),1)

    def test_ep_estimator_and_resolution_uncertainty(self):
        require_scipy(self)
        from maf.observables.flow_resolution import full_resolution
        a=np.array([[.8,.7],[.9,.8],[.7,.6],[.85,.9]])
        b=np.ones_like(a)
        point,error,cov=jackknife_derived(a,b,flow_event_plane.finalize)
        manual=[]
        for i in range(4):
            sample=np.delete(a,i,axis=0).mean(axis=0)
            manual.append(sample[0]/full_resolution(sample[1]))
        manual=np.asarray(manual)
        self.assertAlmostEqual(error[0],np.sqrt(3/4*np.sum((manual-manual.mean())**2)))
        self.assertGreater(error[0],0)

    def test_saving_and_merging_all_flow_methods(self):
        require_scipy(self)
        rng=np.random.default_rng(7)
        for method in ('event_plane','cumulant_2','cumulant_4'):
            obs=Observable('v','v2',bins=[-1,1],flow={'method':method,'reference':{}})
            left=Result([obs],'same','{}',sources={'left':{'metadata':{}}})
            right=Result([copy.deepcopy(obs)],'same','{}',sources={'right':{'metadata':{}}})
            direct=Result([copy.deepcopy(obs)],'same','{}')
            for i in range(8):
                event=event_with_angles(rng.normal(0,.2,6+i%3))
                values=[obs.evaluate(event)]
                direct.add(str(i),'event',event,values)
                (left if i<4 else right).add(str(i),'event',event,values)
            left.merge(right)
            for field in ('value','stat_error','covariance'):
                np.testing.assert_allclose(left.summary()['observables']['v'][field],direct.summary()['observables']['v'][field],atol=1e-12)
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'result.npz'
                left.save(path)
                self.assertEqual(left.summary(),Result.load(path).summary())

    def test_reference_required_for_nonlinear_methods(self):
        for method in ('event_plane','cumulant_2','cumulant_4'):
            with self.assertRaisesRegex(ValueError,'reference'):
                Observable('v','v2',bins=[-1,1],flow={'method':method})

    def test_plot_refuses_old_pt_and_invariant_pt_mixture(self):
        from maf.plotting.figures import distribution
        integrated={'observables':{'pt':{'definition':{'kind':'pt'},'units':'GeV^-1'}}}
        invariant={'observables':{'pt':{'definition':{'kind':'invariant_pt'},'units':'GeV^-2'}}}
        with self.assertRaisesRegex(ValueError,'different normalizations'):
            distribution([integrated,invariant],'pt')
