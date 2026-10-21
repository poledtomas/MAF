"""The all_charged alias selects by actual electric charge across APIs."""
import copy
import tempfile
from pathlib import Path
import unittest
import numpy as np
from maf import Event,Particles,Observable
from maf.config import expand_observables
from maf.core.selection import Selection
from maf.statistics.aggregation import Result


def particles():
    pt=np.array([.5,.6,.7,.8,.9,.95])
    return Particles(np.hypot(1,pt),pt,np.zeros(6),np.zeros(6),np.ones(6),
                     np.array([2212,-2212,2112,111,1000020040,999999]),np.array([1,-1,0,0,2,-1]))


class AllChargedTests(unittest.TestCase):
    def test_selection_uses_charge_including_antiparticles_and_nuclei(self):
        p=particles()
        for alias in ('all_charged',['all_charged']):
            selection=Selection(pdg=alias,pt=[.55,1])
            np.testing.assert_array_equal(selection.mask(p),[False,True,False,False,True,True])
        np.testing.assert_array_equal(Selection().mask(p),np.ones(6,bool))
        np.testing.assert_array_equal(Selection(charged=True).mask(p),[True,True,False,False,True,True])

    def test_short_config_combines_charges_and_keeps_species_expansion(self):
        for alias in ('all_charged',['all_charged']):
            config={'dndy':{'pdg':alias,'y_bins':[-1,1]}}
            original=copy.deepcopy(config)
            expanded=expand_observables(config)
            self.assertEqual(config,original)
            self.assertEqual(len(expanded),1)
            self.assertEqual(expanded[0]['name'],'dndy_all_charged')
            self.assertEqual(expanded[0]['selection'],{'charged':True})
            obs=Observable(**expanded[0])
            a,b=obs.evaluate(Event(particles(),0))
            np.testing.assert_array_equal(a,[4])
        numeric=expand_observables({'dndy':{'pdg':[211,-211],'y_bins':[-1,1]}})
        self.assertEqual([o['name'] for o in numeric],['dndy_211','dndy_-211'])

    def test_observable_and_flow_reference_normalize_and_round_trip(self):
        flow={'method':'cumulant_4','reference':{'pdg':['all_charged']}}
        obs=Observable('v4','v2',axis='pt',bins=[.2,1],selection={'pdg':['all_charged']},flow=flow)
        self.assertEqual(obs.selection,{'charged':True})
        self.assertEqual(obs.flow['reference'],{'charged':True})
        self.assertEqual(flow['reference'],{'pdg':['all_charged']})
        event=Event(particles(),0)
        a,b=obs.evaluate(event)
        np.testing.assert_allclose(a,[12,24,12,24])
        np.testing.assert_array_equal(a,b)
        np.testing.assert_allclose(obs.estimate(a,b),[1])
        result=Result([obs],'same','{}')
        result.add('event','event',event,[(a,b)])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.npz'
            result.save(path)
            loaded=Result.load(path)
            self.assertEqual(loaded.observables[0].selection,{'charged':True})
            self.assertEqual(loaded.observables[0].flow['reference'],{'charged':True})
            np.testing.assert_allclose(loaded.summary()['observables']['v4']['value'],[1])

    def test_conflicting_or_mixed_selectors_are_rejected(self):
        for selection in ({'pdg':['all_charged'],'charged':False},
                          {'pdg':['all_charged',211]}, {'pdg':['all_charged'],'charged':1}):
            with self.subTest(selection=selection):
                with self.assertRaises(ValueError):Selection(**selection)
                with self.assertRaises(ValueError):Observable('n','yield',selection=selection)
                with self.assertRaises(ValueError):
                    Observable('v','v2',bins=[0,1],flow={'method':'cumulant_4','reference':selection})
        with self.assertRaises(ValueError):Selection(pdg=['all'])


if __name__=='__main__':unittest.main()
