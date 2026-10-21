"""Independent ordered-tuple checks of Bilandzic et al. arXiv:1010.0233v2."""
import itertools
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from maf import Observable
from maf.observables import flow_cumulants
from maf.statistics.aggregation import Result
from test_reference_physics import event_with_angles


def explicit_tuples(phi, coordinate, reference, poi, bins, harmonic, four=True):
    """Definition (5), (6), (20), (21): loop over distinct particle indices."""
    refs=np.flatnonzero(reference).tolist()
    nbin=len(bins)-1
    a=np.zeros(2+2*nbin if four else 1+nbin)
    b=np.zeros_like(a)
    offset=2 if four else 1
    for order,feature in ((2,0),(4,1)) if four else ((2,0),):
        def correlation(indices):
            half=order//2
            return np.cos(harmonic*(sum(phi[i] for i in indices[:half])-sum(phi[i] for i in indices[half:])))
        for indices in itertools.permutations(refs,order):
            a[feature]+=correlation(indices);b[feature]+=1
        for k,(low,high) in enumerate(zip(bins[:-1],bins[1:])):
            feature=offset+k+(nbin if order==4 else 0)
            particles=np.flatnonzero(poi & np.isfinite(coordinate) & (coordinate>=low) & (coordinate<high))
            for first in particles:
                for rest in itertools.permutations(refs,order-1):
                    if first not in rest:
                        a[feature]+=correlation((first,)+rest);b[feature]+=1
    return a,b


def manual_v4(a,b):
    """Equations (12), (34), (35) after summing independently enumerated tuples."""
    ref2,ref4=a[0]/b[0],a[1]/b[1]
    nbin=(len(a)-2)//2
    diff2=a[2:2+nbin]/b[2:2+nbin]
    diff4=a[2+nbin:]/b[2+nbin:]
    c4=ref4-2*ref2**2
    return -(diff4-2*diff2*ref2)/(-c4)**.75


class DirectCumulantTests(unittest.TestCase):
    def test_all_small_reference_and_poi_masks(self):
        rng=np.random.default_rng(10100233)
        for size in range(5):
            phi=rng.uniform(-np.pi,np.pi,size)
            coordinate=np.linspace(0,1,size,endpoint=False)
            for refbits,poibits in itertools.product(range(1<<size),repeat=2):
                ref=np.array([bool(refbits & (1<<i)) for i in range(size)],dtype=bool)
                poi=np.array([bool(poibits & (1<<i)) for i in range(size)],dtype=bool)
                for harmonic in (1,2):
                    expected=explicit_tuples(phi,coordinate,ref,poi,[0,.5,1],harmonic)
                    actual=flow_cumulants.calculate(phi,coordinate,ref,poi,[0,.5,1],harmonic,True)
                    np.testing.assert_allclose(actual[0],expected[0],atol=2e-12)
                    np.testing.assert_array_equal(actual[1],expected[1])

    def test_random_partial_disjoint_and_nested_selections(self):
        rng=np.random.default_rng(32035)
        for size in (5,6,7,8):
            for trial in range(8):
                phi=rng.uniform(-np.pi,np.pi,size)
                coordinate=rng.uniform(0,2,size)
                reference=rng.random(size)<.65
                poi=rng.random(size)<.5
                for harmonic in (1,2,3):
                    for four in (False,True):
                        expected=explicit_tuples(phi,coordinate,reference,poi,[0,1,2],harmonic,four)
                        actual=flow_cumulants.calculate(phi,coordinate,reference,poi,[0,1,2],harmonic,four)
                        np.testing.assert_allclose(actual[0],expected[0],atol=1e-10)
                        np.testing.assert_array_equal(actual[1],expected[1])

    def test_minimum_reference_multiplicities_for_disjoint_poi(self):
        for count in range(5):
            size=count+1
            reference=np.arange(size)<count
            poi=~reference
            a,b=flow_cumulants.calculate(np.zeros(size),np.zeros(size),reference,poi,[-1,1],2,True)
            self.assertEqual(b[2],count)
            self.assertEqual(b[3],count*(count-1)*(count-2) if count>=3 else 0)
            np.testing.assert_allclose(a,b)
            if count<4:self.assertEqual(b[1],0)
            if count==1:self.assertEqual(b[0],0)

    def test_previous_version_checkpoint_stays_readable(self):
        obs=Observable('v4','v2',bins=[-1,1],flow={'method':'cumulant_4','reference':{}})
        result=Result([obs],'legacy','{}')
        event=event_with_angles([0,.1,.2,.3,.4])
        result.add('event','event',event,[obs.evaluate(event)])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'old.npz'
            result.save(path)
            with np.load(path,allow_pickle=False) as archive:
                arrays={key:archive[key].copy() for key in archive.files}
            metadata=json.loads(str(arrays['metadata']))
            metadata['version']='0.4.0'
            arrays['metadata']=json.dumps(metadata)
            np.savez_compressed(path,**arrays)
            loaded=Result.load(path)
            np.testing.assert_allclose(loaded.summary()['observables']['v4']['value'],
                                       result.summary()['observables']['v4']['value'])

    def test_ensemble_tuple_weights_differential_flow_and_jackknife(self):
        rng=np.random.default_rng(44)
        obs=Observable('proton_v2_4','v2',axis='pt',bins=[0,1,2],selection={'pdg':[2212]},
                       flow={'method':'cumulant_4','reference':{'pdg':[211]}})
        direct=Result([obs],'same','{}')
        left=Result([obs],'same','{}',sources={'left':{'metadata':{}}})
        right=Result([obs],'same','{}',sources={'right':{'metadata':{}}})
        raw=[]
        for i,(nref,npoi) in enumerate([(4,2),(7,3),(5,2),(6,4),(8,2)]):
            phi=rng.normal(0,.1+.08*i,nref+npoi)
            pt=np.r_[np.full(nref,.7),np.resize([.5,1.5],npoi)]
            pdg=np.r_[np.full(nref,211),np.full(npoi,2212)]
            event=event_with_angles(phi,pt=pt,pdg=pdg)
            values=obs.evaluate(event)
            raw.append(explicit_tuples(phi,pt,pdg==211,pdg==2212,[0,1,2],2))
            np.testing.assert_allclose(values[0],raw[-1][0],atol=1e-10)
            np.testing.assert_array_equal(values[1],raw[-1][1])
            direct.add(str(i),'event',event,[values])
            (left if i<2 else right).add(str(i),'event',event,[values])
        a=np.array([v[0] for v in raw]);b=np.array([v[1] for v in raw])
        expected=manual_v4(a.sum(axis=0),b.sum(axis=0))
        summary=direct.summary()['observables'][obs.name]
        np.testing.assert_allclose(summary['value'],expected,atol=1e-12)
        eventwise=np.array([manual_v4(x,y) for x,y in raw]).mean(axis=0)
        self.assertGreater(np.max(np.abs(expected-eventwise)),1e-3)
        replicas=np.array([manual_v4(a.sum(axis=0)-x,b.sum(axis=0)-y) for x,y in raw])
        delta=replicas-replicas.mean(axis=0)
        covariance=(len(raw)-1)/len(raw)*delta.T@delta
        np.testing.assert_allclose(summary['covariance'],covariance,atol=1e-12)
        left.merge(right)
        np.testing.assert_allclose(left.summary()['observables'][obs.name]['value'],expected,atol=1e-12)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.npz'
            left.save(path)
            np.testing.assert_allclose(Result.load(path).summary()['observables'][obs.name]['value'],expected,atol=1e-12)


if __name__=='__main__':unittest.main()
