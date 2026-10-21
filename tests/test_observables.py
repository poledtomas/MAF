import unittest
import numpy as np
from maf import Observable
from helpers import event


class ObservableTests(unittest.TestCase):
    def test_yield_and_absent_species(self):
        e=event(2)
        self.assertEqual(Observable("p","yield",selection={"pdg":[2212]}).evaluate(e)[0][0],2)
        self.assertEqual(Observable("pbar","yield",selection={"pdg":[-2212]}).evaluate(e)[0][0],0)

    def test_charged_distribution(self):
        e=event(3,pdg=[2212,2112,211],charge=[1,0,1])
        obs=Observable("ch","dndeta",bins=[-1,1])
        np.testing.assert_array_equal(obs.evaluate(e)[0],[2])
        np.testing.assert_allclose(obs.normalization(),[.5])

    def test_nonuniform_bins(self):
        obs=Observable("y","dndy",bins=[-2,0,1])
        np.testing.assert_allclose(obs.normalization(),[.5,1])
        np.testing.assert_array_equal(obs.evaluate(event())[0],[0,1])

    def test_invariant_annulus_including_zero(self):
        obs=Observable("pt","invariant_pt",bins=[0,1,3],selection={"y":[-.5,.5]})
        np.testing.assert_allclose(obs.normalization(),[1/np.pi,1/(8*np.pi)])
        np.testing.assert_array_equal(obs.evaluate(event(pt=1))[0],[0,1])

    def test_invariant_requires_y(self):
        with self.assertRaises(ValueError):
            Observable("pt","invariant_pt",bins=[0,1],selection={"eta":[-.5,.5]})

    def test_midrapidity_density(self):
        obs=Observable("n","yield",selection={"y":[-.25,.25]},density="y")
        self.assertEqual(obs.evaluate(event(2))[0][0]*obs.normalization()[0],4)

    def test_net_proton(self):
        obs=Observable("net","net_proton",bins=[-1,1])
        self.assertEqual(obs.evaluate(event(3,pdg=[2212,-2212,2212]))[0][0],1)

    def test_mean_pt_and_empty(self):
        obs=Observable("mean","mean_pt")
        a,b=obs.evaluate(event(2,pt=[1,3]))
        self.assertEqual(a[0]/b[0],2)
        a,b=obs.evaluate(event(0))
        self.assertEqual((a[0],b[0]),(0,0))

    def test_mt(self):
        obs=Observable("mt","mean_mt_minus_m0",mass=1,selection={"pdg":[2212]})
        a,b=obs.evaluate(event(pt=1,mass=1))
        self.assertAlmostEqual(a[0]/b[0],np.sqrt(2)-1)
        spectrum=Observable("mtd","mt_minus_m0",mass=1,bins=[0,.5,1],selection={"pdg":[2212]})
        np.testing.assert_array_equal(spectrum.evaluate(event())[0],[1,0])

    def test_ratio(self):
        obs=Observable("ratio","ratio",selection={"pdg":[2212]},denominator_pdg=[-2212])
        a,b=obs.evaluate(event(3,pdg=[2212,2212,-2212]))
        self.assertEqual(a[0]/b[0],2)

    def test_bin_validation(self):
        for bins in ([1,0],[0,0],[0,float("nan")]):
            with self.assertRaises(ValueError):
                Observable("bad","pt",bins=bins)
        with self.assertRaises(ValueError):
            Observable("bad","pt",bins=[0,2],selection={"pt":[0,1]})

    def test_half_open_final_edge(self):
        obs=Observable("pt","pt",bins=[0,1,2])
        np.testing.assert_array_equal(obs.evaluate(event(3,pt=[0,1,2]))[0],[1,1])

    def test_scalar_bins_and_invalid_selections_rejected(self):
        for definition in ({"name":"bad","kind":"mean_pt","bins":[0,1,2]},
                           {"name":"bad","kind":"yield","selection":{"pdg":"2212"}},
                           {"name":"bad","kind":"yield","selection":{"charged":"false"}}):
            with self.assertRaises(ValueError):
                Observable(**definition)

    def test_invariant_mt_annulus(self):
        obs=Observable("mt","invariant_mt",mass=1,bins=[1,2,3],selection={"pdg":[2212],"y":[-1,1]})
        np.testing.assert_allclose(obs.normalization(),[1/(6*np.pi),1/(10*np.pi)])
        np.testing.assert_array_equal(obs.evaluate(event(pt=1,mass=1))[0],[1,0])
