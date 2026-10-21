import tempfile
import unittest
import warnings
from pathlib import Path
import numpy as np
from maf import OscarReader, Observable
from maf.core.selection import Selection
from maf.io.oscar import OscarFormatError, OscarWarning
from maf.io.metadata import parse_name, grouping_key
from maf.io.discovery import discover
from helpers import HEADER, row, block, write_run, event


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.path=Path(self.tmp.name)/"particle_lists.oscar.unfinished"
        self.addCleanup(self.tmp.cleanup)

    def read(self,text,**kwargs):
        self.path.write_text(text)
        return list(OscarReader(self.path,**kwargs))

    def test_complete_and_empty(self):
        events=self.read(HEADER+block(0,[row()])+block(1,[]))
        self.assertEqual([len(e.particles) for e in events],[1,0])
        self.assertEqual(events[0].particles.ids.tolist(),[0])

    def test_unknown_impact_retains_nonempty_and_empty_events(self):
        text=HEADER+(block(0,[row()])+block(1,[])).replace("impact 1.0","impact -1.000")+block(2,[row()])
        self.path.write_text(text)
        for mode in ("strict", "warn"):
            with self.subTest(errors=mode), warnings.catch_warnings(record=True) as messages:
                warnings.simplefilter("always")
                reader=OscarReader(self.path,errors=mode)
                events=list(reader)
            self.assertEqual(messages,[])
            self.assertEqual([len(e.particles) for e in events],[1,0,1])
            for e in events[:2]:
                self.assertIsNone(e.metadata["impact"])
                self.assertEqual(e.metadata["impact_raw"],-1.)
                self.assertIn("-1.000",e.metadata["end_fields"])
            self.assertEqual(events[2].metadata["impact"],1.)
            self.assertEqual(reader.diagnostics["events"],3)
            self.assertEqual(reader.diagnostics["unknown_impact"],2)
            self.assertEqual(reader.diagnostics["discarded"],0)
            self.assertEqual(reader.diagnostics["issues"],0)

    def test_missing_impact_is_unknown(self):
        self.path.write_text(HEADER+block(0,[]).replace(" impact 1.0",""))
        reader=OscarReader(self.path)
        event=list(reader)[0]
        self.assertIsNone(event.metadata["impact"])
        self.assertIsNone(event.metadata["impact_raw"])
        self.assertEqual(reader.diagnostics["unknown_impact"],1)

    def test_other_invalid_impact_values_are_rejected(self):
        for value in ("-2", "-0.5", "nan", "inf", "-inf"):
            with self.subTest(impact=value):
                text=HEADER+block(0,[row()]).replace("impact 1.0",f"impact {value}")
                with self.assertRaisesRegex(OscarFormatError,"Invalid event impact parameter"):
                    self.read(text)
                self.path.write_text(text+block(1,[]))
                reader=OscarReader(self.path,errors="warn")
                with self.assertWarns(OscarWarning):
                    events=list(reader)
                self.assertEqual([e.event_id for e in events],[1])
                self.assertEqual(reader.diagnostics["discarded"],1)
                self.assertEqual(reader.diagnostics["unknown_impact"],0)

    def test_last_snapshot_only(self):
        text=HEADER+block(0,[row()],end=False)+block(0,[row(),row(ident=1)])
        events=self.read(text)
        self.assertEqual(len(events),1)
        self.assertEqual(len(events[0].particles),2)

    def test_in_block_ignored(self):
        text=HEADER+"# event 0 in 1\n"+row()+block(0,[])
        self.assertEqual(len(self.read(text)[0].particles),0)

    def test_no_end_marker(self):
        with self.assertRaises(OscarFormatError):
            self.read(HEADER+block(0,[row()],end=False))
        with self.assertWarns(OscarWarning):
            events=self.read(HEADER+block(0,[row()],end=False),errors="warn")
        self.assertEqual(events,[])

    def test_truncated_count(self):
        with self.assertRaises(OscarFormatError):
            self.read(HEADER+"# event 0 out 2\n"+row())

    def test_discard_whole_malformed_event(self):
        bad=row().replace("2212","nan")
        with warnings.catch_warnings(record=True) as messages:
            events=self.read(HEADER+block(0,[bad])+block(1,[row()]),errors="warn")
        self.assertTrue(messages)
        self.assertEqual([e.event_id for e in events],[1])

    def test_shell_validation(self):
        bad="0 0 0 0 1 0.1 1 0 0 2212 0 1\n"
        with self.assertRaises(OscarFormatError):
            self.read(HEADER+block(0,[bad]))

    def test_reordered_columns(self):
        columns=HEADER.splitlines()[0].split()[2:]
        fields=row().split()
        order=list(reversed(range(len(columns))))
        text="#!OSCAR2013Extended particle_lists "+" ".join(columns[i] for i in order)+"\n"
        events=self.read(text+block(0,[" ".join(fields[i] for i in order)+"\n"]))
        self.assertEqual(events[0].particles.pdg.tolist(),[2212])
        np.testing.assert_allclose(events[0].particles.px,[1.])

    def test_headerless_refused(self):
        with self.assertRaises(OscarFormatError):
            self.read(block(0,[row()]))

    def test_bad_units_refused(self):
        with self.assertRaises(OscarFormatError):
            self.read(HEADER.replace("GeV","MeV")+block(0,[]))

    def test_duplicate_event_refused(self):
        with self.assertRaises(OscarFormatError):
            self.read(HEADER+block(0,[])+block(0,[]))

    def test_weight_and_ensemble(self):
        events=self.read(HEADER+block(0,[],ensemble=2,weight=2.5))
        self.assertEqual(events[0].ensemble_id,2)
        self.assertEqual(events[0].weight,2.5)
        with self.assertRaises(OscarFormatError):
            self.read(HEADER+block(0,[],weight=-1))


class MetadataTests(unittest.TestCase):
    def test_legacy_names(self):
        for name,energy,run in (("4.3_AuAu_run__67",4.3,67),("7.7_AuAU_run_100",7.7,100),("7.7AuAutf0p2_eth0p4_ec0p3339_run95_scratch",7.7,95)):
            meta=parse_name(name)
            self.assertEqual(meta["energy"],{"value":energy,"unit":"GeV"})
            self.assertEqual(meta["run_id"],run)
            self.assertEqual(meta["collision_system"],"Au+Au")
            self.assertNotIn("parameters",meta)

    def test_arbitrary_parameter(self):
        rules=[{"pattern":r"visc(?P<v>\d+p\d+)","fields":{"v":{"key":"parameters.viscosity","type":"float"}}}]
        self.assertEqual(parse_name("7.7AuAu_visc0p12",rules)["parameters"]["viscosity"],.12)

    def test_missing_groups_never_mix(self):
        self.assertNotEqual(grouping_key({},["parameters.tf"],"a"),grouping_key({},["parameters.tf"],"b"))
        self.assertNotEqual(grouping_key({},["parameters.tf"],"a"),grouping_key({"parameters":{"tf":.2}},["parameters.tf"],"b"))

    def test_discovery_manifest_and_finished_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=write_run(tmp,"7.7_AuAU_run_100",[block(0,[])])
            (path/"particle_lists.oscar.unfinished").write_text("bad")
            manifest={"defaults":{"tags":["x"]},"runs":{path.name:{"seed":2,"parameters":{"x":.5}}}}
            runs=discover(tmp,manifest=manifest)
            self.assertEqual(len(runs),1)
            self.assertEqual(runs[0].oscar.name,"particle_lists.oscar")
            self.assertEqual(runs[0].metadata["seed"],2)


class KinematicsTests(unittest.TestCase):
    def test_y_eta_and_rest(self):
        p=event(pt=1,pz=2,mass=1).particles
        self.assertAlmostEqual(p.rapidity[0],np.arctanh(2/np.sqrt(6)))
        self.assertAlmostEqual(p.eta[0],np.arcsinh(2))
        self.assertNotAlmostEqual(p.eta[0],p.rapidity[0])
        rest=event(pt=0,pz=0,mass=1).particles
        self.assertEqual(rest.rapidity[0],0)
        self.assertTrue(np.isnan(rest.eta[0]))
        self.assertTrue(Selection(y=[-.1,.1]).mask(rest)[0])

    def test_rapidity_boost(self):
        p=event().particles
        cut=Selection(frame="cm",input_frame="lab",cm_rapidity=.3)
        self.assertAlmostEqual(cut.rapidity(p)[0],-.3)
        with self.assertRaises(ValueError):
            Selection(frame="cm",input_frame="lab")
        with self.assertRaises(ValueError):
            Selection(frame="cm",input_frame="lab",cm_rapidity=.3,eta=[-1,1])

    def test_charge_and_pdg(self):
        p=event(3,pdg=[2212,2112,-211],charge=[1,0,-1]).particles
        np.testing.assert_array_equal(Selection(charged=True).mask(p),[True,False,True])
        np.testing.assert_array_equal(Selection(pdg=[2112]).mask(p),[False,True,False])

    def test_y_eta_distinct_acceptance(self):
        p=event(pt=.1,pz=.2,mass=1).particles
        self.assertTrue(Selection(y=[-.5,.5]).mask(p)[0])
        self.assertFalse(Selection(eta=[-.5,.5]).mask(p)[0])

    def test_undocumented_particle_correction_refused(self):
        with self.assertRaises(ValueError):
            Selection(particle_definition="primary")
