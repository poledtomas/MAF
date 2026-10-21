import copy
import tempfile
import unittest
import warnings
from pathlib import Path
import numpy as np
from maf import Observable
from maf.analysis import analyze_run, analyze, merge_results
from maf.io.discovery import discover
from maf.statistics.aggregation import Result
from maf.statistics.errors import jackknife, bootstrap
from helpers import HEADER, block, row, write_run


CONFIG={"group_by":["energy","collision_system"],"reader":{"errors":"strict"},
        "observables":[{"name":"n","kind":"yield","selection":{"pdg":[2212]}}]}


class StatisticsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def run_result(self,counts,*,name="7.7_AuAu_run_1",sampling=None,config=None):
        path=write_run(self.root,name,[block(i,[row(ident=j) for j in range(count)]) for i,count in enumerate(counts)])
        run=next(r for r in discover(self.root) if r.path==path)
        run.metadata["sampling"]=sampling or {"mode":"event"}
        return analyze_run(run,config or CONFIG)

    def test_event_mean_and_standard_error(self):
        result=self.run_result([0,2,4])
        output=result.summary()
        self.assertEqual(output["events"],3)
        self.assertAlmostEqual(output["observables"]["n"]["value"][0],2)
        self.assertAlmostEqual(output["observables"]["n"]["stat_error"][0],np.std([0,2,4],ddof=1)/np.sqrt(3))

    def test_empty_events_in_denominator(self):
        result=self.run_result([0,2])
        self.assertEqual(result.summary()["observables"]["n"]["value"],[1])

    def test_unknown_impact_events_enter_normalization(self):
        blocks=[block(0,[row(),row(ident=1)]),block(1,[])]
        write_run(self.root,"7.7_AuAu_run_1",[b.replace("impact 1.0","impact -1.000") for b in blocks])
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter("always")
            result=analyze_run(run,CONFIG)
        self.assertEqual(messages,[])
        summary=result.summary()
        self.assertEqual(summary["events"],2)
        self.assertEqual(summary["observables"]["n"]["value"],[1])
        source=next(iter(summary["sources"].values()))
        self.assertEqual(source["diagnostics"]["unknown_impact"],2)
        self.assertEqual(source["diagnostics"]["discarded"],0)
        self.assertEqual(source["diagnostics"]["issues"],0)
        self.assertTrue(any("no impact parameter" in m for m in summary["messages"]))
        self.assertFalse(any("Malformed" in m for m in summary["messages"]))
        result.save(self.root/"result.npz")
        self.assertEqual(Result.load(self.root/"result.npz").summary(),summary)

    def test_impact_cut_requires_known_impact(self):
        write_run(self.root,"7.7_AuAu_run_1",[block(0,[row()]).replace("impact 1.0","impact -1.000")])
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        config=copy.deepcopy(CONFIG)
        config["event_selection"]={"impact":[0,2]}
        with self.assertRaisesRegex(ValueError,"Missing event field required for selection: impact"):
            analyze_run(run,config)

    def test_hydro_oversamples_equal_hydro_weight(self):
        sampling={"mode":"hydro","hydro_ids":{"0:0":"a","1:0":"a","2:0":"b"}}
        result=self.run_result([1,3,10],sampling=sampling)
        summary=result.summary()
        self.assertEqual(summary["events"],3)
        self.assertEqual(summary["independent_units"],2)
        self.assertEqual(summary["observables"]["n"]["value"],[6])
        self.assertEqual(summary["observables"]["n"]["stat_error"],[4])

    def test_hydro_missing_id_refused(self):
        with self.assertRaisesRegex(ValueError,"Missing hydro ID"):
            self.run_result([1],sampling={"mode":"hydro"})

    def test_unknown_independence_warns(self):
        with self.assertWarns(UserWarning):
            result=self.run_result([1,3],sampling={"mode":"run"})
        summary=result.summary()
        self.assertEqual(summary["independent_units"],1)
        self.assertIsNone(summary["observables"]["n"]["stat_error"][0])

    def test_merged_equals_direct(self):
        first=self.run_result([0,2],name="7.7_AuAu_run_1")
        second=self.run_result([4,6],name="7.7_AuAu_run_2")
        merged=merge_results([first,second])[0]
        direct=self.run_result([0,2,4,6],name="7.7_AuAu_run_3")
        for field in ("value", "stat_error", "covariance"):
            np.testing.assert_allclose(merged.summary()["observables"]["n"][field],direct.summary()["observables"]["n"][field],rtol=1e-14)

    def test_duplicate_sources_rejected(self):
        result=self.run_result([1,2])
        with self.assertRaisesRegex(ValueError,"Duplicate source"):
            merge_results([result,copy.deepcopy(result)])

    def test_hydro_shared_across_runs(self):
        a=self.run_result([1,3],name="7.7_AuAu_run_1",sampling={"mode":"hydro","namespace":"production","hydro_ids":{"0:0":"a","1:0":"b"}})
        b=self.run_result([5],name="7.7_AuAu_run_2",sampling={"mode":"hydro","namespace":"production","hydro_ids":{"0:0":"a"}})
        merged=merge_results([a,b])[0]
        self.assertEqual(len(merged.units),2)
        self.assertEqual(merged.summary()["observables"]["n"]["value"],[3])

    def test_saved_round_trip(self):
        result=self.run_result([0,2,4])
        path=self.root/"result.npz"
        result.save(path)
        loaded=Result.load(path)
        self.assertEqual(result.summary(),loaded.summary())
        result.export(self.root/"export")
        self.assertTrue((self.root/"export"/"n.dat").exists())

    def test_weighted_events(self):
        path=write_run(self.root,"7.7_AuAu_run_1",[block(0,[row()],weight=1),block(1,[row(),row(ident=1),row(ident=2)],weight=3)])
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        result=analyze_run(run,CONFIG)
        self.assertEqual(result.summary()["observables"]["n"]["value"],[2.5])

    def test_ratio_covariance(self):
        a=np.array([[1,1],[3,3],[5,5]],float)
        b=np.ones_like(a)
        value,error,cov=jackknife(a,b)
        self.assertAlmostEqual(cov[0,1],error[0]**2)
        zero=jackknife([[0],[1]],[[0],[1]])
        self.assertTrue(np.isnan(zero[1][0]))

    def test_bootstrap_deterministic(self):
        a=np.arange(4)[:,None]
        b=np.ones_like(a)
        np.testing.assert_array_equal(bootstrap(a,b,samples=20,seed=4),bootstrap(a,b,samples=20,seed=4))

    def test_cache_invalidates_on_metadata(self):
        self.run_result([1,2])
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        cache=self.root/"cache"
        analyze_run(run,CONFIG,cache=cache)
        checkpoint=next(cache.glob("*.npz"))
        before=checkpoint.stat().st_mtime_ns
        analyze_run(run,CONFIG,cache=cache)
        self.assertEqual(checkpoint.stat().st_mtime_ns,before)
        run.metadata["seed"]=99
        analyze_run(run,CONFIG,cache=cache)
        self.assertEqual(len(list(cache.glob("*.npz"))),2)

    def test_cache_invalidates_on_software_version(self):
        from unittest.mock import patch
        self.run_result([1,2])
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        cache=self.root/"cache"
        with patch("maf.analysis.__version__","0.3.0"):
            previous=analyze_run(run,CONFIG,cache=cache)
        current=analyze_run(run,CONFIG,cache=cache)
        self.assertNotEqual(previous.signature,current.signature)
        self.assertEqual(len(list(cache.glob("*.npz"))),2)

    def test_parallel_and_serial(self):
        self.run_result([0,2],name="7.7_AuAu_run_1")
        self.run_result([4,6],name="7.7_AuAu_run_2")
        manifest={"defaults":{"sampling":{"mode":"event"}}}
        a=analyze(self.root,CONFIG,self.root/"serial",workers=1,manifest=manifest)[0]
        b=analyze(self.root,CONFIG,self.root/"parallel",workers=2,manifest=manifest)[0]
        self.assertEqual(a.summary(),b.summary())

    def test_dry_run_no_output(self):
        self.run_result([1])
        output=self.root/"dry"
        analyze(self.root,CONFIG,output,dry_run=True)
        self.assertFalse(output.exists())

    def test_bad_run_reported(self):
        self.run_result([1,2])
        path=self.root/"bad"
        path.mkdir()
        (path/"particle_lists.oscar").write_text("bad")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results=analyze(self.root,CONFIG,self.root/"out",manifest={"defaults":{"sampling":{"mode":"event"}}})
        self.assertTrue(results)
        self.assertIn("bad",(self.root/"out"/"failures.json").read_text())

    def test_event_selection_and_limit(self):
        config=copy.deepcopy(CONFIG)
        config["event_selection"]={"impact":[0,2],"scattering_projectile_target":True}
        config["max_events"]=2
        path=write_run(self.root,"7.7_AuAu_run_1",[block(0,[row()]).replace("impact 1.0","impact 3.0"),block(1,[row(),row(ident=1)]),block(2,[])])
        # helpers omit the scattering flag; an explicit requested missing field must fail.
        run=discover(self.root)[0]
        run.metadata["sampling"]={"mode":"event"}
        with self.assertRaisesRegex(ValueError,"Missing event field"):
            analyze_run(run,config)
        config["event_selection"].pop("scattering_projectile_target")
        result=analyze_run(run,config)
        self.assertEqual(result.summary()["events"],1)
        self.assertEqual(result.summary()["observables"]["n"]["value"],[2])
        source=next(iter(result.sources.values()))
        self.assertEqual(source["diagnostics"]["events"],2)
        self.assertEqual(source["selected_events"],1)

    def test_different_sampling_modes_rejected(self):
        first=self.run_result([1,2],name="7.7_AuAu_run_1")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            second=self.run_result([3,4],name="7.7_AuAu_run_2",sampling={"mode":"run"})
        with self.assertRaisesRegex(ValueError,"sampling modes"):
            merge_results([first,second])

    def test_config_typo_rejected(self):
        config=dict(CONFIG,group_by="energy")
        with self.assertRaisesRegex(ValueError,"group_by"):
            self.run_result([1],config=config)
        with self.assertRaisesRegex(ValueError,"Unknown configuration"):
            self.run_result([1],config=dict(CONFIG,typo=3))

    def test_fewer_units_message_removed_after_merge(self):
        first=self.run_result([1],name="7.7_AuAu_run_1")
        second=self.run_result([3],name="7.7_AuAu_run_2")
        summary=merge_results([first,second])[0].summary()
        self.assertFalse(any("Fewer than two" in m for m in summary["messages"]))

    def test_merge_cli_preserves_distinct_configurations(self):
        from maf.cli import main
        first=self.run_result([1,3],name="7.7_AuAu_run_1")
        config=copy.deepcopy(CONFIG)
        config["observables"][0]["scale"]=2
        second=self.run_result([1,3],name="7.7_AuAu_run_2",config=config)
        first.save(self.root/"first.npz")
        second.save(self.root/"second.npz")
        output=self.root/"merged"
        status=main(["merge",str(self.root/"first.npz"),str(self.root/"second.npz"),"--output",str(output)])
        self.assertEqual(status,0)
        self.assertEqual(len(list(output.glob("results/group-*/accumulators.npz"))),2)


class RunPatternMetadataTests(unittest.TestCase):
    def test_pattern_survives_cached_analysis_npz_summary_and_dat(self):
        from unittest.mock import patch
        from maf.cli import plot_summaries
        import json
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            write_run(root/'runs','7.7_AuAu_run_1',[block(0,[row()]),block(1,[])])
            output=root/'output'
            kwargs=dict(root=root/'runs',config=CONFIG,output=output,run_pattern='7.7_AuAu_run_*',
                        manifest={'defaults':{'sampling':{'mode':'event'}}})
            analyze(**kwargs)
            with patch('maf.io.oscar.OscarReader.__iter__',side_effect=AssertionError('cached events must not be read')):
                result=analyze(**kwargs)[0]
            self.assertEqual(result.run_pattern,kwargs['run_pattern'])
            index=json.loads((output/'index.json').read_text())
            loaded=Result.load(index['groups'][0])
            self.assertEqual(loaded.run_pattern,kwargs['run_pattern'])
            summary=plot_summaries([output])[0]
            self.assertEqual(summary['run_pattern'],kwargs['run_pattern'])
            dat=output/index['summaries'][0]
            self.assertEqual(plot_summaries([dat.parent/'n.dat'])[0]['run_pattern'],kwargs['run_pattern'])
