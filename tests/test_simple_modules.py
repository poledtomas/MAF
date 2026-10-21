"""Analytic flow fixtures, short configurations and main-program integration."""
import contextlib
import copy
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from maf import Event, Particles, Observable
from maf.analysis import analyze_run, merge_results
from maf.config import expand_observables, load_job
from maf.io.discovery import discover
from maf.main import main
from maf.observables import dndy, dndeta, pt_spectra, mt_spectra, mean_pt, mean_mt, net_proton, ratios, yields, v1, v2
from maf.statistics.aggregation import Result
from helpers import block, row, write_run


def angular_event(phi, pt=1.0):
    phi = np.asarray(phi, float)
    pt = np.broadcast_to(pt, phi.shape)
    zero = np.zeros(len(phi))
    return Event(Particles(np.hypot(1,pt),pt*np.cos(phi),pt*np.sin(phi),zero,np.ones(len(phi)),
                           np.full(len(phi),2212),np.ones(len(phi),int)),0)


class PlainModulesTests(unittest.TestCase):
    def test_standalone_counts_and_normalization(self):
        for module in (dndy,dndeta,pt_spectra):
            a,b=module.calculate([0,.5,1,3],[0,1,3])
            np.testing.assert_array_equal(a,[2,1])
            self.assertEqual(b,1)
            np.testing.assert_allclose(module.normalization([0,1,3]),[1,.5])
        np.testing.assert_allclose(pt_spectra.normalization([0,1,3],invariant=True),[1/np.pi,1/(8*np.pi)])

    def test_standalone_scalar_and_mt_modules(self):
        self.assertEqual(yields.calculate(0),([0],1))
        self.assertEqual(ratios.calculate(2,4),([2],4))
        self.assertEqual(mean_pt.calculate([1,3]),([4],2))
        np.testing.assert_allclose(mean_mt.calculate([1],1)[0],[np.sqrt(2)-1])
        np.testing.assert_array_equal(mt_spectra.calculate([0,1,2],1,[0,.5,2],subtract_mass=True)[0],[2,1])
        np.testing.assert_array_equal(net_proton.calculate([0,0,0],[2212,-2212,211],[-1,1])[0],[0])

    def test_v1_and_v2_known_angles(self):
        phi=np.array([0,np.pi/2,np.pi,3*np.pi/2])
        for module,harmonic in ((v1,1),(v2,2)):
            a,b=module.calculate(np.cos(phi),np.sin(phi),np.zeros(4),[-2,0,2],0)
            np.testing.assert_allclose(a,[0,0],atol=1e-15)
            np.testing.assert_array_equal(b,[0,4])
            a,b=module.calculate([1],[0],[0],[-1,1],0)
            np.testing.assert_allclose(a,[1])
            np.testing.assert_array_equal(b,[1])

    def test_plane_rotation_and_pt_zero(self):
        for module in (v1,v2):
            a,b=module.calculate([0,0],[1,0],[0,0],[-2,2],np.pi/2)
            np.testing.assert_allclose(a,[1])
            np.testing.assert_array_equal(b,[1])
            for angle in (None,float('nan'),"zero",True):
                with self.assertRaises(ValueError):
                    module.calculate([1],[0],[0],[-1,1],angle)

    def test_flow_is_not_a_density(self):
        for kind in ('v1','v2'):
            obs=Observable(kind,kind,bins=[-2,2],flow={"method":"reaction_plane","angle":0})
            a,b=obs.evaluate(angular_event([0,0]))
            np.testing.assert_allclose(a/b*obs.normalization(),[1])
            missing=Observable(kind,kind,bins=[-1,1],flow={"method":"reaction_plane","metadata_key":"plane"})
            with self.assertRaisesRegex(ValueError,"Missing reaction-plane"):
                missing.evaluate(angular_event([0]))
            event=angular_event([np.pi/2])
            event.metadata['plane']=np.pi/2
            np.testing.assert_allclose(missing.evaluate(event)[0],[1])

    def test_undefined_and_unsupported_flow_rejected(self):
        for flow in (None,{"method":"event_plane","angle":0},{"method":"reaction_plane"},
                     {"method":"reaction_plane","angle":0,"metadata_key":"plane"}):
            with self.assertRaises(ValueError):
                Observable('v','v2',bins=[-1,1],flow=flow)
        with self.assertRaises(ValueError):
            Observable('v','v1',bins=[-1,1],axis='wrong',flow={"method":"reaction_plane","angle":0})


class SimpleConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        write_run(self.root/'runs','7.7_AuAu_run_1',[block(0,[row()]),block(1,[])])
        self.config={"input":"runs","output":"out","sampling":{"mode":"event"},
                     "observables":{"dndy":{"pdg":[2212,-2212],"y_bins":[-1,0,1]},
                                    "pt_spectra":{"pdg":[2212],"y_range":[-.5,.5],"pt_bins":[0,1,2]},
                                    "v1":{"enabled":False},"v2":{"enabled":False}}}

    def path(self,document=None):
        path=self.root/'analysis.json'
        path.write_text(json.dumps(document or self.config))
        return path

    def test_species_expansion_and_disabled_modules(self):
        job=load_job(self.path())
        self.assertEqual(job['root'],self.root/'runs')
        self.assertEqual(job['output'],self.root/'out')
        definitions=job['config']['observables']
        self.assertEqual([d['name'] for d in definitions],['dndy_2212','dndy_-2212','pt_spectra_2212'])
        self.assertEqual(definitions[0]['selection']['pdg'],[2212])
        self.assertEqual(definitions[1]['selection']['pdg'],[-2212])
        self.assertEqual(job['manifest']['defaults']['sampling'],{'mode':'event'})

    def test_legacy_explicit_list_unchanged(self):
        old=[{'name':'both','kind':'yield','selection':{'pdg':[211,-211]}}]
        self.assertEqual(expand_observables(old),old)

    def test_short_flow_and_invariant_configuration(self):
        options={'pt_spectra':{'pdg':[2212],'invariant':True,'pt_bins':[0,1,2],'y_range':[-.5,.5]},
                 'v2':{'pdg':[2212],'axis':'pt','pt_bins':[0,1,2],'flow':{'method':'reaction_plane','angle':0}}}
        expanded=expand_observables(options)
        self.assertEqual(expanded[0]['kind'],'invariant_pt')
        self.assertEqual(expanded[1]['axis'],'pt')
        for definition in expanded:
            Observable(**definition)

    def test_ambiguous_bins_boolean_and_unknown_options_rejected(self):
        for options in ({'dndy':{'enabled':'yes'}},{'dndy':{'pt_bins':[0,1]}},
                        {'dndy':{'y_bins':[-1,1],'bins':[-1,1]}},
                        {'dndy':{'pdg':2212,'y_bins':[-1,1]}}):
            with self.assertRaises(ValueError):
                expand_observables(options)
        document=copy.deepcopy(self.config)
        document['observables']['dndy']['typo']=1
        with self.assertRaises(ValueError):
            load_job(self.path(document))

    def test_yaml_template_parses(self):
        try:
            import yaml
        except ImportError:
            self.skipTest('PyYAML is optional')
        from maf.io.config import load_config
        template=Path(__file__).resolve().parents[1]/'analysis.yaml'
        expanded=expand_observables(load_config(template)['observables'])
        for definition in expanded:
            Observable(**definition)

    def test_main_dry_run_and_single_stream(self):
        path=self.path()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config',str(path),'--dry-run']),0)
        self.assertFalse((self.root/'out').exists())
        from maf.io.oscar import OscarReader
        original=OscarReader.__iter__
        reads=[]
        def track(reader):
            reads.append(reader.path)
            yield from original(reader)
        with patch.object(OscarReader,'__iter__',track), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config',str(path)]),0)
        self.assertEqual(len(reads),1)
        summary=json.loads(next((self.root/'out').rglob('summary.json')).read_text())
        self.assertEqual(summary['events'],2)
        self.assertEqual(summary['observables']['dndy_2212']['value'],[0,.5])
        self.assertEqual(summary['observables']['dndy_-2212']['value'],[0,0])

    def test_legacy_command_dispatch(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['discover',str(self.root/'runs')]),0)

    def test_bad_config_fails_before_writing(self):
        document=copy.deepcopy(self.config)
        document['observables']['v2']={'enabled':True,'y_bins':[-1,1]}
        with self.assertLogs(level='ERROR'),contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config',str(self.path(document))]),1)
        self.assertFalse((self.root/'out').exists())

    def test_metadata_plane_uncertainties_and_merging(self):
        config={'group_by':['energy','collision_system'],'observables':[{'name':'v','kind':'v1','bins':[-2,2],
                'flow':{'method':'reaction_plane','metadata_key':'plane'}}]}
        run=discover(self.root/'runs')[0]
        # One particle at phi=0 in each event; the known reaction plane rotates.
        run.oscar.write_text(run.oscar.read_text().replace(block(1,[]),block(1,[row()])))
        run.metadata['sampling']={'mode':'event','events':{'0:0':{'metadata':{'plane':0}},'1:0':{'metadata':{'plane':np.pi}}}}
        result=analyze_run(run,config)
        summary=result.summary()['observables']['v']
        np.testing.assert_allclose(summary['value'],[0],atol=1e-15)
        np.testing.assert_allclose(summary['stat_error'],[1])
        result.save(self.root/'flow.npz')
        self.assertEqual(Result.load(self.root/'flow.npz').summary(),result.summary())
        # Merge another independent file without changing the mean or weighting.
        write_run(self.root/'runs','7.7_AuAu_run_2',[block(0,[row()])])
        other=discover(self.root/'runs')[1]
        other.metadata['sampling']={'mode':'event','events':{'0:0':{'metadata':{'plane':0}}}}
        merged=merge_results([result,analyze_run(other,config)])[0]
        np.testing.assert_allclose(merged.summary()['observables']['v']['value'],[1/3])

    def test_directory_glob_selects_only_requested_configuration(self):
        prefix='3.9_clasicAuAutf0p6_eth0p5_ec0p25_'
        write_run(self.root/'runs',prefix+'run1_scratch',[block(0,[row()])])
        write_run(self.root/'runs',prefix+'run2_scratch',[block(0,[])])
        write_run(self.root/'runs','3.9_clasicAuAutf0p4_eth0p5_ec0p25_run1',[block(0,[row()])])
        document=copy.deepcopy(self.config)
        document['run_pattern']=prefix+'*'
        job=load_job(self.path(document))
        from maf.analysis import analyze
        runs=analyze(**job,dry_run=True)
        self.assertEqual(len(runs),2)
        for run in runs:
            self.assertTrue(run.path.name.startswith(prefix))
            self.assertEqual(run.metadata['energy'],{'value':3.9,'unit':'GeV'})
            self.assertEqual(run.metadata['collision_system'],'Au+Au')
            self.assertEqual(run.metadata['directory_label'],'clasic')
            self.assertNotIn('parameters',run.metadata)

    def test_all_selected_files_form_one_event_ensemble(self):
        prefix='3.9_clasicAuAutf0p6_eth0p5_ec0p25_'
        # Unequal file sizes expose accidental equal-file weighting. IDs restart per file.
        write_run(self.root/'runs',prefix+'run1_scratch',[block(0,[row(ident=j) for j in range(4)])])
        write_run(self.root/'runs',prefix+'run2_scratch',[block(0,[]),block(1,[row()]),block(2,[row()])])
        write_run(self.root/'runs','3.9_clasicAuAutf0p4_eth0p5_ec0p25_run1',[block(0,[row()])])
        document=copy.deepcopy(self.config)
        document['run_pattern']=prefix+'*'
        document['group_by']=[]
        document['observables']['pt_spectra']['invariant']=True
        path=self.path(document)
        job=load_job(path)
        direct_path=write_run(self.root/'concatenated','same-events',[
            block(0,[row(ident=j) for j in range(4)]),block(1,[]),block(2,[row()]),block(3,[row()])])
        direct_run=discover(direct_path,manifest=job['manifest'])[0]
        direct=analyze_run(direct_run,job['config']).summary()
        from maf.io.oscar import OscarReader
        for workers in (1,2):
            with self.subTest(workers=workers),contextlib.redirect_stdout(io.StringIO()) as output:
                document['output']=f'ensemble-{workers}'
                path=self.path(document)
                self.assertEqual(main(['--config',str(path),'--workers',str(workers)]),0)
                directory=self.root/document['output']
                index=json.loads((directory/'index.json').read_text())
                self.assertEqual(index['failures'],[])
                self.assertEqual(len(index['groups']),1)
                self.assertEqual(Path(index['groups'][0]).parent,directory/'results')
                self.assertFalse(list(directory.glob('group-*')))
                self.assertFalse(list((directory/'results').glob('*.csv')))
                result=Result.load(index['groups'][0])
                summary=result.summary()
                self.assertEqual(summary['group'],{})
                self.assertEqual(summary['events'],4)
                self.assertEqual(summary['independent_units'],4)
                self.assertEqual(summary['sampling_modes'],['event'])
                self.assertEqual(len(summary['sources']),2)
                self.assertEqual(sorted(s['selected_events'] for s in summary['sources'].values()),[1,3])
                self.assertEqual(summary['observables']['dndy_2212']['value'],[0,1.5])
                np.testing.assert_allclose(summary['observables']['dndy_2212']['stat_error'],[0,np.sqrt(3)/2])
                np.testing.assert_allclose(summary['observables']['pt_spectra_2212']['value'],[0,.5/np.pi])
                for name in summary['observables']:
                    for field in ('value','stat_error','covariance'):
                        np.testing.assert_allclose(summary['observables'][name][field],direct['observables'][name][field])
                caches={p:p.stat().st_mtime_ns for p in (directory/'cache').glob('*.npz')}
                self.assertEqual(len(caches),2)
                # Resume must use the two raw checkpoints without rereading particles.
                with patch.object(OscarReader,'__iter__',side_effect=AssertionError('Unnecessary OSCAR reread')):
                    self.assertEqual(main(['--config',str(path),'--workers','1']),0)
                self.assertEqual(caches,{p:p.stat().st_mtime_ns for p in caches})
                resumed=Result.load(index['groups'][0]).summary()
                self.assertEqual(summary,resumed)
            self.assertIn('4 selected events from 2 files',output.getvalue())

    def test_common_ensemble_keeps_flow_calculations_per_event(self):
        prefix='3.9_clasicAuAutf0p6_eth0p5_ec0p25_'
        write_run(self.root/'runs',prefix+'run1_scratch',[block(0,[row(ident=j) for j in range(2)])])
        transverse=[]
        for j in range(3):
            fields=row(ident=j).split()
            fields[6],fields[7]='0','1'
            transverse.append(' '.join(fields)+'\n')
        write_run(self.root/'runs',prefix+'run2_scratch',[block(0,transverse),block(1,[])])
        document=copy.deepcopy(self.config)
        document['run_pattern']=prefix+'*'
        document['group_by']=[]
        for kind in ('v1','v2'):
            document['observables'][kind]={'enabled':True,'pdg':[2212],'y_bins':[-1,1],
                'flow':{'method':'reaction_plane','angle':0.,'weighting':'event'}}
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config',str(self.path(document))]),0)
        index=json.loads((self.root/'out'/'index.json').read_text())
        self.assertEqual(len(index['groups']),1)
        summary=Result.load(index['groups'][0]).summary()
        self.assertEqual(summary['events'],3)
        self.assertEqual(summary['independent_units'],3)
        np.testing.assert_allclose(summary['observables']['v1_2212']['value'],[.5],atol=1e-15)
        np.testing.assert_allclose(summary['observables']['v2_2212']['value'],[0.],atol=1e-15)
        self.assertEqual(summary['observables']['dndy_2212']['value'],[0,5/3])

    def test_bad_directory_glob_rejected(self):
        document=copy.deepcopy(self.config)
        document['run_pattern']=[]
        with self.assertRaisesRegex(ValueError,'run_pattern'):
            load_job(self.path(document))
