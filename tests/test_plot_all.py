"""Saved-result batch plotting: current groups, axes and no event rereads."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from maf.analysis import analyze_run
from maf.io.discovery import discover
from maf.main import main
from helpers import block, row, write_run


class BatchPlotTests(unittest.TestCase):
    def setUp(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('Matplotlib is optional')
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        self.addCleanup(plt.close,'all')
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        path=write_run(self.root/'runs','7.7_AuAu_run_1',[
            block(0,[row(pt=.5)]),block(1,[row(pt=1.2),row(pt=.7,ident=1)]),block(2,[])])
        run=discover(path)[0]
        run.metadata['sampling']={'mode':'event'}
        self.definitions=[
            {'name':'spectrum','kind':'invariant_pt','bins':[0,1,2],'selection':{'pdg':[2212],'y':[-.5,.5]}},
            {'name':'rapidity','kind':'dndy','bins':[-1,0,1]},
            {'name':'charged_eta','kind':'dndeta','bins':[-1,0,1]},
            {'name':'flow','kind':'v2','bins':[-1,1],'flow':{'method':'reaction_plane','angle':0,'weighting':'event'}},
            {'name':'empty_spectrum','kind':'pt','bins':[0,1,2],'selection':{'pdg':[-2212]}}]
        result=analyze_run(run,{'group_by':[],'observables':self.definitions})
        self.results=self.root/'results_7.7_AuAu'
        self.group=self.results/'group-current'
        result.export(self.group)
        # A valid summary is sufficient; plotting must not open this checkpoint.
        (self.group/'accumulators.npz').write_bytes(b'not a readable NPZ')
        (self.results/'index.json').write_text(json.dumps({'groups':[str(self.group/'accumulators.npz')],'failures':[]}))
        stale=self.results/'group-stale'
        stale.mkdir()
        (stale/'summary.json').write_text('stale groups must not be read')
        self.config=self.root/'analysis.json'
        self.config.write_text(json.dumps({'input':'runs','output':'results_{run_pattern}',
            'run_pattern':'7.7_AuAu_*','group_by':[],'sampling':{'mode':'event'},'observables':self.definitions}))

    def test_all_from_config_uses_saved_current_summaries_and_correct_axes(self):
        import matplotlib.pyplot as plt
        from maf.plotting import figures
        original=figures.distribution
        scales={}
        labels={}
        legends={}
        def track(results,observable,**kwargs):
            fig,ax=original(results,observable,**kwargs)
            scales[observable]=ax.get_yscale()
            labels[observable]=(ax.get_xlabel(),ax.get_ylabel())
            legends[observable]=ax.get_legend_handles_labels()[1]
            return fig,ax
        with patch.object(figures,'distribution',side_effect=track), \
             patch('maf.cli.Result.load',side_effect=AssertionError('NPZ should not be loaded')), \
             patch('maf.io.oscar.OscarReader.__iter__',side_effect=AssertionError('OSCAR should not be read')), \
             contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(main(['plot','--all','--config',str(self.config)]),0)
        expected={d['name'] for d in self.definitions}
        self.assertEqual({p.stem for p in (self.results/'plots').glob('*.png')},expected)
        self.assertEqual(scales,{'spectrum':'log','rapidity':'linear','charged_eta':'linear','flow':'linear','empty_spectrum':'linear'})
        self.assertIn('Wrote 5 plots',output.getvalue())
        self.assertEqual(legends,{name:['7.7_AuAu_*'] for name in expected})
        self.assertEqual(labels['charged_eta'],(r'$\eta$',r'$\frac{\mathrm{d}N_{\mathrm{ch}}}{\mathrm{d}\eta}$'))
        self.assertEqual(labels['flow'][1],r'$v_{2}\{\mathrm{RP}\}$')
        self.assertEqual(labels['spectrum'][0],r'$p_{T}$ [GeV]')
        for path in (self.results/'plots').glob('*.png'):
            self.assertGreater(path.stat().st_size,1000)
        self.assertEqual(plt.get_fignums(),[])

    def test_pdf_from_directory_uses_relocated_index_and_explicit_destination(self):
        (self.results/'index.json').write_text(json.dumps({'groups':['/old/location/group-current/accumulators.npz']}))
        output=self.root/'pdf'
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['plot',str(self.results),'--all','--format','pdf','--output',str(output)]),0)
        self.assertEqual(len(list(output.glob('*.pdf'))),len(self.definitions))
        for path in output.glob('*.pdf'):
            self.assertTrue(path.read_bytes().startswith(b'%PDF'))

    def test_multiple_groups_are_separate_and_single_plot_still_works(self):
        second=self.results/'group-second'
        second.mkdir()
        (second/'summary.json').write_text((self.group/'summary.json').read_text())
        (self.results/'index.json').write_text(json.dumps({'groups':[
            str(self.group/'accumulators.npz'),str(second/'accumulators.npz')]}))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['plot',str(self.results/'index.json'),'--all']),0)
        for name in ('group-001','group-002'):
            self.assertEqual(len(list((self.results/'plots'/name).glob('*.png'))),len(self.definitions))
        single=self.root/'one.png'
        self.assertEqual(main(['plot',str(self.group/'summary.json'),'--observable','rapidity','--output',str(single)]),0)
        self.assertGreater(single.stat().st_size,1000)

    def test_saved_pattern_takes_priority_and_custom_labels_are_respected(self):
        from maf.cli import plot_summaries
        from maf.plotting.figures import distribution
        path=self.group/'summary.json'
        saved=json.loads(path.read_text())
        saved['run_pattern']='3.9_AuAu_*'
        path.write_text(json.dumps(saved))
        summary=plot_summaries([path])[0]
        _,ax=distribution([summary],'flow')
        self.assertEqual(ax.get_legend_handles_labels()[1],['3.9_AuAu_*'])
        _,ax=distribution([summary],'flow',labels=['Custom simulation'])
        self.assertEqual(ax.get_legend_handles_labels()[1],['Custom simulation'])
        legacy=plot_summaries([self.results])[0]
        legacy.pop('run_pattern')
        _,ax=distribution([legacy],'flow')
        self.assertEqual(ax.get_legend_handles_labels()[1],['7.7_AuAu'])
