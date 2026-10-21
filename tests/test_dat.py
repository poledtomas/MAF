"""DAT metadata and numeric round trips, graph labels and paired net protons."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from maf.analysis import analyze_run
from maf.cli import plot_summaries
from maf.io.dat import read_dat
from maf.io.discovery import discover
from helpers import block, row, write_run


class DatTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        def particles(n,m):
            return [row(pdg=2212,ident=i) for i in range(n)]+[row(pdg=-2212,charge=-1,ident=n+i) for i in range(m)]
        path=write_run(self.root/'runs','7.7_AuAu_run_1',[
            block(0,particles(2,1)),block(1,particles(4,2)),block(2,[])])
        run=discover(path)[0]
        run.metadata['sampling']={'mode':'event'}
        self.result=analyze_run(run,{'group_by':[],'observables':[
            {'name':'spectrum','kind':'pt','bins':[0,1,2],'selection':{'pdg':[2212]}},
            {'name':'net_proton_dndy','kind':'net_proton','bins':[-1,1]},
            {'name':'absent_mean','kind':'mean_pt','selection':{'pdg':[321]}}]})
        self.result.run_pattern='7.7_AuAu_run_*'
        self.output=self.root/'output'
        self.directory=self.output/'results'
        self.result.export(self.directory)

    def test_numeric_round_trip_and_quantity_unit_headers(self):
        summary=self.result.summary()
        self.assertEqual(len(list(self.directory.glob('*.dat'))),3)
        for name,data in summary['observables'].items():
            path=self.directory/f'{name}.dat'
            headers=path.read_text()
            for field in ('x_quantity','x_unit','y_quantity','y_unit','columns','run_pattern'):
                self.assertIn(f'# {field}:',headers)
            saved_name,saved=read_dat(path)
            self.assertEqual(saved_name,name)
            self.assertEqual(saved['axes'],data['axes'])
            self.assertEqual(saved['units'],data['units'])
            self.assertEqual(saved['run_pattern'],'7.7_AuAu_run_*')
            for field in ('value','stat_error'):
                np.testing.assert_allclose(saved[field],np.asarray(data[field],float),equal_nan=True)
            self.assertEqual(np.loadtxt(path,ndmin=2).shape[1],5)
        _,spectrum=read_dat(self.directory/'spectrum.dat')
        self.assertEqual(spectrum['axes']['x'],{'quantity':'p_T','unit':'GeV'})
        self.assertEqual(spectrum['axes']['y'],{'quantity':'dN/dp_T','unit':'GeV^-1'})
        _,absent=read_dat(self.directory/'absent_mean.dat')
        self.assertTrue(np.isnan(absent['value'][0]))

    def test_net_protons_use_paired_event_difference_and_empty_event(self):
        # Event net counts are 1, 2, 0; width=2. Paired covariance must be retained.
        summary=self.result.summary()
        net=summary['observables']['net_proton_dndy']
        self.assertEqual(summary['events'],3)
        np.testing.assert_allclose(net['value'],[.5])
        np.testing.assert_allclose(net['stat_error'],[1/(2*np.sqrt(3))])
        self.assertEqual(net['axes']['y'],{'quantity':'d(N_p-N_pbar)/dy','unit':'1'})

    def test_plot_reads_numeric_values_and_labels_from_dat_header(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('Optional Matplotlib unavailable')
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        self.addCleanup(plt.close,'all')
        path=self.directory/'net_proton_dndy.dat'
        text=path.read_text().replace('# x_quantity: y','# x_quantity: rapidity y')
        text=text.replace('# y_quantity: d(N_p-N_pbar)/dy','# y_quantity: net-proton density')
        lines=text.splitlines()
        fields=lines[-1].split()
        fields[3]='42'
        lines[-1]='\t'.join(fields)
        path.write_text('\n'.join(lines)+'\n')
        (self.output/'index.json').write_text(json.dumps({
            'groups':[str(self.directory/'accumulators.npz')], 'summaries':['results/summary.json']}))
        summary=plot_summaries([self.output])[0]
        from maf.plotting.figures import distribution
        fig,ax=distribution([summary],'net_proton_dndy')
        self.assertEqual(ax.get_xlabel(),'rapidity y')
        self.assertEqual(ax.get_ylabel(),'net-proton density')
        np.testing.assert_allclose(ax.lines[0].get_ydata(),[42])
        fig,spectrum_ax=distribution([summary],'spectrum')
        self.assertEqual(spectrum_ax.get_xlabel(),r'$p_{T}$ [GeV]')
        self.assertEqual(spectrum_ax.get_ylabel(),r'$\frac{\mathrm{d}N}{\mathrm{d}p_{T}}$ [$\mathrm{GeV}^{-1}$]')
        direct=plot_summaries([path])[0]
        self.assertEqual(direct['run_pattern'],'7.7_AuAu_run_*')
        self.assertEqual(direct['observables']['net_proton_dndy']['axes'],summary['observables']['net_proton_dndy']['axes'])

    def test_different_coordinate_units_cannot_be_overlaid(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('Optional Matplotlib unavailable')
        matplotlib.use('Agg')
        import copy
        from maf.plotting.figures import distribution
        first=self.result.summary()
        second=copy.deepcopy(first)
        second['observables']['spectrum']['axes']['x']['unit']='MeV'
        with self.assertRaisesRegex(ValueError,'axis units'):
            distribution([first,second],'spectrum')

    def test_missing_units_and_stale_dat_are_rejected(self):
        path=self.directory/'spectrum.dat'
        original=path.read_text()
        path.write_text(original.replace('# y_unit: GeV^-1\n',''))
        with self.assertRaisesRegex(ValueError,'Missing DAT'):
            read_dat(path)
        path.write_text(original.replace(f'# analysis_signature: {self.result.signature}', '# analysis_signature: old'))
        with self.assertRaisesRegex(ValueError,'does not belong'):
            plot_summaries([self.directory/'summary.json'])

    def test_export_replaces_its_csv_and_preserves_other_files(self):
        (self.directory/'spectrum.csv').write_text('previous exported CSV')
        (self.directory/'external.csv').write_text('external data')
        self.result.export(self.directory)
        self.assertFalse((self.directory/'spectrum.csv').exists())
        self.assertTrue((self.directory/'spectrum.dat').exists())
        self.assertEqual((self.directory/'external.csv').read_text(),'external data')
