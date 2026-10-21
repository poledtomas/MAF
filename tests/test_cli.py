"""CLI compatibility and progress reporting across fresh, cached and failed inputs."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from maf import __version__
from maf.cli import main as cli_main
from maf.main import main
from maf.statistics.aggregation import Result
from helpers import block, row, write_run


class CommandLineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        write_run(self.root / 'runs', '7.7_AuAu_run_1', [block(0, [row()]), block(1, [])])
        write_run(self.root / 'runs', '7.7_AuAu_run_2', [block(0, [row(), row(ident=1)])])
        self.config = {
            'input': 'runs', 'output': 'output', 'run_pattern': '7.7_AuAu_run_*',
            'group_by': [], 'sampling': {'mode': 'event'}, 'reader': {'errors': 'strict'},
            'observables': {'yields': {'pdg': [2212]}},
        }
        self.config_path = self.root / 'analysis.json'
        self.config_path.write_text(json.dumps(self.config))

    def invoke(self, arguments):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertLogs(level='INFO') as logs:
            status = main(arguments)
        return status, output.getvalue(), '\n'.join(logs.output)

    def test_help_lists_commands_and_version_needs_no_configuration(self):
        for entry in (main, cli_main):
            for arguments in (['--help'], ['--version']):
                with self.subTest(entry=entry.__module__, arguments=arguments):
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as result:
                        entry(arguments)
                    self.assertEqual(result.exception.code, 0)
                    if arguments == ['--version']:
                        self.assertEqual(output.getvalue().strip(), f'maf {__version__}')
                    else:
                        for command in ('run', 'discover', 'analyze', 'merge', 'plot', 'compare'):
                            self.assertIn(command, output.getvalue())
                        self.assertIn('--config analysis.yaml', output.getvalue())
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as result:
            main(['plot', '--help'])
        self.assertEqual(result.exception.code, 0)
        self.assertIn('--all', output.getvalue())

    def test_explicit_and_implicit_run_keep_dry_run_json_and_do_not_write(self):
        outputs = []
        for prefix in ([], ['run'], ['--verbose'], ['--verbose', 'run'], ['run', '--verbose']):
            status, output, logs = self.invoke(prefix + ['--config', str(self.config_path), '--dry-run'])
            self.assertEqual(status, 0)
            outputs.append(json.loads(output))
            self.assertIn('Found 2 files', logs)
            self.assertFalse((self.root / 'output').exists())
        self.assertTrue(all(value == outputs[0] for value in outputs))
        self.assertEqual(outputs[0]['calculations'], ['yields_2212'])
        self.assertEqual(len(outputs[0]['runs']), 2)

    def test_serial_and_parallel_progress_count_cache_hits_and_forced_work(self):
        for workers in (1, 2):
            with self.subTest(workers=workers):
                self.config['output'] = f'output-{workers}'
                self.config_path.write_text(json.dumps(self.config))
                args = ['--config', str(self.config_path), '--workers', str(workers)]
                status, output, logs = self.invoke(args)
                self.assertEqual(status, 0)
                self.assertIn('3 selected events from 2 files', output)
                self.assertIn('Progress: 2/2 files (100.0%)', logs)
                self.assertIn('analyzed=2 cached=0 failed=0', logs)
                self.assertIn('elapsed=', logs)
                directory = self.root / self.config['output']
                before = Result.load(directory / 'results/accumulators.npz').summary()
                checkpoints = {p: p.stat().st_mtime_ns for p in (directory / 'cache').glob('*.npz')}
                # On a resumed serial run, even entering the reader is a regression.
                with patch('maf.io.oscar.OscarReader.__iter__', side_effect=AssertionError('cached input was read')):
                    status, _, logs = self.invoke(['run', '--config', str(self.config_path)])
                self.assertEqual(status, 0)
                self.assertIn('analyzed=0 cached=2 failed=0', logs)
                self.assertEqual(checkpoints, {p: p.stat().st_mtime_ns for p in checkpoints})
                status, _, logs = self.invoke(args)
                self.assertEqual(status, 0)
                self.assertIn('analyzed=0 cached=2 failed=0', logs)
                self.assertEqual(before, Result.load(directory / 'results/accumulators.npz').summary())
                status, _, logs = self.invoke(args + ['--no-resume'])
                self.assertEqual(status, 0)
                self.assertIn('analyzed=2 cached=0 failed=0', logs)
                self.assertEqual(before, Result.load(directory / 'results/accumulators.npz').summary())

    def test_partial_failure_is_counted_and_returns_nonzero_with_usable_results(self):
        write_run(self.root / 'runs', '7.7_AuAu_run_bad', [])
        status, _, logs = self.invoke(['--config', str(self.config_path)])
        self.assertEqual(status, 2)
        self.assertIn('Progress: 3/3 files (100.0%)', logs)
        self.assertIn('analyzed=2 cached=0 failed=1', logs)
        directory = self.root / 'output'
        failures = json.loads((directory / 'failures.json').read_text())
        self.assertEqual(len(failures), 1)
        self.assertIn('run_bad', failures[0]['path'])
        self.assertEqual(Result.load(directory / 'results/accumulators.npz').summary()['events'], 3)

    def test_fail_fast_reports_stopped_progress(self):
        for source in (self.root / 'runs').glob('*/particle_lists.oscar'):
            source.write_text('not an OSCAR file\n')
        status, _, logs = self.invoke(['--config', str(self.config_path), '--fail-fast'])
        self.assertEqual(status, 1)
        self.assertIn('Stopped: 1/2 files (50.0%)', logs)
        self.assertIn('analyzed=0 cached=0 failed=1', logs)
        self.assertFalse((self.root / 'output/results').exists())


if __name__ == '__main__':
    unittest.main()
