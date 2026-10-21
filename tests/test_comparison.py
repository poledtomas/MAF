import copy
import tempfile
import unittest
from pathlib import Path
import numpy as np
from maf.analysis import analyze_run
from maf.io.discovery import discover
from maf.experiment.datasets import Dataset,load_csv
from maf.experiment.comparison import compare,rebin,chi_square,compatibility
from helpers import write_run,block,row


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        write_run(self.root,"7.7_AuAu_run_1",[block(0,[row(pt=.5)]),block(1,[row(pt=1.5),row(pt=.5,ident=1)])])
        run=discover(self.root)[0]
        run.metadata.update({"sampling":{"mode":"event"},"centrality":{"range":[0,5],"unit":"percent"},"parameters":{"x":.2}})
        self.definition={"name":"pt","kind":"pt","bins":[0,1,2],"selection":{"pdg":[2212],"frame":"cm","input_frame":"cm","particle_definition":"inclusive_final"},"units":"GeV^-1"}
        self.result=analyze_run(run,{"group_by":["energy","collision_system"],"observables":[self.definition]})
        self.metadata={"experiment":"synthetic validation","reference":"analytic test fixture, not experimental measurements", "collision_system":"Au+Au","energy":{"value":7.7,"unit":"GeV"},"centrality":run.metadata["centrality"],"definition":copy.deepcopy(self.definition)}
        self.dataset=Dataset(self.metadata,[0,1],[1,2],[1,.5],stat=[.1,.1])

    def test_compatible_comparison(self):
        self.assertEqual(compatibility(self.result,"pt",self.dataset),[])
        output=compare(self.result,"pt",self.dataset)
        np.testing.assert_allclose(output["ratio"],[1,1])

    def test_missing_centrality_rejected(self):
        for source in self.result.sources.values():
            source["metadata"].pop("centrality")
        with self.assertRaisesRegex(ValueError,"centrality"):
            compare(self.result,"pt",self.dataset)

    def test_particle_definition_and_acceptance_rejected(self):
        self.dataset.metadata["definition"]["selection"]["y"]=[-.5,.5]
        with self.assertRaisesRegex(ValueError,"acceptance"):
            compare(self.result,"pt",self.dataset)

    def test_rebin_integral_and_error(self):
        value,error,cov=rebin(self.result,"pt",[0,2])
        np.testing.assert_allclose(value,[.75])
        np.testing.assert_allclose(error,[.25])
        np.testing.assert_allclose(cov,[[.0625]])
        with self.assertRaises(ValueError):
            rebin(self.result,"pt",[0,.5,2])

    def test_csv_asymmetric(self):
        path=self.root/"data.csv"
        path.write_text("low,high,value,stat_low,stat_high,systematic\n0,1,1,.1,.2,.3\n1,2,.5,.05,.1,.1\n")
        dataset=load_csv(path,self.metadata)
        np.testing.assert_allclose(dataset.stat,[[.1,.2],[.05,.1]])

    def test_zero_data_ratio_undefined(self):
        self.dataset.value[0]=0
        self.assertTrue(np.isnan(compare(self.result,"pt",self.dataset)["ratio"][0]))

    def test_chi_square(self):
        self.assertAlmostEqual(chi_square([1,2],[0,0],np.eye(2))["chi2"],5)
        with self.assertRaises(ValueError):
            chi_square([1,2],[0,0],np.zeros((2,2)))

    def test_plots_from_saved_results(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Optional Matplotlib unavailable")
        matplotlib.use("Agg")
        from maf.plotting.figures import distribution,comparison
        figures=[distribution([self.result],"pt")[0],comparison(self.result,"pt",self.dataset)[0]]
        for i,fig in enumerate(figures):
            path=self.root/f"figure-{i}.png"
            fig.savefig(path)
            self.assertGreater(path.stat().st_size,1000)
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_hepdata_percent_asymmetric_and_unknown_labels(self):
        try:
            import yaml
        except ImportError:
            self.skipTest("Optional PyYAML unavailable")
        from maf.experiment.datasets import load_hepdata
        path=self.root/"table.yaml"
        table={"independent_variables":[{"values":[{"low":0,"high":1},{"low":1,"high":2}]}],
               "dependent_variables":[{"values":[{"value":1,"errors":[{"label":"stat","symerror":"10%"},{"label":"sys","asymerror":{"minus":-.2,"plus":.3}}]},
                                                  {"value":.5,"errors":[{"label":"stat","symerror":"10%"},{"label":"sys","asymerror":{"minus":-.1,"plus":.2}}]}]}]}
        path.write_text(yaml.safe_dump(table))
        dataset=load_hepdata(path,self.metadata)
        np.testing.assert_allclose(dataset.stat,[[.1,.1],[.05,.05]])
        np.testing.assert_allclose(dataset.systematic,[[.2,.3],[.1,.2]])
        table["dependent_variables"][0]["values"][0]["errors"][0]["label"]="total"
        path.write_text(yaml.safe_dump(table))
        with self.assertRaisesRegex(ValueError,"Ambiguous"):
            load_hepdata(path,self.metadata)

    def test_equivalent_transformed_frame_compatible(self):
        # The data are specified in CM; simulation lab momenta can be transformed to CM.
        from maf.experiment.comparison import definition_key
        from maf.statistics.aggregation import observable_config
        model=copy.deepcopy(self.definition)
        model["selection"].update({"frame":"cm","input_frame":"lab","cm_rapidity":.3})
        self.assertEqual(definition_key(model),definition_key(self.definition))

    def test_unknown_experimental_frame_rejected(self):
        metadata=copy.deepcopy(self.metadata)
        metadata["definition"]["selection"].update({"frame":"input","input_frame":"unknown"})
        with self.assertRaisesRegex(ValueError,"known lab or cm"):
            Dataset(metadata,[0,1],[1,2],[1,.5])
