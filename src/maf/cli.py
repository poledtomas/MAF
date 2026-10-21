"""Command-line interface; all operations are also ordinary Python functions."""
import argparse
import json
import logging
from pathlib import Path
import sys
from ._version import __version__
from .analysis import analyze, merge_results, write_result_groups
from .io.config import load_config
from .io.dat import load_summary_dat, read_dat
from .io.discovery import discover
from .statistics.aggregation import Result, json_clean


COMMANDS = ("run", "discover", "analyze", "merge", "plot", "compare")


def parser():
    root = argparse.ArgumentParser(
        prog="maf",
        description="MAF: Modular Analysis Framework for SMASH-vHLLE simulations.",
        epilog="Quick start: maf --config analysis.yaml; maf plot --all --config analysis.yaml",
    )
    root.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    root.add_argument("--verbose", action="store_true", help="Enable debug logging")
    sub = root.add_subparsers(dest="command", required=True, title="commands")

    def command(name, help):
        child = sub.add_parser(name, help=help, description=help)
        child.add_argument("--verbose", action="store_true", default=argparse.SUPPRESS,
                           help="Enable debug logging")
        return child

    p = command("run", "Run enabled observables from a configuration (default with --config)")
    p.add_argument("--config", "-c", required=True, help="YAML, JSON or TOML configuration")
    p.add_argument("--dry-run", action="store_true", help="Validate settings and discover runs without writing")
    p.add_argument("--workers", type=int, help="Override the worker count from the configuration")
    p.add_argument("--no-resume", action="store_true", help="Recompute cached runs")
    p.add_argument("--fail-fast", action="store_true", help="Stop at the first failed run")

    for name, help in (("discover", "List OSCAR inputs and their simulation metadata"),
                       ("analyze", "Analyze an explicit input root using an expanded configuration")):
        p = command(name, help)
        p.add_argument("root", help="Input directory or OSCAR file")
        p.add_argument("--config", required=name == "analyze", help="Expanded analysis or discovery configuration")
        p.add_argument("--manifest", help="Input metadata manifest")
        if name == "analyze":
            p.add_argument("--output", required=True, help="Results directory")
            p.add_argument("--workers", type=int, default=1, help="Number of concurrent input files")
            p.add_argument("--no-resume", action="store_true", help="Recompute cached runs")
            p.add_argument("--dry-run", action="store_true", help="Discover inputs without writing")
            p.add_argument("--fail-fast", action="store_true", help="Stop at the first failed run")
    p = command("merge", "Merge saved contributions from distinct sources")
    p.add_argument("inputs", nargs="+", help="Explicit per-run or disjoint aggregate NPZ files")
    p.add_argument("--output", required=True, help="Results directory")
    p = command("plot", "Plot observables from saved results without reading events")
    p.add_argument("inputs", nargs="*", help="Results directory, index.json, DAT, NPZ or summary.json files")
    p.add_argument("--config", "-c", help="Read the results directory and style from an analysis configuration")
    choice = p.add_mutually_exclusive_group(required=True)
    choice.add_argument("--observable", help="Name of one saved observable")
    choice.add_argument("--all", dest="all_observables", action="store_true", help="Plot every saved observable")
    p.add_argument("--output", help="Figure filename, or directory with --all (default: results/plots)")
    p.add_argument("--format", choices=("png", "pdf"), default="png", help="Image format for --all")
    p.add_argument("--logy", action="store_true", help="Force logarithmic y; --all uses it automatically for spectra")
    p = command("compare", "Compare a saved observable with an experimental table")
    p.add_argument("result", help="Saved NPZ result")
    p.add_argument("--observable", required=True, help="Name of the saved observable")
    p.add_argument("--experiment", required=True, help="CSV or HEPData YAML table")
    p.add_argument("--metadata", required=True, help="Experimental metadata JSON/YAML")
    p.add_argument("--output", required=True, help="Numerical comparison JSON")
    p.add_argument("--figure", help="Optional comparison figure filename")
    return root


def _normalize_arguments(arguments):
    """Keep `maf --config ...` as shorthand for the explicit run command."""
    arguments = list(arguments)
    offset = 0
    while offset < len(arguments) and arguments[offset] == "--verbose":
        offset += 1
    if offset < len(arguments) and arguments[offset] not in (*COMMANDS, "-h", "--help", "--version"):
        arguments.insert(offset, "run")
    return arguments


def plot_summaries(inputs):
    """Use current index entries and saved summaries, without rereading events."""
    summaries=[]
    for value in inputs:
        path=Path(value)
        paths=[path]
        run_pattern=None
        if path.is_dir() or path.name=="index.json":
            directory=path if path.is_dir() else path.parent
            index=load_config(directory/"index.json")
            run_pattern=index.get("run_pattern")
            if "summaries" in index:
                paths=[directory/entry for entry in index["summaries"]]
            else:
                # Older results had direct group directories, including moved copies.
                paths=[directory/Path(entry).parent.name/"summary.json" for entry in index["groups"]]
        for source in paths:
            if source.suffix==".npz":
                summary=Result.load(source).summary()
            elif source.suffix==".dat":
                name,data=read_dat(source)
                summary={"group":{},"run_pattern":data.get("run_pattern"),"observables":{name:data}}
            else:
                summary=load_summary_dat(source)
            if not summary.get("run_pattern") and run_pattern:
                summary["run_pattern"]=run_pattern
            # Legacy outputs predate saved run patterns; use their directory label.
            for parent in list(source.parents)[:3]:
                if parent.name.startswith("results_"):
                    summary["simulation_label"]=parent.name[len("results_"):]
                    break
            summaries.append(summary)
    if not summaries:
        raise ValueError("No saved result groups to plot")
    return summaries


def main(argv=None):
    arguments = sys.argv[1:] if argv is None else argv
    args=parser().parse_args(_normalize_arguments(arguments))
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,format="%(levelname)s: %(message)s")
    try:
        if args.command == "run":
            from .config import load_job
            logging.info("Loading configuration: %s", args.config)
            job = load_job(args.config)
            if args.workers is not None:
                if args.workers < 1:
                    raise ValueError("workers must be positive")
                job["workers"] = args.workers
            names = [definition["name"] for definition in job["config"]["observables"]]
            logging.info("Enabled calculations: %s", ", ".join(names))
            results = analyze(**job, dry_run=args.dry_run, resume=not args.no_resume,
                              fail_fast=args.fail_fast)
            if args.dry_run:
                print(json.dumps({"calculations": names, "runs": [str(run.path) for run in results]}, indent=2))
            else:
                events = sum(unit["events"] for result in results for unit in result.units.values())
                files = sum(len(result.sources) for result in results)
                print(f"Wrote {len(results)} groups ({events} selected events from {files} files) to {job['output']}")
                if json.loads((Path(job["output"])/"failures.json").read_text()):
                    return 2
        elif args.command in ("discover","analyze"):
            config=load_config(args.config) if args.config else {}
            manifest=load_config(args.manifest) if args.manifest else None
            if args.command=="discover":
                runs=discover(args.root,rules=config.get("parsing_rules",()),manifest=manifest,config_fields=config.get("config_fields",()))
                print(json.dumps([{"oscar":str(r.oscar),"metadata":r.metadata} for r in runs],indent=2))
            else:
                results=analyze(args.root,config,args.output,workers=args.workers,resume=not args.no_resume,
                    dry_run=args.dry_run,manifest=manifest,fail_fast=args.fail_fast)
                if args.dry_run:
                    print(json.dumps([{"path":str(r.path),"metadata":r.metadata} for r in results],indent=2))
                else:
                    print(f"Wrote {len(results)} metadata groups to {args.output}")
                    if json.loads((Path(args.output)/"failures.json").read_text()):
                        return 2
        elif args.command=="merge":
            write_result_groups(merge_results([Result.load(p) for p in args.inputs]),args.output)
        elif args.command=="plot":
            from .plotting.figures import distribution,all_distributions,_plt
            inputs=args.inputs
            rc_params={}
            run_pattern=None
            if args.config:
                if inputs:
                    raise ValueError("Choose --config or explicit plot inputs")
                from .config import load_job,load_plot_style
                job=load_job(args.config)
                inputs=[job["output"]]
                run_pattern=job["run_pattern"]
                rc_params=load_plot_style(args.config)
            if not inputs:
                raise ValueError("Plot requires --config or saved result inputs")
            if not args.all_observables and not args.output:
                raise ValueError("Plotting one observable requires --output with a figure filename")
            results=plot_summaries(inputs)
            for result in results:
                if not result.get("run_pattern") and run_pattern:
                    result["run_pattern"]=run_pattern
            with _plt().rc_context(rc_params):
                if args.all_observables:
                    output=Path(args.output) if args.output else (Path(inputs[0]) if Path(inputs[0]).is_dir() else Path(inputs[0]).parent)/"plots"
                    count=0
                    for i,result in enumerate(results):
                        directory=output if len(results)==1 else output/f"group-{i+1:03d}"
                        count+=len(all_distributions(result,directory,format=args.format,logy=args.logy))
                    print(f"Wrote {count} plots to {output}")
                else:
                    fig,_=distribution(results,args.observable,logy=args.logy)
                    try:
                        fig.savefig(args.output,metadata={"Creator":"maf"})
                    finally:
                        _plt().close(fig)
        elif args.command=="compare":
            from .experiment.datasets import load_csv,load_hepdata
            from .experiment.comparison import compare
            result=Result.load(args.result)
            loader=load_csv if Path(args.experiment).suffix==".csv" else load_hepdata
            dataset=loader(args.experiment,args.metadata)
            comparison=compare(result,args.observable,dataset)
            Path(args.output).write_text(json.dumps(json_clean(comparison),indent=2,allow_nan=False)+"\n")
            if args.figure:
                from .plotting.figures import comparison as plot_comparison
                fig,_=plot_comparison(result,args.observable,dataset)
                fig.savefig(args.figure)
    except (ValueError,TypeError,KeyError,OSError,ImportError) as exc:
        logging.error("%s",exc)
        return 1
    return 0
