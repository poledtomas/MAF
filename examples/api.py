"""Run from the command line: python examples/api.py ROOT CONFIG OUTPUT."""
import sys
from maf.analysis import analyze
from maf.io.config import load_config
from maf.plotting.figures import distribution

if __name__ == "__main__":
    results = analyze(sys.argv[1], load_config(sys.argv[2]), sys.argv[3], workers=1)
    fig, ax = distribution(results, "proton_dndy")
    fig.savefig("proton_dndy.pdf")
