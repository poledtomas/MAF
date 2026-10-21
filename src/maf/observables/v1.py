"""v1 relative to a known reaction plane, with its original function interface."""
from .flow_reaction_plane import calculate as calculate_harmonic


def calculate(px, py, coordinate, bins, plane_angle, weighting="particle"):
    return calculate_harmonic(px, py, coordinate, bins, plane_angle, 1, weighting)
