"""Excitatory / inhibitory class from a transcriptomic type label."""
import re

# Mouse (Tasic et al. 2018) subclasses plus the human "Inh" prefix.
INHIBITORY_PREFIXES = {"lamp5", "sncg", "vip", "sst", "pvalb", "serpinf1", "meis2", "inh"}


def classify_ttype(label: str | None) -> str | None:
    """Return "E", "I", or None for labels that are missing or not neuronal."""
    if not isinstance(label, str) or not label.strip():
        return None
    first = label.strip().split()[0]
    if first.lower() in INHIBITORY_PREFIXES:
        return "I"
    if first == "Exc" or re.match(r"L\d", first):
        return "E"
    return None
