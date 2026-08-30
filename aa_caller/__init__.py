"""Light wrapper to expose the RT variant caller as a package."""

from .cli import main
from .container import SamContainer
from .reference import FullReference
from .runner import VariantCallResult, call_variants, call_variants_from_args, runner_cli
from .sam import SamEntry
from .validators import parse_amplicons

__all__ = [
    "call_variants",
    "call_variants_from_args",
    "runner_cli",
    "VariantCallResult",
    "main",
    "FullReference",
    "SamContainer",
    "SamEntry",
    "parse_amplicons",
]
