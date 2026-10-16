"""Backward-compatibility shim — re-exports everything from the split modules.

Existing code that does ``from aa_caller.app import X`` continues to work.
"""

from .cli import main  # noqa: F401
from .constants import (  # noqa: F401
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_PROTEIN,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
)
from .container import SamContainer  # noqa: F401
from .genetic_code import GENETIC_CODE, codon_to_aminoacid  # noqa: F401
from .models import Amplicon, Protein, Qual, Variant  # noqa: F401
from .reference import (  # noqa: F401
    FullReference,
    parse_protein_spec,
    resolve_protein_targets,
    resolve_proteins,
)
from .sam import SamEntry  # noqa: F401
from .validators import (  # noqa: F401
    parse_amplicons,
    validate_amplicon_file,
    validate_reference_file,
    validate_sam_file,
)
