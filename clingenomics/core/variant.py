"""Core genomic variant representation.

This is the atom the whole platform revolves around. Kept deliberately thin:
it describes *what* the variant is and where it sits, not *how* it is interpreted.
Interpretation lives in the germline/ and somatic/ engines and hangs off Evidence.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, field_validator


class Assembly(str, Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class VariantType(str, Enum):
    SNV = "SNV"
    INSERTION = "insertion"
    DELETION = "deletion"
    INDEL = "indel"
    MNV = "MNV"
    CNV = "CNV"
    SV = "SV"


class MolecularConsequence(str, Enum):
    """SO-term-aligned consequences, restricted to the ones interpretation cares about."""

    FRAMESHIFT = "frameshift_variant"
    STOP_GAINED = "stop_gained"
    STOP_LOST = "stop_lost"
    START_LOST = "start_lost"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_DONOR = "splice_donor_variant"
    SPLICE_REGION = "splice_region_variant"
    MISSENSE = "missense_variant"
    INFRAME_INDEL = "inframe_indel"
    SYNONYMOUS = "synonymous_variant"
    INTRONIC = "intron_variant"
    UTR = "UTR_variant"
    OTHER = "other"

    @property
    def is_predicted_lof(self) -> bool:
        """Null-type consequences that put PVS1 on the table (pre-decision-tree)."""
        return self in {
            MolecularConsequence.FRAMESHIFT,
            MolecularConsequence.STOP_GAINED,
            MolecularConsequence.SPLICE_ACCEPTOR,
            MolecularConsequence.SPLICE_DONOR,
            MolecularConsequence.START_LOST,
        }


class GenomicVariant(BaseModel):
    """A single normalized variant (VCF-style, left-aligned, one ALT per record)."""

    assembly: Assembly
    chrom: str
    pos: int = Field(gt=0, description="1-based position")
    ref: str
    alt: str
    variant_type: VariantType

    gene_symbol: Optional[str] = None
    transcript: Optional[str] = Field(None, description="e.g. NM_000059.4 / ENST...")
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: Optional[MolecularConsequence] = None

    # population / precomputed annotation carried alongside the variant
    gnomad_af: Optional[float] = Field(None, ge=0.0, le=1.0)
    spliceai_ds_max: Optional[float] = Field(
        None, ge=0.0, le=1.0, description="Max SpliceAI delta score across AG/AL/DG/DL"
    )

    @field_validator("chrom")
    @classmethod
    def _normalize_chrom(cls, v: str) -> str:
        return v[3:] if v.lower().startswith("chr") else v

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    def __str__(self) -> str:  # pragma: no cover - convenience
        g = self.gene_symbol or "?"
        return f"{g} {self.hgvs_c or self.key} ({self.assembly.value})"
