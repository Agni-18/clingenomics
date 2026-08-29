"""Evidence primitives shared across the whole platform.

Design decision (the load-bearing one): we separate a raw *Evidence* observation
from an interpreted *CriterionCall*. An RNA splicing outlier is Evidence; it
*emits* a PS3 CriterionCall at some computed strength. This indirection is what
lets the points-based ACMG system apply any criterion at a modulated strength,
which is exactly what RNA / functional data needs.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, model_validator


class Direction(str, Enum):
    """Whether evidence pushes toward pathogenic or benign."""

    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class Strength(Enum):
    """Evidence strengths and their point weight (magnitude only; sign from Direction).

    Point weights follow the ClinGen/Tavtigian Bayesian points framework.
    STANDALONE is BA1's slot; it carries the VeryStrong magnitude (8) but the
    classifier treats it as an automatic-benign override (see classifier).

    Members are kept distinct (STANDALONE vs VERY_STRONG both weigh 8 but must
    not alias each other), so the weight lives on `.points`, not the value.
    """

    STANDALONE = "standalone"
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"

    @property
    def points(self) -> int:
        return {
            Strength.STANDALONE: 8,
            Strength.VERY_STRONG: 8,
            Strength.STRONG: 4,
            Strength.MODERATE: 2,
            Strength.SUPPORTING: 1,
        }[self]

    @property
    def label(self) -> str:
        return {
            Strength.STANDALONE: "Stand-alone",
            Strength.VERY_STRONG: "Very Strong",
            Strength.STRONG: "Strong",
            Strength.MODERATE: "Moderate",
            Strength.SUPPORTING: "Supporting",
        }[self]


class EvidenceSource(str, Enum):
    POPULATION = "population"
    COMPUTATIONAL = "computational"
    FUNCTIONAL_RNA = "functional_rna"
    FUNCTIONAL_OTHER = "functional_other"
    SEGREGATION = "segregation"
    DE_NOVO = "de_novo"
    ALLELIC = "allelic"
    CLINICAL = "clinical"
    OTHER = "other"


class CriterionCall(BaseModel):
    """One ACMG/AMP criterion applied to a variant, at a (possibly modulated) strength.

    `code` is the ACMG label (PVS1, PS3, PM2, BS1, ...). `direction` and the
    criterion's *default* strength come from the registry, but `strength` here is
    the strength actually applied, so a PS3 can legitimately be recorded as
    Strong, Moderate or Supporting.
    """

    code: str
    direction: Direction
    strength: Strength
    source: EvidenceSource = EvidenceSource.OTHER
    rationale: str = ""
    # optional back-pointer to the raw evidence that produced this call
    evidence_id: Optional[str] = None

    @model_validator(mode="after")
    def _no_standalone_pathogenic(self) -> "CriterionCall":
        if self.strength is Strength.STANDALONE and self.direction is Direction.PATHOGENIC:
            raise ValueError("STANDALONE strength is only valid for benign (BA1)")
        return self

    @property
    def points(self) -> int:
        """Signed point contribution: + for pathogenic, - for benign."""
        sign = 1 if self.direction is Direction.PATHOGENIC else -1
        return sign * self.strength.points

    def __str__(self) -> str:  # pragma: no cover - convenience
        return f"{self.code} [{self.strength.label}] {self.points:+d}pt"
