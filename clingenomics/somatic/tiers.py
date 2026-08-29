"""Somatic variant tiering — AMP/ASCO/CAP 2017 four-tier system.

Somatic classification is a *different problem* from germline ACMG: it ranks
variants by clinical actionability in a tumor context, not by pathogenicity.
This branch is kept fully separate from the germline engine on purpose; the only
thing they share is the GenomicVariant.

Tier I    Strong clinical significance   (evidence level A or B)
Tier II   Potential clinical significance (evidence level C or D)
Tier III  Unknown clinical significance
Tier IV   Benign / likely benign

OncoKB therapeutic levels can be supplied and are mapped onto AMP evidence
levels for tiering.
"""

from __future__ import annotations

from enum import Enum
from typing import List, Optional

from pydantic import BaseModel


class EvidenceLevel(str, Enum):
    """AMP/ASCO/CAP evidence levels."""

    A = "A"  # FDA-approved / professional-guideline biomarker in this tumor type
    B = "B"  # Well-powered studies with expert consensus
    C = "C"  # Approved for a *different* tumor type, or multiple small studies
    D = "D"  # Preclinical / case reports / plausible biological rationale


class SomaticTier(str, Enum):
    TIER_I = "Tier I - Strong clinical significance"
    TIER_II = "Tier II - Potential clinical significance"
    TIER_III = "Tier III - Unknown clinical significance"
    TIER_IV = "Tier IV - Benign or likely benign"


class OncoKBLevel(str, Enum):
    LEVEL_1 = "1"
    LEVEL_2 = "2"
    LEVEL_3A = "3A"
    LEVEL_3B = "3B"
    LEVEL_4 = "4"
    R1 = "R1"
    R2 = "R2"


_ONCOKB_TO_AMP = {
    OncoKBLevel.LEVEL_1: EvidenceLevel.A,
    OncoKBLevel.LEVEL_2: EvidenceLevel.A,
    OncoKBLevel.R1: EvidenceLevel.A,
    OncoKBLevel.LEVEL_3A: EvidenceLevel.B,
    OncoKBLevel.LEVEL_3B: EvidenceLevel.C,
    OncoKBLevel.LEVEL_4: EvidenceLevel.D,
    OncoKBLevel.R2: EvidenceLevel.C,
}


class ClinicalAssertion(BaseModel):
    """One actionability assertion for a somatic variant in a tumor context."""

    description: str
    evidence_level: EvidenceLevel
    therapy: Optional[str] = None
    tumor_type: Optional[str] = None
    resistance: bool = False


class SomaticResult(BaseModel):
    tier: SomaticTier
    assertions: List[ClinicalAssertion]
    is_likely_benign: bool = False
    rationale: str = ""


def classify_somatic(
    assertions: List[ClinicalAssertion],
    *,
    is_likely_benign: bool = False,
) -> SomaticResult:
    """Assign an AMP/ASCO/CAP tier from a variant's clinical assertions.

    `is_likely_benign` is set upstream (e.g. common germline polymorphism leaked
    into a tumor-only call, or population-frequency filter) and short-circuits to
    Tier IV.
    """
    if is_likely_benign:
        return SomaticResult(
            tier=SomaticTier.TIER_IV,
            assertions=assertions,
            is_likely_benign=True,
            rationale="Benign/likely benign — population frequency or germline filter",
        )

    if not assertions:
        return SomaticResult(
            tier=SomaticTier.TIER_III,
            assertions=[],
            rationale="No clinical actionability evidence found",
        )

    levels = {a.evidence_level for a in assertions}
    if levels & {EvidenceLevel.A, EvidenceLevel.B}:
        tier = SomaticTier.TIER_I
        rationale = "At least one Level A/B actionability assertion"
    elif levels & {EvidenceLevel.C, EvidenceLevel.D}:
        tier = SomaticTier.TIER_II
        rationale = "Level C/D actionability only"
    else:  # pragma: no cover - defensive
        tier = SomaticTier.TIER_III
        rationale = "Evidence present but unclassifiable"

    return SomaticResult(tier=tier, assertions=assertions, rationale=rationale)


def assertion_from_oncokb(
    level: OncoKBLevel,
    *,
    therapy: Optional[str] = None,
    tumor_type: Optional[str] = None,
    description: str = "",
) -> ClinicalAssertion:
    """Convenience: build a ClinicalAssertion from an OncoKB therapeutic level."""
    return ClinicalAssertion(
        description=description or f"OncoKB level {level.value}",
        evidence_level=_ONCOKB_TO_AMP[level],
        therapy=therapy,
        tumor_type=tumor_type,
        resistance=level in {OncoKBLevel.R1, OncoKBLevel.R2},
    )
