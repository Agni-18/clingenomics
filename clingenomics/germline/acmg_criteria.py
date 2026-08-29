"""Registry of ACMG/AMP 2015 criteria with their direction and default strength.

The points-based framework lets any of these be applied at a *modulated* strength
(e.g. PVS1_Strong, PM2_Supporting, PS3_Moderate). `make_call` builds a
CriterionCall, defaulting to the criterion's canonical strength but accepting an
override, while enforcing that the override keeps the criterion's direction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from ..core.evidence import CriterionCall, Direction, EvidenceSource, Strength


@dataclass(frozen=True)
class CriterionSpec:
    code: str
    direction: Direction
    default_strength: Strength
    description: str
    deprecated: bool = False


def _p(code, strength, desc, deprecated=False):
    return CriterionSpec(code, Direction.PATHOGENIC, strength, desc, deprecated)


def _b(code, strength, desc, deprecated=False):
    return CriterionSpec(code, Direction.BENIGN, strength, desc, deprecated)


# Canonical ACMG/AMP 2015 criteria. Descriptions are abbreviated on purpose.
ACMG_CRITERIA: Dict[str, CriterionSpec] = {
    s.code: s
    for s in [
        # --- Pathogenic ---
        _p("PVS1", Strength.VERY_STRONG, "Null variant where LoF is a known disease mechanism"),
        _p("PS1", Strength.STRONG, "Same amino acid change as an established pathogenic variant"),
        _p("PS2", Strength.STRONG, "De novo (confirmed) in a patient, no family history"),
        _p("PS3", Strength.STRONG, "Well-established functional studies show damaging effect"),
        _p("PS4", Strength.STRONG, "Prevalence in affecteds significantly increased vs controls"),
        _p("PM1", Strength.MODERATE, "Located in a mutational hotspot / critical functional domain"),
        _p("PM2", Strength.MODERATE, "Absent/rare in population databases"),
        _p("PM3", Strength.MODERATE, "For recessive: detected in trans with a pathogenic variant"),
        _p("PM4", Strength.MODERATE, "Protein length change (in-frame indel / stop-loss)"),
        _p("PM5", Strength.MODERATE, "Novel missense at a residue with a known pathogenic change"),
        _p("PM6", Strength.MODERATE, "Assumed de novo (parentage not confirmed)"),
        _p("PP1", Strength.SUPPORTING, "Co-segregation with disease in affected family members"),
        _p("PP2", Strength.SUPPORTING, "Missense in a gene with low benign missense rate"),
        _p("PP3", Strength.SUPPORTING, "Multiple in-silico lines support a deleterious effect"),
        _p("PP4", Strength.SUPPORTING, "Phenotype/family history highly specific for the gene"),
        _p("PP5", Strength.SUPPORTING, "Reputable source reports pathogenic", deprecated=True),
        # --- Benign ---
        _b("BA1", Strength.STANDALONE, "Allele frequency too high for the disorder (stand-alone)"),
        _b("BS1", Strength.STRONG, "Allele frequency greater than expected for the disorder"),
        _b("BS2", Strength.STRONG, "Observed in healthy adults (full-penetrance disorder)"),
        _b("BS3", Strength.STRONG, "Well-established functional studies show no damaging effect"),
        _b("BS4", Strength.STRONG, "Lack of segregation in affected family members"),
        _b("BP1", Strength.SUPPORTING, "Missense in a gene where only truncating cause disease"),
        _b("BP2", Strength.SUPPORTING, "In trans/cis with a pathogenic variant (dominant)"),
        _b("BP3", Strength.SUPPORTING, "In-frame indel in a repetitive region without known function"),
        _b("BP4", Strength.SUPPORTING, "Multiple in-silico lines support no impact"),
        _b("BP5", Strength.SUPPORTING, "Found in a case with an alternate molecular cause"),
        _b("BP6", Strength.SUPPORTING, "Reputable source reports benign", deprecated=True),
        _b("BP7", Strength.SUPPORTING, "Synonymous with no predicted splice impact"),
    ]
}


def make_call(
    code: str,
    *,
    strength: Optional[Strength] = None,
    rationale: str = "",
    source: Optional[EvidenceSource] = None,
    evidence_id: Optional[str] = None,
) -> CriterionCall:
    """Build a CriterionCall for an ACMG code, optionally at a modulated strength."""
    spec = ACMG_CRITERIA.get(code)
    if spec is None:
        raise KeyError(f"Unknown ACMG criterion: {code!r}")

    applied = strength or spec.default_strength

    # A modulated strength must not flip a benign criterion into the pathogenic
    # magnitude space via STANDALONE, and vice versa.
    if spec.direction is Direction.PATHOGENIC and applied is Strength.STANDALONE:
        raise ValueError(f"{code}: STANDALONE is reserved for BA1")

    return CriterionCall(
        code=code,
        direction=spec.direction,
        strength=applied,
        source=source or EvidenceSource.OTHER,
        rationale=rationale or spec.description,
        evidence_id=evidence_id,
    )
