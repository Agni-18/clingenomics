"""Points-based ACMG/AMP germline classifier (ClinGen/Tavtigian Bayesian framework).

Sum signed points across all applied criteria, then map to a 5-tier call:

    Pathogenic          total >= 10
    Likely Pathogenic    6 <= total <= 9
    Uncertain (VUS)      0 <= total <= 5
    Likely Benign       -6 <= total <= -1
    Benign              total <= -7

Special handling:
  * BA1 (stand-alone benign) forces Benign, per ClinGen, unless explicitly
    exempted by the lab (recorded, not silently ignored).
  * Simultaneous strong pathogenic and strong benign evidence is surfaced as a
    conflict flag; the numeric call still returns, but the report must show it.
"""

from __future__ import annotations

from enum import Enum
from typing import List

from pydantic import BaseModel

from ..core.evidence import CriterionCall, Direction, Strength


class ACMGClassification(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    UNCERTAIN = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


def _tier_from_points(points: int) -> ACMGClassification:
    if points >= 10:
        return ACMGClassification.PATHOGENIC
    if points >= 6:
        return ACMGClassification.LIKELY_PATHOGENIC
    if points >= 0:
        return ACMGClassification.UNCERTAIN
    if points >= -6:
        return ACMGClassification.LIKELY_BENIGN
    return ACMGClassification.BENIGN


class ClassificationResult(BaseModel):
    classification: ACMGClassification
    points: int
    applied: List[CriterionCall]
    flags: List[str] = []

    @property
    def codes(self) -> List[str]:
        return [c.code for c in self.applied]

    def summary(self) -> str:
        codes = ", ".join(str(c) for c in self.applied) or "no criteria"
        return f"{self.classification.value} ({self.points:+d} pts) — {codes}"


def classify(
    calls: List[CriterionCall],
    *,
    honor_ba1_override: bool = True,
) -> ClassificationResult:
    """Classify a variant from its applied criteria."""
    flags: List[str] = []

    # de-duplicate by code, keeping the strongest-magnitude call per code
    by_code: dict[str, CriterionCall] = {}
    for c in calls:
        prev = by_code.get(c.code)
        if prev is None or abs(c.points) > abs(prev.points):
            by_code[c.code] = c
    applied = list(by_code.values())

    total = sum(c.points for c in applied)

    has_path_strong = any(
        c.direction is Direction.PATHOGENIC and c.strength.points >= Strength.STRONG.points
        for c in applied
    )
    has_benign_strong = any(
        c.direction is Direction.BENIGN
        and c.strength is not Strength.STANDALONE
        and c.strength.points >= Strength.STRONG.points
        for c in applied
    )
    if has_path_strong and has_benign_strong:
        flags.append("CONFLICTING: strong pathogenic and strong benign evidence both present")

    ba1 = any(c.code == "BA1" for c in applied)
    if ba1 and honor_ba1_override:
        flags.append("BA1 stand-alone benign applied — classification forced to Benign")
        return ClassificationResult(
            classification=ACMGClassification.BENIGN,
            points=total,
            applied=applied,
            flags=flags,
        )

    return ClassificationResult(
        classification=_tier_from_points(total),
        points=total,
        applied=applied,
        flags=flags,
    )
