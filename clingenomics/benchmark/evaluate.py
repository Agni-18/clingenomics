"""ClinVar concordance benchmark.

Runs the engine over ClinVar-classified variants and compares each engine call to
ClinVar's, producing:

  * strict concordance    — exact 5-tier match
  * clinical concordance  — same direction (P↔LP agree, B↔LB agree, VUS↔VUS)
  * contradiction rate    — opposite direction (the safety-critical metric)
  * a mismatch taxonomy   — undercall / overcall / contradiction
  * a confusion matrix    — ClinVar tier × engine tier

Why the taxonomy matters more than the raw match rate: this engine deliberately
lands at VUS when a VCF row lacks the case-level data (segregation, functional
studies) ClinVar's experts had. An engine=VUS / ClinVar=Pathogenic disagreement
is an *undercall* (expected, safe), not the same as calling the opposite
direction. Contradictions should be ~0; if not, that's a real bug.

Input is ClinVar's own VCF format (CLNSIG / CLNREVSTAT INFO fields) joined with
annotation (gnomAD/REVEL/SpliceAI), so a real annotated ClinVar release VCF is a
drop-in — see load_clinvar_records().
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from ..annotation.features import AnnotationFeatures
from ..annotation.proposer import propose
from ..annotation.vcf import VcfFieldMap, read_vcf
from ..core.variant import GenomicVariant
from ..germline.classifier import ACMGClassification, classify

# --- ClinVar CLNSIG -> our 5-tier classification --------------------------

_CLNSIG_MAP = {
    "pathogenic": ACMGClassification.PATHOGENIC,
    "pathogenic/likely_pathogenic": ACMGClassification.PATHOGENIC,
    "likely_pathogenic": ACMGClassification.LIKELY_PATHOGENIC,
    "uncertain_significance": ACMGClassification.UNCERTAIN,
    "likely_benign": ACMGClassification.LIKELY_BENIGN,
    "benign/likely_benign": ACMGClassification.BENIGN,
    "benign": ACMGClassification.BENIGN,
}

# CLNREVSTAT -> ClinVar star rating
_REVSTAT_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
}


def parse_clnsig(clnsig: Optional[str]) -> Optional[ACMGClassification]:
    """Map a ClinVar CLNSIG string to a 5-tier classification, or None to skip."""
    if not clnsig:
        return None
    key = clnsig.strip().lower().replace(" ", "_")
    return _CLNSIG_MAP.get(key)


def clnrevstat_stars(revstat: Optional[str]) -> int:
    if not revstat:
        return 0
    return _REVSTAT_STARS.get(revstat.strip().lower(), 0)


# --- concordance categories ------------------------------------------------

class Category(str, Enum):
    EXACT = "exact"                      # same 5-tier
    CLINICAL_CONCORDANT = "concordant"   # same direction, different certainty (P vs LP)
    UNDERCALL = "undercall"              # engine more conservative (engine=VUS vs ClinVar call)
    OVERCALL = "overcall"                # engine more aggressive (ClinVar=VUS vs engine call)
    CONTRADICTION = "contradiction"      # opposite direction — investigate


_P_SIDE = {ACMGClassification.PATHOGENIC, ACMGClassification.LIKELY_PATHOGENIC}
_B_SIDE = {ACMGClassification.BENIGN, ACMGClassification.LIKELY_BENIGN}


def _bucket(t: ACMGClassification) -> str:
    if t in _P_SIDE:
        return "P"
    if t in _B_SIDE:
        return "B"
    return "N"  # neutral / VUS


def categorize(clinvar: ACMGClassification, engine: ACMGClassification) -> Category:
    if clinvar is engine:
        return Category.EXACT
    cb, eb = _bucket(clinvar), _bucket(engine)
    if cb == eb:                       # same non-neutral direction, e.g. P vs LP
        return Category.CLINICAL_CONCORDANT
    if {cb, eb} == {"P", "B"}:         # opposite directions
        return Category.CONTRADICTION
    if eb == "N":                      # engine VUS, ClinVar had a directional call
        return Category.UNDERCALL
    return Category.OVERCALL           # ClinVar VUS, engine made a directional call


# --- records & evaluation --------------------------------------------------

@dataclass
class ClinVarRecord:
    variant: GenomicVariant
    features: AnnotationFeatures
    clinvar_tier: ACMGClassification
    stars: int


@dataclass
class BenchmarkRow:
    key: str
    gene: Optional[str]
    clinvar_tier: ACMGClassification
    stars: int
    engine_tier: ACMGClassification
    engine_points: int
    category: Category


class BenchmarkResult:
    def __init__(self, rows: List[BenchmarkRow]):
        self.rows = rows

    @property
    def n(self) -> int:
        return len(self.rows)

    def _count(self, cat: Category) -> int:
        return sum(1 for r in self.rows if r.category is cat)

    @property
    def strict_concordance(self) -> float:
        return self._count(Category.EXACT) / self.n if self.n else 0.0

    @property
    def clinical_concordance(self) -> float:
        good = self._count(Category.EXACT) + self._count(Category.CLINICAL_CONCORDANT)
        return good / self.n if self.n else 0.0

    @property
    def contradiction_rate(self) -> float:
        return self._count(Category.CONTRADICTION) / self.n if self.n else 0.0

    def category_counts(self) -> Dict[str, int]:
        return {c.value: self._count(c) for c in Category}

    def confusion_matrix(self) -> Dict[str, Dict[str, int]]:
        tiers = list(ACMGClassification)
        m = {c.value: {e.value: 0 for e in tiers} for c in tiers}
        for r in self.rows:
            m[r.clinvar_tier.value][r.engine_tier.value] += 1
        return m

    def contradictions(self) -> List[BenchmarkRow]:
        return [r for r in self.rows if r.category is Category.CONTRADICTION]

    def summary_text(self) -> str:
        c = self.category_counts()
        lines = [
            f"ClinVar benchmark — {self.n} variants",
            f"  strict concordance (exact 5-tier): {self.strict_concordance:6.1%}",
            f"  clinical concordance (direction):  {self.clinical_concordance:6.1%}",
            f"  contradiction rate (opposite dir): {self.contradiction_rate:6.1%}",
            "  mismatch taxonomy: "
            f"exact={c['exact']} concordant={c['concordant']} "
            f"undercall={c['undercall']} overcall={c['overcall']} "
            f"contradiction={c['contradiction']}",
        ]
        return "\n".join(lines)


def load_clinvar_records(
    path: str | Path,
    *,
    min_stars: int = 2,
    field_map: Optional[VcfFieldMap] = None,
) -> Iterator[ClinVarRecord]:
    """Read a ClinVar-format annotated VCF, yielding records above a star threshold."""
    for variant, features in read_vcf(path, field_map=field_map):
        tier = parse_clnsig(features.clinvar_sig)
        if tier is None:
            continue  # drug_response, conflicting, not_provided, etc. — skip
        stars = clnrevstat_stars(features.clinvar_review_status)
        if stars < min_stars:
            continue
        yield ClinVarRecord(variant, features, tier, stars)


def evaluate(records: List[ClinVarRecord]) -> BenchmarkResult:
    """Classify each record with the engine and score it against ClinVar."""
    rows: List[BenchmarkRow] = []
    for rec in records:
        ps = propose(rec.variant, rec.features)
        result = classify(ps.accepted_calls())
        rows.append(BenchmarkRow(
            key=rec.variant.key,
            gene=rec.variant.gene_symbol,
            clinvar_tier=rec.clinvar_tier,
            stars=rec.stars,
            engine_tier=result.classification,
            engine_points=result.points,
            category=categorize(rec.clinvar_tier, result.classification),
        ))
    return BenchmarkResult(rows)


def run_benchmark(path: str | Path, *, min_stars: int = 2) -> BenchmarkResult:
    """Convenience: load + evaluate in one call."""
    return evaluate(list(load_clinvar_records(path, min_stars=min_stars)))
