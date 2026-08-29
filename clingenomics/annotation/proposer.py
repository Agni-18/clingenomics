"""Turn annotation features into *draft* ACMG criteria for human review.

Nothing here is auto-applied to a classification. `propose()` returns a
`ProposalSet` of criteria the analyst can accept/reject, plus `flags` for things
that need a human decision and must NOT be mechanically accepted (PVS1).

Covered, defensibly automatable:
  PM2 / BS1 / BA1   from population frequency
  PP3 / BP4         from REVEL (missense) and SpliceAI (splicing)
  BP7               synonymous with no predicted splice impact

Deliberately NOT auto-applied:
  PVS1              emitted as a flag ("run the decision tree") — mechanism,
                    NMD escape and biologically-relevant transcript are human calls
  PS1/PM5/PM3/PP1…  need case-level or curated data outside a single VCF row
"""

from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel

from ..core.evidence import CriterionCall, Direction, EvidenceSource, Strength
from ..core.variant import GenomicVariant, MolecularConsequence
from ..germline.acmg_criteria import make_call
from ..germline.pvs1 import evaluate_pvs1, pvs1_input_from_variant
from .features import AnnotationFeatures, Thresholds


class Proposal(BaseModel):
    """A draft criterion awaiting sign-off, with the numbers that justify it."""

    code: str
    direction: Direction
    strength: Strength
    rationale: str
    evidence: Dict[str, float] = {}
    review_required: bool = True

    def to_call(self) -> CriterionCall:
        """Convert an accepted proposal into a CriterionCall for the classifier."""
        return make_call(
            self.code,
            strength=self.strength,
            rationale=self.rationale,
            source=EvidenceSource.POPULATION
            if self.code in {"PM2", "BS1", "BA1"}
            else EvidenceSource.COMPUTATIONAL,
        )


class ProposalSet(BaseModel):
    proposals: List[Proposal] = []
    flags: List[str] = []

    def accepted_calls(self) -> List[CriterionCall]:
        """Convenience for demos/tests: treat every proposal as accepted."""
        return [p.to_call() for p in self.proposals]


def _frequency_proposals(feat: AnnotationFeatures, th: Thresholds) -> List[Proposal]:
    af = feat.frequency_for_acmg
    if af is None:
        # absent from gnomAD entirely — PM2 supporting
        return [Proposal(
            code="PM2", direction=Direction.PATHOGENIC, strength=Strength.SUPPORTING,
            rationale="Absent from gnomAD", evidence={"popmax_af": 0.0},
        )]
    if af >= th.ba1_af:
        return [Proposal(
            code="BA1", direction=Direction.BENIGN, strength=Strength.STANDALONE,
            rationale=f"Popmax AF {af:.4f} ≥ {th.ba1_af} (stand-alone benign)",
            evidence={"popmax_af": af},
        )]
    if af >= th.bs1_af:
        return [Proposal(
            code="BS1", direction=Direction.BENIGN, strength=Strength.STRONG,
            rationale=f"Popmax AF {af:.4f} ≥ {th.bs1_af} — set a disease-specific threshold",
            evidence={"popmax_af": af},
        )]
    if af <= th.pm2_af:
        return [Proposal(
            code="PM2", direction=Direction.PATHOGENIC, strength=Strength.SUPPORTING,
            rationale=f"Popmax AF {af:.6f} ≤ {th.pm2_af} (ultra-rare)",
            evidence={"popmax_af": af},
        )]
    return []


def _revel_pp3(r: float, th: Thresholds) -> Optional[Strength]:
    if r >= th.revel_pp3_strong:
        return Strength.STRONG
    if r >= th.revel_pp3_moderate:
        return Strength.MODERATE
    if r >= th.revel_pp3_supporting:
        return Strength.SUPPORTING
    return None


def _revel_bp4(r: float, th: Thresholds) -> Optional[Strength]:
    if r <= th.revel_bp4_very_strong:
        return Strength.VERY_STRONG
    if r <= th.revel_bp4_strong:
        return Strength.STRONG
    if r <= th.revel_bp4_moderate:
        return Strength.MODERATE
    if r <= th.revel_bp4_supporting:
        return Strength.SUPPORTING
    return None


def _insilico_proposals(
    variant: GenomicVariant, feat: AnnotationFeatures, th: Thresholds
) -> List[Proposal]:
    """Combined protein (REVEL) + splice (SpliceAI) in-silico evidence.

    Implements the ClinGen combining rules:
      * PP3 is counted ONCE at the strongest applicable strength across protein
        and splice predictors (never both) — Pejaver/Walker double-count guard.
      * For a canonical splice consequence the splice signal is owned by the
        PVS1 pathway, so SpliceAI PP3 is suppressed here.
      * BP4 (benign) requires NO pathogenic signal from either axis.
    """
    props: List[Proposal] = []
    cons = variant.consequence
    ds = feat.spliceai_ds_max
    revel = feat.revel
    is_canonical_splice = cons in (
        MolecularConsequence.SPLICE_DONOR, MolecularConsequence.SPLICE_ACCEPTOR
    )

    # --- pathogenic signals ---
    protein_pp3 = _revel_pp3(revel, th) if (cons is MolecularConsequence.MISSENSE and revel is not None) else None
    splice_pathogenic = (
        ds is not None and ds >= th.splice_pp3_supporting and not is_canonical_splice
    )

    # single PP3 at the strongest strength (protein vs splice-supporting)
    candidates = []
    if protein_pp3 is not None:
        candidates.append((protein_pp3, f"REVEL {revel:.3f}"))
    if splice_pathogenic:
        candidates.append((Strength.SUPPORTING, f"SpliceAI Δ {ds:.2f}"))
    if candidates:
        strength, why = max(candidates, key=lambda c: c[0].points)
        props.append(Proposal(
            code="PP3", direction=Direction.PATHOGENIC, strength=strength,
            rationale=f"{why} — deleterious in-silico support (single PP3, max strength)",
            evidence={"revel": revel or 0.0, "spliceai_ds_max": ds or 0.0},
        ))
        return props  # a pathogenic signal precludes BP4

    # --- benign signals (only when no pathogenic in-silico signal) ---
    # protein-benign for missense, blocked if any positive splice signal exists
    if cons is MolecularConsequence.MISSENSE and revel is not None:
        no_splice_signal = ds is None or ds < th.splice_pp3_supporting
        bp4 = _revel_bp4(revel, th) if no_splice_signal else None
        if bp4 is not None:
            props.append(Proposal(
                code="BP4", direction=Direction.BENIGN, strength=bp4,
                rationale=f"REVEL {revel:.3f} — benign in-silico support",
                evidence={"revel": revel},
            ))
    # splice-benign for variants where splicing is the open question (intronic)
    elif cons is MolecularConsequence.INTRONIC and ds is not None and ds <= th.splice_bp4_moderate:
        props.append(Proposal(
            code="BP4", direction=Direction.BENIGN, strength=Strength.MODERATE,
            rationale=f"SpliceAI Δ {ds:.2f} ≤ {th.splice_bp4_moderate} — no predicted splice impact",
            evidence={"spliceai_ds_max": ds},
        ))
    return props


def propose(
    variant: GenomicVariant,
    features: AnnotationFeatures,
    thresholds: Optional[Thresholds] = None,
) -> ProposalSet:
    """Produce draft criteria + flags for one annotated variant."""
    th = thresholds or Thresholds()
    result = ProposalSet()

    result.proposals.extend(_frequency_proposals(features, th))

    cons = variant.consequence
    ds = features.spliceai_ds_max

    # combined protein (REVEL) + splice (SpliceAI) in-silico, with double-count guard
    result.proposals.extend(_insilico_proposals(variant, features, th))

    # BP7: synonymous with no predicted splice impact
    if cons is MolecularConsequence.SYNONYMOUS and (ds is None or ds < th.splice_bp7_max):
        result.proposals.append(Proposal(
            code="BP7", direction=Direction.BENIGN, strength=Strength.SUPPORTING,
            rationale="Synonymous, no predicted splice impact",
            evidence={"spliceai_ds_max": ds or 0.0},
        ))

    # PVS1 — run the decision tree. Still a review-required draft, but now
    # strength-calibrated rather than a bare flag. Mechanism is assumed here
    # (per project setting) and that assumption is always surfaced as a flag.
    if cons is not None and cons.is_predicted_lof:
        pvs1 = evaluate_pvs1(
            pvs1_input_from_variant(variant, assume_lof_mechanism=True, spliceai_ds_max=ds)
        )
        result.flags.extend(pvs1.flags)
        result.flags.append("PVS1 path: " + " → ".join(pvs1.path))
        core = pvs1.core_strength()
        if pvs1.applicable and core is not None:
            result.proposals.append(Proposal(
                code="PVS1", direction=Direction.PATHOGENIC, strength=core,
                rationale=f"PVS1 [{core.label}] via decision tree — REVIEW (mechanism assumed)",
                evidence={"spliceai_ds_max": ds or 0.0},
            ))

    return result
