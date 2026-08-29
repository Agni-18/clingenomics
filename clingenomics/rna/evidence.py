"""RNA-as-evidence: turn RNA-seq observations into ACMG criterion calls.

This is the integration point that makes the platform "multi-modal". RNA data
enters the germline framework mostly through PS3 / BS3 (functional) and through
the PVS1 splicing decision tree, applied at *modulated* strengths per the
ClinGen SVI RNA recommendations (Walker et al. 2023).

The strength logic below is intentionally explicit and conservative; it is the
part a lab will most want to tune / validate against its own controls, so it is
kept in one place rather than scattered through the caller adapters (FRASER,
OUTRIDER, ASE callers, arriba/STAR-Fusion).
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field

from ..core.evidence import CriterionCall, EvidenceSource, Strength
from ..germline.acmg_criteria import make_call


class RNAObservation(str, Enum):
    ABERRANT_SPLICING = "aberrant_splicing"          # FRASER outlier consistent w/ prediction
    NORMAL_SPLICING = "normal_splicing"              # RNA shows only canonical transcript
    EXPRESSION_UNDER = "expression_underexpression"  # OUTRIDER under-expression (NMD-consistent)
    ASE_MONOALLELIC = "ase_monoallelic"              # allele-specific: variant allele lost
    ASE_BIALLELIC = "ase_biallelic"                  # both alleles expressed (argues vs NMD/LoF)
    FUSION = "fusion"                                # arriba / STAR-Fusion gene fusion


class RNAEvidence(BaseModel):
    """A single RNA-seq observation attached to a variant."""

    observation: RNAObservation
    tissue: str = Field(..., description="Source tissue; affects interpretability")
    well_controlled: bool = Field(
        True, description="Assessed against adequate tissue-matched controls / replicates"
    )
    quantified: bool = Field(
        False, description="Effect is quantified (e.g. % aberrant transcript), not just present"
    )
    # supporting numbers, optional
    padj: Optional[float] = Field(None, ge=0.0, le=1.0, description="Outlier adjusted p-value")
    delta_psi: Optional[float] = Field(None, description="FRASER delta-PSI for splicing outlier")
    evidence_id: Optional[str] = None

    # --- gating context supplied by the variant / annotation ---
    def to_criterion_call(
        self,
        *,
        variant_predicts_splicing: bool = False,
        gene_lof_mechanism: bool = False,
    ) -> Optional[CriterionCall]:
        """Emit the ACMG criterion call this RNA observation supports, or None.

        Args:
            variant_predicts_splicing: variant is a splice-site or SpliceAI-flagged
                candidate (needed to let splicing RNA feed the PVS1 tree vs PS3).
            gene_lof_mechanism: LoF is an established disease mechanism for the gene
                (gates whether NMD-consistent evidence can reach PVS1 strength).
        """
        obs = self.observation

        if obs is RNAObservation.ABERRANT_SPLICING:
            # A confirmed aberrant-splicing event supports a damaging functional
            # effect (PS3). If the variant was already a predicted splice
            # candidate in a LoF gene, RNA confirmation can instead uplift the
            # PVS1 branch; here we express it as PS3 at a strength scaled by
            # control quality and quantitation, and flag the PVS1 linkage.
            strength = Strength.STRONG if (self.well_controlled and self.quantified) else Strength.MODERATE
            rationale = "RNA-seq confirms predicted aberrant splicing"
            if variant_predicts_splicing and gene_lof_mechanism:
                rationale += " (supports PVS1 splicing branch; recorded as PS3 functional)"
            return make_call(
                "PS3", strength=strength, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale=rationale, evidence_id=self.evidence_id,
            )

        if obs is RNAObservation.NORMAL_SPLICING:
            # For a variant predicted to disrupt splicing, RNA showing only the
            # normal transcript is functional evidence *against* an effect (BS3).
            if not variant_predicts_splicing:
                return None
            strength = Strength.STRONG if self.well_controlled else Strength.SUPPORTING
            return make_call(
                "BS3", strength=strength, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale="RNA-seq shows only canonical splicing despite prediction",
                evidence_id=self.evidence_id,
            )

        if obs is RNAObservation.EXPRESSION_UNDER:
            # Significant under-expression consistent with NMD supports a
            # loss-of-function effect (PS3). Reaches Strong only when the gene's
            # mechanism is LoF and the outlier is well controlled.
            if self.well_controlled and gene_lof_mechanism:
                strength = Strength.STRONG
            elif self.well_controlled:
                strength = Strength.MODERATE
            else:
                strength = Strength.SUPPORTING
            return make_call(
                "PS3", strength=strength, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale="Expression outlier: under-expression consistent with NMD",
                evidence_id=self.evidence_id,
            )

        if obs is RNAObservation.ASE_MONOALLELIC:
            # Loss of the variant-bearing allele's transcript (monoallelic
            # expression of the reference) is consistent with NMD → PS3.
            strength = Strength.MODERATE if self.well_controlled else Strength.SUPPORTING
            return make_call(
                "PS3", strength=strength, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale="Allele-specific expression consistent with NMD of variant allele",
                evidence_id=self.evidence_id,
            )

        if obs is RNAObservation.ASE_BIALLELIC:
            # Both alleles expressed argues against NMD-mediated LoF for a
            # variant predicted to trigger it → mild benign functional support.
            return make_call(
                "BS3", strength=Strength.SUPPORTING, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale="Biallelic expression argues against predicted NMD",
                evidence_id=self.evidence_id,
            )

        if obs is RNAObservation.FUSION:
            # Fusions are primarily a somatic/structural finding; in the germline
            # framework a validated activating/disrupting fusion is functional
            # evidence. Left as PS3-supporting pending gene-specific curation.
            return make_call(
                "PS3", strength=Strength.SUPPORTING, source=EvidenceSource.FUNCTIONAL_RNA,
                rationale="RNA fusion detected (requires gene-specific curation)",
                evidence_id=self.evidence_id,
            )

        return None
