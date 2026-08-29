"""Raw annotation features carried alongside a variant, plus threshold config.

`AnnotationFeatures` holds the interpretation-relevant facts that don't live on
GenomicVariant itself (population sub-fields, in-silico scores, ClinVar). These
are *facts*, not judgements — the proposer turns them into draft criteria in a
separate step so the numbers stay auditable.

`Thresholds` centralises every cut-off the proposer uses. The in-silico REVEL /
SpliceAI cut-offs are derived from the ClinGen SVI calibration work
(Pejaver et al. 2022) but MUST be verified against current guidance before
clinical use — they are here as reviewable defaults, not settled truth.
"""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field


class AnnotationFeatures(BaseModel):
    """Facts extracted from a VCF's INFO field for one variant."""

    gnomad_af: Optional[float] = Field(None, ge=0.0, le=1.0)
    gnomad_popmax_af: Optional[float] = Field(
        None, ge=0.0, le=1.0, description="Grpmax/popmax AF — the number ACMG frequency codes use"
    )
    gnomad_ac: Optional[int] = Field(None, ge=0)

    revel: Optional[float] = Field(None, ge=0.0, le=1.0)
    cadd_phred: Optional[float] = None
    spliceai_ds_max: Optional[float] = Field(None, ge=0.0, le=1.0)

    clinvar_sig: Optional[str] = None
    clinvar_review_status: Optional[str] = Field(
        None, description="ClinVar CLNREVSTAT — used to derive a star rating for benchmarking"
    )

    @property
    def frequency_for_acmg(self) -> Optional[float]:
        """Popmax if available, else global AF — the value frequency codes act on."""
        if self.gnomad_popmax_af is not None:
            return self.gnomad_popmax_af
        return self.gnomad_af


class Thresholds(BaseModel):
    """All proposer cut-offs in one place, so a lab can tune and version them.

    In-silico values are the ClinGen-calibrated thresholds:
      * REVEL  — Pejaver et al. 2022, Am J Hum Genet (PMID 36413997)
      * SpliceAI — Walker et al. 2023, ClinGen SVI Splicing Subgroup (PMID 37352859)
    Gene-specific VCEPs sometimes override these (e.g. simplified ≥0.7/≤0.3, or
    disease-specific frequency cut-offs) — tune per panel where a VCEP exists.
    """

    # frequency (ACMG defaults; gene/disease-specific overrides exist)
    ba1_af: float = 0.05          # stand-alone benign
    bs1_af: float = 0.01          # above disease-max-credible — needs disease-specific tuning
    pm2_af: float = 0.0001        # absent / ultra-rare -> PM2_supporting

    # REVEL (missense, protein-impact) — Pejaver 2022, calibrated to ACMG strengths
    revel_pp3_strong: float = 0.932
    revel_pp3_moderate: float = 0.773
    revel_pp3_supporting: float = 0.644
    revel_bp4_supporting: float = 0.290
    revel_bp4_moderate: float = 0.183
    revel_bp4_strong: float = 0.016
    revel_bp4_very_strong: float = 0.003   # REVEL is the tool that reaches BP4_VeryStrong

    # SpliceAI (splice-impact) — Walker 2023 ClinGen SVI Splicing Subgroup
    splice_pp3_supporting: float = 0.20    # max Δ ≥ 0.20 -> PP3 supporting (spliceogenic)
    splice_bp4_moderate: float = 0.10      # max Δ ≤ 0.10 -> BP4 moderate (non-spliceogenic)
    splice_bp7_max: float = 0.10           # BP7 needs max Δ < 0.10 (and not in a splice motif)
