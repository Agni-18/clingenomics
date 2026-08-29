"""PVS1 decision tree — Abou Tayoun et al. 2018 (ClinGen SVI refinement).

A null-type consequence (nonsense, frameshift, canonical splice, start-loss) is
NOT automatically PVS1 Very Strong. The strength depends on:
  * whether loss-of-function is a disease mechanism for the gene,
  * whether the variant triggers NMD (or escapes it), and
  * for NMD-escaping / in-frame events, whether the affected region is critical.

This module encodes that tree. It NEVER returns Very Strong without an NMD
determination — when exon context is missing it caps the result at a
configurable provisional strength (default Strong) and records exactly what's
missing, so the analyst knows what to supply. Everything it emits is a draft
requiring sign-off.

Key inputs a lab supplies to strengthen the call (from transcript/exon
annotation or the RNA layer): NMD status, reading-frame outcome for splice
variants, and region criticality.
"""

from __future__ import annotations

from enum import Enum
from typing import List, Optional

from pydantic import BaseModel

from ..core.evidence import Strength
from ..core.variant import GenomicVariant, MolecularConsequence, VariantType


class NMD(str, Enum):
    PREDICTED = "predicted"   # PTC triggers nonsense-mediated decay
    ESCAPES = "escapes"       # last exon / last 50nt of penultimate exon / single-exon
    UNKNOWN = "unknown"


class PVS1Strength(str, Enum):
    NA = "not_applicable"
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"

    def to_core(self) -> Optional[Strength]:
        return {
            PVS1Strength.VERY_STRONG: Strength.VERY_STRONG,
            PVS1Strength.STRONG: Strength.STRONG,
            PVS1Strength.MODERATE: Strength.MODERATE,
            PVS1Strength.SUPPORTING: Strength.SUPPORTING,
            PVS1Strength.NA: None,
        }[self]


class PVS1Input(BaseModel):
    """Everything the tree can use. Most fields are Optional — unknowns are handled."""

    consequence: MolecularConsequence
    variant_type: VariantType = VariantType.SNV
    gene_symbol: Optional[str] = None

    lof_is_mechanism: bool = True
    lof_mechanism_assumed: bool = True     # if True we always flag the assumption
    biologically_relevant_transcript: bool = True

    # NMD: give it directly, or supply exon context to derive it
    nmd: NMD = NMD.UNKNOWN
    is_last_exon: Optional[bool] = None
    is_in_last_50nt_penultimate: Optional[bool] = None
    is_single_exon: Optional[bool] = None

    # for canonical splice variants
    preserves_reading_frame: Optional[bool] = None
    causes_exon_skip_or_cryptic: Optional[bool] = None

    # for NMD-escaping / in-frame events
    region_critical_or_gt10pct: Optional[bool] = None

    spliceai_ds_max: Optional[float] = None

    def resolved_nmd(self) -> NMD:
        if self.nmd is not NMD.UNKNOWN:
            return self.nmd
        ctx = (self.is_single_exon, self.is_last_exon, self.is_in_last_50nt_penultimate)
        if all(c is None for c in ctx):
            return NMD.UNKNOWN
        if self.is_single_exon or self.is_last_exon or self.is_in_last_50nt_penultimate:
            return NMD.ESCAPES
        return NMD.PREDICTED


class PVS1Result(BaseModel):
    strength: PVS1Strength
    applicable: bool
    path: List[str] = []
    flags: List[str] = []

    def core_strength(self) -> Optional[Strength]:
        return self.strength.to_core()


class PVS1Config(BaseModel):
    """Tunables for how conservative the tree is when information is missing."""

    provisional_cap_when_nmd_unknown: PVS1Strength = PVS1Strength.STRONG
    provisional_cap_when_criticality_unknown: PVS1Strength = PVS1Strength.MODERATE


def evaluate_pvs1(inp: PVS1Input, config: Optional[PVS1Config] = None) -> PVS1Result:
    """Run the PVS1 decision tree for one variant."""
    cfg = config or PVS1Config()
    path: List[str] = []
    flags: List[str] = []

    if inp.lof_mechanism_assumed:
        flags.append(
            f"LoF mechanism ASSUMED for {inp.gene_symbol or 'gene'} — "
            f"confirm against ClinGen/gene-specific curation"
        )
    if not inp.lof_is_mechanism:
        return PVS1Result(
            strength=PVS1Strength.NA, applicable=False,
            path=["LoF is not an established mechanism for this gene → PVS1 N/A"],
            flags=flags,
        )
    if not inp.biologically_relevant_transcript:
        return PVS1Result(
            strength=PVS1Strength.NA, applicable=False,
            path=["Variant not in a biologically-relevant transcript → PVS1 N/A"],
            flags=flags,
        )

    cons = inp.consequence

    # ---- Nonsense / frameshift branch ------------------------------------
    if cons in (MolecularConsequence.STOP_GAINED, MolecularConsequence.FRAMESHIFT):
        path.append(f"{cons.value}: premature termination codon")
        return _ptc_branch(inp, cfg, path, flags)

    # ---- Canonical splice branch (±1,2) ----------------------------------
    if cons in (MolecularConsequence.SPLICE_DONOR, MolecularConsequence.SPLICE_ACCEPTOR):
        path.append(f"{cons.value}: canonical splice site")
        return _splice_branch(inp, cfg, path, flags)

    # ---- Initiation codon branch -----------------------------------------
    if cons is MolecularConsequence.START_LOST:
        path.append("start_lost: initiation codon")
        flags.append("Check for a downstream alternative start codon / known pathogenic variants")
        return PVS1Result(
            strength=PVS1Strength.MODERATE, applicable=True,
            path=path + ["No NMD applies; PVS1 capped at Moderate for start-loss"],
            flags=flags,
        )

    # ---- Stop-loss: PVS1 not applicable; use PM4 -------------------------
    if cons is MolecularConsequence.STOP_LOST:
        return PVS1Result(
            strength=PVS1Strength.NA, applicable=False,
            path=["stop_lost → PVS1 N/A; consider PM4 (protein elongation)"],
            flags=flags,
        )

    # ---- Whole/partial-gene deletion (CNV/SV) ----------------------------
    if inp.variant_type in (VariantType.CNV, VariantType.SV):
        flags.append("CNV/SV: provide exon span and reading-frame outcome to refine PVS1")
        return PVS1Result(
            strength=cfg.provisional_cap_when_nmd_unknown, applicable=True,
            path=path + ["CNV/SV LoF: provisional pending exon/frame details"],
            flags=flags,
        )

    return PVS1Result(
        strength=PVS1Strength.NA, applicable=False,
        path=[f"{cons.value if cons else 'consequence'} not a PVS1 null type"],
        flags=flags,
    )


def _ptc_branch(inp, cfg, path, flags) -> PVS1Result:
    nmd = inp.resolved_nmd()

    if nmd is NMD.PREDICTED:
        path.append("NMD predicted → biologically-relevant transcript → PVS1 Very Strong")
        return PVS1Result(strength=PVS1Strength.VERY_STRONG, applicable=True, path=path, flags=flags)

    if nmd is NMD.ESCAPES:
        path.append("NMD escaped (last exon / last 50nt penultimate / single-exon)")
        crit = inp.region_critical_or_gt10pct
        if crit is True:
            path.append("Truncates >10% or removes a critical region → PVS1_Strong")
            return PVS1Result(strength=PVS1Strength.STRONG, applicable=True, path=path, flags=flags)
        if crit is False:
            path.append("Region not critical → PVS1_Moderate")
            return PVS1Result(strength=PVS1Strength.MODERATE, applicable=True, path=path, flags=flags)
        flags.append("Region criticality unknown — provide protein-domain / downstream-pathogenic info")
        return PVS1Result(
            strength=cfg.provisional_cap_when_criticality_unknown, applicable=True,
            path=path + ["Criticality unknown → conservative Moderate"], flags=flags,
        )

    # NMD unknown
    flags.append(
        "NMD status not evaluated — provide exon context "
        "(is_last_exon / last-50nt-penultimate / single-exon) to reach Very Strong"
    )
    return PVS1Result(
        strength=cfg.provisional_cap_when_nmd_unknown, applicable=True,
        path=path + [f"NMD unknown → capped at {cfg.provisional_cap_when_nmd_unknown.value}"],
        flags=flags,
    )


def _splice_branch(inp, cfg, path, flags) -> PVS1Result:
    frame = inp.preserves_reading_frame
    nmd = inp.resolved_nmd()

    # frame-disrupting splice that triggers NMD behaves like a PTC
    if frame is False and nmd is NMD.PREDICTED:
        path.append("Frame-disrupting + NMD predicted → PVS1 Very Strong")
        return PVS1Result(strength=PVS1Strength.VERY_STRONG, applicable=True, path=path, flags=flags)

    if frame is True or nmd is NMD.ESCAPES:
        path.append("In-frame skip or NMD-escaping splice outcome")
        crit = inp.region_critical_or_gt10pct
        if crit is True:
            return PVS1Result(strength=PVS1Strength.STRONG, applicable=True,
                              path=path + ["Critical region removed → PVS1_Strong"], flags=flags)
        if crit is False:
            return PVS1Result(strength=PVS1Strength.MODERATE, applicable=True,
                              path=path + ["Non-critical region → PVS1_Moderate"], flags=flags)
        flags.append("Region criticality unknown for in-frame splice outcome")
        return PVS1Result(strength=cfg.provisional_cap_when_criticality_unknown, applicable=True,
                          path=path + ["Criticality unknown → conservative Moderate"], flags=flags)

    # splicing outcome unknown
    flags.append(
        "Splicing outcome not evaluated — provide reading-frame / NMD result "
        "(SpliceAI frame analysis or RNA splicing evidence) to refine PVS1"
    )
    if inp.spliceai_ds_max is not None:
        path.append(f"SpliceAI Δ {inp.spliceai_ds_max:.2f} supports a splice effect")
    return PVS1Result(
        strength=cfg.provisional_cap_when_nmd_unknown, applicable=True,
        path=path + [f"Outcome unknown → capped at {cfg.provisional_cap_when_nmd_unknown.value}"],
        flags=flags,
    )


def pvs1_input_from_variant(
    variant: GenomicVariant,
    *,
    assume_lof_mechanism: bool = True,
    spliceai_ds_max: Optional[float] = None,
) -> PVS1Input:
    """Build a PVS1Input from a variant, inferring only what a VCF row supports."""
    return PVS1Input(
        consequence=variant.consequence or MolecularConsequence.OTHER,
        variant_type=variant.variant_type,
        gene_symbol=variant.gene_symbol,
        lof_is_mechanism=assume_lof_mechanism,
        lof_mechanism_assumed=assume_lof_mechanism,
        spliceai_ds_max=spliceai_ds_max if spliceai_ds_max is not None else variant.spliceai_ds_max,
    )
