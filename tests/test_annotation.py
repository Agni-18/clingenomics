"""Tests for the annotation adapter: VCF reader + proposer, driven by sample.vcf."""

from pathlib import Path

from clingenomics.annotation.features import AnnotationFeatures, Thresholds
from clingenomics.annotation.proposer import propose
from clingenomics.annotation.vcf import VcfFieldMap, read_vcf
from clingenomics.core.evidence import Strength
from clingenomics.core.variant import MolecularConsequence
from clingenomics.germline.classifier import ACMGClassification, classify

SAMPLE = Path(__file__).resolve().parents[1] / "data" / "sample.vcf"


def _load():
    return list(read_vcf(SAMPLE))


# --- reader ----------------------------------------------------------------

def test_reader_parses_all_records():
    records = _load()
    assert len(records) == 6


def test_reader_extracts_gene_and_consequence():
    v0, _ = _load()[0]
    assert v0.gene_symbol == "BRCA2"
    assert v0.consequence is MolecularConsequence.SPLICE_DONOR


def test_reader_computes_spliceai_max():
    v0, f0 = _load()[0]
    assert abs((f0.spliceai_ds_max or 0) - 0.91) < 1e-6  # max of 0.02/0.00/0.91/0.05


def test_reader_uses_grpmax_for_frequency():
    _, f_common = _load()[2]  # MSH2 common
    assert f_common.frequency_for_acmg == 0.134


# --- proposer --------------------------------------------------------------

def test_absent_variant_gets_pm2():
    v, f = _load()[0]  # BRCA2, AC=0, no AF fields -> absent
    ps = propose(v, f)
    assert any(p.code == "PM2" for p in ps.proposals)


def test_lof_variant_gets_pvs1_proposal_with_review_flags():
    v, f = _load()[0]  # splice_donor -> PVS1 tree runs
    ps = propose(v, f)
    pvs1 = [p for p in ps.proposals if p.code == "PVS1"]
    assert len(pvs1) == 1
    assert pvs1[0].review_required is True
    # mechanism-assumed + a decision path are always surfaced
    assert any("mechanism ASSUMED" in flag for flag in ps.flags)
    assert any(flag.startswith("PVS1 path:") for flag in ps.flags)
    # from a bare VCF row (no exon context) it must NOT reach Very Strong
    assert pvs1[0].strength is not Strength.VERY_STRONG


def test_canonical_splice_pp3_suppressed_by_pvs1_double_count_guard():
    # BRCA2 splice_donor: SpliceAI is high, but PVS1 owns the splice signal,
    # so SpliceAI PP3 must NOT also be proposed (no double counting).
    v, f = _load()[0]
    ps = propose(v, f)
    assert not any(p.code == "PP3" for p in ps.proposals)
    assert any(p.code == "PVS1" for p in ps.proposals)


def test_revel_very_low_reaches_bp4_very_strong():
    from clingenomics.annotation.features import AnnotationFeatures
    from clingenomics.core.variant import Assembly, GenomicVariant, VariantType
    v = GenomicVariant(assembly=Assembly.GRCH38, chrom="1", pos=100, ref="A", alt="G",
                       variant_type=VariantType.SNV, consequence=MolecularConsequence.MISSENSE)
    ps = propose(v, AnnotationFeatures(revel=0.002, gnomad_popmax_af=0.0009))
    bp4 = next(p for p in ps.proposals if p.code == "BP4")
    assert bp4.strength is Strength.VERY_STRONG


def test_high_revel_rare_missense_proposes_pp3_and_pm2():
    v, f = _load()[1]  # BRCA1, REVEL 0.962, ultra-rare
    ps = propose(v, f)
    codes = {p.code for p in ps.proposals}
    assert "PM2" in codes
    pp3 = next(p for p in ps.proposals if p.code == "PP3")
    assert pp3.strength is Strength.STRONG


def test_common_variant_proposes_ba1():
    v, f = _load()[2]  # MSH2, grpmax 0.134
    ps = propose(v, f)
    assert any(p.code == "BA1" for p in ps.proposals)


def test_low_revel_moderate_freq_proposes_bp4_and_bs1():
    v, f = _load()[3]  # MLH1, REVEL 0.081, grpmax 0.024
    ps = propose(v, f)
    codes = {p.code for p in ps.proposals}
    assert "BP4" in codes
    assert "BS1" in codes


def test_synonymous_low_splice_proposes_bp7():
    v, f = _load()[4]  # CFTR synonymous, low SpliceAI
    ps = propose(v, f)
    assert any(p.code == "BP7" for p in ps.proposals)


# --- end to end: proposals flow into a draft classification ----------------

def test_brca1_missense_classifies_from_accepted_proposals():
    v, f = _load()[1]  # BRCA1: PM2_supporting (1) + PP3_strong (4) = 5 -> VUS draft
    ps = propose(v, f)
    result = classify(ps.accepted_calls())
    assert result.points == 5
    assert result.classification is ACMGClassification.UNCERTAIN


def test_common_variant_classifies_benign_via_ba1():
    v, f = _load()[2]
    ps = propose(v, f)
    result = classify(ps.accepted_calls())
    assert result.classification is ACMGClassification.BENIGN


def test_custom_field_map_is_respected():
    # point the reader at a non-existent AF tag; frequency should read as absent
    fm = VcfFieldMap(gnomad_popmax_af="NOT_A_TAG", gnomad_af="ALSO_MISSING")
    _, f = next(read_vcf(SAMPLE, field_map=fm))
    assert f.frequency_for_acmg is None
