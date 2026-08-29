"""Tests for the PVS1 decision tree (Abou Tayoun et al. 2018)."""

from clingenomics.core.evidence import Strength
from clingenomics.core.variant import MolecularConsequence, VariantType
from clingenomics.germline.pvs1 import (
    NMD,
    PVS1Config,
    PVS1Input,
    PVS1Strength,
    evaluate_pvs1,
)


def _ptc(**kw):
    base = dict(consequence=MolecularConsequence.STOP_GAINED, variant_type=VariantType.SNV)
    base.update(kw)
    return PVS1Input(**base)


def test_ptc_with_nmd_is_very_strong():
    r = evaluate_pvs1(_ptc(nmd=NMD.PREDICTED))
    assert r.strength is PVS1Strength.VERY_STRONG
    assert r.core_strength() is Strength.VERY_STRONG


def test_ptc_escaping_nmd_critical_region_is_strong():
    r = evaluate_pvs1(_ptc(nmd=NMD.ESCAPES, region_critical_or_gt10pct=True))
    assert r.strength is PVS1Strength.STRONG


def test_ptc_escaping_nmd_noncritical_is_moderate():
    r = evaluate_pvs1(_ptc(nmd=NMD.ESCAPES, region_critical_or_gt10pct=False))
    assert r.strength is PVS1Strength.MODERATE


def test_ptc_unknown_nmd_is_capped_and_flagged():
    r = evaluate_pvs1(_ptc(nmd=NMD.UNKNOWN))
    assert r.strength is PVS1Strength.STRONG  # default provisional cap
    assert r.strength is not PVS1Strength.VERY_STRONG
    assert any("NMD status not evaluated" in fl for fl in r.flags)


def test_nmd_derived_from_exon_context():
    # last exon -> escapes NMD -> not Very Strong
    r = evaluate_pvs1(_ptc(is_last_exon=True, region_critical_or_gt10pct=False))
    assert r.strength is PVS1Strength.MODERATE
    # middle exon -> predicted NMD -> Very Strong
    r2 = evaluate_pvs1(_ptc(is_last_exon=False, is_single_exon=False, is_in_last_50nt_penultimate=False))
    assert r2.strength is PVS1Strength.VERY_STRONG


def test_lof_not_mechanism_is_not_applicable():
    r = evaluate_pvs1(_ptc(nmd=NMD.PREDICTED, lof_is_mechanism=False, lof_mechanism_assumed=False))
    assert r.applicable is False
    assert r.strength is PVS1Strength.NA


def test_mechanism_assumption_is_flagged():
    r = evaluate_pvs1(_ptc(nmd=NMD.PREDICTED, lof_mechanism_assumed=True))
    assert any("ASSUMED" in fl for fl in r.flags)


def test_start_lost_is_moderate():
    r = evaluate_pvs1(PVS1Input(consequence=MolecularConsequence.START_LOST))
    assert r.strength is PVS1Strength.MODERATE


def test_stop_lost_is_not_applicable():
    r = evaluate_pvs1(PVS1Input(consequence=MolecularConsequence.STOP_LOST))
    assert r.applicable is False


def test_splice_frame_disrupting_nmd_is_very_strong():
    r = evaluate_pvs1(PVS1Input(
        consequence=MolecularConsequence.SPLICE_DONOR,
        preserves_reading_frame=False, nmd=NMD.PREDICTED,
    ))
    assert r.strength is PVS1Strength.VERY_STRONG


def test_splice_unknown_outcome_is_capped():
    r = evaluate_pvs1(PVS1Input(
        consequence=MolecularConsequence.SPLICE_DONOR, spliceai_ds_max=0.91,
    ))
    assert r.strength is PVS1Strength.STRONG
    assert any("Splicing outcome not evaluated" in fl for fl in r.flags)


def test_conservative_config_lowers_provisional_cap():
    cfg = PVS1Config(provisional_cap_when_nmd_unknown=PVS1Strength.MODERATE)
    r = evaluate_pvs1(_ptc(nmd=NMD.UNKNOWN), config=cfg)
    assert r.strength is PVS1Strength.MODERATE
