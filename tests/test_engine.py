"""Tests for the interpretation engine: points classifier, RNA mapping, somatic tiers."""

import pytest

from clingenomics.core.evidence import Direction, Strength
from clingenomics.germline.acmg_criteria import make_call
from clingenomics.germline.classifier import ACMGClassification, classify
from clingenomics.rna.evidence import RNAEvidence, RNAObservation
from clingenomics.somatic.tiers import (
    ClinicalAssertion,
    EvidenceLevel,
    OncoKBLevel,
    SomaticTier,
    assertion_from_oncokb,
    classify_somatic,
)


# --- points-based ACMG classifier -----------------------------------------

def test_pvs1_alone_is_likely_pathogenic():
    r = classify([make_call("PVS1")])
    assert r.points == 8
    assert r.classification is ACMGClassification.LIKELY_PATHOGENIC


def test_pvs1_plus_pm2_reaches_pathogenic():
    r = classify([make_call("PVS1"), make_call("PM2")])
    assert r.points == 10
    assert r.classification is ACMGClassification.PATHOGENIC


def test_modulated_strength_changes_points():
    # PVS1 applied at Strong (splicing decision tree) + PP3 supporting = 4 + 1
    r = classify([make_call("PVS1", strength=Strength.STRONG), make_call("PP3")])
    assert r.points == 5
    assert r.classification is ACMGClassification.UNCERTAIN


def test_balanced_evidence_is_vus():
    r = classify([make_call("PM2", strength=Strength.SUPPORTING), make_call("BP4")])
    assert r.points == 0
    assert r.classification is ACMGClassification.UNCERTAIN


def test_two_benign_strong_is_benign():
    r = classify([make_call("BS1"), make_call("BS2")])
    assert r.points == -8
    assert r.classification is ACMGClassification.BENIGN


def test_ba1_forces_benign_regardless_of_points():
    # even with strong pathogenic evidence present, BA1 override wins
    r = classify([make_call("BA1"), make_call("PVS1")])
    assert r.classification is ACMGClassification.BENIGN
    assert any("BA1" in f for f in r.flags)


def test_conflict_is_flagged():
    r = classify([make_call("PVS1"), make_call("BS3")])
    assert any("CONFLICTING" in f for f in r.flags)


def test_duplicate_code_keeps_strongest():
    r = classify([make_call("PS3", strength=Strength.SUPPORTING), make_call("PS3", strength=Strength.STRONG)])
    assert r.points == 4


def test_standalone_pathogenic_is_rejected():
    with pytest.raises(ValueError):
        make_call("PVS1", strength=Strength.STANDALONE)


# --- RNA-as-evidence mapping ----------------------------------------------

def test_aberrant_splicing_maps_to_ps3_strong_when_quantified():
    ev = RNAEvidence(
        observation=RNAObservation.ABERRANT_SPLICING,
        tissue="fibroblast", well_controlled=True, quantified=True,
    )
    call = ev.to_criterion_call(variant_predicts_splicing=True, gene_lof_mechanism=True)
    assert call.code == "PS3"
    assert call.direction is Direction.PATHOGENIC
    assert call.strength is Strength.STRONG


def test_normal_splicing_needs_a_splicing_prediction_to_apply_bs3():
    ev = RNAEvidence(observation=RNAObservation.NORMAL_SPLICING, tissue="blood")
    assert ev.to_criterion_call(variant_predicts_splicing=False) is None
    call = ev.to_criterion_call(variant_predicts_splicing=True)
    assert call.code == "BS3"


def test_underexpression_reaches_strong_only_in_lof_gene():
    ev = RNAEvidence(observation=RNAObservation.EXPRESSION_UNDER, tissue="fibroblast", well_controlled=True)
    assert ev.to_criterion_call(gene_lof_mechanism=True).strength is Strength.STRONG
    assert ev.to_criterion_call(gene_lof_mechanism=False).strength is Strength.MODERATE


def test_rna_evidence_flows_into_final_classification():
    # a predicted splice variant, absent from gnomAD (PM2), with RNA confirmation
    rna = RNAEvidence(
        observation=RNAObservation.ABERRANT_SPLICING,
        tissue="fibroblast", well_controlled=True, quantified=True,
    )
    rna_call = rna.to_criterion_call(variant_predicts_splicing=True, gene_lof_mechanism=True)
    result = classify([
        make_call("PVS1", strength=Strength.STRONG),  # splicing branch, RNA-supported
        make_call("PM2"),
        rna_call,                                       # PS3 Strong
    ])
    # 4 (PVS1_Strong) + 2 (PM2) + 4 (PS3) = 10
    assert result.points == 10
    assert result.classification is ACMGClassification.PATHOGENIC


# --- somatic tiers ---------------------------------------------------------

def test_level_a_assertion_is_tier_i():
    a = ClinicalAssertion(description="EGFR L858R -> osimertinib", evidence_level=EvidenceLevel.A)
    assert classify_somatic([a]).tier is SomaticTier.TIER_I


def test_only_level_d_is_tier_ii():
    a = ClinicalAssertion(description="preclinical", evidence_level=EvidenceLevel.D)
    assert classify_somatic([a]).tier is SomaticTier.TIER_II


def test_no_assertions_is_tier_iii():
    assert classify_somatic([]).tier is SomaticTier.TIER_III


def test_likely_benign_short_circuits_to_tier_iv():
    a = ClinicalAssertion(description="x", evidence_level=EvidenceLevel.A)
    assert classify_somatic([a], is_likely_benign=True).tier is SomaticTier.TIER_IV


def test_oncokb_level_1_maps_to_tier_i():
    a = assertion_from_oncokb(OncoKBLevel.LEVEL_1, therapy="drug", tumor_type="NSCLC")
    assert a.evidence_level is EvidenceLevel.A
    assert classify_somatic([a]).tier is SomaticTier.TIER_I


def test_oncokb_resistance_flag():
    a = assertion_from_oncokb(OncoKBLevel.R1)
    assert a.resistance is True
