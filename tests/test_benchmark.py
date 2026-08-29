"""Tests for the ClinVar benchmark harness (mechanics, not a validation claim)."""

from pathlib import Path

from clingenomics.benchmark.evaluate import (
    Category,
    categorize,
    clnrevstat_stars,
    load_clinvar_records,
    parse_clnsig,
    run_benchmark,
)
from clingenomics.germline.classifier import ACMGClassification as C

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "clinvar_benchmark_synthetic.vcf"


# --- CLNSIG / review-status parsing ---------------------------------------

def test_parse_clnsig_variants():
    assert parse_clnsig("Pathogenic") is C.PATHOGENIC
    assert parse_clnsig("Likely_benign") is C.LIKELY_BENIGN
    assert parse_clnsig("Pathogenic/Likely_pathogenic") is C.PATHOGENIC
    assert parse_clnsig("drug_response") is None
    assert parse_clnsig(None) is None


def test_review_status_stars():
    assert clnrevstat_stars("reviewed_by_expert_panel") == 3
    assert clnrevstat_stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert clnrevstat_stars("no_assertion_provided") == 0
    assert clnrevstat_stars(None) == 0


# --- categorisation logic --------------------------------------------------

def test_exact_match():
    assert categorize(C.PATHOGENIC, C.PATHOGENIC) is Category.EXACT


def test_same_direction_is_clinical_concordant():
    assert categorize(C.PATHOGENIC, C.LIKELY_PATHOGENIC) is Category.CLINICAL_CONCORDANT
    assert categorize(C.BENIGN, C.LIKELY_BENIGN) is Category.CLINICAL_CONCORDANT


def test_engine_vus_vs_directional_is_undercall():
    assert categorize(C.PATHOGENIC, C.UNCERTAIN) is Category.UNDERCALL
    assert categorize(C.BENIGN, C.UNCERTAIN) is Category.UNDERCALL


def test_clinvar_vus_vs_directional_is_overcall():
    assert categorize(C.UNCERTAIN, C.LIKELY_PATHOGENIC) is Category.OVERCALL


def test_opposite_direction_is_contradiction():
    assert categorize(C.PATHOGENIC, C.BENIGN) is Category.CONTRADICTION
    assert categorize(C.LIKELY_BENIGN, C.LIKELY_PATHOGENIC) is Category.CONTRADICTION


# --- end-to-end over the fixture ------------------------------------------

def test_fixture_loads_and_filters_by_stars():
    recs = list(load_clinvar_records(FIXTURE, min_stars=2))
    assert len(recs) > 0
    assert all(r.stars >= 2 for r in recs)


def test_benchmark_runs_and_metrics_are_bounded():
    result = run_benchmark(FIXTURE, min_stars=2)
    assert result.n > 0
    for v in (result.strict_concordance, result.clinical_concordance, result.contradiction_rate):
        assert 0.0 <= v <= 1.0
    # counts sum to n
    assert sum(result.category_counts().values()) == result.n


def test_no_contradictions_on_curated_fixture():
    # the synthetic set is engineered clean; a contradiction here would be a real bug
    result = run_benchmark(FIXTURE, min_stars=2)
    assert result.contradiction_rate == 0.0


def test_confusion_matrix_totals_match():
    result = run_benchmark(FIXTURE, min_stars=2)
    cm = result.confusion_matrix()
    total = sum(cm[c][e] for c in cm for e in cm[c])
    assert total == result.n
