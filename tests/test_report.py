"""Tests for the HTML report renderer."""

from pathlib import Path

from clingenomics.annotation.vcf import read_vcf
from clingenomics.report.html_report import build_report_for, render_variant_report, write_report
from clingenomics.annotation.proposer import propose
from clingenomics.germline.classifier import classify

SAMPLE = Path(__file__).resolve().parents[1] / "data" / "sample.vcf"


def _first():
    return next(read_vcf(SAMPLE))


def test_report_is_standalone_html():
    v, f = _first()
    html = build_report_for(v, f)
    assert html.startswith("<!doctype html>")
    assert "<style>" in html  # CSS embedded, no external deps
    assert "http://" not in html and "https://" not in html  # fully self-contained


def test_report_shows_classification_and_points():
    v, f = _first()  # BRCA2 -> VUS +5
    ps = propose(v, f)
    result = classify(ps.accepted_calls())
    html = render_variant_report(v, f, ps, result)
    assert result.classification.value in html
    assert f"{result.points:+d}" in html


def test_report_lists_every_proposed_criterion():
    v, f = _first()
    ps = propose(v, f)
    result = classify(ps.accepted_calls())
    html = render_variant_report(v, f, ps, result)
    for p in ps.proposals:
        assert p.code in html


def test_report_marks_draft_and_flags_pvs1_assumption():
    v, f = _first()  # BRCA2 has the PVS1 mechanism-assumed flag
    html = build_report_for(v, f)
    assert "Draft" in html or "DRAFT" in html
    assert "mechanism" in html.lower()


def test_write_report_creates_file(tmp_path):
    v, f = _first()
    html = build_report_for(v, f)
    out = write_report(tmp_path / "r.html", html)
    assert out.exists()
    assert out.read_text(encoding="utf-8").startswith("<!doctype html>")
