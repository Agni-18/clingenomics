"""Standalone HTML clinical report renderer — zero dependencies.

Renders one germline variant interpretation as a self-contained, print-friendly
HTML page: variant identity, the draft classification, a transparent points
ledger (the arithmetic, shown not hidden), the annotation facts, the PVS1
decision path, and the review flags that require human sign-off.

Everything produced is a DRAFT — proposals are review-required and PVS1 assumes
LoF mechanism — so the report says so, loudly and in print.
"""

from __future__ import annotations

import html
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

from ..annotation.features import AnnotationFeatures, Thresholds
from ..annotation.proposer import ProposalSet, propose
from ..core.variant import GenomicVariant
from ..germline.classifier import ACMGClassification, ClassificationResult, classify

# tier -> (background, foreground) for the classification banner
_TIER_COLORS = {
    ACMGClassification.PATHOGENIC: ("#8f1d14", "#ffffff"),
    ACMGClassification.LIKELY_PATHOGENIC: ("#b4530b", "#ffffff"),
    ACMGClassification.UNCERTAIN: ("#6b5d2a", "#ffffff"),
    ACMGClassification.LIKELY_BENIGN: ("#2f6d44", "#ffffff"),
    ACMGClassification.BENIGN: ("#1d5a37", "#ffffff"),
}

_CSS = """
* { box-sizing: border-box; }
body {
  font-family: "Segoe UI", -apple-system, Helvetica, Arial, sans-serif;
  color: #1c2530; background: #eef1f4; margin: 0; padding: 32px;
  line-height: 1.5;
}
.report {
  max-width: 820px; margin: 0 auto; background: #ffffff;
  border: 1px solid #d4dae0; border-radius: 10px; overflow: hidden;
  box-shadow: 0 2px 10px rgba(20,40,60,.06);
}
.draft-strip {
  background: repeating-linear-gradient(45deg,#3a4a5a,#3a4a5a 12px,#33424f 12px,#33424f 24px);
  color: #fff; text-align: center; font-size: 12px; letter-spacing: .12em;
  text-transform: uppercase; padding: 7px; font-weight: 600;
}
.head { padding: 26px 30px 18px; border-bottom: 1px solid #e4e9ee; }
.eyebrow { font-size: 11px; letter-spacing: .16em; text-transform: uppercase; color: #6a7684; margin: 0 0 6px; }
.gene { font-size: 28px; font-weight: 700; margin: 0; letter-spacing: -.01em; }
.hgvs { font-family: "Cascadia Code", "Consolas", monospace; font-size: 15px; color: #33424f; margin: 4px 0 0; }
.meta { margin: 14px 0 0; font-size: 13px; color: #55636f; display: flex; flex-wrap: wrap; gap: 4px 22px; }
.meta span b { color: #1c2530; font-weight: 600; }
.banner { display: flex; align-items: center; justify-content: space-between; padding: 16px 30px; }
.banner .label { font-size: 20px; font-weight: 700; }
.banner .pts { font-size: 15px; font-weight: 600; font-variant-numeric: tabular-nums; opacity: .95; }
.section { padding: 22px 30px; border-top: 1px solid #eef1f4; }
.section h2 { font-size: 12px; letter-spacing: .14em; text-transform: uppercase; color: #6a7684; margin: 0 0 12px; }
table { width: 100%; border-collapse: collapse; font-size: 13.5px; }
th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid #eef1f4; vertical-align: top; }
th { font-size: 11px; letter-spacing: .06em; text-transform: uppercase; color: #7a8794; font-weight: 600; }
td.code { font-weight: 700; font-family: "Cascadia Code","Consolas",monospace; white-space: nowrap; }
td.pts { text-align: right; font-variant-numeric: tabular-nums; font-weight: 700; white-space: nowrap; }
.pts-pos { color: #8f1d14; } .pts-neg { color: #1d5a37; }
tr.total td { border-top: 2px solid #d4dae0; border-bottom: none; font-weight: 700; padding-top: 12px; }
.facts { display: grid; grid-template-columns: repeat(2,1fr); gap: 8px 26px; font-size: 13.5px; }
.facts div { display: flex; justify-content: space-between; border-bottom: 1px dotted #e4e9ee; padding: 5px 0; }
.facts .k { color: #55636f; } .facts .v { font-weight: 600; font-variant-numeric: tabular-nums; }
.path { font-size: 13px; color: #33424f; background: #f6f8fa; border: 1px solid #e4e9ee;
  border-radius: 6px; padding: 12px 14px; line-height: 1.7; }
.flags { list-style: none; margin: 0; padding: 0; }
.flags li { font-size: 13px; padding: 10px 12px 10px 34px; position: relative;
  background: #fff8ec; border: 1px solid #f0e0bf; border-radius: 6px; margin-bottom: 8px; }
.flags li::before { content: "!"; position: absolute; left: 12px; top: 9px; font-weight: 800;
  color: #b4530b; background: #f6e4c2; width: 16px; height: 16px; border-radius: 50%;
  text-align: center; line-height: 16px; font-size: 11px; }
.foot { padding: 18px 30px 24px; font-size: 11.5px; color: #7a8794; border-top: 1px solid #eef1f4; }
.foot code { font-family: "Cascadia Code","Consolas",monospace; }
@media print {
  body { background: #fff; padding: 0; }
  .report { border: none; box-shadow: none; max-width: none; }
}
"""


def _esc(x) -> str:
    return html.escape(str(x)) if x is not None else "—"


def _fmt_af(v: Optional[float]) -> str:
    if v is None:
        return "—"
    if v == 0:
        return "0 (absent)"
    return f"{v:.6f}".rstrip("0").rstrip(".") if v < 0.001 else f"{v:.4f}"


def render_variant_report(
    variant: GenomicVariant,
    features: AnnotationFeatures,
    proposals: ProposalSet,
    classification: ClassificationResult,
) -> str:
    """Return a complete standalone HTML document for one variant."""
    bg, fg = _TIER_COLORS[classification.classification]

    # points ledger rows from the proposals (each is a review-required draft)
    rows = []
    for p in sorted(proposals.proposals, key=lambda x: -x.to_call().points):
        pts = p.to_call().points
        cls = "pts-pos" if pts > 0 else "pts-neg"
        rows.append(
            f"<tr><td class='code'>{_esc(p.code)}</td>"
            f"<td>{_esc(p.strength.label)}</td>"
            f"<td>{_esc(p.rationale)}</td>"
            f"<td class='pts {cls}'>{pts:+d}</td></tr>"
        )
    ledger = "".join(rows) or "<tr><td colspan='4'>No criteria proposed.</td></tr>"

    # separate the PVS1 decision path from the other review flags
    pvs1_path = next((f for f in proposals.flags if f.startswith("PVS1 path:")), None)
    other_flags = [f for f in proposals.flags if not f.startswith("PVS1 path:")]
    flags_html = "".join(f"<li>{_esc(f)}</li>" for f in other_flags)

    path_section = ""
    if pvs1_path:
        crumbs = pvs1_path.replace("PVS1 path:", "").strip()
        path_section = (
            "<div class='section'><h2>PVS1 decision path</h2>"
            f"<div class='path'>{_esc(crumbs)}</div></div>"
        )

    flags_section = ""
    if flags_html:
        flags_section = (
            "<div class='section'><h2>Review required before sign-out</h2>"
            f"<ul class='flags'>{flags_html}</ul></div>"
        )

    conflict = ""
    if classification.flags:
        conflict = "".join(f"<li>{_esc(f)}</li>" for f in classification.flags)
        conflict = f"<div class='section'><h2>Classifier notes</h2><ul class='flags'>{conflict}</ul></div>"

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(variant.gene_symbol or variant.key)} — variant report (DRAFT)</title>
<style>{_CSS}</style></head>
<body><div class="report">
<div class="draft-strip">Draft · research use · not for clinical diagnosis</div>
<div class="head">
  <p class="eyebrow">Germline variant interpretation</p>
  <h1 class="gene">{_esc(variant.gene_symbol or "Unknown gene")}</h1>
  <p class="hgvs">{_esc(variant.hgvs_c or variant.key)}{(" · " + _esc(variant.hgvs_p)) if variant.hgvs_p else ""}</p>
  <div class="meta">
    <span><b>Genomic:</b> {_esc(variant.key)}</span>
    <span><b>Assembly:</b> {_esc(variant.assembly.value)}</span>
    <span><b>Consequence:</b> {_esc(variant.consequence.value if variant.consequence else "—")}</span>
    <span><b>Transcript:</b> {_esc(variant.transcript or "—")}</span>
  </div>
</div>
<div class="banner" style="background:{bg};color:{fg};">
  <span class="label">{_esc(classification.classification.value)}</span>
  <span class="pts">Draft · {classification.points:+d} points</span>
</div>
<div class="section">
  <h2>Evidence &amp; points ledger</h2>
  <table>
    <thead><tr><th>Criterion</th><th>Strength</th><th>Rationale</th><th style="text-align:right">Points</th></tr></thead>
    <tbody>
      {ledger}
      <tr class="total"><td colspan="3">Total (all proposals accepted)</td>
        <td class="pts">{classification.points:+d}</td></tr>
    </tbody>
  </table>
</div>
<div class="section">
  <h2>Annotation</h2>
  <div class="facts">
    <div><span class="k">gnomAD popmax AF</span><span class="v">{_fmt_af(features.gnomad_popmax_af)}</span></div>
    <div><span class="k">gnomAD global AF</span><span class="v">{_fmt_af(features.gnomad_af)}</span></div>
    <div><span class="k">gnomAD allele count</span><span class="v">{_esc(features.gnomad_ac)}</span></div>
    <div><span class="k">REVEL</span><span class="v">{_esc(features.revel)}</span></div>
    <div><span class="k">SpliceAI Δ max</span><span class="v">{_esc(features.spliceai_ds_max)}</span></div>
    <div><span class="k">ClinVar</span><span class="v">{_esc(features.clinvar_sig)}</span></div>
  </div>
</div>
{path_section}
{flags_section}
{conflict}
<div class="foot">
  Generated {ts}. All criteria are machine-proposed drafts requiring analyst review; PVS1
  assumes loss-of-function is the disease mechanism unless overridden. In-silico strengths use
  ClinGen-calibrated thresholds (REVEL: Pejaver 2022, <code>PMID 36413997</code>; SpliceAI:
  Walker 2023, <code>PMID 37352859</code>). Not a validated clinical device.
</div>
</div></body></html>"""


def build_report_for(
    variant: GenomicVariant,
    features: AnnotationFeatures,
    thresholds: Optional[Thresholds] = None,
) -> str:
    """Convenience: propose → classify (all accepted) → render, in one call."""
    ps = propose(variant, features, thresholds)
    result = classify(ps.accepted_calls())
    return render_variant_report(variant, features, ps, result)


def write_report(path: str | Path, html_text: str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(html_text, encoding="utf-8")
    return p
