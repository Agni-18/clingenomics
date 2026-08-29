"""Run: python -m clingenomics.report.demo

Reads the bundled sample VCF and writes one HTML report per variant into
./reports/, plus an index.html linking them. Open reports/index.html in a browser.
"""

from __future__ import annotations

import html
from pathlib import Path

from ..annotation.proposer import propose
from ..annotation.vcf import read_vcf
from ..germline.classifier import classify
from .html_report import render_variant_report, write_report

SAMPLE = Path(__file__).resolve().parents[2] / "data" / "sample.vcf"
OUTDIR = Path("reports")


def _slug(variant) -> str:
    g = (variant.gene_symbol or "variant").replace("/", "_")
    return f"{g}_{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}.html"


def main() -> None:
    index_rows = []
    for variant, features in read_vcf(SAMPLE):
        ps = propose(variant, features)
        result = classify(ps.accepted_calls())
        report = render_variant_report(variant, features, ps, result)
        fname = _slug(variant)
        write_report(OUTDIR / fname, report)
        index_rows.append(
            f"<tr><td><a href='{html.escape(fname)}'>{html.escape(variant.gene_symbol or variant.key)}</a></td>"
            f"<td><code>{html.escape(variant.hgvs_c or variant.key)}</code></td>"
            f"<td>{html.escape(result.classification.value)}</td>"
            f"<td style='text-align:right'>{result.points:+d}</td></tr>"
        )
        print(f"  wrote {OUTDIR / fname}  ({result.classification.value}, {result.points:+d})")

    index = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Variant reports (DRAFT)</title>
<style>
 body{{font-family:"Segoe UI",Arial,sans-serif;color:#1c2530;background:#eef1f4;padding:40px;}}
 .wrap{{max-width:760px;margin:0 auto;background:#fff;border:1px solid #d4dae0;border-radius:10px;padding:28px 32px;}}
 h1{{font-size:22px;margin:0 0 4px;}} p.sub{{color:#6a7684;margin:0 0 20px;font-size:14px;}}
 table{{width:100%;border-collapse:collapse;font-size:14px;}}
 th,td{{text-align:left;padding:9px 10px;border-bottom:1px solid #eef1f4;}}
 th{{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:#7a8794;}}
 a{{color:#1a5091;font-weight:600;text-decoration:none;}} a:hover{{text-decoration:underline;}}
 code{{font-family:Consolas,monospace;font-size:13px;}}
</style></head><body><div class="wrap">
<h1>Germline variant reports</h1>
<p class="sub">Draft · research use · {len(index_rows)} variants from sample.vcf</p>
<table><thead><tr><th>Gene</th><th>Variant</th><th>Draft classification</th><th style="text-align:right">Points</th></tr></thead>
<tbody>{''.join(index_rows)}</tbody></table>
</div></body></html>"""
    write_report(OUTDIR / "index.html", index)
    print(f"\nOpen {OUTDIR / 'index.html'} in your browser.")


if __name__ == "__main__":
    main()
