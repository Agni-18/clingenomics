"""Run: python -m clingenomics.benchmark.demo

Runs the ClinVar concordance benchmark over the bundled SYNTHETIC fixture and
prints metrics, the mismatch taxonomy, and a confusion matrix.

NOTE: the bundled data is synthetic — these numbers exercise the harness, they
are NOT a real validation result. For a real number, pass an annotated ClinVar
release VCF:  python -m clingenomics.benchmark.demo path/to/clinvar_annotated.vcf
"""

from __future__ import annotations

import sys
from pathlib import Path

from ..germline.classifier import ACMGClassification
from .evaluate import run_benchmark

DEFAULT = Path(__file__).resolve().parents[2] / "data" / "clinvar_benchmark_synthetic.vcf"

_SHORT = {
    ACMGClassification.PATHOGENIC: "P",
    ACMGClassification.LIKELY_PATHOGENIC: "LP",
    ACMGClassification.UNCERTAIN: "VUS",
    ACMGClassification.LIKELY_BENIGN: "LB",
    ACMGClassification.BENIGN: "B",
}


def _print_confusion(result) -> None:
    tiers = list(ACMGClassification)
    cm = result.confusion_matrix()
    labels = [_SHORT[t] for t in tiers]
    print("\n  confusion matrix (rows = ClinVar, cols = engine)")
    print("            " + "".join(f"{l:>5}" for l in labels))
    for t in tiers:
        row = cm[t.value]
        print(f"    {_SHORT[t]:>7}   " + "".join(f"{row[e.value]:>5}" for e in tiers))


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    is_default = path == DEFAULT

    result = run_benchmark(path)
    if is_default:
        print("=== SYNTHETIC fixture — demonstrates the harness, NOT a real result ===\n")
    print(result.summary_text())
    _print_confusion(result)

    contradictions = result.contradictions()
    print(f"\n  contradictions (opposite direction): {len(contradictions)}")
    for r in contradictions:
        print(f"    {r.gene} {r.key}: ClinVar={_SHORT[r.clinvar_tier]} "
              f"engine={_SHORT[r.engine_tier]} ({r.engine_points:+d})")

    # show the undercalls too — the expected, safe disagreements
    under = [r for r in result.rows if r.category.value == "undercall"]
    if under:
        print(f"\n  undercalls (engine conservative — expected): {len(under)}")
        for r in under[:8]:
            print(f"    {r.gene} {r.key}: ClinVar={_SHORT[r.clinvar_tier]} "
                  f"engine={_SHORT[r.engine_tier]} ({r.engine_points:+d})")


if __name__ == "__main__":
    main()
