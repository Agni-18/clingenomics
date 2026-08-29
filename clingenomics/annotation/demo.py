"""Run: python -m clingenomics.annotation.demo — pipeline over the bundled sample VCF."""
from pathlib import Path
from .vcf import read_vcf
from .proposer import propose
from ..germline.classifier import classify

SAMPLE = Path(__file__).resolve().parents[2] / "data" / "sample.vcf"


def main() -> None:
    for v, f in read_vcf(SAMPLE):
        ps = propose(v, f)
        draft = classify(ps.accepted_calls())
        print(f"{v.gene_symbol or '?':6} {v.hgvs_c or v.key:16} "
              f"{v.consequence.value if v.consequence else '-':24}")
        for p in ps.proposals:
            print(f"     proposed: {p.code:4} [{p.strength.label}]  {p.rationale}")
        for fl in ps.flags:
            print(f"     FLAG    : {fl}")
        print(f"     => DRAFT (all accepted): {draft.classification.value} "
              f"({draft.points:+d} pts)\n")


if __name__ == "__main__":
    main()
