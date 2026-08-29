"""Pure-Python VCF reader (no cyvcf2 / pysam — installs everywhere, incl. Windows).

Yields (GenomicVariant, AnnotationFeatures) pairs. Field names vary across
annotation stacks, so every INFO tag is configurable via `VcfFieldMap`; the
defaults match a common VEP + SpliceAI + gnomAD layout (and the bundled
data/sample.vcf).

Scope note: this reads one ALT per record. Decompose and normalise upstream
(`bcftools norm -m- -f ref.fa`) before feeding real multiallelic VCFs; that's a
solved problem and not worth re-implementing here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from ..core.variant import Assembly, GenomicVariant, MolecularConsequence, VariantType
from .features import AnnotationFeatures


@dataclass
class VcfFieldMap:
    """INFO tag names to read. Override to match your annotation pipeline."""

    gnomad_af: str = "gnomAD_AF"
    gnomad_popmax_af: str = "gnomAD_AF_grpmax"
    gnomad_ac: str = "gnomAD_AC"
    spliceai: str = "SpliceAI"          # packed: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...
    revel: str = "REVEL"
    cadd_phred: str = "CADD_PHRED"
    csq: str = "CSQ"                    # VEP packed field
    clinvar_sig: str = "CLNSIG"
    clinvar_review_status: str = "CLNREVSTAT"


# VEP SO-term -> our consequence enum (most-severe term wins; list is not exhaustive)
_VEP_CONSEQUENCE = {
    "frameshift_variant": MolecularConsequence.FRAMESHIFT,
    "stop_gained": MolecularConsequence.STOP_GAINED,
    "stop_lost": MolecularConsequence.STOP_LOST,
    "start_lost": MolecularConsequence.START_LOST,
    "splice_acceptor_variant": MolecularConsequence.SPLICE_ACCEPTOR,
    "splice_donor_variant": MolecularConsequence.SPLICE_DONOR,
    "splice_region_variant": MolecularConsequence.SPLICE_REGION,
    "missense_variant": MolecularConsequence.MISSENSE,
    "inframe_insertion": MolecularConsequence.INFRAME_INDEL,
    "inframe_deletion": MolecularConsequence.INFRAME_INDEL,
    "synonymous_variant": MolecularConsequence.SYNONYMOUS,
    "intron_variant": MolecularConsequence.INTRONIC,
    "5_prime_UTR_variant": MolecularConsequence.UTR,
    "3_prime_UTR_variant": MolecularConsequence.UTR,
}


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info == ".":
        return out
    for field_ in info.split(";"):
        if "=" in field_:
            k, v = field_.split("=", 1)
            out[k] = v
        else:
            out[field_] = "true"
    return out


def _csq_format(header_lines: List[str], csq_tag: str) -> Optional[List[str]]:
    """Extract the pipe-delimited CSQ subfield order from the VCF header."""
    marker = f"ID={csq_tag},"
    for line in header_lines:
        if line.startswith("##INFO=") and marker in line and "Format:" in line:
            fmt = line.split("Format:")[1].strip().rstrip('">').strip()
            return fmt.split("|")
    return None


def _spliceai_ds_max(raw: str) -> Optional[float]:
    """Max delta score across AG/AL/DG/DL, over all gene entries in the field."""
    best: Optional[float] = None
    for entry in raw.split(","):
        parts = entry.split("|")
        if len(parts) < 6:
            continue
        for ds in parts[2:6]:  # DS_AG, DS_AL, DS_DG, DS_DL
            try:
                val = float(ds)
            except ValueError:
                continue
            if best is None or val > best:
                best = val
    return best


def _infer_type(ref: str, alt: str) -> VariantType:
    if len(ref) == 1 and len(alt) == 1:
        return VariantType.SNV
    if len(ref) < len(alt):
        return VariantType.INSERTION
    if len(ref) > len(alt):
        return VariantType.DELETION
    return VariantType.MNV


def _first_float(info: Dict[str, str], key: str) -> Optional[float]:
    v = info.get(key)
    if v is None:
        return None
    try:
        return float(v.split(",")[0])
    except ValueError:
        return None


def read_vcf(
    path: str | Path,
    *,
    field_map: Optional[VcfFieldMap] = None,
    assembly: Assembly = Assembly.GRCH38,
) -> Iterator[Tuple[GenomicVariant, AnnotationFeatures]]:
    """Stream (variant, features) pairs from a VCF file."""
    fm = field_map or VcfFieldMap()
    header_lines: List[str] = []
    csq_fields: Optional[List[str]] = None

    with open(path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                csq_fields = _csq_format(header_lines, fm.csq)
                continue
            if not line.strip():
                continue

            cols = line.split("\t")
            chrom, pos, _id, ref, alts = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = _parse_info(cols[7]) if len(cols) > 7 else {}

            for alt in alts.split(","):
                gene = hgvs_c = hgvs_p = None
                consequence: Optional[MolecularConsequence] = None

                if csq_fields and fm.csq in info:
                    # take the first transcript block; pick-canonical is a later refinement
                    block = info[fm.csq].split(",")[0].split("|")
                    csq = dict(zip(csq_fields, block))
                    terms = csq.get("Consequence", "").split("&")
                    for t in terms:
                        if t in _VEP_CONSEQUENCE:
                            consequence = _VEP_CONSEQUENCE[t]
                            break
                    gene = csq.get("SYMBOL") or None
                    hgvs_c = csq.get("HGVSc") or None
                    hgvs_p = csq.get("HGVSp") or None

                spliceai_max = (
                    _spliceai_ds_max(info[fm.spliceai]) if fm.spliceai in info else None
                )
                popmax = _first_float(info, fm.gnomad_popmax_af)
                gaf = _first_float(info, fm.gnomad_af)

                variant = GenomicVariant(
                    assembly=assembly,
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref,
                    alt=alt,
                    variant_type=_infer_type(ref, alt),
                    gene_symbol=gene,
                    hgvs_c=hgvs_c,
                    hgvs_p=hgvs_p,
                    consequence=consequence,
                    gnomad_af=gaf,
                    spliceai_ds_max=spliceai_max,
                )

                ac = info.get(fm.gnomad_ac)
                features = AnnotationFeatures(
                    gnomad_af=gaf,
                    gnomad_popmax_af=popmax,
                    gnomad_ac=int(ac.split(",")[0]) if ac and ac.split(",")[0].isdigit() else None,
                    revel=_first_float(info, fm.revel),
                    cadd_phred=_first_float(info, fm.cadd_phred),
                    spliceai_ds_max=spliceai_max,
                    clinvar_sig=info.get(fm.clinvar_sig),
                    clinvar_review_status=info.get(fm.clinvar_review_status),
                )
                yield variant, features
