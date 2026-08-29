# Clinical Genomics Platform

Multi-modal germline + somatic variant interpretation: DNA variant calling →
annotation → **points-based ACMG** interpretation with **RNA-as-evidence**
(splicing, NMD/ASE, expression outliers, fusions) → **somatic tiering** →
clinical report.

## Where the value lives

You don't reimplement variant callers or RNA outlier tools — you orchestrate
them and own the **interpretation engine**. That engine is what this repo starts
with, because it's self-contained, testable without external data, and it's the
part that encodes clinical judgement.

## Architecture (layers)

```
                +-------------------------------------------------------+
 Orchestration  |  Nextflow / Snakemake                                 |
   (external)    |  DNA: GATK / DeepVariant / DRAGEN  (SNV/indel/CNV/SV) |
                 |  RNA: DROP (FRASER, OUTRIDER), arriba / STAR-Fusion   |
                +-----------------------------|-------------------------+
                                              v
                +-------------------------------------------------------+
 Annotation     |  VEP / SnpEff · gnomAD · ClinVar · SpliceAI · OncoKB  |
   (adapters)    |  -> raw features attached to GenomicVariant          |
                +-----------------------------|-------------------------+
                                              v
                +=======================================================+
 Interpretation ||  core/      GenomicVariant, Evidence, CriterionCall  ||   <-- built
   ENGINE       ||  germline/  ACMG criteria registry + points classifier||   <-- built
   (this repo)   ||  rna/       RNA observations -> ACMG calls (PS3/BS3)  ||   <-- built
                 ||  somatic/   AMP/ASCO/CAP 2017 tiers + OncoKB mapping  ||   <-- built
                +=======================================================+
                                              v
                +-------------------------------------------------------+
 Report         |  report/  Jinja2 -> HTML/PDF clinical report          |   <-- next
                +-------------------------------------------------------+
```

## The core design decision

Raw **`Evidence`** is separated from an interpreted **`CriterionCall`**. An RNA
splicing outlier is evidence; it *emits* a `PS3` call at a computed strength.
This indirection is what makes the points-based system pay off: any criterion
can fire at a modulated strength (VeryStrong=8, Strong=4, Moderate=2,
Supporting=1), which is exactly what functional/RNA data needs.

Germline classification thresholds (ClinGen/Tavtigian points):

| Total points | Classification        |
|-------------:|-----------------------|
| ≥ 10         | Pathogenic            |
| 6 … 9        | Likely Pathogenic     |
| 0 … 5        | Uncertain (VUS)       |
| −1 … −6      | Likely Benign         |
| ≤ −7         | Benign                |

`BA1` (stand-alone benign) forces Benign; conflicting strong P/B evidence is
flagged rather than silently averaged.

## RNA → ACMG mapping (rna/evidence.py)

| RNA observation                 | Emits        | Notes                                   |
|---------------------------------|--------------|-----------------------------------------|
| Aberrant splicing (FRASER)      | PS3          | Strong if quantified + controlled; feeds PVS1 splicing branch |
| Normal splicing despite prediction | BS3       | only when the variant predicted splicing |
| Under-expression (OUTRIDER)     | PS3          | Strong only in a LoF-mechanism gene     |
| ASE monoallelic (variant lost)  | PS3          | NMD-consistent                          |
| ASE biallelic                   | BS3          | argues against predicted NMD            |
| Fusion (arriba/STAR-Fusion)     | PS3 (supp.)  | needs gene-specific curation            |

The strength logic is deliberately centralised in one file — it's the part a lab
will tune against its own controls, and it should be validated against the
current ClinGen SVI RNA recommendations before clinical use.

## Somatic branch (somatic/tiers.py)

Fully separate from germline ACMG. AMP/ASCO/CAP 2017 four-tier by actionability
(Tier I–IV), with an OncoKB level → AMP evidence level mapping.

## Run

```bash
pip install pydantic pytest
python -m pytest -q          # 19 tests
```

## Roadmap / not yet built

1. **Annotation adapters** — DONE (v0.2): pure-Python VCF reader
   (`annotation/vcf.py`) + feature→criteria proposer (`annotation/proposer.py`).
   Proposes PM2/BS1/BA1 (frequency), PP3/BP4 (REVEL, SpliceAI), BP7; flags PVS1
   candidates for the decision tree instead of auto-applying. Run the demo:
   `python -m clingenomics.annotation.demo` (or see the snippet in step 4 below).
2. **PVS1 decision tree** — DONE (v0.3): Abou Tayoun et al. 2018 tree in
   `germline/pvs1.py`. Calibrates PVS1 strength by NMD / region criticality;
   never returns Very Strong without an NMD determination; caps + flags on
   missing info. Wired into the proposer (review-required, mechanism-assumed).
3. **Orchestration** — Nextflow/Snakemake wrapping the callers + DROP.
4. **DROP / caller adapters** — FRASER/OUTRIDER/arriba outputs → `RNAEvidence`.
5. **Report layer** — DONE (v0.4): zero-dependency HTML report in
   `report/html_report.py` (points ledger, PVS1 path, review flags, DRAFT
   banner, calibrated-threshold citations). Run `python -m clingenomics.report.demo`
   → writes reports/*.html + index.html. (PDF export is a later add-on.)
6. **Gene/disease context** — inheritance, mechanism, PM3 phasing, sign-out.

## In-silico thresholds (calibrated)

PP3/BP4 strengths use ClinGen-calibrated cut-offs, editable in
`annotation/features.py::Thresholds`:
  * REVEL — Pejaver et al. 2022 (PMID 36413997): PP3 ≥0.644/0.773/0.932
    (supp/mod/strong); BP4 ≤0.290/0.183/0.016/0.003 (supp/mod/strong/vstrong).
  * SpliceAI — Walker et al. 2023 (PMID 37352859): PP3 Δ≥0.20 (supporting);
    BP4 Δ≤0.10 (moderate); BP7 Δ<0.10.
PP3 is counted once at the strongest applicable strength across protein/splice
predictors, and is suppressed on canonical splice variants (PVS1 owns that
signal) to avoid double-counting. Gene-specific VCEPs may override these.

## Important caveat

The point weights and thresholds are stable, but the specific ACMG/RNA/somatic
*rule wording and strengths* are refined periodically by ClinGen SVI. Validate
against current guidance before any clinical use. This is decision-support
scaffolding, not a validated clinical device.
