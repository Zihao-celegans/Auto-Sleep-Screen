# SKAT MMP Pipeline — full data flow, input to output

This document traces every file in the SKAT analysis from the raw source data to
the final ranked results, and explains how each step connects. It exists so the
analysis is fully reproducible (senior-editor data-availability requirement).

---

## Overview

```
 MMP.vcf (6.4 GB, raw)                  combined_phenotype.csv
        │                                        │
        │  scripts/make_gene_variants.py         │  (per-strain SIS phenotype)
        ▼                                        │
 gene_variants.txt ──┐                           │
                     │  run_pipeline.R Step 1     │
                     ▼                            │
                 MMP.SSID (gene → variant map)    │
                                                  │
 MMP.bed / MMP.bim / MMP.fam ─────────────────────┤  run_pipeline.R Step 2
        (PLINK genotypes)                         ▼
                                          MMP.fam (phenotype column updated)
                     │                            │
                     └──────────┬─────────────────┘
                                ▼  run_pipeline.R Steps 3–4 (SKAT)
                     ┌──────────┴───────────┐
                     ▼                      ▼
        SKAT_all-pvals.results   SKAT_all_reduced_940.results
        SKAT_all-qvals.results   (≥4-variant subset; final results)
```

---

## Source data

| File | What it is | Where it comes from |
|---|---|---|
| `MMP.vcf` | Full Million Mutation Project VCF: ~841k variant records × 2,007 strain genotype columns (~6.4 GB). Reference WS225. | Thompson et al. 2013 ([PMC3787271](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787271/)); WormBase FTP. **Not deposited in this repo** (too large; cite to its permanent home). |
| `combined_phenotype.csv` | `Strain,Phenotype` — the normalized continuous SIS phenotype for each screened strain (939 strains). | This study's behavioral screen. |

## Deposited inputs (`inputs/inputs.zip`)

| File | What it is | Produced by |
|---|---|---|
| `gene_variants.txt` | Header-stripped, coding-filtered variant list: 8 VCF info columns for 191,938 protein-altering variants (no genotype columns). | `scripts/make_gene_variants.py` (from `MMP.vcf`) |
| `MMP.bed` / `MMP.bim` / `MMP.fam` | PLINK binary genotypes for the screened strains. `MMP.fam` column 6 holds the phenotype. | `plink --vcf MMP.vcf --keep <strains> --make-bed` (see README) |

---

## Step 0 — Build `gene_variants.txt` from the VCF  (`scripts/make_gene_variants.py`)

Turns the 6.4 GB raw VCF into the 191,938-variant coding list.

**Keeps a variant if it is protein-altering and maps to a real gene:**
- **SNV** — nonsynonymous coding change (`AAC=X->Y`, X≠Y: missense or stop). Synonymous and non-coding SNVs dropped.
- **Indel** — any induced small indel (`INDEL=` present).
- **Structural variant** — carries a coding annotation (`SVTYPE=` and `CODING=`).
- **Gene requirement** — `SN=` must be a real gene (not a chromosome name, empty/NA, or a non-standard `E_…` / `cTel…` model).

**Drops:** intergenic, intronic, synonymous, UTR-only, ncRNA, and unassigned variants.
Then keeps only the 8 fixed VCF columns (drops the ~2,007 genotype columns). No strain subsetting.

```bash
python3 scripts/make_gene_variants.py --vcf MMP.vcf --out inputs/gene_variants.txt \
        --verify inputs/gene_variants.txt   # optional byte-for-byte check
```

Verified: regenerates the deposited `gene_variants.txt` **exactly** (191,938 variants, identical content).

> **Note:** This coding filter is why some annotated genes have no gene set and are
> never tested — e.g. *flp-13* (F33D4.3), whose only MMP variants are two introns and
> one non-coding exon variant, none of which is a nonsynonymous coding change.

---

## Steps 1–4 — SKAT association  (`scripts/run_pipeline.R`)

1. **Build the SSID** — parse `gene_variants.txt`, extract each variant's gene from the
   `SN=`/`CODING=` tag → `inputs/MMP.SSID` (a gene→variant map defining 19,749 gene sets).
2. **Attach phenotypes** — `left_join` `combined_phenotype.csv` onto `MMP.fam` (keyed on
   strain, preserving BED order) → updates `MMP.fam` column 6.
3. **Full SKAT run** — `Generate_SSD_SetID` + `SKAT.SSD.All` over all gene sets →
   `outputs/SKAT_all-pvals.results`; add FDR q-values, sort → `SKAT_all-qvals.results`.
4. **Reduced run** — restrict to genes with ≥4 annotated variants, re-run →
   `outputs/SKAT_all_reduced_940.results` (**the final results**).

```bash
unzip inputs/inputs.zip -d inputs/
Rscript scripts/run_pipeline.R
```

Intermediates written to `inputs/` (regenerated each run, git-ignored):
`MMP.SSID`, `MMP-reduced.SSID`, `MMP.SSD`, `MMP.info`, `MMP-reduced.SSD`, `MMP-reduced.info`.

---

## Outputs (`outputs/`) — see `outputs/README.md` for column-level detail

| File | Contents |
|---|---|
| `SKAT_all-pvals.results` | P-value per gene set, alphabetical. 19,749 defined / 18,070 testable / 1,679 untested. |
| `SKAT_all-qvals.results` | Same, sorted by P-value, with FDR q-values. |
| `SKAT_all_reduced_940.results` | Genes with ≥4 variants, sorted by P-value — **final results**. |
| `validation_genes.tsv` | Known sleep genes with real ranks/percentiles (validation). |

### Gene-count reconciliation (every number that appears in the paper)

| Number | Definition |
|---|---|
| 19,749 | gene sets *defined* (≥1 coding variant in `gene_variants.txt`) |
| 18,070 | *testable* (`N.Marker.Test > 0`) — the actual ranked list |
| 1,679 | *untested* (`N.Marker.Test = 0`, P = 1) |
| 15,786 | reduced file (`N.Marker.All ≥ 4`) |
| 6,663 | genes with `N.Marker.Test ≥ 5` (the abstract figure) |
| 191,938 | coding variants in `gene_variants.txt` |
