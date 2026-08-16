# SKAT MMP Analysis

This runs a SKAT (Sequence Kernel Association Test) analysis to quantify the association between individual genes in
Million Mutation Project (MMP) *C. elegans* strains and a continuous stress-induced sleep (SIS) phenotype. It walks
through every file, from the raw MMP data to the final ranked results, so the whole thing can be
reproduced.

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
        SKAT_all-pvals.results   SKAT_filtered_Na5.results
        SKAT_all-qvals.results   (≥5-allele subset; final results)
```

---

## What each file and step does (quick reference)

**Scripts**

| File | What it does |
|---|---|
| `make_gene_variants.py` | Filters the raw MMP VCF down to [non-synonymous](https://en.wikipedia.org/wiki/Non-synonymous_substitution) variants. |
| `run_pipeline.R` | Builds the gene-to-variant map, attaches phenotypes, and runs SKAT. |
| `Make_SSID_file.R` | Standalone helper to build the gene-to-variant map from the command line. |

The following scripts are kept in `scripts/debug/` for reference and are not part of the main pipeline:

| File | What it does |
|---|---|
| `SKAT_Run1.R` | Older step-by-step version of the prep (SSID + phenotype setup). |
| `SKAT_prelim.R` | Interactive version of the SKAT run. |

**Inputs**

The ready-to-use inputs are bundled in `inputs/inputs.zip`, which contains `gene_variants.txt`,
`MMP.bed`, `MMP.bim`, and `MMP.fam`. Unzip it into `inputs/` before running the pipeline
(`combined_phenotype.csv` is already unzipped in `inputs/`).

| File | What it is |
|---|---|
| `MMP.vcf` | The raw Million Mutation Project variant data (not included, too large). Source: [MMP data portal](https://genome.sfu.ca/mmp/). |
| `gene_variants.txt` | The filtered variant list, containing only the [non-synonymous](https://en.wikipedia.org/wiki/Non-synonymous_substitution) variants. |
| `MMP.bed` / `.bim` / `.fam` | Strain genotypes in PLINK format and the FAM also holds the phenotype. |
| `combined_phenotype.csv` | The sleep phenotype Z-scores for each screened MMP strain. |
| `MMP.SSID` | A map telling SKAT which variants belong to which gene (built during the run). |

**Outputs**

| File | What it is |
|---|---|
| `SKAT_all-pvals.results` | A p-value for every gene after B correction, lower p implies stronger association with changes in phenoytpe |
| `SKAT_all-qvals.results` | The sorted version of SKAT_all-pvals.results with an additional column showing q-values after FDR control. |
| `SKAT_filtered_Na5.results` | The ranked list limited to genes with at least 5 tested non-synonymous alleles (6,663 genes) — the final list used for candidate selection. |
| `validation_genes.tsv` | The 15 known sleep genes used to validate the ranking, with their rank and percentile. |

**Steps (in order)**

1. **Filter variants** — `make_gene_variants.py` turns the raw VCF into `gene_variants.txt`.
2. **Map genes to variants** — `run_pipeline.R` builds `MMP.SSID`.
3. **Attach phenotypes** — match each strain's sleep phenotype to its genotype.
4. **Run SKAT** — test every gene, write the P-value results.
5. **Refine** — re-rank the genes with ≥5 tested non-synonymous alleles to get the final ranked list.

---

## Dependencies

**R packages:**
```r
install.packages(c("SKAT", "fdrtool", "dplyr", "plyr", "stringr"))
```

**PLINK** (only needed if rebuilding `MMP.bed` from the VCF):
[www.cog-genomics.org/plink](https://www.cog-genomics.org/plink/)

---

## Running the pipeline

The ready-to-use inputs are in `inputs.zip`, so most runs are just two commands:

```bash
unzip inputs/inputs.zip -d inputs/
Rscript scripts/run_pipeline.R
```

`run_pipeline.R` does four things in order:

| Step | What it does | Produces |
|---|---|---|
| 1 | Build the gene-to-variant map (defines 19,749 genes) | `MMP.SSID` |
| 2 | Match each strain's phenotype to its genotype | updated `MMP.fam` |
| 3 | Run SKAT on every gene, then add FDR q-values and sort | `SKAT_all-pvals.results`, `SKAT_all-qvals.results` |
| 4 | Re-rank the genes that have ≥5 tested non-synonymous alleles | `SKAT_filtered_Na5.results` (final) |

Temporary files (`MMP.SSID`, `MMP.SSD`, etc.) are written to `inputs/` and regenerated each run.

---

## Source data

| File | What it is | Where it's from |
|---|---|---|
| `MMP.vcf` | The full MMP VCF: ~841k variants across 2,007 strains (~6.4 GB). | Thompson et al. 2013 ([PMC3787271](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787271/)); download from the [MMP data portal](https://genome.sfu.ca/mmp/). **Not included here** — too large. |
| `combined_phenotype.csv` | The sleep phenotype for each of the 940 screened strains (`Strain,Phenotype`). | This study's behavioral screen. |

---

## Building the variant list (`gene_variants.txt`)

We keep a variant only if it changes a protein **and** belongs to an annotated gene:

- **SNVs** (single-nucleotide variants)
- **Indels** (insertions/deletions)
- **Structural variants**

Everything else is dropped: synonymous, intronic, intergenic, UTR, and non-coding-RNA variants. This leaves 191,938 variants (from ~841k in the raw VCF).

```bash
python3 scripts/make_gene_variants.py --vcf MMP.vcf --out inputs/gene_variants.txt \
        --verify inputs/gene_variants.txt   # optional: check against a reference
```

This reproduces the deposited `gene_variants.txt` exactly.

> A plain `grep` for coding tags does not reproduce this file. It keeps the synonymous and
> intronic variants too. Use the script.
>
> This filter is also why some genes never appear in the results. For example, *flp-13* is
> represented by three variants, but all are in the non-coding region, so they're dropped here.

---

## Rebuilding the genotype files (`MMP.bed` / `.bim` / `.fam`)

Only needed if you want to regenerate the PLINK files from the raw VCF. Build the strain-keep
list, then convert:

```bash
awk -F',' 'NR>1 {print $1, $1}' inputs/combined_phenotype.csv > keep_strains.txt

plink --vcf MMP.vcf \
      --keep keep_strains.txt \
      --double-id \
      --allow-extra-chr \
      --allow-no-sex \
      --make-bed \
      --out inputs/MMP
```

> `--allow-extra-chr` is needed because *C. elegans* uses Roman-numeral chromosome names (I–V, X).
> Strains not in the VCF are dropped silently. The strain order in `MMP.fam` must match `MMP.bed`;
> `run_pipeline.R` keeps them aligned when it updates phenotypes, so don't reorder the FAM by hand.

---

## Outputs

See [`outputs/README.md`](outputs/README.md) for the column details.

### What the gene counts mean

| Number | Meaning |
|---|---|
| 19,749 | genes defined (had ≥1 non-synonymous variant) |
| 18,070 | tested and ranked |
| 1,679 | had no usable marker, left unranked (P = 1) |
| 6,663 | genes with ≥5 tested non-synonymous alleles (the final filtered file) |
| 191,938 | variants in `gene_variants.txt` |
