# SKAT MMP Analysis — data flow from input to output

This runs a SKAT (Sequence Kernel Association Test) analysis linking gene-level variants in the
Million Mutation Project (MMP) *C. elegans* strains to a continuous sleep phenotype. It walks
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
        SKAT_all-pvals.results   SKAT_all_reduced_940.results
        SKAT_all-qvals.results   (≥4-variant subset; final results)
```

---

## What each file and step does (quick reference)

**Scripts**

| File | What it does |
|---|---|
| `make_gene_variants.py` | Filters the raw MMP VCF down to protein-altering variants. |
| `run_pipeline.R` | Builds the gene-to-variant map, attaches phenotypes, and runs SKAT. |
| `SKAT_Run1.R` | Older step-by-step version of the prep (SSID + phenotype setup). |
| `Make_SSID_file.R` | Standalone helper to build the gene-to-variant map from the command line. |
| `SKAT_prelim.R` | Interactive version of the SKAT run. |

**Inputs**

| File | What it is |
|---|---|
| `MMP.vcf` | The raw Million Mutation Project variant data (not included, too large). |
| `gene_variants.txt` | The filtered variant list. |
| `MMP.bed` / `.bim` / `.fam` | Strain genotypes in PLINK format; the FAM also holds the phenotype. |
| `combined_phenotype.csv` | The sleep phenotype for each screened strain. |
| `MMP.SSID` | A map telling SKAT which variants belong to which gene (built during the run). |

**Outputs**

| File | What it is |
|---|---|
| `SKAT_all-pvals.results` | A P-value for every gene. |
| `SKAT_all-qvals.results` | The same, sorted, with FDR q-values added. |
| `SKAT_all_reduced_940.results` | Results limited to genes with enough variants — the final list. |
| `validation_genes.tsv` | Known sleep genes and where they landed, used to check the method. |

**Steps (in order)**

1. **Filter variants** — `make_gene_variants.py` turns the raw VCF into `gene_variants.txt`.
2. **Map genes to variants** — `run_pipeline.R` builds `MMP.SSID`.
3. **Attach phenotypes** — match each strain's sleep phenotype to its genotype.
4. **Run SKAT** — test every gene, write the P-value results.
5. **Refine** — re-run on genes with enough variants to get the final ranked list.

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
| 4 | Re-run on genes with ≥4 variants | `SKAT_all_reduced_940.results` (final) |

Temporary files (`MMP.SSID`, `MMP.SSD`, etc.) are written to `inputs/` and regenerated each run.

---

## Source data

| File | What it is | Where it's from |
|---|---|---|
| `MMP.vcf` | The full MMP VCF: ~841k variants across 2,007 strains (~6.4 GB). | Thompson et al. 2013 ([PMC3787271](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787271/)); WormBase. **Not included here** — too large. |
| `combined_phenotype.csv` | The sleep phenotype for each of the 939 screened strains (`Strain,Phenotype`). | This study's behavioral screen. |

---

## Building the variant list (`gene_variants.txt`)

We keep a variant only if it changes a protein **and** belongs to a real gene:

- **SNVs** — kept if nonsynonymous (missense or stop); synonymous ones are dropped.
- **Indels** — kept (small induced insertions/deletions).
- **Structural variants** — kept if they have a coding effect.

Everything else is dropped: synonymous, intronic, intergenic, UTR, and non-coding-RNA variants,
plus anything not assigned to a named gene. This leaves 191,938 variants (from ~841k in the raw VCF).

```bash
python3 scripts/make_gene_variants.py --vcf MMP.vcf --out inputs/gene_variants.txt \
        --verify inputs/gene_variants.txt   # optional: check against a reference
```

This reproduces the deposited `gene_variants.txt` exactly.

> A plain `grep` for coding tags does not reproduce this file. It keeps the synonymous and
> intronic variants too. Use the script.
>
> This filter is also why some genes never appear in the results. For example, *flp-13* is
> mutated in the MMP, but all of its variants are non-coding, so they're dropped here.

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

| File | Contents |
|---|---|
| `SKAT_all-pvals.results` | A P-value for every gene (19,749 total; 18,070 tested). |
| `SKAT_all-qvals.results` | Same, sorted by P-value, with FDR q-values. |
| `SKAT_all_reduced_940.results` | Genes with ≥4 variants — the final results. |
| `validation_genes.tsv` | Known sleep genes with their ranks (used to check the method). |

### What the gene counts mean

| Number | Meaning |
|---|---|
| 19,749 | genes defined (had ≥1 protein-altering variant) |
| 18,070 | tested and ranked |
| 1,679 | had no usable marker, left unranked (P = 1) |
| 15,786 | genes with ≥4 variants (the reduced file) |
| 191,938 | variants in `gene_variants.txt` |
