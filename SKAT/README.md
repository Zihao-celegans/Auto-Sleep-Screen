# SKAT MMP Analysis

This runs a SKAT (Sequence Kernel Association Test) analysis to quantify the association between individual genes in
Million Mutation Project (MMP) *C. elegans* strains and a continuous stress-induced sleep (SIS) phenotype. It walks
through every file, from the raw MMP data to the final ranked results, so the whole thing can be
reproduced.

---

## 1. Overview

This pipeline takes two inputs.

### 1.1 Inputs

* The MMP variant data (`MMP.vcf`) — the mutations carried by each strain.
* The sleep phenotype for each strain (`combined_phenotype.csv`).

### 1.2 What it does

**What it does** — it filters the variants down to non-synonymous variants, groups them
by gene, matches each strain's genotype to its phenotype, and runs SKAT to score every gene.

### 1.3 Outputs

* `SKAT_all-pvals.results` / `SKAT_all-qvals.results` — a p-value and q-value for every gene.
* `SKAT_filtered_Na5.results` — the final ranked list, limited to genes with enough alleles to
  test reliably (≥5 tested non-synonymous alleles).

The rest of this README lists each file and the exact steps.

---

## 2. What each file and step does

### 2.1 Scripts

| File                    | What it does                                                                                               |
| ----------------------- | ---------------------------------------------------------------------------------------------------------- |
| `mmp_txt_to_vcf.py`     | Converts the MMP tabular download into `MMP.vcf` (see "Building `MMP.vcf`").                               |
| `make_gene_variants.py` | Filters the raw MMP VCF down to [non-synonymous](#45-building-the-variant-list-gene_variantstxt) variants. |
| `run_pipeline.R`        | Builds the gene-to-variant map, attaches phenotypes, and runs SKAT.                                        |
| `Make_SSID_file.R`      | Standalone helper to build the gene-to-variant map from the command line.                                  |

The following scripts are kept in `scripts/debug/` for reference and are not part of the main pipeline:

| File            | What it does                                                     |
| --------------- | ---------------------------------------------------------------- |
| `SKAT_Run1.R`   | Older step-by-step version of the prep (SSID + phenotype setup). |
| `SKAT_prelim.R` | Interactive version of the SKAT run.                             |

### 2.2 Inputs

The ready-to-use inputs are bundled in `inputs/inputs.zip`, which contains `gene_variants.txt`,
`MMP.bed`, `MMP.bim`, and `MMP.fam`. Unzip it into `inputs/` before running the pipeline
(`combined_phenotype.csv` is already unzipped in `inputs/`).

| File                        | What it is                                                                                                                                                                                              |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MMP.vcf`                   | The raw Million Mutation Project variant data (not included, too large). Built from the [MMP data portal](https://genome.sfu.ca/mmp/) table via `scripts/mmp_txt_to_vcf.py` — see "Building `MMP.vcf`". |
| `gene_variants.txt`         | The filtered variant list, containing only the [non-synonymous](#45-building-the-variant-list-gene_variantstxt) variants.                                                                               |
| `MMP.bed` / `.bim` / `.fam` | Strain genotypes in PLINK format and the FAM also holds the phenotype.                                                                                                                                  |
| `combined_phenotype.csv`    | The sleep phenotype Z-scores for each screened MMP strain.                                                                                                                                              |
| `MMP.SSID`                  | A map telling SKAT which variants belong to which gene (built during the run).                                                                                                                          |

### 2.3 Outputs

| File                        | What it is                                                                                                                                  |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `SKAT_all-pvals.results`    | A p-value for every gene after Bonferroni correction, a lower p implies stronger association with changes in phenoytpes                     |
| `SKAT_all-qvals.results`    | The sorted version of SKAT_all-pvals.results with an additional column showing q-values after FDR control.                                  |
| `SKAT_filtered_Na5.results` | The ranked list limited to genes with at least 5 tested non-synonymous alleles (6,663 genes) — the final list used for candidate selection. |
| `validation_genes.tsv`      | The 15 known sleep genes used to validate the ranking, with their rank and percentile.                                                      |

### 2.4 Steps (in order)

1. **Filter variants** — `make_gene_variants.py` turns the raw VCF into `gene_variants.txt`.
2. **Map genes to variants** — `run_pipeline.R` builds `MMP.SSID`.
3. **Attach phenotypes** — match each strain's sleep phenotype to its genotype.
4. **Run SKAT** — test every gene, write the P-value results.
5. **Refine** — re-rank the genes with ≥5 tested non-synonymous alleles to get the final ranked list.

---

## 3. Dependencies

### 3.1 R packages

```r
install.packages(c("SKAT", "fdrtool", "dplyr", "plyr", "stringr"))
```

### 3.2 PLINK

**PLINK** (only needed if rebuilding `MMP.bed` from the VCF):
[www.cog-genomics.org/plink](https://www.cog-genomics.org/plink/)

---

## 4. Running the pipeline

Users can run these two commands to recreate the results, but if they want they can look at the steps below for the end-to-end creation.

### 4.1 Recreating the results

The ready-to-use inputs are in `inputs.zip`, so most runs are just two commands:

```bash
unzip inputs/inputs.zip -d inputs/
Rscript scripts/run_pipeline.R
```

`run_pipeline.R` does four things in order:

| Step | What it does                                                 | Produces                                           |
| ---- | ------------------------------------------------------------ | -------------------------------------------------- |
| 1    | Build the gene-to-variant map (defines 19,749 genes)         | `MMP.SSID`                                         |
| 2    | Match each strain's phenotype to its genotype                | updated `MMP.fam`                                  |
| 3    | Run SKAT on every gene, then add FDR q-values and sort       | `SKAT_all-pvals.results`, `SKAT_all-qvals.results` |
| 4    | Re-rank the genes that have ≥5 tested non-synonymous alleles | `SKAT_filtered_Na5.results` (final)                |

Temporary files (`MMP.SSID`, `MMP.SSD`, etc.) are written to `inputs/` and regenerated each run.

---

### 4.2 Source data

| File                     | What it is                                                                     | Where it's from                                                                                                                                                                                              |
| ------------------------ | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `MMP.vcf`                | The full MMP VCF: ~841k variants across 2,007 strains (~6.4 GB).               | Built from the MMP variant table (see below), or the original combined-calls VCF (Thompson et al. 2013, [PMC3787271](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787271/)). **Not included here** — too large. |
| `combined_phenotype.csv` | The sleep phenotype for each of the 940 screened strains (`Strain,Phenotype`). | This study's behavioral screen.                                                                                                                                                                              |

---

### 4.3 Building `MMP.vcf` from the MMP download

The MMP data portal (https://genome.sfu.ca/mmp/) distributes the variant calls as a tab-separated table (one row per strain × variant), not as a VCF. Download "Results from mutagenized strains" (file size: 76 MB, if it comes as `.xlsx`, export the first sheet to a `.txt`/`.tsv`), then convert it to a VCF:

```bash
python3 scripts/mmp_txt_to_vcf.py --txt mmp_mut_strains_data.txt --out MMP.vcf
```

This produces the multi-sample `MMP.vcf` used by the steps below (857,114 variants × 2,007
strains). The INFO field is written with the `GF=`, `SN=`, `PN=`, `AAC=`, and `INDEL=` tags
that `make_gene_variants.py` reads.

---

### 4.4 Building the variant list (`gene_variants.txt`)

We keep a variant only if it is non-synonymous:

* **SNVs** (single-nucleotide variants)
* **Indels** (insertions/deletions)
* **Structural variants**

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

### 4.5 Rebuilding the genotype files (`MMP.bed` / `.bim` / `.fam`)

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

## 5. Outputs

See [`outputs/README.md`](outputs/README.md) for the column details.

### 5.1 What the gene counts mean

| Number  | Meaning                                                                                                   |
| ------- | --------------------------------------------------------------------------------------------------------- |
| 19,749  | genes represented by ≥1 non-synonymous alleles in the MMP strains we screened                             |
| 18,070  | tested and ranked                                                                                         |
| 1,679   | genes not represented by any non-synonymous alleles in the MMP strains we screened, left unranked (P = 1) |
| 6,663   | genes with ≥5 tested non-synonymous alleles (the final filtered file)                                     |
| 191,938 | variants in `gene_variants.txt`                                                                           |
