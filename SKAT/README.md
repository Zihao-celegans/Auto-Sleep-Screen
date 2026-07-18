# Sequence Kernel Association Test (SKAT)

Sequence Kernel Association Test (SKAT) pipeline for gene-level association analysis of stress-induced sleep (SIS) phenotype in the Million Mutation Project (MMP) *C. elegans* strains.

---

## Folder Structure

```
SKAT/
├── README.md
├── inputs/
│   ├── inputs.zip               # MMP.bed, MMP.bim, MMP.fam, gene_variants.txt
│   └── combined_phenotype.csv   # Phenotype data: Strain, Phenotype (normalized)
├── scripts/
│   ├── run_pipeline.R           # Main script — run this
│   ├── SKAT_Run1.R              # Step-by-step prep (SSID + FAM update)
│   ├── Make_SSID_file.R         # Command-line SSID generator
│   └── SKAT_prelim.R            # Interactive SKAT analysis
└── outputs/
    ├── SKAT_all-pvals.results         # P-values for all 19,749 defined gene sets (18,070 testable)
    ├── SKAT_all-qvals.results         # P-values + FDR q-values, sorted
    └── SKAT_all_reduced_940.results   # Final variants
```

---

## Dependencies

**R packages:**
```r
install.packages(c("SKAT", "fdrtool", "dplyr", "plyr", "stringr"))
```

**PLINK** (only needed if regenerating `MMP.bed` from a VCF):  
Download from [www.cog-genomics.org/plink](https://www.cog-genomics.org/plink/)

---

## Running the Pipeline

**Step 1** — Unzip the input files:
```bash
unzip inputs/inputs.zip -d inputs/
```

**Step 2** — Run from the `SKAT_MMP/` root directory:
```bash
Rscript scripts/run_pipeline.R
```

This will:
1. Read `inputs/gene_variants.txt` → generate `inputs/MMP.SSID`
2. Join `inputs/combined_phenotype.csv` onto `inputs/MMP.fam` (preserving strain order) → update phenotype column in place
3. Run SKAT on all gene sets → `outputs/SKAT_all-pvals.results` and `outputs/SKAT_all-qvals.results`
4. Filter genes → run SKAT again → `outputs/SKAT_all_reduced_940.results`

Intermediate files (`MMP.SSID`, `MMP-reduced.SSID`, `MMP.SSD`, `MMP.info`) are written to `inputs/`.

---

## Regenerating `MMP.bed` from a VCF

The source VCF comes from the **Million Mutation Project** (Thompson et al. 2013):
- Paper: [PMC3787271](https://pmc.ncbi.nlm.nih.gov/articles/PMC3787271/)
- Data: [ftp.wormbase.org](ftp://ftp.wormbase.org/pub/wormbase/) → navigate to `species/c_elegans/` → variations folder

If you need to rebuild the PLINK binary files from the source VCF:

**Step 1** — Prepare a keep file listing the strains you want (one per line, `FID IID` format, where FID = IID = strain name):
```bash
awk -F',' 'NR>1 {print $1, $1}' inputs/combined_phenotype.csv > keep_strains.txt
```

**Step 2** — Convert VCF to PLINK BED, keeping only those strains:
```bash
plink --vcf MMP.vcf \
      --keep keep_strains.txt \
      --double-id \
      --allow-extra-chr \
      --allow-no-sex \
      --make-bed \
      --out inputs/MMP
```

**Step 3** — Re-run the pipeline as normal.

---

## Preparing `gene_variants.txt` from a VCF

`gene_variants.txt` is a header-stripped, **coding / protein-altering** subset of the VCF
(191,938 variants; the per-strain genotype columns are dropped). A variant is kept if it is:

- a **nonsynonymous** SNV (`AAC=X->Y`, X≠Y — missense or stop), **or**
- an **induced indel** (`INDEL=` present), **or**
- a **structural variant** with a coding annotation (`SVTYPE=` and `CODING=`),

and maps to a **real gene** (`SN=` is a gene, not a chromosome name / empty / `E_…` / `cTel…`).
Intergenic, intronic, synonymous, UTR-only, ncRNA, and unassigned variants are dropped.

Generate it with the provided script (regenerates the deposited file byte-for-byte):

```bash
python3 scripts/make_gene_variants.py --vcf MMP.vcf --out inputs/gene_variants.txt \
        --verify inputs/gene_variants.txt   # optional: check against a reference
```

> A plain `grep -E "SN=|CODING="` does **not** reproduce `gene_variants.txt` — it keeps
> ~838k lines (all annotated variants, including synonymous and intronic). The coding /
> nonsynonymous filter above is required. See `PIPELINE.md` for the full data flow.

---

## Input File Formats

| File | Format | Notes |
|---|---|---|
| `MMP.bed` / `MMP.bim` / `MMP.fam` | Binary PLINK | FAM col 6 = phenotype; strain order must match BED |
| `gene_variants.txt` | Tab-separated VCF body (no header), coding/nonsynonymous only | Built by `scripts/make_gene_variants.py`; see "Preparing gene_variants.txt" |
| `combined_phenotype.csv` | CSV with header `Strain,Phenotype` | Normalized continuous trait values |

**Critical:** The strain order in `MMP.fam` must match the sample encoding in `MMP.bed`. Never reorder the FAM independently — `run_pipeline.R` uses `left_join` (keyed on FAM) to guarantee this when updating phenotypes.

---

## Output Files

See [`outputs/README.md`](outputs/README.md) for full per-column descriptions and the
gene-count breakdown.

| File | Contents |
|---|---|
| `SKAT_all-pvals.results` | SetID, P.value, N.Marker.All, N.Marker.Test for all 19,749 defined gene sets (of which **18,070 are testable** — see note) |
| `SKAT_all-qvals.results` | Same + FDR Q.value (fdrtool), sorted by P.value |
| `SKAT_all_reduced_940.results` | Genes with ≥4 variants, sorted by P.value — **use this for final results** |

> **Note on gene counts:** 19,749 gene sets are *defined*, but only **18,070** contain at
> least one testable (polymorphic) variant and receive a real P-value. The remaining
> **1,679** have no testable marker (`N.Marker.Test = 0`, `P = 1`) and are unranked. Genes
> absent entirely (e.g. *flp-13* / F33D4.3) carry no qualifying coding variant in the
> analyzed MMP strains. See [`outputs/README.md`](outputs/README.md).

