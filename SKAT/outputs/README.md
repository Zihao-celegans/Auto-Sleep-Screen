# Outputs: File and column descriptions

SKAT gene-level association results for the MMP SIS screen. Files are whitespace-delimited
with a header row.

## Gene counts

We defined 19,749 gene sets - one per gene that has at least one protein-altering variant
in our filtered variant list (see below).

Of these, 18,070 were tested and given a P-value. 

| Count | Meaning |
|---|---|
| 19,749 | gene sets defined |
| 18,070 | tested (`N.Marker.Test > 0`) — the ranked list |
| 1,679 | not tested (`N.Marker.Test = 0`, P = 1) |

> In short, the results rank the 18,070 tested genes. The 1,679 untested genes (P = 1) sit
> at the bottom as placeholders. The list is not a ranking of all 19,749 genes.

## Why some genes are missing entirely

Before running SKAT, we kept only these variants: nonsynonymous SNVs, coding indels, and coding structural variants. We removed everything else:
synonymous, intronic, intergenic, UTR, and non-coding-RNA variants (see the top-level `README.md`
and `scripts/make_gene_variants.py`). Genes mutated only by non-coding or silent variants drop out. For example, *flp-13* (F33D4.3) is mutated in the
MMP. It has three variants but all three are non-coding, so our filter removed them and no gene set was built for it. 

## Files

### `SKAT_all-pvals.results`
All 19,749 gene sets

| Column | Description |
|---|---|
| `SetID` | Gene sequence name (e.g. `ZK1067.1`) |
| `P.value` | SKAT P-value (1 if untested) |
| `N.Marker.All` | Variants assigned to the gene |
| `N.Marker.Test` | Variants actually used by SKAT (0 = untested) |

### `SKAT_all-qvals.results`
Same rows, sorted by P-value, plus an FDR `Q.value` (`fdrtool`). The 1,679 untested genes (P = 1)
sort to the bottom; the real ranking is the tested genes above them.

### `SKAT_filtered_Na5.results`
The final list: genes with at least 5 tested non-synonymous alleles (`N.Marker.Test >= 5`),
sorted by P-value. This is 6,663 genes and is the list used for candidate selection. We chose
the threshold of 5 by analyzing the median (PMR) and mean percentile rankings of known sleep
genes across allele-count cutoffs.



## Reproducing
Run `Rscript scripts/run_pipeline.R` from the repo root. See the top-level `README.md`.
