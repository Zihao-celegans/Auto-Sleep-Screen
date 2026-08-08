# Outputs — file and column descriptions

SKAT gene-level association results for the MMP SIS screen. Files are whitespace-delimited
with a header row.

## Gene counts (read this first)

We defined **19,749 gene sets** — one per gene that has at least one protein-altering variant
in our filtered variant list (see below).

Of these, **18,070 were tested** and given a P-value. 

| Count | Meaning |
|---|---|
| 19,749 | gene sets defined |
| 18,070 | tested (`N.Marker.Test > 0`) — the ranked list |
| 1,679 | not tested (`N.Marker.Test = 0`, P = 1) |

> **In short:** the results rank the 18,070 tested genes. The 1,679 untested genes (P = 1) sit
> at the bottom as placeholders — the list is not a ranking of all 19,749 genes.

## Why some genes are missing entirely

Before running SKAT, we kept only these variants: nonsynonymous SNVs, coding indels, and coding structural variants. We removed everything else:
synonymous, intronic, intergenic, UTR, and non-coding-RNA variants (see the top-level `README.md`
and `scripts/make_gene_variants.py`). Genes mutated only
by non-coding or silent variants drop out. For example, *flp-13* (F33D4.3) is mutated in the
MMP. It has three variants but all three are non-coding, so our filter removed them and no gene set was built for it. 

## Files

### `SKAT_all-pvals.results`
All 19,749 gene sets

| Column | Description |
|---|---|
| `SetID` | Gene sequence name (e.g. `ZK1067.1`) |
| `P.value` | SKAT P-value (1 if untested) |
| `N.Marker.All` | Variants assigned to the gene |
| `N.Marker.Test` | Variants actually used by SKAT (**0 = untested**) |

### `SKAT_all-qvals.results`
Same rows, sorted by P-value, plus an FDR `Q.value` (`fdrtool`). The 1,679 untested genes (P = 1)
sort to the bottom; the real ranking is the tested genes above them.

### `SKAT_all_reduced_940.results`
Genes with ≥4 annotated variants, re-run and sorted by P-value with fresh q-values. This is the filter found in original SKAT analysis. We found our own threshold by analyzing the median and mean percentile rankings.



## Reproducing
Run `Rscript scripts/run_pipeline.R` from the repo root. See the top-level `README.md`.
