# Outputs — file and column descriptions

This folder contains the SKAT gene-level association results for the MMP SIS screen.
Each file is whitespace-delimited with a header row.

## Gene counts (read this first)

The gene-set definition (`MMP.SSID`) defines **19,749 gene sets** — one per gene that
has at least one annotated coding variant in the source VCF.

Not every defined gene set can be tested. SKAT drops variants that are monomorphic
across the analyzed strains, so a gene set can end up with zero testable markers:

| Category | Count | Meaning |
|---|---|---|
| Gene sets defined | 19,749 | one per gene with ≥1 annotated coding variant |
| Gene sets **tested** (`N.Marker.Test` > 0) | 18,070 | received a real SKAT P-value and are ranked |
| Gene sets **untested** (`N.Marker.Test` = 0) | 1,679 | no polymorphic marker in the analyzed strains; assigned P = 1, **not ranked** |

Two further filtered counts appear in the paper/pipeline:

| Number | Definition |
|---|---|
| 15,786 | `SKAT_all_reduced_940.results` — genes with `N.Marker.All ≥ 4` (≥4 *annotated* variants) |
| 6,663 | genes with `N.Marker.Test ≥ 5` (≥5 *testable/polymorphic* variants) |

> The results are therefore a **ranked list of the 18,070 tested gene sets**, padded with
> 1,679 untested gene sets carrying a placeholder P = 1. They are not a ranking of all
> 19,749 genes.

**Genes absent from these files** (e.g. *flp-13* / F33D4.3) carry **no qualifying coding
variant** in the MMP strains analyzed under the `SN=` / `CODING=` INFO filter, so no gene
set was defined for them. This is an inherent limitation of association testing on a fixed
mutation panel: SKAT can only test genes that the panel actually mutates.

## Files

### `SKAT_all-pvals.results`
All 19,749 gene sets in **alphabetical SetID order** (as returned by `SKAT.SSD.All`).

| Column | Description |
|---|---|
| `SetID` | Gene model / sequence name (e.g. `ZK1067.1`) |
| `P.value` | SKAT association P-value (1 for untested sets) |
| `N.Marker.All` | Variants assigned to the gene set |
| `N.Marker.Test` | Variants actually used (polymorphic in analyzed strains); **0 = untested** |

### `SKAT_all-qvals.results`
Same rows as above **sorted ascending by P.value**, with an added FDR column.

| Column | Description |
|---|---|
| `SetID`, `P.value`, `N.Marker.All`, `N.Marker.Test` | as above |
| `Q.value` | FDR q-value from `fdrtool` (`cutoff.method = "fndr"`) |

Because the 1,679 untested sets all have P = 1, they sort to the bottom; the meaningful
ranking is over the 18,070 tested sets above them.

### `SKAT_all_reduced_940.results`
Restricted to gene sets whose gene had **≥ 4 annotated variants** in the SSID, then
re-run and sorted by P.value with fresh FDR q-values. **This is the file used for the
final results.** (`940` refers to the 939-strain phenotype panel used.)

Note: this file still contains 343 rows with `N.Marker.Test = 0` — genes that met the
≥4-annotated-variant threshold but whose variants were monomorphic in the analyzed
strains. Treat only rows with `N.Marker.Test > 0` as tested.

## Reproducing
See the top-level `README.md`. Run `Rscript scripts/run_pipeline.R` from the repo root.
