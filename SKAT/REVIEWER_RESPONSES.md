# Responses to reviewer comments — SKAT analysis & data deposit

This document answers every reviewer comment that concerns the SKAT analysis, the
gene counts, and the GitHub/data deposit. All numbers are verified against the
deposited files. Manuscript-prose-only comments (wording, citations, figures) are
tracked separately in the author response letter.

---

## Reviewer 1, comment 2 — describe the output files; reconcile gene counts

**Comment:** The output files need description; the gene counts (19,749 vs. genes with
P=1/zero markers) are unclear; some genes (e.g. flp-13) are missing; the "6,663" in the
abstract is unexplained.

**Response — all addressed:**

- **File documentation added:** `outputs/README.md` now describes every file and column.
  The full data flow is documented in `PIPELINE.md`.
- **Gene counts reconciled.** Of 19,749 gene sets *defined*, only **18,070 are testable**
  (`N.Marker.Test > 0`) and receive a real P-value; the other **1,679** have no testable
  marker (`N.Marker.Test = 0`, P = 1) and are unranked. The results are therefore a ranked
  list of 18,070 tested gene sets, not 19,749. The manuscript text has been corrected
  accordingly.
- **The "6,663" is explained:** it is the number of gene sets with **≥5 testable variants**
  (`N.Marker.Test ≥ 5`) — verified exactly (1,704 with exactly 5 + 4,959 with more). Full
  reconciliation:

  | Number | Definition |
  |---|---|
  | 19,749 | gene sets defined |
  | 18,070 | testable (ranked list) |
  | 1,679 | untested (P=1) |
  | 15,786 | reduced file (≥4 annotated variants) |
  | 6,663 | abstract figure (≥5 testable variants) |

- **Missing genes (flp-13) explained.** *flp-13* (F33D4.3) is absent because it has no
  **nonsynonymous coding** variant in the MMP. Its three MMP variants are two introns and
  one non-coding exon variant, all removed by the coding filter (see below). This is an
  inherent limitation of association testing on a fixed mutation panel: SKAT can only test
  genes the panel mutates in a protein-altering way. This caveat is now stated in
  `outputs/README.md`.

---

## Reviewer 1, comment 3 — validation transparency (the "first 53%" claim)

**Comment:** The list of known sleep genes is a small subset and never enumerated; the
claim that all known regulators rank in the first 53% is contradicted by let-23 (ZK1067.1,
rank 16,016/19,749) and npr-1 (C39E6.6, rank 13,445/19,749).

**Response — the claim was corrected.** The reviewer is right; the "53%" statement is not
supported. `outputs/validation_genes.tsv` now lists the known sleep genes explicitly with
their WormBase/sequence IDs, P-values, `N.Marker.Test`, rank, and percentile. The data show
**partial enrichment**, not uniform enrichment:

| Gene | SetID | P | Percentile (of 18,070 tested) |
|---|---|---|---|
| egl-4 | F55A8.2 | 0.015 | 3.5% |
| daf-16 | R13H8.1 | 0.017 | 3.8% |
| aptf-1 | F21A9.2 | 0.063 | 8.6% |
| npr-1 | C39E6.6 | 0.74 | 74% |
| unc-31 | ZK897.1 | 0.88 | 87% |
| let-23 | ZK1067.1 | 0.89 | 89% |
| flp-13 | F33D4.3 | — | absent (no coding variant) |

The revised text reports this honestly: several core sleep genes (egl-4, daf-16, aptf-1)
enrich strongly, while others do not — consistent with the MMP capturing loss-of-function-
biased coding alleles and with UV-SIS depending on only a subset of general sleep machinery.
*(The canonical known-gene list and the gene→SetID mappings should be finalized against
WormBase before submission.)*

---

## Reviewer 1, comment 4 — justify the 83-candidate cutoff; supplemental table

**Comment:** More detail is needed on how 83 candidates were chosen; provide a supplemental
table showing the top of the ranked list with the 83 candidates indicated.

**Response — pending data.** A supplemental table will show the top of
`SKAT_all_reduced_940.results` with the 83 tested candidates flagged and the selection rule
stated explicitly. *(This requires the author-held list of the 83 tested candidates; it can
be built directly against the deposited ranked list once that list is provided.)*

---

## Reviewer 2, comment 1 — deposit per-strain baseline & post-UV quiescence scores

**Comment:** The baseline and post-UV quiescence scores for each MMP strain should be
included (e.g. in Table S1).

**Response — pending data.** `combined_phenotype.csv` currently holds only the normalized
composite phenotype. The raw baseline and post-UV component scores will be added as columns
to the strain table. *(Requires the author-held raw scores.)*

---

## Reviewer 2, comment 2 — document robot vs. manual strains

**Comment:** Document which strains were screened by robot (189) vs. manually (752).

**Response — pending data.** A `method` column (automated/manual) will be added to the
strain table. *(Requires the author-held strain-to-method assignment.)*

---

## Senior editor — reproducibility, permanent repository, data availability

**Comment:** All computational methods must be documented and available at a permanent
repository (Zenodo/FigShare); follow the Screen Report and G3 data-availability guidance.

**Response:**

- **Reproducibility gap closed.** The previously documented command for building
  `gene_variants.txt` (`grep -E "SN=|CODING="`) did **not** reproduce the deposited file
  (it yields 838,661 lines, not 191,938). The real filter — keep nonsynonymous coding SNVs,
  induced indels, and coding structural variants, on real genes — is now implemented and
  validated in `scripts/make_gene_variants.py` (regenerates the file byte-for-byte) and
  documented in `PIPELINE.md`. The README has been corrected.
- **Full pipeline documented:** `PIPELINE.md` traces every file from `MMP.vcf` to the final
  results, with exact commands.
- **Permanent DOI:** the repository will be archived on Zenodo (linked GitHub release) and
  the DOI cited in Data Availability. *(One-time setup; pending.)*
- **Source VCF:** `MMP.vcf` is not redeposited (size); it is cited to Thompson et al. 2013 /
  WormBase, per that project's requested citation.

---

## Status summary

| Item | Status |
|---|---|
| Output-file documentation (`outputs/README.md`) | ✅ done |
| Full pipeline documentation (`PIPELINE.md`) | ✅ done |
| Gene-count reconciliation incl. 6,663 | ✅ done |
| flp-13 / missing-gene explanation | ✅ done |
| Validation table replacing "53%" (`validation_genes.tsv`) | ✅ done |
| Reproducible filter script (`make_gene_variants.py`) | ✅ done & verified |
| README filter description corrected | ✅ done |
| 83-candidate supplemental table | ⏳ needs author list |
| Per-strain baseline/post-UV scores + method column | ⏳ needs author raw data |
| Zenodo DOI | ⏳ one-time setup |
