# SIS/SKAT paper — complete summary of changes to make

All values verified against the deposited `~/SKAT_MMP/outputs/SKAT_all-pvals.results` using the
**canonical SetIDs** (from the MMP VCF / WormBase WS225 — the mapping the pipeline actually used).
Threshold N_a is applied on **N.Marker.Test** (alleles *tested*), matching the paper's `N_total` row.

---

## A. Table 1 — REPLACE the PMR and N_kg rows (known-gene list now 17 genes)

**Decision applied:** let-23 (ZK1067.1) and npr-1 (C39E6.6) are added to the known-gene list
(15 → **17 genes**), since they are established sleep/behavior regulators. This changes PMR and
N_kg (both genes rank low, so PMR rises), and updates the Table 1 legend to "17 known genes."

**Old (paper, 15 genes):**

| N_a | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| PMR | 19.26 | 15.85 | 13.06 | 13.46 | 12.76 | 16.97 | 20.10 | 20.50 | 22.44 |
| N_kg | 15 | 13 | 11 | 9 | 7 | 5 | 4 | 3 | 3 |
| N_total | 18070 | 14902 | 11736 | 8835 | 6663 | 4959 | 3745 | 2797 | 2161 |

**New (corrected, 17 genes — use these):**

| N_a | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| **PMR** | **24.90** | **23.87** | **20.74** | **21.35** | **24.25** | **24.26** | **35.65** | **35.17** | **25.30** |
| **N_kg** | **17** | **14** | **10** | **8** | **7** | **7** | **6** | **6** | **5** |
| **N_total** | **18070** | **14903** | **11737** | **8836** | **6663** | **4959** | **3745** | **2798** | **2162** |

- **PMR row**: fully replaced. Higher than the 15-gene version because let-23 (88.6%) and
  npr-1 (74.4%) both rank low; the median still stays ~21–25%, well under the 50% null.
- **N_kg row**: replaced (now starts at 17).
- **N_total row**: unchanged from the data (a few cells differ by 1 from the old paper values —
  a tie/boundary detail; use these).
- **Optimal N_a**: still **3** (min PMR = 20.74%).
- **Legend**: change "known genes" count to **17**; add let-23 and npr-1 (with citations) to the
  list of named known genes in the results text.

---

## B. Optimal N_a — CHANGE 5 → 3

- **Old:** "The optimized value for N_a = 5."
- **DECISION (committed):** move the optimized threshold to **N_a = 3** (min PMR = 20.74% over
  the 17-gene list). Change "N_a = 5" → "N_a = 3" everywhere.
- **Cascade — the headline gene count changes: 6,663 → 11,737.** Both prose mentions of "6,663
  genes … at least five non-synonymous alleles" become **"11,737 genes … at least three
  non-synonymous alleles."** Two places to edit:
  1. **Abstract (lines ~17-18):** "…to prioritize **6,663** genes…" → "…to prioritize **11,737**
     genes…"
  2. **Results (lines ~84-85):** "a ranking for **6,663** genes that have **at least five**
     non-synonymous alleles tested" → "a ranking for **11,737** genes that have **at least three**
     non-synonymous alleles tested"
  - (The number **6,663 still appears once**, correctly, inside the Table 1 N_total row at the
    N_a = 5 column — leave that as-is.)

---

## C. Validation / "53%" claim — REWORD (now over the 17-gene list)

- **Old:** "all the known sleep regulators rank among the first 53% of the list."
- **Fact (17 genes, correct SetIDs, N_a = 1):** with let-23 and npr-1 now included, **the "first
  53%" claim is false** — 4 of 17 exceed 53%: siss-1 (65.6%), grk-2 (66.1%), npr-1 (74.4%),
  let-23 (88.6%). But the majority still rank high: 12 of 17 within the top ~25%.
- **New wording (recommended):** drop "all … within 53%" entirely and report the enrichment
  statistic instead: *"The known regulators are significantly enriched toward the top of the
  ranked list (median percentile 20.7% at the optimized threshold N_a = 3; permutation p = [recompute]
  vs. a 50% null). A minority — including let-23 and npr-1 — rank low; these are expected misses
  for an MMP/SKAT approach (see below)."*

## C3. The two specific sentences you asked about — recomputed values

Both sentences change. Here are the actual recomputed numbers (same permutation test the paper
describes: median percentile of the known-gene ranks vs. a 50% null, 10,000 trials).

### Sentence 1 — the p-value (line ~254-255)
> Old: "…median percentile … equals 50% … we obtain a **p value of 0.0034** (N = 10,000
> simulation trials)…"

**New (17-gene list):**

| N_a | median percentile | **p (enrichment, one-sided)** |
|---|---|---|
| 1 (full list) | 24.9% | **0.012** |
| **3 (optimum)** | **20.7%** | **0.012** |
| 5 | 24.3% | 0.059 *(not significant)* |

→ Replace **0.0034** with **p ≈ 0.012** (at N_a = 1 or the N_a = 3 optimum). Still supports
enrichment. **Note:** at N_a = 5 the test is no longer significant (p ≈ 0.06) — another reason
not to keep N_a = 5 as the headline threshold.
*(For reference, the 15-gene list without let-23/npr-1 gives p ≈ 0.009 at N_a = 1, 0.016 at
N_a = 3. Adding the two low-ranked genes weakens significance only slightly.)*

### Sentence 2 — the N_a value and "top 1,000" (line ~260-261)
> Old: "The optimized value for **N_a = 5**. With N_a = 5, all the known genes cluster within the
> **top 1,000** in our list."

**New:**
- Optimized value → **N_a = 3** (min PMR = 20.74%).
- "top 1,000" → **false and must be removed/replaced.** Worst known-gene rank is 10,425 (of
  11,737) at N_a = 3, and 5,924 (of 6,663) even at N_a = 5. There is no threshold at which all
  known genes fall in the top 1,000.
- Suggested replacement: *"The optimized value is N_a = 3, at which the known regulators rank
  with a median percentile of 20.7% (permutation p ≈ 0.012)."* — i.e. drop the "top 1,000"
  claim and cite the median percentile instead.

---

## C2. let-23 and npr-1 — now IN the table + a response-letter note

Per your decision, these are **added to the known-gene list** (reflected in the 17-gene Table 1
above). Verified against the deposited data — **the reviewer's numbers are exactly right**:

| gene | SetID | P | tested alleles | rank / 19,749 | percentile |
|---|---|---|---|---|---|
| let-23 / EGFR | ZK1067.1 | 0.888 | 10 | **16,016** | 88.6% |
| npr-1 | C39E6.6 | 0.739 | 2 | **13,445** | 74.4% |

**Key point:** these two genes are **NOT in the paper's 15-gene validation list** (Table 1) —
the reviewer pulled them independently. So they do **not** change Table 1, PMR, or N_kg. They
need a **direct answer in the response letter**, framed honestly:

- **let-23** has 10 tested alleles (well-powered) yet ranks at 88.6% → this is a genuine
  *absence of association signal*, not low power. Biologically defensible: let-23/EGFR is an
  essential gene; loss-of-function is lethal, so the MMP carries only mild/hypomorphic coding
  variants that need not perturb the EGF→SIS signaling relevant here. SKAT tests random coding
  alleles, which won't recapitulate the specific signaling change.
- **npr-1** has only 2 tested alleles → underpowered, and its known behavioral effect is driven
  by a *specific natural variant* (215F/215V), not random coding mutations — so a low rank is
  expected and uninformative.

**Suggested response wording:** *"We agree that let-23 and npr-1 rank low. Neither was among the
literature-curated regulators we used for validation (Table 1). Both are special cases for an
MMP/SKAT approach: let-23/EGFR is essential, so the panel contains only viable hypomorphic
alleles, and npr-1's characterized effect derives from a specific natural polymorphism rather
than the random coding alleles SKAT aggregates. This illustrates a general limitation — SKAT on
a loss-of-function-biased panel is expected to miss genes whose SIS role depends on gain-of-
function or specific non-null alleles — which we now state explicitly."*

---

## D. Numbers that STAY (verified correct — do NOT change)

| Number | Where | Status |
|---|---|---|
| **19,749** genes / ranked list | abstract, results, methods | ✓ correct |
| **18,070** testable | (add this clarification) | ✓ |
| **96%** of genome | methods | ✓ |
| **4.24 ± 4.63** alleles/gene | methods | ✓ exact (mean±SD of N.Marker.Test over 19,749) |
| **941 = 189 robot + 752 manual** | methods | ✓ arithmetic |
| **6,663** (= ≥5 tested alleles) | abstract/results | ✓ correct *as a count* — but see B re: whether it stays the headline |

---

## E. SetIDs — for Table 1 / validation supplement (verified, from WS225)

ceh-17 → **D1007.1** · ceh-14 → **F46C8.5** · aptf-1 → **K06A1.1** · kin-29 → **F58H12.1** ·
hda-4 → **C10E2.3** · rom-4 → **Y116A8C.16** · adm-4 → **ZK154.7** · frm-10 → **F25H9.5** ·
siss-1 → **F45D3.2** · unc-108 → **F53F10.4** · eel-1 → **Y67D8C.5** · goa-1 → **C26C6.2** ·
npr-4 → **C16D6.2** · flp-11 → **K02G10.4** · grk-2 → **W02B3.2**

Full per-gene ranks/percentiles: `~/SKAT_MMP/outputs/validation_genes.tsv` (updated).

---

## F. Author-held items still needed (not derivable from the SKAT data)

These are unaffected by the SetID work but must be settled for the revision:

- Candidate screen counts: **83 candidate genes / 106 strains / 3 validated (strd-1, cla-1,
  egl-8)** and the abstract's **"19 additional genes"** — update if new genes were discovered.
- Per-strain robot/manual assignment (Reviewer 2 #2) and baseline/post-UV scores (Reviewer 2 #1).

---

## Decision needed from you

1. **Headline threshold:** keep 6,663 (≥5) as the abstract/results figure, or switch to the new
   optimum N_a = 3 → 11,737 genes? (Affects B and every dependent mention.)
2. **"53%" rewording:** which option in C.
3. Then I finalize Table 1, recompute the enrichment p-value, and (if you want) apply the edits
   into a copy of the manuscript.
