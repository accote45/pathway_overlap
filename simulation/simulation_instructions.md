# GSR simulation analysis — build & run instructions (Tier 2, MAGMA only)

**For:** Claude Code
**Goal:** Run a ground-truth simulation that demonstrates the Gene Swap Randomization (GSR) null improves GWAS pathway prioritization *accuracy* — not just concordance with imperfect external resources. This addresses Reviewer 1's request for simulation validation (R1 Comment 3) and the underlying accuracy question (R1 Comment 2).

Read this whole document before writing code. Where it says "reuse the existing pipeline," inspect the repo first and reuse rather than reimplement.

---

## 0. The one claim this simulation must establish

> Under a pathway architecture where the ground truth is known, standard and pathway-size-preserving (PS) enrichment falsely flag null pathways that are enriched **only** because they share signal-carrying multi-pathway ("hub") genes with truly-enriched pathways. GSR removes these false positives while retaining the truly-enriched pathways, and therefore ranks pathways closer to the planted truth.

Every output should serve this claim. The critical comparison is **"contaminated null" pathways** (null biology, but containing signal-carrying hub genes) across the three methods.

This is the controlled, ground-truth analog of Figures 2 and 5 in the manuscript.

---

## 1. Assets to locate first (do not rebuild these)

Inspect the existing repository and environment and record the exact paths in a `config.yaml` before doing anything else:

1. **Existing GSR pipeline:** `https://github.com/accote45/pathway_overlap` (Nextflow). Identify:
   - the step that generates the 1,000 GSR-randomized pathway databases (BiRewire),
   - the step that generates the 1,000 PS (pathway-size-only) randomized databases,
   - the MAGMA gene-set enrichment step,
   - the step that computes **empirical p-values** and **standardized effect sizes (SES)** from the null databases.
   You will feed simulated inputs into these existing steps. Do not reimplement BiRewire, MAGMA calls, or the empirical-p/SES math.
2. **MAGMA binary** (v1.10, matching the paper).
3. **Reference panel** used in the paper: the 1000 Genomes European PLINK fileset (`g1000_eur`). Used both for MAGMA gene analysis and for LD when simulating summary statistics.
4. **MAGMA gene location file** and **SNP location file** for the build matching the reference panel. The real gene IDs and coordinates in this file are what make the simulation realistic.
5. **BiRewire R package** (already a dependency of the pipeline).

If any asset is missing, stop and report it rather than substituting.

---

## 2. Pathway architecture (the simulated gene-by-pathway matrix)

Build a synthetic pathway database with fully known structure. Use **real gene IDs** sampled from the MAGMA gene-location file so that SNPs map to genes and real LD/gene-gene correlation is preserved.

**Fixed design (defaults):**
- **100 pathways, 50 genes each** → 5,000 gene-slots.
- **30 truly-enriched pathways**, in 3 effect tiers of 10 each: **low / medium / high**.
- **70 truly-null pathways** (no planted biology of their own).

**Overlap (multi-pathway / "hub" genes).** Overlap is the fraction of gene-slots occupied by hub genes that sit in multiple pathways. Core conditions:
- **0% overlap:** every slot is a unique gene (5,000 distinct genes, each in exactly 1 pathway). No hubs. This is the essential negative control — GSR and PS **must be near-identical here**.
- **5% overlap:** 250 hub-slots. Default hub design: **25 hub genes, each placed in 10 pathways** (25 × 10 = 250). Remaining 4,750 slots are unique genes. Total distinct genes ≈ 4,775.
- Optional extra conditions if time permits: **10% and 20%** overlap (scale hub count accordingly).

**Where hubs go (this creates the contamination):** place each hub gene into **at least one high-tier enriched pathway** and **several truly-null pathways**. A null pathway that receives ≥1 signal-carrying hub becomes a **contaminated null**. This yields three pathway classes for analysis:
- **Truly enriched** (30) → should be flagged (measures power).
- **Contaminated null** (subset of the 70 nulls containing ≥1 signal hub) → should **not** be flagged (the key FPR test).
- **Clean null** (nulls with no hub genes) → baseline FPR.

**Realism guards (per the manuscript's own caution about not doing anything that wouldn't occur in real data):**
- Keep hub multiplicity modest (in ~5–10% of pathways), proportionally comparable to real MSigDB multi-pathway genes rather than extreme.
- Sample unique genes spread across the genome (random draw from the gene-loc file is fine); avoid accidental genomic clustering.
- Optionally run one robustness condition with **variable pathway sizes** drawn to resemble the MSigDB size distribution (right-skewed, median ~37) instead of fixed 50, to mirror Figure 4's size-dependence. Keep fixed-50 as the primary design for clarity.

**Output format:** write the simulated pathway database in the **exact format the existing pipeline's MSigDB matrix uses** (same gene-set file layout MAGMA `--set-annot` / the pipeline expects), so the GSR and PS randomization steps run unchanged. Also save a **ground-truth table**: `pathway_id, class (enriched_low/enriched_med/enriched_high/contaminated_null/clean_null), n_signal_hubs, planted_effect`.

---

## 3. Signal model — per-tier effect-size shift

Each **gene** has a single association status (a gene is one gene, even if it sits in many pathways).

- **Null genes:** no signal.
- **Associated genes:** genes belonging to an enriched pathway, plus the signal-carrying hub genes. Each associated gene gets **1 causal SNP** (its lead/most-central SNP in its ±35 kb window) assigned a mean Z-shift **μ** determined by tier.

**Per-tier mean shift (μ), at the causal SNP, in Z units:**

| Tier | Role | μ (medium condition) |
|---|---|---|
| low | enriched low pathways | μ_low |
| medium | enriched medium pathways | μ_med |
| high | enriched high pathways | μ_high |
| hub | signal-carrying multi-pathway genes | μ_hub (set = μ_high) |

**Effect-size axis (the small / medium / large conditions):** scale all μ by a global factor `s ∈ {s_small, s_med, s_large}`. Do **not** hard-code magic numbers — instead **calibrate** the μ's once (in the 0%-overlap condition) to hit these operating points, then reuse them across all conditions:
- clean-null FPR ≈ 0.05,
- high-tier power ≈ 0.8,
- low-tier power ≈ 0.2.

This guarantees results are interpretable and avoids saturation (all-significant) or floor (all-null) regimes. Report the calibrated μ values in the log.

The **gradient** low→med→high gives a *known ranking* of true enrichment, which the ranking-accuracy metric (Section 6) uses.

---

## 4. Simulating GWAS summary statistics (MVN-under-LD)

Do **not** simulate genotypes or phenotypes. Simulate per-SNP Z-scores directly, region by region, using real LD from the reference panel.

For each gene region (SNPs within the gene's ±35 kb window, matching the paper):
1. Extract the SNP list and compute the **LD matrix R** for those SNPs from `g1000_eur` (PLINK `--r square`). Regularize to PSD: `R ← R + εI` (small ε, e.g. 1e-3) and/or shrink toward identity if not positive-definite.
2. Build a **mean vector m**: `m[causal SNP] = μ` (0 for all other SNPs in the region; 0 for all SNPs in null genes).
3. Sample `Z ~ MVN(m, R)` via Cholesky (`R = L Lᵀ`, `Z = m + L·z`, `z ~ N(0, I)`).
4. Convert to p-values: `P = 2·Φ(−|Z|)`.

Assemble all regions into one genome-wide summary-statistics file with columns MAGMA needs (`SNP`, `P`), and set a fixed large sample size (e.g. `N = 200000`) held constant across all conditions.

**Fast fallback for a first smoke-test only:** draw null SNP Z's i.i.d. `N(0,1)` (skip the LD covariance) to confirm the plumbing end-to-end, then switch to MVN-under-LD for all reported results. Note in the log which was used — independent draws slightly miscalibrate the gene-level null and should not appear in final figures.

---

## 5. Running MAGMA + GSR (per replicate)

Gene analysis depends only on the summary stats, so **run it once per replicate**; the enrichment step is then repeated over the null databases (the pipeline already does this in bulk).

```
# Step 1 — annotate (once; depends only on gene-loc + snp-loc, so can be cached across reps)
magma --annotate window=35,35 \
      --snp-loc <snp_loc> --gene-loc <gene_loc> \
      --out sim_annot

# Step 2 — gene analysis (once per replicate; depends on simulated sumstats)
magma --bfile <g1000_eur> \
      --pval sim_sumstats.txt N=200000 \
      --gene-annot sim_annot.genes.annot \
      --out sim_genes

# Step 3 — gene-set enrichment (reuse existing pipeline step)
#   run on: (a) the real simulated pathway database  -> original p-values
#           (b) the 1,000 GSR-randomized databases    -> GSR empirical p + SES
#           (c) the 1,000 PS-randomized databases      -> PS empirical p + SES
#   MAGMA takes many sets per --set-annot file, so batch them.
```

Then use the **existing pipeline code** to compute, per pathway:
- **Original** p-value (from the real simulated database),
- **PS** empirical p-value + SES,
- **GSR** empirical p-value + SES.

Keep the number of null databases at **1,000** to match the paper (drop to 500 only for a first pass; note it if so).

---

## 6. Metrics (per replicate, then aggregated over replicates)

For each method ∈ {Original, PS, GSR}:

1. **Contaminated-null FPR (the money metric).** Proportion of *contaminated-null* pathways called significant (empirical p < 0.05), as a function of overlap level. Expectation: Original and PS climb with overlap; GSR stays ≈ 0.05.
2. **Clean-null FPR.** Same for clean nulls — should be ≈ 0.05 for all methods and all overlap levels (calibration check).
3. **Power.** Proportion of truly-enriched pathways called significant, broken out by tier (low/med/high). GSR should **retain** power (not over-correct).
4. **Ranking accuracy vs. planted truth** (this is the direct "accuracy" evidence R1 asked for):
   - **AUPRC / AUROC** for discriminating truly-enriched vs. all-null pathways, per method, per overlap level.
   - **Spearman** between each method's ranking statistic and the planted enrichment magnitude (encode null=0, low<med<high).
   - **Top-K precision** (K = 30, the number of truly-enriched pathways).
   GSR ≥ Original/PS, with the gap widening as overlap increases.
5. **Mechanistic close-up.** Regression of pathway −log10(p) on the **number of signal-carrying hub genes** the pathway contains. Slope should be **positive for Original/PS** and **≈ 0 for GSR**.

**Controls that must be reported:**
- **0%-overlap control:** GSR ≈ PS ≈ Original on every metric (no hubs → nothing to correct).
- **Null-hub control:** rerun the 5% architecture with hub genes assigned **no signal**. Contaminated-null FPR should stay ≈ 0.05 for all methods, and GSR ≈ Original. This isolates that the inflation is specifically due to *signal-carrying* multi-pathway genes (mirrors Fig 5).

**Replicates:** 100 per condition (reduce to 50 for a first pass). Report mean ± 95% CI (or boxplots across replicates).

---

## 7. Figures to produce

1. **Contaminated-null FPR vs. overlap level**, one line per method. *(Main-text candidate — the controlled analog of Fig 2.)*
2. **Power by effect tier**, grouped bar/box per method (shows GSR keeps power).
3. **Ranking accuracy (AUPRC or Spearman-to-truth) vs. overlap**, per method (shows the accuracy gain, widening with overlap).
4. **−log10(p) vs. number of signal-hub genes**, per method, with fitted slope (the mechanism).
5. **Null-hub control panel:** contaminated-null FPR flat across methods.

Save all underlying numbers as tidy CSVs alongside the figures.

---

## 8. Sanity checks (run and log before trusting results)

- **Simulator behaves:** correlation between planted gene-level signal and MAGMA gene-level Z is strongly positive; enriched pathways outrank nulls under the Original statistic in the 0%-overlap condition.
- **Calibration:** clean-null FPR ≈ 0.05 across methods and conditions.
- **Randomization is well-posed on the simulated matrix (ties to R2 Comment 3):** run the BiRewire convergence diagnostic (`birewire.analysis.bipartite`) on the simulated gene-by-pathway matrix and confirm the Jaccard index drops from 1 and plateaus at/ before the analytic swap bound. Confirm `birewire.similarity()` between two randomized databases is low and comparable to original-vs-randomized. This shows the constrained randomization space is non-trivial for the simulated architecture.

---

## 9. Suggested repository layout

```
gsr_sim/
  config.yaml                 # all discovered paths + fixed parameters
  01_build_architecture.R/py  # pathway matrix + ground-truth table (Sec 2)
  02_assign_signal.R/py       # gene effect assignment, per-tier μ (Sec 3)
  03_simulate_sumstats.py     # MVN-under-LD Z simulation (Sec 4)
  04_run_magma.sh             # annotate + gene analysis (Sec 5)
  05_run_enrichment.*         # calls existing pipeline GSR/PS + empirical-p/SES (Sec 5)
  06_metrics.py               # FPR / power / ranking / mechanism (Sec 6)
  07_figures.py               # figures 1–5 (Sec 7)
  run_all.sh                  # loops conditions × replicates, parallelizable
  outputs/
    <condition>/<rep>/...     # intermediate
    summary/                  # aggregated CSVs + figures
```

Parameterize everything (overlap level, effect scale `s`, replicate seed) via `config.yaml` + CLI so conditions can be swept and parallelized. **Set and log a seed per replicate** for reproducibility.

---

## 10. Recommended run order (time-boxed)

Run the **minimum viable cell first**, confirm the money plot appears, then expand:

1. **MVP:** overlap ∈ {0%, 5%}, medium effect size only, plus the **null-hub control**, 50 replicates. This alone yields Figure 1 (FPR vs overlap), the 0%-overlap control, and the null-hub control — enough to state the core claim.
2. **Add** the effect-size axis (small/large) and the ranking-accuracy metric/figure.
3. **Add** 10%/20% overlap and the variable-pathway-size robustness condition only if time remains.

---

## 11. How each output maps back to the reviewers (for the response letter)

- **R1 Comment 3 (run a simulation):** the whole analysis, with known architecture.
- **R1 Comment 2 (is the re-ranking a real accuracy improvement?):** Section 6.4 ranking accuracy vs. planted truth — direct accuracy, not proxy concordance.
- **Manuscript Fig 2 / Fig 5 analog:** Section 6.1 contaminated-null FPR and the null-hub control (inflation only when multi-pathway genes carry signal).
- **R2 Comment 3 (randomization diagnostics):** Section 8 BiRewire diagnostics on the simulated matrix.

Report results factually. If GSR does **not** outperform on some metric/condition, say so — the honest division of labor is that the simulation isolates the overlap mechanism against ground truth, while the real-data analysis covers LD realism, additional tools (PascalX, PRSet, GSA-MiXeR), and biological plausibility.