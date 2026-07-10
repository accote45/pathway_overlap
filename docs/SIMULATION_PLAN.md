# Simulation Validation Plan

A build specification for a small simulation pipeline that validates the
gene-set randomization (GSR) adjustment on data with **known pathway
architecture**. This document is the reference we will script from.

---

## 1. Purpose

Reviewer comment:

> "Some simulation validation would be both fairly conventional and improve the
> characterization of the method across a well-defined set of pathways (with
> 'known' architecture in the data)."

We will build pathways where we *know the truth* — which pathways are genuinely
enriched, which are not, and how they overlap — then show the method behaves
correctly. Because we control the architecture, we can measure what real traits
never let us measure: false-positive rate, power, calibration, and whether the
method's ranking matches the true ranking.

## 2. What the method claims (what the simulation must prove)

The contribution is the GSR adjustment, specifically that:

- **BiReWire** (null preserving pathway size **and** each gene's pathway-
  membership count) is better calibrated than
- **KeepPathSize** (null preserving size only), and than
- the **raw asymptotic MAGMA p-value**,

**because BiReWire controls for gene overlap between pathways** (shared genes /
hub genes appearing in many sets).

The simulation must demonstrate, on known-truth data:

1. Gene overlap creates *false* enrichment in innocent ("bystander") pathways.
2. Raw-p and KeepPathSize are fooled by it; BiReWire is not.
3. BiReWire still detects genuinely enriched pathways (it does not win merely by
   being conservative).

**Clean structural prediction we expect to confirm:** with *no* overlap, every
gene sits in exactly one pathway, so BiReWire and KeepPathSize are nearly
identical and both behave well. With shared hubs, hub genes have membership
count > 1, the two nulls diverge, and BiReWire is the one that stays calibrated.
That divergence *is* the result.

## 3. Design in one paragraph

Simulate a **quantitative trait** on a real UKB genotype subset. Place causal
variants inside the genes of pathways we designate as "enriched," with per-gene
effect sizes graded across strength tiers. Run a standard GWAS, then run the
**existing MAGMA + GSR pipeline unchanged**. Because we know which pathways carry
signal and which only *borrow* it through shared genes, every performance metric
has a ground truth. MAGMA only; we argue the GSR adjustment is tool-agnostic
because it wraps any per-pathway statistic.

## 4. Key efficiency insight (drives the whole pipeline shape)

Pathway **overlap is purely a gene-membership (GMT) construct.** The genotypes,
the simulated phenotype, the GWAS, and the per-gene Z-scores are **identical**
whether pathways overlap or not — only *which pathways a gene is listed in*
changes. Therefore:

- Simulate phenotype **once**, run GWAS **once per replicate**, run the MAGMA
  **gene step once per replicate.**
- Reuse those gene-level results for **both** overlap settings (disjoint and
  shared-hub); only the cheap **gene-set step + null comparison + scoring**
  differ.

So the expensive compute scales with **replicates only**, not replicates ×
overlap settings.

---

## 5. Inputs and prerequisites

| Input | Description | Source |
|---|---|---|
| Genotype subset | ~30,000 unrelated UKB individuals, QC'd, PLINK bfile | UKB; `scripts/ukb_data_preparation/` |
| Gene location file | MAGMA `--gene-loc` (gene, chr, start, end) | existing pipeline input |
| PCs | Top ~10 genotype PCs for GWAS covariates | compute on subset |
| Tools | `gcta64` (phenotype sim), `plink2` (GWAS), `magma` (1.10), `R` (BiReWire, tidyverse, data.table, GSA) | cluster modules |

Reused pipeline components (unchanged):

- `modules/magma/magma.nf` — gene analysis + gene-set steps.
- `scripts/core/generate_birewire_gmts.R`, `generate_keeppathsize_gmts.R` — null
  gene-set generation.
- `scripts/core/calc_empirical.r` — empirical p + standardized effect size.
- `scripts/core/calc_fpr.r` — false-positive rate helper.

New components (small; Section 11):

1. Architecture builder (R).
2. Causal-loci / effect builder (R).
3. Scoring + figures (R).
4. A small Nextflow workflow (`sim.nf`) tying it together.

---

## 6. The synthetic pathway architecture (Generator 1)

**Collection:** 100 pathways × 50 genes. Genes drawn from the real gene-loc file
(so real gene sizes, positions, and LD come along — this keeps MAGMA's size/LD
confounds realistic, which is exactly what the method corrects).

**Gene pool split:**
- **Signal genes** — carry causal variants; assigned a per-gene effect tier.
- **Background genes** — no causal variants; fill pathways to size 50.

**Strength tiers (the ground-truth ranking).** Partition the 100 pathways into
4 tiers with per-gene effect multipliers 0 : 1× : 2× : 3× (matching the
discussion notes "second twice as high, third three times as high"):

| Tier | # pathways (default) | Composition | True enrichment |
|---|---|---|---|
| T0 null | 40 | 50 background genes | none |
| T1 low | 20 | k signal genes (1× effect) + (50−k) background | weak |
| T2 med | 20 | k signal genes (2× effect) + (50−k) background | moderate |
| T3 high | 20 | k signal genes (3× effect) + (50−k) background | strong |

Default `k = 10` signal genes per enriched pathway (tunable). Each pathway gets a
**true enrichment score** = sum of its genes' per-gene heritability
contributions; this continuous score is what we rank-correlate against the
method's output (rank-recovery metric).

**Two overlap settings (same enriched/null structure in both):**

- **Disjoint (0% overlap).** Every one of the 5,000 gene slots is a unique gene;
  no gene appears in two pathways. Membership count = 1 for all genes.
- **Shared-hub (~5% overlap, adversarial).** Convert a block of the T0 null
  pathways into **bystanders**: replace some of their background genes with
  **signal genes borrowed from enriched (T2/T3) pathways** (the "hubs"). A
  bystander thus shares high-signal genes with a truly-enriched pathway but has
  **no independent signal of its own.** Ground truth: **bystanders are null.**

Overlap fraction is set so shared genes ≈ 5% of memberships (tunable). The hub
genes now have membership count > 1 — the case where BiReWire and KeepPathSize
diverge.

**Outputs of Generator 1:**
- `arch_disjoint.gmt`, `arch_sharedhub.gmt` — the two pathway collections.
- `causal_genes.tsv` — signal genes with their per-gene effect tier (shared by
  both settings; identical signal).
- `manifest.tsv` — per pathway: `pathway_id, tier, is_enriched, is_bystander,
  n_signal_genes, n_shared_hub_genes, true_enrichment_score`.

---

## 7. Causal variants and phenotype (Generator 2 + GCTA)

**Pick causal SNPs.** For each signal gene, select causal SNPs from the bfile
within the gene body ± window (default 35 kb, matching the MAGMA annotation
window). Default: a fixed number of causal SNPs per gene (e.g. 1–5), or all SNPs
in the gene — a parameter.

**Assign effect weights.** Each causal SNP gets an effect weight equal to its
pathway's tier multiplier (1× / 2× / 3×). To keep the ground truth tied to
**genes, not gene size**, normalize so each gene contributes the same target h²
within its tier regardless of how many SNPs it has (divide weight by the gene's
causal-SNP count). *(Alternative, as a robustness check: skip normalization so
larger genes carry more signal — this deliberately stresses the size confound.
Parameter `size_confound = off|on`.)*

**Write GCTA causal-loci file** `causal_loci.txt`: two columns `SNP_id effect`.

**Simulate phenotype** with GCTA (fresh noise per replicate, all replicates in
one call):

```
gcta64 --bfile <subset> \
       --simu-qt \
       --simu-causal-loci causal_loci.txt \
       --simu-hsq <H2>          # default 0.4 total trait heritability
       --simu-rep <R>           # default 20 (tunable to 50)
       --out sim
```

Produces `sim.phen` with one column per replicate.

---

## 8. GWAS (per replicate)

Standard linear association on the same subset (which also serves as MAGMA's LD
reference):

```
plink2 --bfile <subset> \
       --glm hide-covar \
       --pheno sim.phen --pheno-col-nums <rep> \
       --covar pcs.txt --covar-col-nums 1-10 \
       --out gwas_rep<rep>
```

Emit a MAGMA-ready summary file with columns `SNP  P` and a fixed sample size
`N = <subset size>` passed to MAGMA as `N=<n>`.

## 9. MAGMA + GSR (reuse existing pipeline)

Per replicate (shared across both overlap settings):

1. **Gene analysis** — `magma --bfile <subset> --pval gwas use=SNP,P N=<n>
   --gene-annot <annot>` → `.genes.raw` (per-gene Z + gene-gene LD).

Then per replicate **× overlap setting**:

2. **Null gene sets** — run `generate_birewire_gmts.R` and
   `generate_keeppathsize_gmts.R` on that setting's architecture GMT to produce
   `N_null` (default 1000) randomized collections. *(These depend only on the
   GMT, not on the phenotype, so they can be generated once per overlap setting
   and reused across replicates.)*
3. **Gene-set test** — MAGMA gene-set analysis on the real architecture GMT and
   on each null GMT.
4. **Empirical stats** — `calc_empirical.r` → per pathway `empirical_pval`,
   `std_effect_size` under BiReWire and KeepPathSize; keep raw MAGMA `P` for the
   third comparison arm.

---

## 10. Scoring and outputs

For each **replicate × overlap setting × method** (`birewire`, `keeppathsize`,
`raw_p`), join results to `manifest.tsv` and compute:

| Metric | Definition | Answers |
|---|---|---|
| **Bystander FPR** | fraction of bystander pathways with significant enrichment (e.g. emp-p < 0.05, enriched direction) | the headline: does overlap cause false positives? |
| **Null FPR / calibration** | fraction of true-null (T0, non-bystander) pathways flagged; QQ / uniformity of emp-p | Type-I error control |
| **Power / TPR** | fraction of truly-enriched pathways flagged, by tier | sensitivity, and that BiReWire isn't just conservative |
| **Rank recovery** | Spearman(`true_enrichment_score`, `std_effect_size`) across all pathways | mirrors the real-data OpenTargets validation, with known truth |
| **AUC** | enriched vs null discrimination using emp-p / std-effect as score | overall separation |

Aggregate across replicates (mean ± CI). **Headline figures:**

1. **Bystander FPR vs overlap setting, by method** — expect all near-nominal
   under disjoint; under shared-hub, raw-p and KeepPathSize inflate, BiReWire
   stays ≈ 5%. *(The money figure.)*
2. **Power vs tier, by method** — BiReWire tracks the others (not conservative).
3. **Emp-p QQ under the null, by method** — calibration.
4. **Rank recovery + AUC** — bar/table by method × overlap.

Because the whole simulation has clean ground truth, this **replaces the separate
"calibration cell"** from the earlier draft — FPR and calibration now fall
directly out of the null pathways.

## 11. Experiment matrix

Expensive axis (phenotype → GWAS → gene analysis): **replicates only.**

| Axis | Levels | Where it acts |
|---|---|---|
| Replicates | 20 (default), up to 50 | phenotype noise → full recompute |
| Overlap | disjoint, shared-hub | GMT + gene-set step only (cheap) |
| Strength tier | T0/T1/T2/T3 | inside one phenotype (per-gene effects) |
| Effect size (optional 2nd sim) | small / med / large total-h² | rerun GCTA with different `--simu-hsq` |
| Pathway size (optional, last) | 50 vs larger | alternate architecture GMT |

Default first pass: 20 replicates × 2 overlap settings, one h². Add effect-size
and size sweeps only if the first pass looks good.

---

## 12. Proposed repository layout

```
scripts/simulation/
  build_architecture.R      # Generator 1: GMTs + manifest + causal_genes
  build_causal_loci.R       # Generator 2: causal_loci.txt for GCTA
  score_simulation.R        # metrics + figures
sim.nf                      # small Nextflow workflow (reuses modules/magma, core scripts)
docs/SIMULATION_PLAN.md     # this document
```

**`sim.nf` process graph:**

```
build_architecture ─┬─> arch_disjoint.gmt ──> null_sets(disjoint)
                    ├─> arch_sharedhub.gmt ─> null_sets(sharedhub)
                    ├─> causal_genes.tsv ──> build_causal_loci ──> causal_loci.txt
                    └─> manifest.tsv ─────────────────────────────┐
                                                                  │
causal_loci ──> simulate_phenotype (GCTA, R reps) ──> per rep:   │
                 gwas ──> magma_gene_analysis ──┐                 │
                                               (fan out × overlap × method)
                 magma_geneset(real+null) ──> calc_empirical ──> score_simulation
                                                                  ▲
                                                    manifest ─────┘
```

## 13. Parameters (defaults)

| Param | Default | Notes |
|---|---|---|
| `n_pathways` | 100 | |
| `genes_per_pathway` | 50 | |
| `tier_sizes` | 40/20/20/20 (T0/T1/T2/T3) | |
| `tier_effect_mult` | 0/1/2/3 | |
| `signal_genes_per_pathway` (k) | 10 | |
| `overlap_fraction` | 0.05 | shared-hub setting |
| `n_bystander_pathways` | e.g. 10–20 (from T0 block) | |
| `causal_snps_per_gene` | 1–5 (param) | |
| `annotation_window_kb` | 35 | matches MAGMA |
| `trait_h2` | 0.4 | GCTA `--simu-hsq` |
| `n_replicates` (R) | 20 | up to 50 |
| `n_subjects` | ~30,000 unrelated | |
| `n_pcs` | 10 | GWAS covariates |
| `n_null_sets` | 1000 | BiReWire / KeepPathSize |
| `size_confound` | off | gene-size normalization toggle |

## 14. Realism guardrails / caveats to state in the paper

- Real genotypes, real LD, real gene boundaries and sizes — nothing structurally
  synthetic about the data-generating genome.
- Quantitative trait (not case/control): standard for method validation, maximal
  power per N. No trait ascertainment needed; every subject gets a simulated
  phenotype.
- Include PCs as GWAS covariates to mirror real analyses.
- Assign genes to pathways without regard to genomic location, and check that
  same-pathway and cross-pathway genes are not concentrated in shared LD blocks
  (avoid spurious cross-pathway correlation); MAGMA's gene-gene correlation
  handles residual LD.
- Uniform ~5% shared-hub overlap is simpler than real MSigDB's heavy-tailed
  overlap; we acknowledge this and keep two overlap settings for speed (a
  heavy-tailed scenario is a possible extension, not in the first pass).

## 15. How this complements the existing paper

The paper's current FPR — randomizing genes while retaining real GWAS
associations — stays as the **real-data** false-positive result. This simulation
is the **independent, known-architecture** complement: its ground truth
(causal variants) is generated separately from the gene-shuffling null, so it is
not circular, and it exercises the full SNP → gene → pathway path.

## 16. Build order (milestones)

1. **`build_architecture.R`** — GMTs + manifest. Zero heavy compute; unblocks
   everything and lets us eyeball the design. *(first deliverable)*
2. **`build_causal_loci.R`** + a single GCTA/GWAS/MAGMA-gene run on one replicate
   — confirm signal flows end-to-end and gene Z-scores track the tiers.
3. **Wire `sim.nf`** with the reused MAGMA + null + empirical steps; run the
   disjoint setting, few replicates — sanity check (BiReWire ≈ KeepPathSize).
4. **Add shared-hub setting** — confirm the divergence and the bystander FPR gap.
5. **`score_simulation.R`** — metrics + the four headline figures.
6. Scale replicates to 20–50; optional effect-size / pathway-size sweeps.

## 17. Open decisions before scripting

- Confirm `n_subjects` (~30k) and that an unrelated QC'd UKB subset + PCs are
  available, or whether we generate them first.
- Confirm `causal_snps_per_gene` policy (fixed count vs all SNPs in gene).
- Confirm number of bystander pathways and the exact `overlap_fraction`.
- Confirm replicate count for the first pass (20 proposed).
</content>
