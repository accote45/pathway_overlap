# GSR simulation sub-pipeline

Ground-truth simulation for the revision: does the Gene Swap Randomization (GSR)
null improve pathway prioritisation **accuracy**, not just concordance with
imperfect external resources? Addresses reviewer R1 comments 2 and 3 and R2
comment 3. Specification: [simulation_instructions.md](simulation_instructions.md).

## It does not touch the main pipeline

This directory is a self-contained Nextflow pipeline with its own `main.nf`,
`nextflow.config` and `modules/`. It **reads** the parent repository's
`scripts/core/*.R` by absolute path and modifies nothing outside `simulation/`.
Nothing was cloned or duplicated.

The three statistical steps the reviewers care about are the study's own code,
called unmodified:

| Step | Code actually executed |
|---|---|
| GSR (BiRewire) null databases | `../scripts/core/generate_birewire_gmts.R` |
| PS (pathway-size) null databases | `../scripts/core/generate_keeppathsize_gmts.R` |
| Empirical p-value + standardised effect size | `../scripts/core/calc_empirical.r` |

## Running it

```bash
cd simulation

./run_all.sh calibrate     # 1. solve for mu (Sec. 3) - REQUIRED FIRST
#            paste results/calibration/calibration_mu.env into nextflow.config
./run_all.sh sweep         # 2. the MVP sweep
```

The sweep refuses to start if `mu_low/mu_med/mu_high` are unset — Sec. 3 forbids
hard-coded effect sizes, so they must come from the calibration pilot.

`./run_all.sh smoke` runs a tiny independent-Z version first if you want to check
the plumbing before committing cluster time. Its results are **not** valid for
figures (no LD).

## What it does

1. **Architecture** (`01`) — 100 pathways × 50 genes over real gene IDs from the
   MAGMA gene-location file. 30 truly enriched (10 low / 10 med / 10 high),
   35 contaminated nulls (contain signal-carrying hub genes), 35 clean nulls.
   Hubs sit in ≥1 high-tier pathway plus several nulls. Redrawn every replicate.
2. **Reference** (`02`) — one-time gene/SNP index and eligible-gene pool.
3. **Summary statistics** (`03`) — `Z ~ MVN(m, R)` per LD block, with `R` read
   straight out of the `g1000_eur` panel. `m` is zero except at each causal
   gene's most-central SNP. No genotypes or phenotypes are simulated.
4. **Enrichment** — MAGMA gene analysis once per replicate, then the real
   database plus 1,000 GSR and 1,000 PS null databases.
5. **Metrics** (`06`) and **figures** (`08`) — Sec. 6 and Sec. 7.

## Two implementation decisions worth knowing

**Batched MAGMA gene-set calls.** The parent pipeline runs one MAGMA job per null
database. At conditions × replicates × 1,000 nulls that is ~10⁵–10⁶ LSF jobs.
Because MAGMA scores each set in `--set-annot` independently, `04` concatenates
100 null databases into one call with prefixed set names and `05` splits the
output back into per-database `.gsa.out` files with the original names — so
`calc_empirical.r` runs unchanged. The `verify_batching` process asserts on real
MAGMA output that batched and standalone results are identical, and fails the run
if they are not.

**Null databases are per-architecture, not per-batch.** They depend only on the
gene-by-pathway matrix, so they are generated once per (condition, replicate) and
reused across all batches.

## Layout

```
main.nf                       SWEEP + CALIBRATE, selected by `--stage`
nextflow.config               all paths and parameters (single source of truth)
conf/extended.config          the 10% / 20% overlap conditions
run_all.sh                    smoke | calibrate | sweep | extended
modules/sim_setup.nf          reference, annotation, architecture, sumstats
modules/sim_enrichment.nf     MAGMA, randomization, empirical p, batch check
modules/sim_analysis.nf       metrics, aggregation, figures, diagnostics
scripts/sim_common.R          shared helpers incl. a dependency-free .bed reader
scripts/01..10_*.R            the numbered stages of the specification
results/
  reference/                  gene pool, SNP loc, annotation
  architectures/<cond>/<rep>/ ground truth + diagnostics
  sumstats/<cond>/<rep>/      per-replicate simulation diagnostics
  empirical/                  calc_empirical.r output
  metrics/                    per-replicate metrics
  summary/                    aggregated CSVs, sanity_checks.txt, figures/
  verification/               batching equivalence report
  birewire_diagnostics/       Sec. 8 randomization diagnostics
  calibration/                the fitted mu values
```

## Nextflow version

Developed and verified against **Nextflow 26.04.3** (the version on Minerva) with
its strict parser. `nextflow lint main.nf modules/ nextflow.config` is clean —
no errors, no deprecation warnings. Also runs on 21.10.

Stage selection uses `--stage calibrate` rather than `-entry`, because the strict
parser requires the entry workflow to be the anonymous one.

## Verification status

Verified off-cluster on a synthetic PLINK panel built for the purpose:

- `.bed` reader reproduces a known genotype matrix exactly, including mid-file seeks.
- Simulator: mean Z at causal SNPs matches the planted μ (standardised deviations
  mean −0.003, sd 1.019 over 136 replicates, t-test p = 0.98); non-causal
  var(Z) = 0.9995; empirical Z correlation vs. panel LD has slope 0.973 and RMSE
  0.083 against a sampling-noise floor of 0.086.
- Batch → split round-trip reproduces MAGMA output exactly, and
  `calc_empirical.r` consumes the split files unmodified.
- AUROC / AUPRC / top-K helpers match hand-computed values including ties.
- Calibration inverts a known power curve to within 1%.
- Full Nextflow DAG runs to completion in stub mode on Nextflow 26 with correct
  task counts, for both `--stage sweep` and `--stage calibrate`;
  `prepare_reference`, `build_architecture` and `simulate_sumstats` also run for
  real end-to-end, including with numeric parameters supplied on the command line
  (which arrive as Strings and must be coerced before use in a range).

**Not yet exercised** (no MAGMA / BiRewire on the development machine): the MAGMA
invocations themselves, and the two parent randomization scripts. Those are the
first things to watch on the cluster; `verify_batching` is the guard for the
riskiest of them.
