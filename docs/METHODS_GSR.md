# Methods — Gene-Set Randomization (GSR) Adjustment

This is the core method of the pipeline. Standard pathway-enrichment tools
(MAGMA, PRSet, PascalX, GSA-MiXeR) report a per-pathway statistic, but those
statistics are confounded by **pathway size** and **pathway overlap** (genes
shared across many pathways, and hub genes that appear in many sets). The GSR
adjustment recalibrates each pathway's result against an empirical null built
from randomized gene sets that preserve those nuisance properties.

Implementation: [`../scripts/core/calc_empirical.r`](../scripts/core/calc_empirical.r),
run by the `calc_empirical_pvalues` Nextflow process.

---

## 1. Build the null gene sets

For the real gene-set collection (MSigDB C2), generate `N` randomized collections
(study default `N = 1000`) under two null models:

- **BiReWire** — degree-preserving rewiring of the gene × pathway bipartite graph.
  Preserves **both** each pathway's size **and** each gene's membership frequency
  (how many pathways it belongs to). This controls for hub genes / pathway overlap.
- **KeepPathSize** — randomly re-samples genes for each pathway, preserving **only**
  pathway size. This controls for size alone.

Comparing the two isolates the contribution of gene-overlap structure beyond size.

## 2. Run the same enrichment tool on real and null sets

Each tool is run once on the real GMT and once per randomized GMT, producing, for
every pathway, the tool's native statistic:

| Tool | p-value column | effect column | notes |
|------|----------------|---------------|-------|
| MAGMA | `P` | `BETA` | competitive gene-set test |
| PRSet | `P` | `Coefficient` | also keeps `Competitive.P` |
| PascalX | `Pvalue` | — | p-value only, no effect size |
| GSA-MiXeR | — | `enrich` | enrichment only; selected via `loglike_aic` |

Tool→column mapping is the `tool_config` list at the top of `calc_empirical.r`.

## 3. Compute the two adjusted statistics

For each pathway, let the real statistic be compared against that **same pathway's**
distribution across the `N` null replicates.

### Empirical p-value (rank-based)
```
empirical_pval = (n_more_extreme + 1) / (n_total_random + 1)
```
where `n_more_extreme` = number of null replicates of that pathway whose p-value is
≤ the real p-value. The `+1` numerator/denominator is the standard Davison–Hinkley
finite-sampling correction (no pathway can reach an empirical p of exactly 0).

- Computed for MAGMA, PRSet, PascalX (`calc_pvalue = TRUE`).
- **Not** computed for GSA-MiXeR (no native p-value); it uses the effect-size
  statistic below only.

### Standardized effect size (z-like)
```
std_effect_size = (beta_real - mean(beta_null)) / sd(beta_null)
```
i.e. how many null standard deviations the real effect sits above the null mean for
that pathway. `NA` when no effect column exists (PascalX) or when `sd(beta_null) = 0`.

## 4. Why use these instead of the raw tool output

- The raw `p_value` is asymptotic and can be **anti-conservative** for large or
  highly-overlapping pathways.
- `empirical_pval` and `std_effect_size` are calibrated against pathways with matched
  size/overlap, so cross-pathway ranking is fair.
- **Use `empirical_pval` (and `std_effect_size`) for ranking and significance — not
  the raw `p_value`,** which is retained only for reference.

## 5. Output

One file per trait × tool × randomization method:
```
{trait}_{tool}_{birewire|keeppathsize}_empirical_pvalues.txt
```
Key columns: `pathway_name`, `empirical_pval`, `std_effect_size`,
`mean_beta_perm`, `sd_beta_perm`, `n_more_extreme`, `n_total_random`, plus the raw
`p_value`/effect for reference. These files are the input to every validation
analysis (see [VALIDATION.md](VALIDATION.md)).

## 6. Inputs the script expects (positional args)

```
calc_empirical.r <trait> <tool_method> <real_results_file> <random_dir> [gmt_path]
```
- `<tool_method>` e.g. `magma_birewire`; the base tool is parsed before the `_`.
- `<random_dir>` directory holding the `N` null result files for that trait/method.
- `[gmt_path]` **required for GSA-MiXeR** (restricts results to pathways defined in
  the GMT); ignored for other tools.
