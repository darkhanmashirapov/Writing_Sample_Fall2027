# 4. ARDL Bounds Testing: Pre-War Relationships and Long-Run Dynamics

## 4.1 Motivation and Methodology

This section estimates ARDL bounds-testing models for Kazakhstan's and
Russia's corporate credit-risk spreads against domestic liquidity and global
macro-financial determinants, testing for a long-run (cointegrating)
relationship and estimating the corresponding error-correction model (ECM).
This directly addresses H1 (pre-war similarity in spread determinants) and
establishes the baseline against which the Markov-switching (Section 5) and
VAR/IRF (Section 6) analyses test for post-war structural change.

**Why ARDL bounds testing fits this data.** As established in the
stationarity analysis (Section 3), Kazakhstan's credit spreads are integrated
of order 1, I(1), while Russia's are integrated of order 0, I(0) — a mixed
order of integration across the two countries' dependent variables. The ARDL
bounds-testing approach (Pesaran, Shin & Smith, 2001) is specifically
designed for this setting: it accommodates a mix of I(0) and I(1) regressors
without requiring a common order of integration across all series, provided
none are I(2) (already confirmed and resolved in Section 2). This avoids the
need to pre-test for cointegration using methods such as Engle-Granger or
Johansen, both of which require all series to share the same order of
integration — a requirement this dataset does not satisfy.

Each model is estimated in its Unrestricted Error Correction (UECM)
representation, from which the bounds test is computed. Under Case 3
(unrestricted constant, no trend), the test compares an F-statistic against
lower (I(0)) and upper (I(1)) critical value bounds: exceeding the upper
bound supports cointegration; falling below the lower bound supports no
cointegration; falling between the two bounds is inconclusive.

## 4.2 Variable Specification

**Dependent variables.** Credit-risk spreads (corporate yield minus
government yield) at three maturities per country: `kz_spread_1y_m`,
`kz_spread_5y_m`, `kz_spread_10y_m`, and the Russian equivalents.

**Regressor sets.** Three specifications were tested, escalating in
complexity:

| Specification | Kazakhstan regressors | Russia regressors |
|---|---|---|
| A (minimal) | TONIA, Brent oil, VIX | RUONIA, Brent oil, VIX |
| B (mixed) | TONIA, Brent oil, VIX (unchanged) | Full set: RUONIA, oil, VIX, inflation, GDP, external debt/GDP, reserves/M2, volatility |
| C (full) | Full set: TONIA, oil, VIX, inflation, GDP, external debt/GDP, reserves/M2, volatility | Full set (identical to B) |

All models enforce a minimum lag order of 1 for every regressor (a
requirement of the UECM bounds-testing procedure) and use `lags=1` for the
dependent variable's own autoregressive term.

**Sample splits.** Two sample windows were used: (i) pre-war only (before
February 2022), and (ii) the full sample with a binary war dummy
(`war_dummy` = 1 from February 2022 onward) included as an additional
regressor.

**A binding data constraint.** Kazakhstan's government bond yield series
begins only in late 2019, giving a pre-war sample of just 27 observations,
compared to 109 for Russia (whose spread series begins 2013). Even the
minimal specification (4 regressors including the maturity-matched US
Treasury yield, later dropped — see Section 4.3) requires roughly 8
parameters once lag terms are counted; the full specification requires
approximately 18. Against 27 observations, the fuller specifications
approach or exceed the classical rule of thumb of ~10 observations per
estimated parameter. This constraint shapes the interpretation of every
Kazakhstan result in this section and is not a peripheral caveat.

## 4.3 Specification A: Minimal Model, Pre-War

**Table 4.1 — Specification A results, pre-war sample**

| Country | Maturity | n | EC term (p-value) | Bounds F-stat | Verdict |
|---|---|---|---|---|---|
| KZ | 1yr | 27 | -0.072 (0.568) | 2.17 | No cointegration |
| KZ | 5yr | 27 | -0.265 (0.068) | 1.11 | No cointegration |
| KZ | 10yr | 27 | -0.302 (0.024) | 1.99 | No cointegration |
| RU | 1yr | 109 | -0.318 (<0.001) | 6.33 | Cointegration |
| RU | 5yr | 109 | -0.328 (<0.001) | 8.78 | Cointegration |
| RU | 10yr | 109 | -0.328 (<0.001) | 7.63 | Cointegration |

*Note: the US Treasury yield was dropped from the Kazakhstan 1yr model by
automatic lag-order selection in initial testing, plausibly reflecting
collinearity with VIX/oil given the small sample; it was excluded from all
subsequent Specification A runs across maturities for consistency.*

**Findings.** All three Russian maturities show strong, statistically
confident cointegration, with bounds F-statistics ranging from 6.33 to 8.78
and error-correction coefficients tightly clustered around -0.32 to -0.33
(p<0.001 throughout). All three Kazakhstan maturities fail to reject the
null of no cointegration under the joint bounds test.

A subtlety is worth flagging: Kazakhstan's individual error-correction
coefficients approach or reach conventional significance as maturity
increases (10yr: p=0.024, technically significant on its own), yet the joint
bounds F-test still indicates no cointegration even at 10yr. This is not a
contradiction — the bounds test is a joint test across all lagged-level
terms and is the authoritative verdict in the Pesaran-Shin-Smith framework;
a single coefficient's own significance does not override it.

The pattern splits cleanly by *country*, not by maturity: Russia passes at
all three maturities, Kazakhstan fails at all three. A genuine economic
difference in underlying dynamics would more plausibly produce
maturity-dependent variation within each country; a uniform pass/fail split
that tracks sample size (109 vs. 27) almost exactly is the signature expected
if the deciding factor is statistical power rather than a true structural
difference between the two countries' pre-war credit markets.

## 4.4 Specification B: Minimal Kazakhstan, Full Russia, Pre-War

**Table 4.2 — Specification B results, pre-war sample**

| Country | Maturity | n | EC term (p-value) | Bounds F-stat | Verdict |
|---|---|---|---|---|---|
| KZ | 1yr | 27 | -0.072 (0.568) | 2.17 | No cointegration |
| RU | 1yr | 74 | -0.503 (<0.001) | 5.31 | Cointegration |
| KZ | 5yr | 27 | -0.265 (0.068) | 1.11 | No cointegration |
| RU | 5yr | 74 | -0.418 (<0.001) | 2.83 | Inconclusive |
| KZ | 10yr | 27 | -0.302 (0.024) | 1.99 | No cointegration |
| RU | 10yr | 74 | -0.306 (0.008) | 1.51 | No cointegration |

**Findings.** Adding Russia's full regressor set (inflation, GDP, external
debt/GDP, reserves/M2, volatility) reduces its usable pre-war sample from 109
to 74 observations, since several added variables — particularly the
annualized fiscal ratios — have shorter coverage than RUONIA/oil/VIX alone.

Russia's cointegration result proves sensitive to this specification change:
it holds cleanly only at 1yr; the 5yr result weakens to inconclusive, and the
10yr result weakens to no cointegration — the same verdict as Kazakhstan.
This mirrors the same statistical-power mechanism documented for Kazakhstan
throughout this section, now visible in Russia's results once its effective
sample size and parameter count move closer to Kazakhstan's constrained
regime. The strength of Russia's pre-war cointegration finding should
therefore be reported as specification-dependent rather than uniformly
robust; Specification A provides the strongest, most stable evidence across
all three maturities.

## 4.5 Specification C: Full Model, Both Countries, Pre-War

**Table 4.3 — Specification C results, pre-war sample**

| Country | Maturity | n | EC term (p-value) | Bounds F-stat | Verdict |
|---|---|---|---|---|---|
| KZ | 1yr | 27 | -0.669 (0.149) | 1.18 | No cointegration |
| RU | 1yr | 74 | -0.503 (<0.001) | 5.31 | Cointegration |
| KZ | 5yr | 27 | -0.091 (0.804) | 1.50 | No cointegration |
| RU | 5yr | 74 | -0.418 (<0.001) | 2.83 | Inconclusive |
| KZ | 10yr | 27 | -0.018 (0.941) | 2.22 | Inconclusive |
| RU | 10yr | 74 | -0.306 (0.008) | 1.51 | No cointegration |

*Russia's results are identical to Specification B, since its full regressor
set was already in use there.*

**Findings.** Kazakhstan's full specification did not fail outright, but its
estimates become progressively less meaningful as maturity increases: the
error-correction coefficient's p-value rises from 0.149 (1yr) to 0.804 (5yr)
to 0.941 (10yr) — by 10yr the coefficient carries essentially no statistical
information. This is the empirical realization of the small-sample concern
raised in Section 4.2: with 8 regressors and enforced lags estimated on only
27 observations, the model technically converges but produces estimates too
imprecise to interpret, rather than genuine evidence of no relationship.

**Synthesis across Specifications A–C.** Kazakhstan fails to show
cointegration under every specification tested, with results becoming less
reliable — not more informative — as regressors are added, given its fixed
27-observation constraint. Russia shows robust cointegration only under the
minimal specification across all maturities; its result weakens
progressively at longer maturities once the full regressor set is
introduced. Taken together, these results support reporting Specification A
as the primary specification for both countries, with Specifications B and C
retained as documented robustness checks rather than alternative preferred
models.

## 4.6 Full-Sample Estimation with a War Dummy

To test whether the long-run relationship extends across the full sample
once the war's potential effect on the equilibrium level is accounted for,
Specification A was re-estimated on the full 1987–2026 sample (constrained by
each variable's own coverage) with a binary war dummy added as a regressor.

**Table 4.4 — Specification A, full sample with war dummy**

| Country | Maturity | n | EC term (p) | war\_dummy.L1 (p) | Verdict |
|---|---|---|---|---|---|
| KZ | 1yr | 82 | -0.200 (0.001) | -0.033 (0.925) | Cointegration |
| KZ | 5yr | 82 | -0.115 (0.039) | — | No cointegration |
| KZ | 10yr | 82 | -0.149 (0.010) | — | No cointegration |
| RU | 1yr | 164 | -0.334 (<0.001) | -0.074 (0.498) | Cointegration |
| RU | 5yr | 164 | -0.201 (<0.001) | -0.238 (0.075) | Cointegration |
| RU | 10yr | 164 | -0.164 (<0.001) | -0.256 (0.088) | Inconclusive |

**Key finding.** Kazakhstan's 1yr spread shows confirmed cointegration here,
in contrast to every pre-war-only specification (Sections 4.3–4.5), where it
consistently failed. With the sample roughly tripling (27 → 82 observations)
via inclusion of the post-war period, this supports the interpretation that
Kazakhstan's pre-war-only null results reflected insufficient statistical
power rather than a genuine absence of long-run relationship. The
maturity-dependent pattern persists, however: Kazakhstan's 5yr and 10yr
spreads still fail to show cointegration even with the larger full sample,
echoing the term-structure heterogeneity identified in the exploratory
analysis (Section 3), where oil-decoupling and domestic-liquidity effects
were found to strengthen specifically with maturity.

**An important methodological correction.** Across every model in Table 4.4,
the war dummy itself is statistically insignificant (all p ≥ 0.075). This
means the cointegration detected in these full-sample models should **not**
be interpreted as evidence that the war produced a detectable shift in the
long-run equilibrium level — the dummy carries no independent explanatory
power in this specification. Rather, both countries exhibit a genuine
long-run relationship between spread and domestic-liquidity/oil/VIX across
the *entire* sample period, largely independent of whether the war dummy is
included.

This result is best read as evidence that **a simple constant-shift dummy is
not an adequate tool for detecting the war's structural effect** within the
ARDL bounds-testing framework — a finding that directly motivates the
subsequent Markov-switching analysis (Section 5, which detects regime
changes endogenously rather than assuming a known break date) and the
VAR/IRF analysis (Section 6, which tests for changes in dynamic transmission
rather than a simple level shift) as necessary complements, not redundant
robustness checks.

Consequently, the pre-war-only results in Table 4.1 remain the most direct
and interpretable test of H1's "pre-war similarity" claim: Russia shows
robust pre-war cointegration at all three maturities; Kazakhstan's sample is
too short to test this claim with confidence in isolation. The full-sample
results in Table 4.4 are a separate finding — a long-run relationship exists
across the full period for both countries — and should not be conflated with
a war-driven shift specifically.

## 4.7 Robustness Check: Expanded Specification for Kazakhstan's 5yr and 10yr Spreads

Given Kazakhstan's persistent failure to show cointegration at longer
maturities, four candidate variables — external debt/GDP, reserves/M2,
inflation, and GDP — were each added individually (not jointly, to avoid
repeating the overfitting problem documented in Section 4.5) to the
full-sample, war-dummy specification.

**Table 4.5 — Expanded specification test, Kazakhstan 5yr/10yr**

| Maturity | Added variable | n | EC term (p-value) | Bounds F-stat | Verdict |
|---|---|---|---|---|---|
| 5yr | External debt/GDP | 77 | -0.108 (0.079) | 1.22 | No cointegration |
| 5yr | Reserves/M2 | 80 | -0.137 (0.020) | 1.30 | No cointegration |
| 5yr | Inflation | 81 | -0.127 (0.073) | 1.60 | No cointegration |
| 5yr | GDP | 77 | -0.088 (0.144) | 1.37 | No cointegration |
| 10yr | External debt/GDP | 77 | -0.171 (0.006) | 2.33 | No cointegration |
| 10yr | Reserves/M2 | 80 | -0.165 (0.020) | 2.05 | No cointegration |
| 10yr | Inflation | 81 | -0.132 (0.047) | 1.64 | No cointegration |
| 10yr | GDP | 77 | -0.129 (0.042) | 2.11 | No cointegration |

**Findings.** None of the four candidate variables recovers cointegration at
either maturity, even though several individual error-correction
coefficients reach conventional significance on their own (e.g., 10yr with
external debt/GDP, p=0.006). The joint bounds F-test remains decisively
below the lower critical bound in every case, consistent with the
coefficient-versus-joint-test distinction noted in Section 4.3.

This is treated as a genuine, informative null result rather than an
unresolved gap. Kazakhstan's 5yr and 10yr spreads do not exhibit a detectable
long-run relationship with domestic liquidity, global risk factors, or any
of the fiscal/macro variables tested here, even after accounting for the war
via a level-shift dummy. Combined with the exploratory finding that
oil-decoupling and liquidity-sensitivity effects were strongest specifically
at the long end (Section 3), this suggests Kazakhstan's long-maturity spread
may be governed by maturity-specific factors — market liquidity conditions,
issuance-specific dynamics, or investor composition — that fall outside the
macro-financial variable set examined in this paper. This is documented as a
legitimate scope limitation rather than a modeling failure to be pursued
further.

## 4.8 Interpreting the Long-Run Coefficients: Russia's Pre-War Baseline

Given Specification A's pre-war results (Table 4.1) are the most robust in
this analysis, the underlying coefficients are interpreted in full as the
baseline against which post-war transmission changes (Sections 5–6) are
compared.

**Table 4.6 — Russia, Specification A, pre-war: long-run coefficients**

| Variable | 1yr coef. (p) | 5yr coef. (p) | 10yr coef. (p) |
|---|---|---|---|
| Error-correction (own lag) | -0.318 (<0.001) | -0.324 (<0.001) | -0.329 (<0.001) |
| RUONIA | 0.027 (0.079) | 0.061 (<0.001) | 0.099 (<0.001) |
| Brent oil | 0.0061 (0.003) | 0.0024 (0.183) | 0.0012 (0.478) |
| VIX | 0.0204 (0.006) | 0.0111 (0.120) | 0.0078 (0.222) |

**The error-correction speed is remarkably stable across maturities**
(-0.318 to -0.329, all p<0.001), implying that roughly 32% of any deviation
from long-run equilibrium is corrected each month, for a half-life of
approximately 1.8 months, regardless of maturity. This indicates one
coherent equilibrium-reversion mechanism operating similarly across
Russia's credit curve pre-war.

**Domestic liquidity conditions (RUONIA) strengthen with maturity** — the
long-run coefficient roughly quadruples from 0.027 at 1yr (only marginally
significant) to 0.099 at 10yr (highly significant), with precision
increasing alongside magnitude. This is consistent with an economic
interpretation in which short-term interbank rate movements are often
transient (reflecting liquidity-management operations) and are therefore
only partially priced into the shortest-maturity spread, while persistent
liquidity conditions are more fully reflected in the compensation demanded
for longer-duration credit risk.

**Global risk factors (oil, VIX) show the opposite maturity pattern** —
both significant at 1yr (oil: p=0.003; VIX: p=0.006) but fading to
insignificance at 5yr and 10yr. This suggests Russia's short-end spread was,
pre-war, more exposed to immediate global risk-sentiment shocks, while the
long end was more structurally anchored to domestic monetary conditions.

**Oil's positive sign at 1yr** (higher oil prices associated with a *higher*
spread) is notable, since it runs counter to the "oil exporter, higher oil
→ lower credit risk" relationship documented for Kazakhstan. Given this
effect fades to insignificance at longer maturities, it is best interpreted
as a short-end-specific phenomenon — plausibly linked to capital-flow or
ruble-volatility dynamics that coincided with oil-price movements within
this particular pre-war sample — rather than a robust, curve-wide
relationship between oil prices and Russian credit risk.

**Overall interpretation.** Russia's pre-war credit market exhibits a
maturity-differentiated transmission mechanism: short-term spreads are
more sentiment- and global-shock-driven, while long-term spreads are more
structurally anchored to domestic monetary conditions. This coherent
baseline is the reference point against which the post-war VAR/IRF analysis
(Section 6) should assess whether — and how — this transmission mechanism
changed.

## 4.9 Summary and Implications for Subsequent Analysis

1. **Method justification.** ARDL bounds testing was the appropriate choice
   given the confirmed mixed order of integration (Kazakhstan I(1), Russia
   I(0)) across the two countries' dependent variables, a setting in which
   Johansen or Engle-Granger cointegration testing cannot be validly applied.

2. **Specification choice.** The minimal specification (Section 4.3) is the
   most defensible primary result for both countries; fuller specifications
   (Sections 4.4–4.5) demonstrate the overfitting risk inherent in adding
   regressors against Kazakhstan's small pre-war sample, and reveal that
   Russia's own cointegration finding is not fully robust to specification
   choice at longer maturities.

3. **Pre-war cointegration (H1).** Russia shows robust, consistent
   cointegration across all three maturities pre-war. Kazakhstan's 27-
   observation pre-war sample does not permit an equally confident test —
   documented as a data limitation rather than evidence against a
   relationship, a conclusion strengthened by Kazakhstan's 1yr spread
   achieving cointegration once tested on its full sample (Section 4.6).

4. **The war dummy does not drive cointegration.** A corrected, fair
   comparison (Section 4.6) shows the war dummy is statistically
   insignificant across every model tested. A simple level-shift dummy is
   not an adequate tool for capturing the war's structural effect in this
   framework — motivating the Markov-switching and VAR/IRF analyses that
   follow.

5. **Kazakhstan's maturity-specific limitation.** The 5yr and 10yr spreads
   fail to show cointegration under every specification and every candidate
   additional variable tested (Section 4.7), suggesting their post-war
   dynamics are governed by factors outside this paper's macro-financial
   variable set — a legitimate scope limitation to state directly.

6. **Carrying forward.** Russia's stable pre-war error-correction speed and
   maturity-differentiated sensitivity to domestic versus global factors
   (Section 4.8) form the baseline against which Section 6's post-war
   VAR/IRF results should be compared. Section 5's Markov-switching model
   should be evaluated for whether it detects a regime break near February
   2022, given the war dummy's failure to capture a level shift here.
