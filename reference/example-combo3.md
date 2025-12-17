# Three-drug combination example

Example using a combination of two experimental drugs, with EXNEX and
differential discounting.

## Details

This dataset involves a hypothetical dose-escalation study of
combination therapy with three treatment components. From two previous
studies `HistAgent1` and `HistAgent2`, historical data is available on
each of the treatments as single-agents, as well as two of the two-way
combinations. However, due to a difference in treatment schedule between
the `Combo` study and the historical studies, a stratification (through
`stratum_id`) is made between the groups to allow differential
discounting of the alternate-schedule data. The association is as below.

|                |                   |
|----------------|-------------------|
| group_id (j):  | stratum_id (s_j): |
| Combo (1)      | BID (1)           |
| HistAgent1 (2) | QD (2)            |
| HistAgent2 (3) | QD (2)            |

For additional robustness, EXNEX priors are used for all group-level
treatment components while not for the interaction parameters. This is
to limit the amount of borrowing in case of significant heterogeneity
across groups.

The complete model is as follows. As a function of doses
\\d_1,d_2,d_3\\, the DLT rate in group \\j\\ is, for \\j = 1,\ldots,3\\,
\$\$\text{logit}\\ \pi_j(d_1,d_2,d_3) = \text{logit}\Bigl( 1 -
\prod\_{i=1}^3 (1-\pi\_{ij}(d_i))\Bigr) +
\eta\_{j}^{(12)}\frac{d_1}{d_1^\*}\frac{d_2}{d_2^\*} +
\eta\_{j}^{(13)}\frac{d_1}{d_1^\*}\frac{d_3}{d_3^\*} +
\eta\_{j}^{(23)}\frac{d_2}{d_2^\*}\frac{d_3}{d_3^\*} +
\eta\_{j}^{(123)}\frac{d_1}{d_1^\*}\frac{d_2}{d_2^\*}\frac{d_3}{d_3^\*}.\$\$

In group \\j\\ each treatment component \\i\\ toxicity is modeled with
logistic regression, \$\$\text{logit}\\ \pi\_{ij}(d_i) = \log\\
\alpha\_{ij} + \beta\_{ij} \\ \log\\ \Bigl(\frac{d_i}{d_i^\*}\Bigr).\$\$
The intercept and log-slope parameters \\\boldsymbol\theta\_{ij} =
(\log\\ \alpha\_{ij}, \log\\ \beta\_{ij})\\ are are given an EXNEX prior
\$\$\boldsymbol\theta\_{ij} \sim p\_{ij} \text{BVN}(\boldsymbol\mu_i,
\boldsymbol\Sigma\_{ij}) + (1-p\_{ij}) \text{BVN}(\boldsymbol m\_{ij},
\boldsymbol S\_{ij}),\$\$ where the exchangeability weights are all
\\p\_{ij} = 0.9\\. The NEX parameters are set to \\\boldsymbol m\_{ij} =
(\text{logit}(1/3), \log\\ 1)\\, \\\boldsymbol S\_{ij} =
\text{diag}(2^2, 1^2)\\ for all components \\i=1,2,3\\ and groups \\j =
1,2,3\\, and the EX parameters are modeled hierarchically. The mean of
the exchangeable part has the distribution \$\$\boldsymbol\mu_i =
(\mu\_{\alpha i}, \mu\_{\beta i}) \sim \text{BVN}(\boldsymbol m_i,
\boldsymbol S_i),\$\$ with \\\boldsymbol m_i = (\text{logit}(1/3), \log
1)\\ and \\\boldsymbol S_i = \text{diag}(2^2, 1^2)\\ for each component
\\i = 1,2,3\\. For differentially discounting data from each schedule
(QD and BID), the covariance parameters for the exchangeable part
\$\$\Sigma\_{ij} = \left( \begin{array}{cc} \tau^2\_{\alpha s_j i} &
\rho_i \tau\_{\alpha s_j i} \tau\_{\beta s_j i}\\ \rho_i \tau\_{\alpha
s_j i} \tau\_{\beta s_j i} & \tau^2\_{\beta s_j i} \end{array}
\right).\$\$ are allowed to vary across groups \\j\\ depending on their
mapping to strata \\s(j)\\ as described above. For stratum \\s=1\\
(`BID`, which contains only the group \\j = 1\\ (`Combo`)), the standard
deviations are modeled as \$\$\tau\_{\alpha 1 i} \sim
\text{Log-Normal}(\log\\0.25, (\log 4 / 1.96)^2)\$\$ \$\$\tau\_{\beta 1
i} \sim \text{Log-Normal}(\log\\0.125, (\log 4 / 1.96)^2).\$\$ Whereas
in stratum \\s=2\\ (`QD`, which contains the historical groups \\j=2,3\\
(`HistData1`, `HistData2`)), the standard deviations are
\$\$\tau\_{\alpha 2 i} \sim \text{Log-Normal}(\log\\0.5, (\log 4 /
1.96)^2)\$\$ \$\$\tau\_{\beta 2 i} \sim \text{Log-Normal}(\log\\0.25,
(\log 4 / 1.96)^2).\$\$

For all interaction parameters \\\eta\_{j}^{(12)}\\,
\\\eta\_{j}^{(13)}\\, \\\eta\_{j}^{(23)}\\, and \\\eta\_{j}^{(123)}\\
(\\j = 1,2,3\\), the following prior is assumed: \$\$\eta\_{j}^{(\cdot)}
\sim \text{N}(\mu\_{\eta}^{(\cdot)},{\tau\_{\eta s_j}^{(\cdot)}}^2).\$\$
The exchangeability weights are \\p\_{\eta j}^{(\cdot)} = 0.9\\ for all
parameters with EXNEX. Here, for each \\\mu\_{\eta}^{(12)}\\,
\\\mu\_{\eta}^{(13)}\\, \\\mu\_{\eta}^{(23)}\\, and
\\\mu\_{\eta}^{(123)}\\, we take \$\$\mu\_{\eta}^{(\cdot)} \sim
\text{N}(0, 1/2),\$\$ and for each \\\tau\_{\eta s}^{(12)}\\,
\\\tau\_{\eta s}^{(13)}\\, \\\tau\_{\eta s}^{(23)}\\, and \\\tau\_{\eta
s}^{(123)}\\, \$\$\tau\_{\eta s}^{(\cdot)} \sim
\text{Log-Normal}(\log(0.25), (\log 2 / 1.96)^2),\$\$ for both strata
\\s = 1,2\\.

Below is the syntax for specifying this model in `blrm_exnex`.

## References

Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use
of co-data in clinical trials. *Statistics in Biopharmaceutical
Research*, 8(3), 345-354.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

## example combo3

library(abind)

dref <- c(500, 500, 1000)
num_comp <- 3
num_inter <- choose(3, 2) + 1
num_strata <- nlevels(hist_combo3$stratum_id)
num_groups <- nlevels(hist_combo3$group_id)

blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref[1])) |
      1 + I(log(drug_B / dref[2])) |
      1 + I(log(drug_C / dref[3])) |
      0
      + I(drug_A / dref[1] * drug_B / dref[2])
        + I(drug_A / dref[1] * drug_C / dref[3])
        + I(drug_B / dref[2] * drug_C / dref[3])
        + I(drug_A / dref[1] * drug_B / dref[2] * drug_C / dref[3]) |
      stratum_id / group_id,
  data = hist_combo3,
  prior_EX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(1/3), 0, diag(c(2^2, 1)))), FALSE),
  prior_EX_tau_comp = list(replicate(num_comp,
                                     mixmvnorm(c(1, log(c(0.25, 0.125)),
                                               diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE),
                           replicate(num_comp,
                                     mixmvnorm(c(1, log(2 * c(0.25, 0.125)),
                                               diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE)),
  prior_EX_mu_inter = mixmvnorm(c(1, rep.int(0, num_inter),
                                     diag((rep.int(sqrt(2) / 2, num_inter))^2))),
  prior_EX_tau_inter = replicate(num_strata,
                                 mixmvnorm(c(1, rep.int(log(0.25), num_inter),
                                             diag((rep.int(log(2) / 1.96, num_inter))^2))), FALSE),
  prior_EX_prob_comp = matrix(0.9, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1.0, nrow = num_groups, ncol = num_inter),
  prior_is_EXNEX_comp = rep(TRUE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_tau_dist = 1,
  prior_PD = FALSE
)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
## Recover user set sampling defaults
options(.user_mc_options)
```
