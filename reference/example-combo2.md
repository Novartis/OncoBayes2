# Two-drug combination example

Example using a combination of two experimental drugs.

## Details

The following example is described in the reference Neuenschwander, B.
et al (2016). The data are described in the help page for
`codata_combo2`. In the study `trial_AB`, the risk of DLT was studied as
a function of dose for two drugs, drug A and drug B. Historical
information on the toxicity profiles of these two drugs was available
from single agent trials `trial_A` and `trial_B`. Another study `IIT`
was run concurrently to `trial_AB`, and studies the same combination.

The model described in Neuenschwander, et al (2016) is adapted as
follows. For groups \\j = 1,\ldots, 4\\ representing each of the four
sources of data mentioned above, \$\$\text{logit}\\ \pi\_{1j}(d_1) =
\log\\ \alpha\_{1j} + \beta\_{1j} \\ \log\\
\Bigl(\frac{d_1}{d_1^\*}\Bigr),\$\$ and \$\$\text{logit}\\
\pi\_{2j}(d_2) = \log\\ \alpha\_{2j} + \beta\_{2j} \\ \log\\
\Bigl(\frac{d_2}{d_2^\*}\Bigr),\$\$ are logistic regressions for the
single-agent toxicity of drugs A and B, respectively, when administered
in group \\j\\. Conditional on the regression parameters
\\\boldsymbol\theta\_{1j} = (\log \\ \alpha\_{1j}, \log \\
\beta\_{1j})\\ and \\\boldsymbol\theta\_{2j} = (\log \\ \alpha\_{2j},
\log \\ \beta\_{2j})\\, the toxicity \\\pi\_{j}(d_1, d_2)\\ for the
combination is modeled as the "no-interaction" DLT rate,
\$\$\tilde\pi\_{j}(d_1, d_2) = 1 - (1-\pi\_{1j}(d_1) )(1-
\pi\_{2j}(d_2))\$\$ with a single interaction term added on the log odds
scale, \$\$\text{logit} \\ \pi\_{j}(d_1, d_2) = \text{logit} \\
\tilde\pi\_{j}(d_1, d_2) + \eta_j
\frac{d_1}{d_1^\*}\frac{d_2}{d_2^\*}.\$\$ A hierarchical model across
the four groups \\j\\ allows dose-toxicity information to be shared
through common hyperparameters.

For the component parameters \\\boldsymbol\theta\_{ij}\\,
\$\$\boldsymbol\theta\_{ij} \sim \text{BVN}(\boldsymbol \mu_i,
\boldsymbol\Sigma_i).\$\$ For the mean, a further prior is specified as
\$\$\boldsymbol\mu_i = (\mu\_{\alpha i}, \mu\_{\beta i}) \sim
\text{BVN}(\boldsymbol m_i, \boldsymbol S_i),\$\$ while in the
manuscript the prior \\\boldsymbol m_i = (\text{logit}\\ 0.1, \log 1)\\
and \\\boldsymbol S_i = \text{diag}(3.33^2, 1^2)\\ for each \\i = 1,2\\
is used, we deviate here and use instead \\\boldsymbol m_i =
(\text{logit}\\ 0.2, \log 1)\\ and \\\boldsymbol S_i = \text{diag}(2^2,
1^2)\\. For the standard deviations and correlation parameters in the
covariance matrix, \$\$\boldsymbol\Sigma_i = \left( \begin{array}{cc}
\tau^2\_{\alpha i} & \rho_i \tau\_{\alpha i} \tau\_{\beta i}\\ \rho_i
\tau\_{\alpha i} \tau\_{\beta i} & \tau^2\_{\beta i} \end{array}
\right), \$\$ the specified priors are \\\tau\_{\alpha i} \sim
\text{Log-Normal}(\log\\ 0.25, ((\log 4) / 1.96)^2)\\,

\\\tau\_{\beta i} \sim \text{Log-Normal}(\log\\ 0.125, ((\log 4) /
1.96)^2)\\, and \\\rho_i \sim \text{U}(-1,1)\\ for \\i = 1,2\\.

For the interaction parameters \\\eta_j\\ in each group, the
hierarchical model has \$\$\eta_j \sim \text{N}(\mu\_\eta,
\tau^2\_\eta),\$\$ for \\j = 1,\ldots, 4\\, with \\\mu\_\eta \sim
\text{N}(0, 1.121^2)\\ and \\\tau\_\eta \sim \text{Log-Normal}(\log\\
0.125, ((\log 4) / 1.96)^2).\\

Below is the syntax for specifying this fully exchangeable model in
`blrm_exnex`.

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

dref <- c(6, 960)

num_comp   <- 2 # two investigational drugs
num_inter  <- 1 # one drug-drug interaction needs to be modeled
num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
num_strata <- 1 # no stratification needed

blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref[1])) |
      1 + I(log(drug_B / dref[2])) |
      0 + I(drug_A / dref[1] * drug_B / dref[2]) |
      group_id,
  data = codata_combo2,
  prior_EX_mu_comp  = list(mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
                           mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))),
  prior_EX_tau_comp = list(mixmvnorm(c(1,
                                       log(0.250), log(0.125),
                                       diag(c(log(4)/1.96, log(4)/1.96)^2))),
                           mixmvnorm(c(1,
                                       log(0.250), log(0.125),
                                       diag(c(log(4)/1.96, log(4)/1.96)^2)))),
  prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
  prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
  prior_is_EXNEX_comp = rep(FALSE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
  prior_tau_dist = 1
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
