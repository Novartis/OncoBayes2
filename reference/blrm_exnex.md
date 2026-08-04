# Bayesian Logistic Regression Model for N-compounds with EXNEX

Bayesian Logistic Regression Model (BLRM) for N compounds using
EXchangability and NonEXchangability (EXNEX) modeling.

## Usage

``` r
blrm_exnex(
  formula,
  data,
  prior_EX_mu_comp,
  prior_EX_mu_mean_comp,
  prior_EX_mu_sd_comp,
  prior_EX_tau_comp,
  prior_EX_tau_mean_comp,
  prior_EX_tau_sd_comp,
  prior_EX_corr_eta_comp,
  prior_EX_mu_inter,
  prior_EX_mu_mean_inter,
  prior_EX_mu_sd_inter,
  prior_EX_tau_inter,
  prior_EX_tau_mean_inter,
  prior_EX_tau_sd_inter,
  prior_EX_corr_eta_inter,
  prior_is_EXNEX_inter,
  prior_is_EXNEX_comp,
  prior_EX_prob_comp,
  prior_EX_prob_inter,
  prior_NEX_mu_comp,
  prior_NEX_mu_mean_comp,
  prior_NEX_mu_sd_comp,
  prior_NEX_mu_inter,
  prior_NEX_mu_mean_inter,
  prior_NEX_mu_sd_inter,
  prior_tau_dist,
  sample_map = FALSE,
  iter = getOption("OncoBayes2.MC.iter", 2000),
  warmup = getOption("OncoBayes2.MC.warmup", 1000),
  save_warmup = getOption("OncoBayes2.MC.save_warmup", TRUE),
  thin = getOption("OncoBayes2.MC.thin", 1),
  init = getOption("OncoBayes2.MC.init", 0.5),
  chains = getOption("OncoBayes2.MC.chains", 4),
  cores = getOption("mc.cores", 1L),
  control = getOption("OncoBayes2.MC.control", list()),
  backend = getOption("OncoBayes2.MC.backend", "rstan"),
  prior_PD = FALSE,
  verbose = FALSE
)

# S3 method for class 'blrmfit'
print(x, ..., prob = 0.95, digits = 2)
```

## Arguments

- formula:

  the model formula describing the linear predictors of the model. The
  lhs of the formula is a two-column matrix which are the number of
  occured events and the number of times no event occured. The rhs of
  the formula defines the linear predictors for the marginal models for
  each drug component, then the interaction model and at last the
  grouping and optional stratum factors of the models. These elements of
  the formula are separated by a vertical bar. The marginal models must
  follow a intercept and slope form while the interaction model must not
  include an interaction term. See the examples below for an example
  instantiation.

- data:

  optional data frame containing the variables of the model. If not
  found in `data`, the variables are taken from `environment(formula)`.

- prior_EX_mu_comp:

  List of bivariate normal mixture priors for intercept and slope
  parameters \\\boldsymbol\mu_i = (\mu\_{\alpha i}, \mu\_{\beta i})\\ of
  each component. In case of a single drug model, then a mixture prior
  is accepted as well.

- prior_EX_mu_mean_comp, prior_EX_mu_sd_comp:

  **\[deprecated\]** Please use `prior_EX_mu_comp` instead. Mean and sd
  for the prior on the mean parameters \\\boldsymbol\mu_i =
  (\mu\_{\alpha i}, \mu\_{\beta i})\\ of each component. Two column
  matrix (intercept, log-slope) with one row per component.

- prior_EX_tau_comp:

  List of bivariate normal mixture priors for heterogeniety parameter
  \\\boldsymbol\tau\_{si} = (\tau\_{\alpha s i}, \tau\_{\beta s i})\\ of
  each stratum. If no differential discounting is required (i.e. if
  there is only one stratum \\s = 1\\), then it suffices to provide a
  bivariate normal mixture prior instead of a list with just one
  element.

- prior_EX_tau_mean_comp, prior_EX_tau_sd_comp:

  **\[deprecated\]** Please use `prior_EX_tau_comp` instead. Prior mean
  and sd for heterogeniety parameter \\\boldsymbol\tau\_{si} =
  (\tau\_{\alpha s i}, \tau\_{\beta s i})\\ of each stratum. If no
  differential discounting is required (i.e. if there is only one
  stratum \\s = 1\\), then it is a two-column matrix (intercept,
  log-slope) with one row per component. Otherwise it is a
  three-dimensional array whose first dimension indexes the strata,
  second dimension indexes the components, and third dimension of length
  two for (intercept, log-slope).

- prior_EX_corr_eta_comp:

  Prior LKJ correlation parameter for each component given as numeric
  vector. If missing, then a 1 is assumed corresponding to a marginal
  uniform prior of the correlation.

- prior_EX_mu_inter:

  Multivariate normal mixture prior for interaction parameter vector
  \\\boldsymbol{\mu\_{\eta}}\\. Dimension must correspond to the number
  of interactions.

- prior_EX_mu_mean_inter, prior_EX_mu_sd_inter:

  **\[deprecated\]** Please use `prior_EX_mu_inter` instead. Prior mean
  and sd for population mean parameters \\\mu\_{\eta k}\\ of each
  interaction parameter. Vector of length equal to the number of
  interactions.

- prior_EX_tau_inter:

  List of multivariate normal mixture priors for heterogeniety
  interaction parameter vector \\\boldsymbol{\tau\_{\eta s}}\\ of each
  stratum. If no differential discounting is required (i.e. if there is
  only one stratum \\s = 1\\), then it suffices to provide a mixture
  prior instead of a list with just one element.

- prior_EX_tau_mean_inter, prior_EX_tau_sd_inter:

  **\[deprecated\]** Please use `prior_EX_tau_inter` instead. Prior mean
  and sd for heterogeniety parameter \\\tau\_{\eta s k}\\ of each
  stratum. Matrix with one column per interaction and one row per
  stratum.

- prior_EX_corr_eta_inter:

  Prior LKJ correlation parameter for interaction given as numeric. If
  missing, then a 1 is assumed corresponding to a marginal uniform prior
  of the correlations.

- prior_is_EXNEX_inter:

  Defines if non-exchangability is admitted for a given interaction
  parameter. Logical vector of length equal to the number of
  interactions. If missing `FALSE` is assumed for all interactions.

- prior_is_EXNEX_comp:

  Defines if non-exchangability is admitted for a given component.
  Logical vector of length equal to the number of components. If missing
  `TRUE` is assumed for all components.

- prior_EX_prob_comp:

  Prior probability \\p\_{ij}\\ for exchangability of each component per
  group. Matrix with one column per component and one row per group.
  Values must lie in \\\[0-1\]\\ range.

- prior_EX_prob_inter:

  Prior probability \\p\_{\eta k j}\\ for exchangability of each
  interaction per group. Matrix with one column per interaction and one
  row per group. Values must lie in \\\[0-1\]\\ range.

- prior_NEX_mu_comp:

  List of bivariate normal mixture priors \\\boldsymbol m\_{ij}\\ and
  \\\boldsymbol s\_{ij} = \text{diag}(\boldsymbol S\_{ij})\\ of each
  component for non-exchangable case. If missing set to the same prior
  as given for the EX part. It is required that the specification be the
  same across groups j.

- prior_NEX_mu_mean_comp, prior_NEX_mu_sd_comp:

  **\[deprecated\]** Please use `prior_NEX_mu_comp` instead. Prior mean
  \\\boldsymbol m\_{ij}\\ and sd \\\boldsymbol s\_{ij} =
  \text{diag}(\boldsymbol S\_{ij})\\ of each component for
  non-exchangable case. Two column matrix (intercept, log-slope) with
  one row per component. If missing set to the same prior as given for
  the EX part. It is required that the specification be the same across
  groups j.

- prior_NEX_mu_inter:

  Multivariate normal mixture prior (mean \\m\_{\eta k j}\\, sd
  \\s\_{\eta k j}\\ and covariance) for the interaction parameter vector
  for non-exchangable case. Dimension must correspond to the number of
  interactions. If missing set to the same prior as given for the EX
  part.

- prior_NEX_mu_mean_inter, prior_NEX_mu_sd_inter:

  **\[deprecated\]** Please use `prior_NEX_mu_inter` instead. Prior mean
  \\m\_{\eta k j}\\ and sd \\s\_{\eta k j}\\ for each interaction
  parameter for non-exchangable case. Vector of length equal to the
  number of interactions. If missing set to the same prior as given for
  the EX part.

- prior_tau_dist:

  Defines the distribution used for heterogeniety parameters. Choices
  are 0=fixed to it's mean, 1=log-normal, 2=truncated normal or `NULL`
  shutting off the hierarchical structure of the model.

- sample_map:

  Logical flag (defaults to `FALSE`) controlling inclusion of MAP priors
  for each stratum defined as part of the generated posterior. If set to
  `TRUE` then the posterior samples will contain `map_log_beta` and
  `map_eta` variables.

- iter:

  number of iterations (including warmup).

- warmup:

  number of warmup iterations.

- save_warmup:

  save warmup samples (`TRUE` / `FALSE`). Only if set to `TRUE`, then
  all random variables are saved in the posterior. This substantially
  increases the storage needs of the posterior.

- thin:

  period of saving samples.

- init:

  positive number to specify uniform range on unconstrained space for
  random initialization. See
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- chains:

  number of Markov chains.

- cores:

  number of cores for parallel sampling of chains.

- control:

  additional sampler parameters for NuTS algorithm.

- backend:

  sets Stan backend to be used. Possible choices are `"rstan"` (default)
  or `"cmdstanr"`.

- prior_PD:

  Logical flag (defaults to `FALSE`) indicating if to sample the prior
  predictive distribution instead of conditioning on the data.

- verbose:

  Logical flag (defaults to `FALSE`) controlling if additional output
  like stan progress is reported.

- x:

  `blrmfit` object to print

- ...:

  not used in this function

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- digits:

  number of digits to show

## Value

The function returns a S3 object of type `blrmfit`.

## Details

`blrm_exnex` is a flexible function for Bayesian meta-analytic modeling
of binomial count data. In particular, it is designed to model counts of
the number of observed dose limiting toxicities (DLTs) by dose, for
guiding dose-escalation studies in Oncology. To accommodate dose
escalation over more than one agent, the dose may consist of
combinations of study drugs, with any number of treatment components.

In the simplest case, the aim is to model the probability \\\pi\\ that a
patient experiences a DLT, by complementing the binomial likelihood with
a monotone logistic regression

\$\$\text{logit}\\\pi(d) = \log\\\alpha + \beta \\ t(d),\$\$

where \\\beta \> 0\\. Most typically, \\d\\ represents the dose, and
\\t(d)\\ is an appropriate transformation, such as \\t(d) = \log (d \big
/ d^\*)\\. A joint prior on \\\boldsymbol \theta = (\log\\\alpha,
\log\\\beta)\\ completes the model and ensures monotonicity \\\beta \>
0\\.

Many extensions are possible. The function supports general combination
regimens, and also provides framework for Bayesian meta-analysis of
dose-toxicity data from multiple historical and concurrent sources.

For an example of a single-agent trial refer to
[`example-single-agent()`](https://opensource.nibr.com/OncoBayes2/reference/example-single-agent.md).

## Functions

- `print(blrmfit)`: print function.

## Combination of two treatments

For a combination of two treatment components, the basic modeling
framework is that the DLT rate \\\pi(d_1,d_2)\\ is comprised of (1) a
"no-interaction" baseline model \\\tilde \pi(d_1,d_2)\\ driven by the
single-agent toxicity of each component, and (2) optional interaction
terms \\\gamma(d_1,d_2)\\ representing synergy or antagonism between the
drugs. On the log-odds scale,

\$\$\text{logit} \\\pi(d_1,d_2) = \text{logit} \\ \tilde \pi(d_1,d_2) +
\eta \\ \gamma(d_1,d_2). \$\$

The "no interaction" part \\\tilde \pi(d_1,d_2)\\ represents the
probability of a DLT triggered by either treatment component acting
*independently*. That is, \$\$ \tilde \pi(d_1,d_2) = 1- (1 -
\pi_1(d_1))(1 - \pi_2(d_2)). \$\$ In simple terms, P(no DLT for
combination) = P(no DLT for drug 1) \* P(no DLT from drug 2). To
complete this part, the treatment components can then be modeled with
monotone logistic regressions as before.

\$\$\text{logit} \\ \pi_i(d_i) = \log\\ \alpha_i + \beta_i \\
t(d_i),\$\$

where \\t(d_i)\\ is a monotone transformation of the doses of the
respective drug component \$i\$, such as \\t(d_i) = \log (d_i \big /
d_i^\*)\\.

The inclusion of an interaction term \\\gamma(d_1,d_2)\\ allows DLT
rates above or below the "no-interaction" rate. The magnitude of the
interaction term may also be made dependent on the doses (or other
covariates) through regression. As an example, one could let

\$\$\gamma(d_1, d_2) = \frac{d_1}{d_1^\*} \frac{d_1}{d_2^\*}.\$\$

The specific functional form is specified in the usual notation for a
design matrix. The interaction model must respect the constraint that
whenever any dose approaches zero, then the interaction term must vanish
as well. Therefore, the interaction model must not include an intercept
term which would violate this consistency requirement. A dual
combination example can be found in
[`example-combo2()`](https://opensource.nibr.com/OncoBayes2/reference/example-combo2.md).

## General combinations

The model is extended to general combination treatments consisting of
\\N\\ components by expressing the probability \\\pi\\ on the logit
scale as

\$\$ \text{logit} \\ \pi(d_1,\ldots,d_N) = \text{logit} \Bigl( 1 -
\prod\_{i = 1}^N ( 1 - \pi_i(d_i) ) \Bigr) + \sum\_{k=1}^K \eta_k \\
\gamma_k(d_1,\ldots,d_N), \$\$

Multiple drug-drug interactions among the \\N\\ components are now
possible, and are represented through the \\K\\ interaction terms
\\\gamma_k\\.

Regression models can be again be specified for each \\\pi_i\\ and
\\\gamma_k\\, such as

\$\$ \text{logit}\\ \pi_i(d_i) = \log\\ \alpha_i + \beta_i \\ t(d_i)
\$\$

Interactions for some subset \\I(k) \subset \\1,\ldots,N \\\\ of the
treatment components can be modeled with regression as well, for example
on products of doses,

\$\$ \gamma_k(d_1,\ldots,d_N) = \prod\_{i \in I(k)}
\frac{d_i}{d_i^\*}.\$\$

For example, \\I(k) = \\1,2,3\\\\ results in the three-way interaction
term

\$\$ \frac{d_1}{d_1^\*} \frac{d_2}{d_2^\*} \frac{d_3}{d_3^\*} \$\$

for drugs 1, 2, and 3.

For a triple combination example please refer to
[`example-combo3()`](https://opensource.nibr.com/OncoBayes2/reference/example-combo3.md).

## Meta-analytic framework

Information on the toxicity of a drug may be available from multiple
studies or sources. Furthermore, one may wish to stratify observations
within a single study (for example into groups of patients corresponding
to different geographic regions, or multiple dosing `dose_info`
corresponding to different schedules).

`blrm_exnex` provides tools for robust Bayesian hierarchical modeling to
jointly model data from multiple sources. An additional index \\j=1,
\ldots, J\\ on the parameters and observations denotes the \\J\\ groups.
The resulting model allows the DLT rate to differ across the groups. The
general \\N\\-component model becomes

\$\$ \text{logit} \\ \pi_j(d_1,\ldots,d_N) = \text{logit} \Bigl( 1 -
\prod\_{i = 1}^N ( 1 - \pi\_{ij}(d_i) ) \Bigr) + \sum\_{k=1}^K
\eta\_{kj} \\ \gamma\_{k}(d_1,\ldots,d_N), \$\$

for groups \\j = 1,\ldots,J\\. The component toxicities \\\pi\_{ij}\\
and interaction terms \\\gamma\_{k}\\ are modelled, as before, through
regression. For example, \\\pi\_{ij}\\ could be a logistic regression on
\\t(d_i) = \log(d_i/d_i^\*)\\ with intercept and log-slope \\\boldsymbol
\theta\_{ij}\\, and \\\gamma\_{k}\\ regressed with coefficient
\\\eta\_{kj}\\ on a product \\\prod\_{i\in I(k)} (d_i/d_i^\*)\\ for some
subset \\I(k)\\ of components.

Thus, for \\j=1,\ldots,J\\, we now have group-specific parameters
\\\boldsymbol\theta\_{ij} = (\log\\ \alpha\_{ij}, \log\\ \beta\_{ij})\\
and \\\boldsymbol\nu\_{j} = (\eta\_{1j}, \ldots, \eta\_{Kj})\\ for each
component \\i=1,\ldots,N\\ and interaction \\k=1,\ldots,K\\.

The structure of the prior on
\\(\boldsymbol\theta\_{i1},\ldots,\boldsymbol\theta\_{iJ})\\ and
\\(\boldsymbol\nu\_{1}, \ldots, \boldsymbol\nu\_{J})\\ determines how
much information will be shared across groups \\j\\. Several modeling
choices are available in the function.

- *EX (Full exchangeability):* One can assume the parameters are
  conditionally exchangeable given hyperparameters

  \$\$\boldsymbol \theta\_{ij} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \theta i}, \boldsymbol \Sigma\_{\boldsymbol \theta
  i} \bigr), \$\$

  independently across groups \\j = 1,\ldots, J\\ and treatment
  components \\i=1,\ldots,N\\. The covariance matrix \\\boldsymbol
  \Sigma\_{\boldsymbol \theta i}\\ captures the patterns of cross-group
  heterogeneity, and is parametrized with standard deviations
  \\\boldsymbol \tau\_{\boldsymbol\theta i} = (\tau\_{\alpha i},
  \tau\_{\beta i})\\ and the correlation \\\rho_i\\. Similarly for the
  interactions, the fully-exchangeable model is

  \$\$\boldsymbol \nu\_{j} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \nu}, \boldsymbol \Sigma\_{\boldsymbol \nu}
  \bigr)\$\$

  for groups \\j = 1,\ldots, J\\ and interactions \\k=1,\ldots,K\\, and
  the prior on the covariance matrix \\\boldsymbol \Sigma\_{\boldsymbol
  \nu}\\ captures the amount of heterogeneity expected in the
  interaction terms a-priori. The covariance is again parametrized with
  standard deviations \\(\tau\_{\eta 1}, \ldots, \tau\_{\eta K})\\ and
  its correlation matrix.

- *Differential discounting:* For one or more of the groups
  \\j=1,\ldots,J\\, larger deviations of \\\boldsymbol\theta\_{ij}\\ may
  be expected from the mean \\\boldsymbol\mu_i\\, or of the interactions
  \\\eta\_{kj}\\ from the mean \\\mu\_{\eta,k}\\. Such differential
  heterogeneity can be modeled by mapping the groups \\j = 1,\ldots,J\\
  to *strata* through \\s_j \in \\1,\ldots,S\\\\, and modifying the
  model specification to \$\$\boldsymbol \theta\_{ij} \sim
  \text{N}\bigl( \boldsymbol \mu\_{\boldsymbol \theta i}, \boldsymbol
  \Sigma\_{\boldsymbol \theta ij} \bigr), \$\$ where \$\$\boldsymbol
  \Sigma\_{\boldsymbol \theta ij} = \left( \begin{array}{cc}
  \tau^2\_{\alpha s_j i} & \rho_i \tau\_{\alpha s_j i} \tau\_{\beta s_j
  i}\\ \rho_i \tau\_{\alpha s_j i} \tau\_{\beta s_j i} & \tau^2\_{\beta
  s_j i} \end{array} \right).\$\$ For the interactions, the model
  becomes \$\$\boldsymbol \nu\_{j} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \nu}, \boldsymbol \Sigma\_{\boldsymbol \nu j}
  \bigr),\$\$ where the covariance matrix \\\boldsymbol
  \Sigma\_{\boldsymbol \nu j}\\ is modelled as stratum specific standard
  deviations \\(\tau\_{\eta 1 s_j}, \ldots, \tau\_{\eta K s_j})\\ and a
  stratum independent correlation matrix. Each stratum \\s=1,\ldots,S\\
  then corresponds to its own set of standard deviations \\\tau\\
  leading to different discounting per stratum. Independent priors are
  specified for the component parameters \\\tau\_{\alpha s i}\\ and
  \\\tau\_{\beta s i}\\ and for the interaction parameters \\\tau\_{\eta
  s k}\\ for each stratum \\s=1,\ldots,S\\. Inference for strata \\s\\
  where the prior is centered on larger values of the \\\tau\\
  parameters will exhibit less shrinkage towards the the means,
  \\\boldsymbol\mu\_{\boldsymbol \theta i}\\ and \\\boldsymbol
  \mu\_{\boldsymbol \nu}\\ respectively.

- *EXNEX (Partial exchangeability):* Another mechansim for increasing
  robustness is to introduce mixture priors for the group-specific
  parameters, where one mixture component is shared across groups, and
  the other is group-specific. The result, known as an
  EXchangeable-NonEXchangeable (EXNEX) type prior, has a form

  \$\$\boldsymbol \theta\_{ij} \sim p\_{\boldsymbol \theta ij}\\
  \text{N}\bigl( \boldsymbol \mu\_{\boldsymbol \theta i}, \boldsymbol
  \Sigma\_{\boldsymbol \theta i} \bigr) +(1-p\_{\boldsymbol \theta
  ij})\\ \text{N}\bigl(\boldsymbol m\_{\boldsymbol \theta ij},
  \boldsymbol S\_{\boldsymbol \theta ij}\bigr)\$\$

  when applied to the treatment-component parameters, and

  \$\$\boldsymbol \nu\_{kj} \sim p\_{\boldsymbol \nu\_{kj}}
  \\\text{N}\bigl(\mu\_{\boldsymbol \nu}, \boldsymbol
  \Sigma\_{\boldsymbol \nu}\bigr)\_k + (1-p\_{\boldsymbol \nu\_{kj}})\\
  \text{N}(m\_{\boldsymbol \nu\_{kj}}, s^2\_{\boldsymbol \nu\_{kj}})\$\$

  when applied to the interaction parameters. The *exchangeability
  weights* \\p\_{\boldsymbol \theta ij}\\ and \\p\_{\boldsymbol
  \nu\_{kj}}\\ are fixed constants in the interval \\\[0,1\]\\ that
  control the degree to which inference for group \\j\\ is informed by
  the exchangeable mixture components. Larger values for the weights
  correspond to greater exchange of information, while smaller values
  increase robustness in case of outlying observations in individual
  groups \\j\\.

## References

Neuenschwander B, Roychoudhury S, Schmidli H (2016). “On the use of
co-data in clinical trials.” *Statistics in Biopharmaceutical Research*,
**8**(3), 345–354.
[doi:10.1080/19466315.2016.1174149](https://doi.org/10.1080/19466315.2016.1174149)
.

Neuenschwander B, Wandel S, Roychoudhury S, Bailey S (2016). “Robust
exchangeability designs for early phase clinical trials with multiple
strata.” *Pharmaceutical Statistics*, **15**(2), 123–134.
[doi:10.1002/pst.1730](https://doi.org/10.1002/pst.1730) .

Neuenschwander B, Branson M, Gsponer T (2008). “Critical aspects of the
Bayesian approach to phase I cancer trials.” *Statistics in Medicine*,
**27**(13), 2420–2439.
[doi:10.1002/sim.3230](https://doi.org/10.1002/sim.3230) .

Neuenschwander B, Matano A, Tang Z, Roychoudhury S, Wandel S, Bailey S
(2014). “A Bayesian Industry Approach to Phase I Combination Trials in
Oncology.” In *Statistical Methods in Drug Combination Studies*, volume
69. CRC Press. [doi:10.1201/b17965-9](https://doi.org/10.1201/b17965-9)
.

## Examples

``` r
.user_mc_options <- options()

# fit an example model. See documentation for "combo3" example
example_model("combo3")
#> Running combo3 example:
#>  ## example combo3
#> 
#>  library(abind)
#> 
#>  dref <- c(500, 500, 1000)
#>  num_comp <- 3
#>  num_inter <- choose(3, 2) + 1
#>  num_strata <- nlevels(hist_combo3$stratum_id)
#>  num_groups <- nlevels(hist_combo3$group_id)
#> 
#>  blrmfit <- blrm_exnex(
#>    cbind(num_toxicities, num_patients - num_toxicities) ~
#>      1 + I(log(drug_A / dref[1])) |
#>        1 + I(log(drug_B / dref[2])) |
#>        1 + I(log(drug_C / dref[3])) |
#>        0
#>        + I(drug_A / dref[1] * drug_B / dref[2])
#>          + I(drug_A / dref[1] * drug_C / dref[3])
#>          + I(drug_B / dref[2] * drug_C / dref[3])
#>          + I(drug_A / dref[1] * drug_B / dref[2] * drug_C / dref[3]) |
#>        stratum_id / group_id,
#>    data = hist_combo3,
#>    prior_EX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(1/3), 0, diag(c(2^2, 1)))), FALSE),
#>    prior_EX_tau_comp = list(replicate(num_comp,
#>                                       mixmvnorm(c(1, log(c(0.25, 0.125)),
#>                                                 diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE),
#>                             replicate(num_comp,
#>                                       mixmvnorm(c(1, log(2 * c(0.25, 0.125)),
#>                                                 diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE)),
#>    prior_EX_mu_inter = mixmvnorm(c(1, rep.int(0, num_inter),
#>                                       diag((rep.int(sqrt(2) / 2, num_inter))^2))),
#>    prior_EX_tau_inter = replicate(num_strata,
#>                                   mixmvnorm(c(1, rep.int(log(0.25), num_inter),
#>                                               diag((rep.int(log(2) / 1.96, num_inter))^2))), FALSE),
#>    prior_EX_prob_comp = matrix(0.9, nrow = num_groups, ncol = num_comp),
#>    prior_EX_prob_inter = matrix(1.0, nrow = num_groups, ncol = num_inter),
#>    prior_is_EXNEX_comp = rep(TRUE, num_comp),
#>    prior_is_EXNEX_inter = rep(FALSE, num_inter),
#>    prior_tau_dist = 1,
#>    prior_PD = FALSE
#>  )

# print a summary of the prior
prior_summary(blrmfit, digits = 3)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Mixture configuration
#> ---------------------
#> EXNEX components : 3 
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) I(log(drug_C/dref[3])) 
#>                      1                      1                      1 
#> 
#> EXNEX interactions: 0 
#> interaction
#>                  I(drug_A/dref[1] * drug_B/dref[2]) 
#>                                                   0 
#>                  I(drug_A/dref[1] * drug_C/dref[3]) 
#>                                                   0 
#>                  I(drug_B/dref[2] * drug_C/dref[3]) 
#>                                                   0 
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 
#>                                                   0 
#> 
#> Prior probability for exchangeability per group
#>             component
#> group        I(log(drug_A/dref[1])) I(log(drug_B/dref[2]))
#>   Combo                         0.9                    0.9
#>   HistAgent1                    0.9                    0.9
#>   HistAgent2                    0.9                    0.9
#>             component
#> group        I(log(drug_C/dref[3]))
#>   Combo                         0.9
#>   HistAgent1                    0.9
#>   HistAgent2                    0.9
#> 
#>             interaction
#> group        I(drug_A/dref[1] * drug_B/dref[2])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_A/dref[1] * drug_C/dref[3])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_B/dref[2] * drug_C/dref[3])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])
#>   Combo                                                        1
#>   HistAgent1                                                   1
#>   HistAgent2                                                   1
#> 
#> EXchangable hyperparameter priors
#> ---------------------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope    rho
#> component              mix                                                                       
#> I(log(drug_A/dref[1])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_B/dref[2])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_C/dref[3])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> 
#> Heterogeneity tau_log_beta (log-normal)
#>                                         prior weight m_tau_intercept m_tau_log_slope s_tau_intercept s_tau_log_slope    rho
#> stratum   component              mix                                                                                       
#> stratum_1 I(log(drug_A/dref[1])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#>           I(log(drug_B/dref[2])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#>           I(log(drug_C/dref[3])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#> stratum_2 I(log(drug_A/dref[1])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#>           I(log(drug_B/dref[2])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#>           I(log(drug_C/dref[3])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#> 
#> Correlation LKJ
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) I(log(drug_C/dref[3])) 
#>                      1                      1                      1 
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior     w  m[1]  m[2]  m[3]  m[4]  s[1]  s[2]  s[3]  s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> mix                                                                                                                    
#> comp1       1.000 0.000 0.000 0.000 0.000 0.707 0.707 0.707 0.707    0.000    0.000    0.000    0.000    0.000    0.000
#> 
#> Heterogeneity tau_eta (log-normal)
#>                 prior      w   m[1]   m[2]   m[3]   m[4]   s[1]   s[2]   s[3]   s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> stratum   mix                                                                                                                             
#> stratum_1 comp1        1.000 -1.386 -1.386 -1.386 -1.386  0.354  0.354  0.354  0.354    0.000    0.000    0.000    0.000    0.000    0.000
#> stratum_2 comp1        1.000 -1.386 -1.386 -1.386 -1.386  0.354  0.354  0.354  0.354    0.000    0.000    0.000    0.000    0.000    0.000
#> 
#> Correlation LKJ
#> interaction 
#>           1 
#> 
#> NonEXchangable priors
#> ---------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope    rho
#> component              mix                                                                       
#> I(log(drug_A/dref[1])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_B/dref[2])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_C/dref[3])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior     w  m[1]  m[2]  m[3]  m[4]  s[1]  s[2]  s[3]  s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> mix                                                                                                                    
#> comp1       1.000 0.000 0.000 0.000 0.000 0.707 0.707 0.707 0.707    0.000    0.000    0.000    0.000    0.000    0.000

# print a summary of the posterior (model parameters)
print(blrmfit)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Number of observations: 18 
#> Number of groups      : 3 
#> Number of strata      : 2 
#> Number of components  : 3 
#> Number of interactions: 4 
#> EXNEX components      : 3 
#> EXNEX interactions    : 0 
#> 
#> Observations per group:
#>        Group n Stratum n_total
#> 1      Combo 3     BID      17
#> 2 HistAgent1 7      QD      52
#> 3 HistAgent2 8      QD      32
#> 
#> Groups per stratum:
#>   Stratum Groups n_total
#> 1     BID      1      17
#> 2      QD      2      84
#> 
#> Component posterior:
#> Population mean posterior mu_log_beta
#> intercept:
#>                        mean se_mean  sd q2.5  q50 q97.5 n_eff Rhat
#> I(log(drug_A/dref[1])) -1.5   0.023 1.1 -3.9 -1.4  0.33  2541    1
#> I(log(drug_B/dref[2])) -3.1   0.018 1.1 -5.3 -3.1 -1.11  3866    1
#> I(log(drug_C/dref[3])) -1.2   0.021 1.4 -3.9 -1.2  1.70  4367    1
#> log-slope:
#>                          mean se_mean   sd q2.5    q50 q97.5 n_eff Rhat
#> I(log(drug_A/dref[1]))  0.153   0.013 0.99 -1.8  0.174   2.0  5748    1
#> I(log(drug_B/dref[2])) -0.072   0.012 0.84 -1.9 -0.021   1.4  5112    1
#> I(log(drug_C/dref[3]))  0.248   0.012 0.89 -1.6  0.295   1.9  5074    1
#> 
#> Population heterogeniety posterior tau_log_beta
#> intercept:
#>                            mean se_mean   sd  q2.5  q50 q97.5 n_eff Rhat
#> BID,I(log(drug_A/dref[1])) 0.31  0.0037 0.25 0.058 0.24  0.96  6711    1
#> QD,I(log(drug_A/dref[1]))  0.61  0.0068 0.44 0.129 0.49  1.77  5562    1
#> BID,I(log(drug_B/dref[2])) 0.31  0.0039 0.24 0.065 0.25  0.93  7030    1
#> QD,I(log(drug_B/dref[2]))  0.63  0.0070 0.47 0.133 0.50  1.83  5915    1
#> BID,I(log(drug_C/dref[3])) 0.32  0.0037 0.25 0.062 0.25  0.95  6826    1
#> QD,I(log(drug_C/dref[3]))  0.62  0.0070 0.48 0.122 0.49  1.87  6779    1
#> log-slope:
#>                            mean se_mean   sd  q2.5  q50 q97.5 n_eff Rhat
#> BID,I(log(drug_A/dref[1])) 0.16  0.0020 0.13 0.031 0.12  0.50  7286    1
#> QD,I(log(drug_A/dref[1]))  0.32  0.0042 0.25 0.060 0.25  1.01  5760    1
#> BID,I(log(drug_B/dref[2])) 0.16  0.0021 0.13 0.031 0.12  0.51  7167    1
#> QD,I(log(drug_B/dref[2]))  0.31  0.0039 0.25 0.059 0.24  0.96  6763    1
#> BID,I(log(drug_C/dref[3])) 0.16  0.0018 0.12 0.032 0.13  0.49  6373    1
#> QD,I(log(drug_C/dref[3]))  0.32  0.0036 0.25 0.061 0.25  0.95  6590    1
#> 
#> Population correlation posterior rho_log_beta
#>                           mean se_mean   sd  q2.5     q50 q97.5 n_eff Rhat
#> I(log(drug_A/dref[1])) -0.0013  0.0070 0.57 -0.94  0.0019  0.94  6276    1
#> I(log(drug_B/dref[2])) -0.0108  0.0076 0.58 -0.95 -0.0385  0.95  5208    1
#> I(log(drug_C/dref[3]))  0.0126  0.0074 0.57 -0.95 -0.0059  0.95  5555    1
#> 
#> Interaction model posterior:
#> Population mean posterior mu_eta
#>                                                       mean se_mean   sd q2.5
#> I(drug_A/dref[1] * drug_B/dref[2])                  -0.498  0.0076 0.51 -1.5
#> I(drug_A/dref[1] * drug_C/dref[3])                   0.163  0.0084 0.72 -1.2
#> I(drug_B/dref[2] * drug_C/dref[3])                   0.074  0.0085 0.69 -1.3
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.071  0.0083 0.69 -1.3
#>                                                        q50 q97.5 n_eff Rhat
#> I(drug_A/dref[1] * drug_B/dref[2])                  -0.489  0.48  4427    1
#> I(drug_A/dref[1] * drug_C/dref[3])                   0.158  1.62  7284    1
#> I(drug_B/dref[2] * drug_C/dref[3])                   0.068  1.47  6766    1
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.074  1.46  7011    1
#> 
#> Population heterogeniety posterior tau_eta
#>                                                         mean se_mean    sd q2.5
#> BID,I(drug_A/dref[1] * drug_B/dref[2])                  0.26  0.0012 0.096 0.12
#> QD,I(drug_A/dref[1] * drug_B/dref[2])                   0.27  0.0013 0.099 0.12
#> BID,I(drug_A/dref[1] * drug_C/dref[3])                  0.26  0.0012 0.098 0.12
#> QD,I(drug_A/dref[1] * drug_C/dref[3])                   0.27  0.0013 0.098 0.13
#> BID,I(drug_B/dref[2] * drug_C/dref[3])                  0.27  0.0013 0.098 0.13
#> QD,I(drug_B/dref[2] * drug_C/dref[3])                   0.26  0.0012 0.094 0.13
#> BID,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 0.26  0.0013 0.096 0.13
#> QD,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.27  0.0012 0.096 0.12
#>                                                          q50 q97.5 n_eff Rhat
#> BID,I(drug_A/dref[1] * drug_B/dref[2])                  0.25  0.49  6547    1
#> QD,I(drug_A/dref[1] * drug_B/dref[2])                   0.25  0.51  6123    1
#> BID,I(drug_A/dref[1] * drug_C/dref[3])                  0.25  0.50  7571    1
#> QD,I(drug_A/dref[1] * drug_C/dref[3])                   0.25  0.50  7552    1
#> BID,I(drug_B/dref[2] * drug_C/dref[3])                  0.25  0.50  8213    1
#> QD,I(drug_B/dref[2] * drug_C/dref[3])                   0.25  0.49  8232    1
#> BID,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 0.25  0.50  6336    1
#> QD,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.25  0.50  7379    1
#> 
#> Population correlation posterior Sigma_corr_eta
#>                                                                                                            mean
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.0000
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.0072
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.0027
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  -0.0079
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.0072
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.0000
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.0072
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  -0.0150
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.0027
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.0072
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.0000
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.0041
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0079
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0150
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0041
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.0000
#>                                                                                                         se_mean
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                        NA
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   5.6e-03
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   5.3e-03
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  5.2e-03
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   5.6e-03
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   1.4e-18
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   7.1e-03
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  6.7e-03
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   5.3e-03
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   7.1e-03
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   1.1e-18
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  8.2e-03
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  5.2e-03
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  6.7e-03
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  8.2e-03
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 2.4e-18
#>                                                                                                              sd
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                   0.0e+00
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   4.4e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   4.4e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  4.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   4.4e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   8.4e-17
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   4.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  4.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   4.4e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   4.5e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   6.7e-17
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  4.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  4.5e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  4.5e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  4.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 6.4e-17
#>                                                                                                          q2.5
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.00
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.81
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.81
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  -0.80
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.81
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.00
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.81
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  -0.82
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.81
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.81
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.00
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.82
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.80
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.82
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.82
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.00
#>                                                                                                             q50
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.0000
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.0088
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.0095
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  -0.0044
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.0088
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.0000
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.0106
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  -0.0208
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.0095
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.0106
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.0000
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.0080
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0044
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0208
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.0080
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.0000
#>                                                                                                         q97.5
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.00
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.81
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.81
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                   0.80
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.81
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.00
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.81
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                   0.81
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.81
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.81
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.00
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                   0.81
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.80
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.81
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.81
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.00
#>                                                                                                         n_eff
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                      NA
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    6164
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    6990
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                   7382
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                    6164
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    4109
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    3997
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                   4575
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                    6990
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    3997
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    3849
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                   2877
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   7382
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   4575
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   2877
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  4175
#>                                                                                                         Rhat
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                     NA
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                      1
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                      1
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                     1
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                      1
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                      1
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                      1
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                     1
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                      1
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                      1
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                      1
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                     1
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     1
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     1
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     1
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])    1

# summary of posterior for DLT rate by dose for observed covariate levels
summ <- summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(hist_combo3, summ))
#>    stratum_id   group_id drug_A drug_B drug_C num_toxicities num_patients
#> 1         BID      Combo    400    800     80              0            4
#> 2         BID      Combo    400    800    160              1            8
#> 3         BID      Combo    400    800    240              2            5
#> 4          QD HistAgent1      0    240      0              0            5
#> 5          QD HistAgent1      0    400      0              0           10
#> 6          QD HistAgent1      0    800      0              1           11
#> 7          QD HistAgent1    240     80      0              0            4
#> 8          QD HistAgent1    400    400      0              1            6
#> 9          QD HistAgent1    400    800      0              0           11
#> 10         QD HistAgent1    400   1000      0              1            5
#> 11         QD HistAgent2      0      0     80              0            3
#> 12         QD HistAgent2      0      0    160              0            3
#> 13         QD HistAgent2      0      0    320              0            6
#> 14         QD HistAgent2      0      0    480              0            3
#> 15         QD HistAgent2      0      0    640              1            5
#> 16         QD HistAgent2    400      0    160              1            6
#> 17         QD HistAgent2    400      0    320              2            5
#> 18         QD HistAgent2    400      0    240              1            1
#>          mean         sd         2.5%         50%      97.5% [0,0.16]
#> 1  0.17351171 0.07605072 5.857934e-02 0.161587553 0.35313861  0.49225
#> 2  0.18912604 0.07947909 6.640411e-02 0.177392047 0.37722395  0.40525
#> 3  0.20707208 0.09031008 6.927814e-02 0.193764819 0.42022014  0.33875
#> 4  0.02298783 0.02298376 3.732154e-04 0.015837196 0.08369477  0.99925
#> 5  0.03506052 0.02857199 2.515864e-03 0.027590078 0.10752114  0.99725
#> 6  0.07728761 0.05555395 9.264175e-03 0.064743297 0.21378562  0.92050
#> 7  0.08005803 0.06669722 3.205476e-03 0.061543854 0.24977519  0.87725
#> 8  0.11792238 0.06086212 2.845231e-02 0.109144238 0.26332406  0.78050
#> 9  0.11004942 0.05428188 3.098655e-02 0.101205956 0.24039318  0.83125
#> 10 0.11201947 0.06571285 2.495420e-02 0.098501395 0.28113515  0.79975
#> 11 0.02006749 0.03109868 8.156197e-08 0.006249682 0.11076434  0.99475
#> 12 0.03189990 0.04015456 6.333634e-06 0.015800587 0.14489272  0.98325
#> 13 0.05857092 0.05630983 4.679099e-04 0.041063062 0.20777150  0.93300
#> 14 0.09339636 0.07555067 4.169594e-03 0.074649583 0.28411075  0.82325
#> 15 0.13997493 0.10372139 1.234474e-02 0.115863145 0.39751696  0.65575
#> 16 0.24936663 0.11041498 7.254350e-02 0.236169556 0.49930278  0.22350
#> 17 0.27610209 0.10839309 1.003781e-01 0.263397368 0.51407497  0.13950
#> 18 0.26198607 0.10853293 8.953842e-02 0.247457524 0.50334952  0.17725
#>    (0.16,0.33] (0.33,1]
#> 1      0.47025  0.03750
#> 2      0.53950  0.05525
#> 3      0.56100  0.10025
#> 4      0.00075  0.00000
#> 5      0.00275  0.00000
#> 6      0.07800  0.00150
#> 7      0.11700  0.00575
#> 8      0.21575  0.00375
#> 9      0.16750  0.00125
#> 10     0.19200  0.00825
#> 11     0.00525  0.00000
#> 12     0.01675  0.00000
#> 13     0.06650  0.00050
#> 14     0.16700  0.00975
#> 15     0.28150  0.06275
#> 16     0.54750  0.22900
#> 17     0.57125  0.28925
#> 18     0.56850  0.25425

# summary of posterior for DLT rate by dose for new set of covariate levels
newdata <- expand.grid(
  stratum_id = "BID", group_id = "Combo",
  drug_A = 400, drug_B = 800, drug_C = c(320, 400, 600, 800),
  stringsAsFactors = FALSE
)
summ_pred <- summary(blrmfit, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(newdata, summ_pred))
#>   stratum_id group_id drug_A drug_B drug_C      mean        sd       2.5%
#> 1        BID    Combo    400    800    320 0.2274812 0.1079270 0.06906200
#> 2        BID    Combo    400    800    400 0.2500150 0.1300658 0.06478604
#> 3        BID    Combo    400    800    600 0.3121192 0.1902993 0.05380448
#> 4        BID    Combo    400    800    800 0.3766590 0.2438401 0.03813018
#>         50%     97.5% [0,0.16] (0.16,0.33] (0.33,1]
#> 1 0.2107673 0.4850876  0.29800     0.53600  0.16600
#> 2 0.2289374 0.5698641  0.27425     0.49225  0.23350
#> 3 0.2746696 0.7549850  0.23600     0.36850  0.39550
#> 4 0.3343236 0.8905154  0.21950     0.27475  0.50575

# update the model after observing additional data
newdata$num_patients <- rep(3, nrow(newdata))
newdata$num_toxicities <- c(0, 1, 2, 2)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
blrmfit_new <- update(blrmfit,
  data = rbind(hist_combo3, newdata) %>%
    arrange(stratum_id, group_id)
)

# updated posterior summary
summ_upd <- summary(blrmfit_new, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(newdata, summ_upd))
#>   stratum_id group_id drug_A drug_B drug_C num_patients num_toxicities
#> 1        BID    Combo    400    800    320            3              0
#> 2        BID    Combo    400    800    400            3              1
#> 3        BID    Combo    400    800    600            3              2
#> 4        BID    Combo    400    800    800            3              2
#>        mean         sd      2.5%       50%     97.5% [0,0.16] (0.16,0.33]
#> 1 0.2707876 0.08036741 0.1329872 0.2643171 0.4419459  0.07025     0.70375
#> 2 0.3139414 0.09110246 0.1570042 0.3071101 0.5083422  0.02975     0.55875
#> 3 0.4383488 0.12884865 0.2010431 0.4322673 0.7080711  0.00600     0.19900
#> 4 0.5649788 0.16690687 0.2343986 0.5715281 0.8712571  0.00525     0.09025
#>   (0.33,1]
#> 1   0.2260
#> 2   0.4115
#> 3   0.7950
#> 4   0.9045
## Recover user set sampling defaults
options(.user_mc_options)
```
