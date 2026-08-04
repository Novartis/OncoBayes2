# Single Agent Example

Example using a single experimental drug.

## Details

The single agent example is described in the reference Neuenschwander,
B. et al (2008). The data are described in the help page for `hist_SA`.
In this case, the data come from only one study, with the treatment
being only single agent. Hence the model specified does not involve a
hierarchical prior for the intercept and log-slope parameters. The model
described in Neuenschwander, et al (2008) is adapted as follows:
\$\$\text{logit}\\ \pi(d) = \log\\ \alpha + \beta \\ \log\\
\Bigl(\frac{d}{d^\*}\Bigr),\$\$ where \\d^\* = 250\\, and the prior for
\\\boldsymbol\theta = (\log\\ \alpha, \log\\ \beta)\\ is
\$\$\boldsymbol\theta \sim \text{N}(\boldsymbol m, \boldsymbol S),\$\$
and \\\boldsymbol m = (\text{logit}\\ 0.5, \log\\ 1)\\ and \\\boldsymbol
S = \text{diag}(2^2, 1^2)\\ are constants.

The above model is non-hierarchical. To disable the hierarchical model
structure of the `blrm_exnex` framework, the user can specify the option
`prior_tau_dist=NULL`. This will internally set all the heterogeniety
parameters (\\\tau^2\_\alpha\\ and \\\tau^2\_\beta\\) to zero.

## References

Neuenschwander B, Branson M, Gsponer T (2008). “Critical aspects of the
Bayesian approach to phase I cancer trials.” *Statistics in Medicine*,
**27**(13), 2420–2439.
[doi:10.1002/sim.3230](https://doi.org/10.1002/sim.3230) .

## Examples

``` r
.user_mc_options <- options()

## Example from Neuenschwander, B., et al. (2009). Stats in Medicine

dref <- 50

## Since there is no prior information the hierarchical model
## is not used in this example by setting tau to (almost) 0.
blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
      1 + log(drug_A / dref) |
      0 |
      group_id,
  data = hist_SA,
  prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
  ## Setting prior_tau_dist=NULL disables the hierarchical prior which is
  ## not required in this example as we analyze a single trial.
  prior_tau_dist = NULL,
  prior_PD = FALSE
)
## Recover user set sampling defaults
options(.user_mc_options)
```
