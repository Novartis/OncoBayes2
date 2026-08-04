# OncoBayes2

Bayesian logistic regression model with optional
EXchangeability-NonEXchangeability parameter modelling for flexible
borrowing from historical or concurrent data-sources. The safety model
can guide dose-escalation decisions for adaptive Oncology phase I
dose-escalation trials which involve an arbitrary number of drugs.

## Global Options

|  |  |  |
|----|----|----|
| Option | Default | Description |
| `OncoBayes2.MC.warmup` | 1000 | MCMC warmup iterations |
| `OncoBayes2.MC.iter` | 2000 | total MCMC iterations |
| `OncoBayes2.MC.save_warmup` | TRUE | save warmup samples |
| `OncoBayes2.MC.chains` | 4 | MCMC chains |
| `OncoBayes2.MC.thin` | 1 | MCMC thinning |
| `OncoBayes2.MC.control` | `list(adapt_delta=0.99,` | sets `control` argument for Stan call |
|  | `stepsize=0.1`) |  |
| `OncoBayes2.MC.backend` | rstan | Backend used to run Stan (`rstan` or `cmdstanr`) |
| `OncoBayes2.abbreviate.min` | 0 | Minimal length of variable names |
|  |  | when abbreviating variable names. |
|  |  | The default 0 disables abbreviation. |

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

Stan Development Team (2024). RStan: the R interface to Stan. R package
version 2.32.6. https://mc-stan.org

## See also

Useful links:

- <https://opensource.nibr.com/OncoBayes2/>

- Report bugs at <https://github.com/Novartis/OncoBayes2/issues>

## Author

**Maintainer**: Sebastian Weber <sebastian.weber@novartis.com>

Authors:

- Sebastian Weber <sebastian.weber@novartis.com>

- Lukas A. Widmer <lukas_andreas.widmer@novartis.com>

- Andrew Bean <andrew.bean@novartis.com>

Other contributors:

- Novartis Pharma AG \[copyright holder\]

- Trustees of Columbia University (R/stanmodels.R, configure,
  configure.win) \[copyright holder\]
