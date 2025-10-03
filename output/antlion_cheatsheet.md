---
output:
  html_document: default
  pdf_document: default
---

# Discrete‑time survival analysis in **brms**
### Take‑aways from an antlion eclosion experiment

---

## 1  Data wrangling
* **Start point:** one row per insect — `id`, `temperature`, `time`, `status` (1 = eclosed, 0 = censored).
* **Expand to person‑period table:**

```r
long <- raw %>% 
  mutate(end    = if_else(is.na(time), exp_len, time),
         status = if_else(is.na(time), 0, 1)) %>% 
  tidyr::uncount(end, .id = "day") %>% 
  mutate(event = as.integer(day == end & status == 1),
         day_z = scale(day)[, 1])
```

---

## 2  Model syntax

```r
formula <- bf(event | trials(1) ~ temperature +
              s(day_z, by = temperature, k = 10) +
              (1 | p | id))
family  = bernoulli(link = "cloglog")
```

---

## 3  Priors and regularisation

| Component                     | Prior specification                 | Why?                               |
|-------------------------------|-------------------------------------|------------------------------------|
| Fixed effects & intercept     | `normal(0, 0.4)`                   | Mild shrinkage                     |
| ID random‑effect SD           | `normal(0, 0.5), lb = 0` (half‑N)  | Avoids σ → 0 funnel                |
| Spline SDs (`sds_*`)          | `exponential(1.5)` (mean ≈ 0.67)   | Allows wiggle, limits over‑fit     |

> Always shrink **SD** parameters, **not** the individual `r_*` weights.

---

## 4  HMC diagnostics

* **Divergences** → geometry; raise `adapt_delta` or tighten SD priors.  
* **Treedepth warnings** → raise `max_treedepth`; unrelated to divergences.

Helper to extract divergent draws (works for *rstan* and *cmdstanr*):

```r
get_divergent <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit$fit)) {
    np <- posterior::nuts_params(fit$fit)
    return(as.vector(np[, , "divergent__"]) == 1)
  }
  np <- bayesplot::nuts_params(fit$fit)
  if ("Parameter" %in% names(np))
    return(np %>% filter(Parameter == "divergent__") %>% pull(Value) == 1)
  np$divergent__ == 1
}
```

---

## 5  Spline tips & prior‑predictive checks

* `k` too low **and** tight SD → almost linear curve; interaction invisible.  
* Target **posterior EDF ≈ 3–6** per smooth by tuning `k` and the SD prior.  
* Use `sample_prior = "only"` to plot **prior hazard curves** and calibrate before fitting to data.

---

## 6  Recommended workflow

1. **Prior‑predictive check** (PPC₀).  
2. Fit with `adapt_delta ≥ 0.95`, `max_treedepth = 12`.  
3. Check divergences, R‑hat, ESS.  
4. Adjust SD priors or `k` until **zero divergences**.  
5. Posterior‑predictive checks.  
6. Compare candidate models with `loo()`.

---

> Stan diagnostics ask *“Did we sample the distribution you wrote?”*  
> Prior‑predictive checks ask *“Did you write the right distribution?”*

*Cheat‑sheet generated from ChatGPT conversation – 20 May 2025*
