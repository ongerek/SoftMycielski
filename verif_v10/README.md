# verif_v10 — Error probability, and figure generation

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(10)` from the repository root.

| module | does |
|---|---|
| `v10_errprob.m` V10.1 | 0–1 loss of the plug-in decision against the Bayes floor, across temperatures, sources and `n`; writes `errprob_data.mat` |
| `v10_maketikz.m` V10.2 | emits three TikZ/pgfplots figure bodies into `../figures/tikz/` |
| `v10_plots_png.m` V10.3 | Octave-rendered PNG versions into `../figures/png/`; skipped automatically when no graphics toolkit is present |
| `err_excess2.m` | helper: error-probability excess for the order-2 source |

**What to look for.** The decision-level advantage of adaptivity is
*larger* than the codelength advantage, and grows faster in `n`. That is
not a restatement of the log-loss result: 0–1 loss rewards only a
correct arg max, which is where the endpoint-pinned bridge concentrates
mass. Note also one place where the ordering does not hold — on the
high-floor source, all subcritical predictors decide almost identically
and the metric cannot discriminate.

**Writes** `errprob_data.mat`, and figures into `../figures/`.
**Runtime** tens of minutes for V10.1; the figure modules are fast.

For the PNG module on a bare container:
`apt-get install gnuplot-nox fonts-freefont-otf` (the fonts matter as
much as the backend).
