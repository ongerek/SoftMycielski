# Soft Mycielski Prediction and the Critical Temperature of Attention

Verification code, figure sources and reproduction instructions for the
manuscript *Soft Mycielski Prediction and the Critical Temperature of
Attention*.

The paper studies the predictor

$$p_\beta(a \mid x_1^n) = (1-\nu)\, \frac{\sum_{m < n} \mathbf{1}[x_{m+1}=a]\, e^{\beta \ell(m,n)}}{\sum_{m < n} e^{\beta \ell(m,n)}} + \frac{\nu}{A},$$

where $\ell(m,n)$ is the length of the agreement between the suffix
ending at position $m$ and the suffix ending at position $n$. This is a
single attention operation whose logits are linear in the match length,
so it is simultaneously an idealised induction head and a
finite-temperature relaxation of the Ehrenfeucht–Mycielski predictor.
The inverse temperature $\beta$ turns out to have a critical value: the
entropy rate $H$ of the source.

Everything in this repository is **verification code**: it recomputes
the numbers that appear in the paper and prints them to the console.
There is no build step and nothing to install beyond Octave.

## Requirements

* GNU Octave 8 or later, **or** MATLAB. No toolboxes are used; all code
  is written to the intersection of the two languages.
* Optional, for the PNG figures only: any Octave graphics toolkit. On a
  bare Debian/Ubuntu container that means
  `apt-get install gnuplot-nox fonts-freefont-otf` — the fonts are
  needed as well as the backend, or Octave's text renderer fails.
* Optional, for the geometry checker: Python 3 and `pdftotext`
  (`poppler-utils`).

Tested with GNU Octave 8.4 on Ubuntu 24.04.

## Quick start

```octave
>> run_all_verifications          % all eleven packages
>> run_all_verifications(1:4)     % a subset
>> run_all_verifications(7)       % a single package
```

or run one package directly:

```sh
cd verif_v5 && octave --no-gui --quiet runme.m
```

Each package is self-contained apart from the shared routines in
`common/`, which every `runme.m` puts on the path automatically, so the
packages work whichever directory you launch them from.

**Runtime.** The full sweep is dominated by `verif_v8` (Monte-Carlo
constants at large horizon), `verif_v9` (beyond-Markov testbeds) and the
temperature sweeps in `verif_v2`/`verif_v3`. Expect on the order of one
to two hours single-threaded. The packages are independent.

## Layout

```
common/                shared generators, match-length and predictor routines
verif_v1 ... verif_v11 the eleven verification packages, each with runme.m
figures/tikz/          figure sources (TikZ/pgfplots, self-contained)
figures/png/           Octave-rendered PNG versions of the data figures
figures/checkfig.py    geometry checker for compiled figures
run_all_verifications.m
```

## What each package verifies

| package | verifies | principal claims |
|---|---|---|
| `verif_v1` | Abel/geometric-mixture identity to machine precision; the zero-temperature limit; order-1 and order-2 temperature sweeps | Proposition 1, Theorem 1 |
| `verif_v2` | the critical-annealing law $\beta^\star(n)\to H^-$; self-critical schedules; the Mycielski-78 streaming predictor; order-weight profiles either side of the boundary | Proposition 2, the annealing and profile tables |
| `verif_v3` | oracle vs. self-critical excess across $n$; the context-restart trie (M78-R) against plain M78 and the exact predictor | the adaptivity table, the M78-R result |
| `verif_v4` | cross-entropy scaling of $\beta^\star(n)$ at three source entropies; the head-weight atom | the scaling and atom tables |
| `verif_v5` | the bridge representation: crossing statistic, endpoint pinning, pinned-atom moments, small-maximum law shape; adaptivity ratio at $n=32000$ | Lemma (crossing), Theorem (bridge) |
| `verif_v6` | the untuned critical temperature $\beta=H$ against the oracle; the rotation identity; arcsine vs. uniform argmax laws | Theorem (fixed upper bound), Lemma (pinned small maximum) |
| `verif_v7` | the second-order misallocation law, by exact conditional KL | Lemma (second-order misallocation) |
| `verif_v8` | the conjectured constants by Monte Carlo in the idealised channel model at horizons unreachable by sequence simulation; Rao–Blackwellised real-system anchors at $n=10^5$ and $2.5\times10^5$ | the constants table, the Conjecture |
| `verif_v9` | hidden-Markov and piecewise-stationary testbeds with Bayes, fixed-share and sigmoid gating | the oracle inequality, the tracking proposition |
| `verif_v10` | error probability (0–1 loss) against the Bayes floor; emits the three data figures in TikZ and PNG | the 0–1 loss columns, Figures (approach, separation, error probability) |
| `verif_v11` | the worked example behind the explanatory figures: match lengths and soft distributions on `abracadabra`, the LZ78 parse and trie, the realised fluctuation path and its bridge | the explanatory figures |

Package-level detail is in each `verif_vN/README.md`.

## Reading the output

The packages print tables to the console; those tables are the
deliverable. Where a script states a prediction before the numbers
(“`[prediction: ...]`”, “`[predicted ...]`”), the prediction was
registered before the measurement was taken, and the printed comparison
is the test. Several scripts deliberately print quantities that
*disagree* with an earlier hypothesis; those disagreements are discussed
in the paper rather than hidden.

Two estimator conventions recur and matter:

* **Rao–Blackwellisation.** Where the true conditional law of the
  testbed is known, excess loss is measured as an exact conditional
  Kullback–Leibler divergence rather than by averaging realised
  log-losses. At large $n$ the excess falls below the sampling noise of
  realised losses, so the raw estimator is not merely noisier but
  unusable; `verif_v7` shows the effect directly.
* **Seeds.** Every script sets its own seeds explicitly, so repeated
  runs reproduce the printed numbers exactly on the same Octave version.

## Figures

`figures/tikz/` holds eight self-contained `tikzpicture` bodies. Three
of them (`fig_approach`, `fig_separation`, `fig_errprob`) are generated
by `verif_v10/v10_maketikz.m` and carry an "auto-generated" header; the
other five are hand-written explanatory diagrams whose content is
verified by `verif_v11`. To use one:

```latex
\usepackage{tikz,pgfplots}
\usetikzlibrary{decorations.pathreplacing,arrows.meta,positioning,calc,fit,shapes.geometric}
\pgfplotsset{compat=1.16}
...
\input{figures/tikz/fig_mechanism.tex}
```

`figures/checkfig.py` reports the extent of a compiled figure and any
overlapping text, which is useful when adapting the figures to a
different column width:

```sh
python3 figures/checkfig.py myfigure.pdf 18.19
```

## Citing

See `CITATION.cff` (GitHub's citation widget reads this) or
`CITATION.bib`. Please cite the paper rather than this repository alone.

## License

MIT — see [LICENSE] 
