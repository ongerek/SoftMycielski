# verif_v9 — Beyond Markov: hidden-Markov and piecewise-stationary sources

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(9)` from the repository root.

| module | does |
|---|---|
| `v9_hmm.m` V9.1 | a two-state hidden-Markov source (infinite observable order) with the exact predictive law from a forward filter; streaming experts combined by Bayes mixture, fixed share, and the closed-form sigmoid gate |
| `v9_piecewise.m` V9.2 | a source whose parameter alternates between two regimes, with a windowed expert and tracking |

**What to look for.** V9.1: the Bayes mixture matches its best expert to
three decimals, so principled gating is free at first order; the
un-fitted sigmoid gate recovers about two thirds of the gap, which is
reported as the imperfection it is. V9.2: global experts pay the
regime-averaging floor, and the fixed-share tracker beats *every*
individual expert within *every* regime by re-selecting after each
switch. The full-horizon Bayes overhead saturates its theoretical bound
almost exactly.

**Runtime** tens of minutes.
