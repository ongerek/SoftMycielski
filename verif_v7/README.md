# verif_v7 — The second-order misallocation law

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(7)` from the repository root.

| module | does |
|---|---|
| `v7_secondorder.m` V7.1 | excess loss from mixing weight `w` of a wrong predictor into the truth, against the predicted quadratic law |

**What to look for.** `excess/w^2` should be flat and `excess/w` should
not. It is, to three decimals, against predictions computed in closed
form beforehand.

This module also carries a methodological lesson recorded in the paper:
the *raw* Monte-Carlo estimator fails here at small `w`, precisely
because the first-order term cancels in conditional mean, so realised
transition noise swamps the `w^2` signal. The script therefore computes
the conditional Kullback–Leibler divergence exactly. The
ill-conditioning is the phenomenon being tested.

**Runtime** seconds.
