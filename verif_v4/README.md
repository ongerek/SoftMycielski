# verif_v4 — Cross-entropy scaling and the head-weight atom

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(4)` from the repository root.

| module | does |
|---|---|
| `v4_scaling.m` V4.1 | the annealing experiment repeated at two further source entropies, reporting the normalised offset and the excess scaled by `log n` |
| `v4_scaling.m` V4.2 | moments of the realised order-0 weight `W_0` at fixed temperature |

**What to look for.** V4.2 discriminates two accounts of the head
channel: the idealised profile predicts `E[W_0^2]/E[W_0]` proportional
to `eps` (so declining toward zero), while the excursion account
predicts it stays constant. It stays constant. V4.1 is the measurement
that *refuted* an earlier conjectured scaling law for the optimal
offset; the low-entropy source has a supercritical optimum throughout
the accessible range, which is discussed rather than smoothed over.

**Runtime** tens of minutes.
