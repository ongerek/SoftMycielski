# verif_v8 — The constants, by Monte Carlo and by measurement

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(8)` from the repository root.

| module | does |
|---|---|
| `v8_constants.m` V8.1 | the idealised channel model simulated directly at large horizon, giving the bridge shape integral and both rate constants |
| `v8_tiltsweep.m` V8.1b | a tilt sweep for the tuned fixed-temperature constant |
| `v8_precision.m` V8.1c | a high-precision block at the largest horizons |
| `v8_bign.m` V8.2 | Rao–Blackwellised real-system anchors at `n = 1e5` and `2.5e5` |

**What to look for.** The model needs only the surprisal walk, not the
sequence, so it runs at horizons corresponding to astronomically long
sequences — which is what makes the constants accessible at all. Two
honest caveats are printed rather than hidden: the bridge shape integral
does *not* converge to its Brownian value (a persistent discrete
correction at the scale that matters), and the real-system constants sit
at roughly twice the bare model values at every horizon reachable by
sequence simulation. That factor is the one question the paper's
conjecture retains.

**Runtime** the longest package, of order tens of minutes to an hour.
