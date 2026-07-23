# verif_v6 — Criticality without tuning, and the argmax laws

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(6)` from the repository root.

| module | does |
|---|---|
| `v6_critical.m` V6.1 | the untuned critical temperature `beta = H` against the oracle `beta*(n)` |
| `v6_critical.m` V6.2 | the rotation (Vervaat-type) identity on realised profiles |
| `v6_critical.m` V6.3 | the end-window argmax probability under both tilts |

**What to look for.** V6.1: the ratio stays bounded, so no tuning is
needed to reach the fixed-temperature rate class. V6.3 is the sharpest
confirmation in the paper — scaled by `sqrt(log n)` the fixed-tilt
end-window probability is flat, and scaled by `log n` the self-tilt one
is flat, matching the arcsine and uniform argmax laws to a few per cent.
V6.2 reports a stable finite-size gap between the two sides of the
rotation identity; it is a finite-size effect, quantified rather than
suppressed.

**Runtime** tens of minutes.
