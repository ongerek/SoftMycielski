# verif_v2 — Critical annealing, self-tuning, and the streaming predictor

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(2)` from the repository root.

| module | does |
|---|---|
| `anneal_test.m` V2.1 | the annealing law: `beta*(n)` approaching `H` from below over a dyadic grid of `n` |
| `selftune_test.m` V2.2 | parameter-free schedules `log t / (L_t + 1)` and a damped variant, against fixed and oracle temperatures on three sources |
| `m78_and_profile_test.m` V2.3/V2.4 | the Mycielski-78 streaming predictor (speed and the context-reset penalty) and order-weight profiles either side of the boundary |

**What to look for.** The self-critical schedules beat the *oracle*
fixed temperature on every source without tuning anything: that is the
observation the paper's average-case theory was written to explain. In
V2.4 the profile is decaying below `H`, flat at `H`, and climbing above
it — the three phases, measured.

**Runtime** tens of minutes (temperature sweeps at several `n`).
