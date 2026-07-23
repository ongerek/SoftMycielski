# verif_v1 — Identity, limit, and the temperature–entropy law

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(1)` from the repository root.

| module | does |
|---|---|
| `runme.m` V1.1 | checks the Abel-summation identity of Proposition 1 on three chains |
| `runme.m` V1.2 | the zero-temperature limit: `\|p_beta - p_inf\|_1` against `e^{-beta}` |
| `runme.m` V1.3 | temperature sweep over six source entropies, with quadratic refinement of the minimiser |
| `runme.m` V1.4 | the match-length law `L_t H / log t -> 1` |
| `order2_test.m` V1.5 | the same sweep on an order-2 source whose order-1 statistics are exactly uniform |
| `abel_check.m` | helper: relative error between the direct sum and the geometric-mixture form |

**What to look for.** V1.1 should report a relative error at machine
precision (order `1e-15`); it is an algebraic identity, so anything
larger is a bug, not noise. V1.3 is the first sign of the phase
boundary: `beta*/H` sits near one across entropies that differ by more
than an order of magnitude. V1.5 is the mechanism test — on a source
where single symbols carry no information, the optimal temperature
shifts markedly upward, which is what the smoothing-phase account
predicts and a "temperature is just a smoothing knob" account does not.

**Writes** `results_v1.txt`. **Runtime** a few minutes.
