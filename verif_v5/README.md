# verif_v5 — Bridge representation and pinned-atom diagnostics

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(5)` from the repository root.

| module | does |
|---|---|
| `v5_bridge.m` V5.1–V5.3 | the crossing statistic `log(t P_L)`; the endpoint statistic under fixed and self-critical tilts; pinned-atom moments; the shape of the small-maximum law |
| `v5_ratio32k.m` V5.4 | the adaptivity ratio at `n = 32000` |

**What to look for.** The endpoint spread is the whole point: under the
self-critical tilt it is flat in `n`, while under a fixed temperature it
grows — that is endpoint pinning, measured. The moment ratio shrinks by
a factor matching the predicted one to within a per cent. The CDF-shape
test (linear against quadratic) is ordering-consistent but has *not*
reached its asymptotic regime at the horizons simulated here, and the
script says so.

**Runtime** tens of minutes.
