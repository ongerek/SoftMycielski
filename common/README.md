# common/

Routines shared by two or more verification packages. Every
`verif_vN/runme.m` adds this directory to the path automatically, so
nothing here needs to be copied or installed.

| file | purpose |
|---|---|
| `gen_markov.m` | order-1 circulant testbed: alphabet size `A`, successor probability `1-delta`, the remaining mass spread uniformly. Returns the sequence, the transition matrix and the closed-form entropy rate. The chain is doubly stochastic, so the stationary law is uniform and single symbols carry no information about the next one — every predictable structure has to come from context. |
| `gen_markov2.m` | the order-2 companion, successor rule `s(i,j) = ((i+j) mod A) + 1`. |
| `matchlens.m` | match lengths `ell(m,t)` of every past position against the suffix ending at `t`. |
| `soft_myc_dist.m` | the soft Mycielski predictive distribution for a given inverse temperature and smoothing weight. |
| `beta_sweep.m` | sequential log-loss over a grid of inverse temperatures, evaluated on a common suffix window. |
| `m78_predict.m` | the Mycielski-78 streaming predictor over an LZ78 trie (plain parse and path-update variants). |

The M78-R predictor with per-order suffix tracking lives in
`verif_v3/m78r_predict.m`, since only that package uses it.

Numerical conventions: seeds are always passed explicitly, logs are
natural (losses are in nats), and no toolbox functions are used.
