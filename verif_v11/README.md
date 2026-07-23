# verif_v11 — The worked example behind the explanatory figures

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(11)` from the repository root.

| module | does |
|---|---|
| `v11_figdata.m` V11.1 | match lengths and continuations on `x = abracadabra`; soft predictive distributions across temperatures including both limits; the LZ78 parse, the trie, and the M78-R suffix nodes; a realised fluctuation path with its chord and bridge |

**What to look for.** Every number printed here appears in a figure, so
the figures illustrate a true computation rather than a plausible
sketch. Two checks are built in: the `beta = 0` row must equal the
empirical marginal exactly, and the `beta -> Inf` row must be the hard
Mycielski rule. The example is also honest about the dictionary's cost —
exact search finds a four-symbol match while at this length the trie
indexes the past only to depth two.

The companion geometry checker for compiled figures is
`../figures/checkfig.py`.

**Runtime** under a second.
