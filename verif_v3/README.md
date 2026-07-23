# verif_v3 — Adaptivity across n, and the context-restart trie

Run with `octave --no-gui --quiet runme.m`, or `run_all_verifications(3)` from the repository root.

| module | does |
|---|---|
| `v3_anneal_self.m` V3.1/V3.2 | oracle vs. self-critical excess across `n`, and the same loss binned by the realised match length |
| `v3_trie.m` V3.3 | plain LZ78-parse M78, the path-update variant, and M78-R against the exact predictor |
| `m78r_predict.m` | the M78-R predictor: one live trie node per context order, advanced by a single child lookup per symbol |

**What to look for.** The binned table in V3.2 localises where
adaptivity pays: the gain is largest at *short* realised matches, which
is the bias-dilution mechanism rather than a tail effect. In V3.3,
M78-R with a self-tuned temperature comes within a few thousandths of a
nat of the entropy rate while remaining fully streaming and
parameter-free — the best figure measured in the paper.

**Runtime** tens of minutes.
