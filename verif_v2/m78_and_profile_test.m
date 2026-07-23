%M78_AND_PROFILE_TEST  Modules V2.3 and V2.4.
%
%   V2.3  Streaming soft Mycielski-78 vs exact soft Mycielski:
%         log-loss penalty of the dictionary restriction and wall-time
%         per prediction.  (LZ78 bound: c(n) <= n/((1-eps)log_A n).)
%   V2.4  Order-weight profiles W_k on REAL counts across the phase
%         boundary beta < H, beta = H, beta > H  (Proposition 2).

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;
[~, ~, H] = gen_markov(10, A, delta);

% ---------------- V2.3 ---------------------------------------------------
n = 8000;  t0 = 4000;  beta = 0.7;  seeds = [1 2];
fprintf('[V2.3] Soft Mycielski-78 vs exact  (delta=%.2f, H=%.4f, beta=%.2f)\n', ...
        delta, H, beta);
Lex = 0; Llz = 0; Lpa = 0; Tex = 0; Tlz = 0; Tpa = 0; NP = 0;
for s = seeds
  x = gen_markov(n, A, delta, 7000 + s);
  % exact predictor over the eval window
  tstart = tic();
  le = 0;
  for t = t0:(n-1)
    ell  = matchlens(x, t);
    p    = soft_myc_dist(ell, x(1:t), beta, A, nu);
    le   = le - log(p(x(t+1)));
  end
  Tex = Tex + toc(tstart) / (n - t0);
  Lex = Lex + le / (n - t0);
  % trie predictors over the full stream (evaluated on the same window)
  [ll, np, tp] = m78_predict(x, A, beta, nu, t0, 'lz');
  Llz = Llz + ll;  Tlz = Tlz + tp;  NP = NP + np;
  [lp, ~,  tp] = m78_predict(x, A, beta, nu, t0, 'path');
  Lpa = Lpa + lp;  Tpa = Tpa + tp;
end
ns = numel(seeds);
fprintf('   %-22s %10s %14s\n', 'predictor', 'loss', 'time/symbol');
fprintf('   %-22s %10.4f %12.2e s\n', 'exact soft Mycielski', Lex/ns, Tex/ns);
fprintf('   %-22s %10.4f %12.2e s\n', 'M78 (lz update)',      Llz/ns, Tlz/ns);
fprintf('   %-22s %10.4f %12.2e s\n', 'M78 (path update)',    Lpa/ns, Tpa/ns);
fprintf('   phrases c(n) = %d,  LZ bound n/log_A(n) = %.0f\n', ...
        round(NP/ns), n / (log(n)/log(A)));
fprintf('   [H = %.4f; trie pays a context-reset penalty -- quantified here]\n\n', H);

% ---------------- V2.4 ---------------------------------------------------
fprintf('[V2.4] Order-weight profiles W_k on real counts (n = 4000)\n');
n = 4000;
x = gen_markov(n, A, delta, 42);
bset = [0.5, H, 1.6];
lbl  = {'beta = 0.5 (< H)', 'beta = H       ', 'beta = 1.6 (> H)'};
% average profiles over the last 200 time steps
Kmax = 0;  prof = zeros(3, 64);
for t = (n-199):n
  ell = matchlens(x, t);
  L   = max(ell);
  Kmax = max(Kmax, L);
  h = zeros(1, L+1);
  for k = 0:L
    h(k+1) = sum(ell == k);
  end
  for b = 1:3
    w = h .* exp(bset(b) * ((0:L) - L));       % overflow-safe
    w = w / sum(w);
    prof(b, 1:L+1) = prof(b, 1:L+1) + w;
  end
end
prof = prof / 200;
fprintf('   %-18s', 'k =');
fprintf(' %6d', 0:Kmax); fprintf('\n');
for b = 1:3
  fprintf('   %-18s', lbl{b});
  fprintf(' %6.3f', prof(b, 1:Kmax+1)); fprintf('\n');
end
fprintf('   [expected: geometric decay / near-flat / concentration at top]\n');
