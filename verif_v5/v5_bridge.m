%V5_BRIDGE  Modules V5.1-V5.3: the bridge representation of the
%           self-critical tilt, tested on the actual predictor.
%
%   Define D_k = kH - log(1/P_k) for the realized length-k suffix
%   (computable exactly on the testbed).  Theory:
%
%   V5.1 (endpoint pinning / crossing).  The Ornstein-Wyner-Ziv
%        crossing gives  log(t P_{L_t}) = O_P(1), and the endpoint
%        statistic  G = log(m_{L_t}/m_0)  of the tilted profile is
%        O(1)-concentrated under beta^self = log t/(L_t+1), while
%        under a fixed beta it has spread sigma_v sqrt(L) ~ 2.7.
%
%   V5.2 (pinned atom).  The head weight W_0 under the self tilt
%        scales like 1/log t (bridge, quadratic small-max law),
%        against 1/sqrt(log t) under fixed beta (free walk, linear).
%
%   V5.3 (small-max law shape).  Empirical CDF of the profile-max
%        gap  M = max_k log(m_k/m_0):  P[M<=s] approx linear in s
%        for fixed beta, approx quadratic (convex) for self tilt.

clear all; more off;
A = 4;  delta = 0.35;  seeds = [1 2 3];
[~, ~, H] = gen_markov(10, A, delta);
sv = sqrt((1-delta)*log(1-delta)^2 + delta*log(delta/(A-1))^2 - H^2);
fprintf('source: delta=%.2f, H=%.4f, sigma_v=%.4f\n\n', delta, H, sv);

ns     = [1000 4000 16000];
bfixs  = [0.714 0.747 0.776];      % oracle beta*(n) from Module V2.1
smax_s = [0.5 1.0 2.0];

for ii = 1:numel(ns)
  n = ns(ii);  bfix = bfixs(ii);  t0 = floor(n/2) + 1;
  % accumulators: [fix, self] x {G-mean, G-sq, W0, W0^2, atom, cdf(3)}
  Gm = zeros(1,2); Gs = zeros(1,2); W1 = zeros(1,2); W2 = zeros(1,2);
  At = zeros(1,2); CD = zeros(2,3); Xm = 0; Xs = 0; cnt = 0;
  for s = seeds
    x = gen_markov(n, A, delta, 41000 + 7*ii + s);
    % per-step surprisal of transitions:  surp(i) = -log p(x(i+1)|x(i))
    succ = (x(2:end) == mod(x(1:end-1), A) + 1);
    surp = -succ*log(1-delta) - (~succ)*log(delta/(A-1));
    S    = [0 cumsum(surp)];               % S(i) = sum surp(1..i-1)... S(j+1)=sum_{i<=j}
    for t = t0:(n-1)
      ell = matchlens(x, t);
      L   = max(ell);
      if L < 2, continue; end
      % log(1/P_L) of realized suffix x(t-L+1..t)
      lp  = log(A) + (S(t) - S(t-L+1));    % pi uniform + transitions in window
      Xm  = Xm + (log(t) - lp);  Xs = Xs + (log(t) - lp)^2;
      h   = histc(ell, 0:L);
      bs  = [bfix, log(t)/(L+1)];
      for j = 1:2
        w  = h .* exp(bs(j) * ((0:L) - L));
        w  = w / w(1);                     % relative to head
        G  = log(w(end));                  % endpoint statistic
        M  = log(max(w));                  % profile-max gap
        W0 = 1 / sum(w);
        Gm(j) = Gm(j) + G;   Gs(j) = Gs(j) + G^2;
        W1(j) = W1(j) + W0;  W2(j) = W2(j) + W0^2;
        At(j) = At(j) + (W0 > 0.3);
        CD(j,:) = CD(j,:) + (M <= smax_s);
      end
      cnt = cnt + 1;
    end
  end
  Gm = Gm/cnt; Gs = sqrt(max(Gs/cnt - Gm.^2, 0));
  W1 = W1/cnt; W2 = W2/cnt; At = At/cnt; CD = CD/cnt;
  Xm = Xm/cnt; Xs = sqrt(max(Xs/cnt - Xm^2, 0));
  fprintf('n = %5d  (log n = %.2f, sv*sqrt(L)~%.2f, beta_fix = %.3f)\n', ...
          n, log(n), sv*sqrt(log(n)/H), bfix);
  fprintf('  crossing  log(t P_L): mean %+.3f  std %.3f   [both O(1)]\n', Xm, Xs);
  fprintf('  %-6s | %8s %8s | %8s %9s %8s | %6s %6s %6s\n', 'tilt', ...
          'G mean', 'G std', 'E[W0]', 'E[W0^2]', 'P[>.3]', ...
          'P(M<.5)', 'P(<1)', 'P(<2)');
  lab = {'fixed', 'self'};
  for j = 1:2
    fprintf('  %-6s | %+8.3f %8.3f | %8.4f %9.4f %8.4f | %6.3f %6.3f %6.3f\n', ...
            lab{j}, Gm(j), Gs(j), W1(j), W2(j), At(j), CD(j,1), CD(j,2), CD(j,3));
  end
  fprintf(['  scaled: E[W0^2]*sqrt(logn) fix = %.3f,  E[W0^2]*logn self = %.3f\n\n'], ...
          W2(1)*sqrt(log(n)), W2(2)*log(n));
end
fprintf(['[predictions: G std ~ O(1) self vs ~sv*sqrt(L) fixed;\n' ...
         ' P(M<=s) ratios P(1)/P(2): ~0.5 fixed (linear), ~0.25 self (quadratic);\n' ...
         ' scaled atom columns roughly constant across n]\n']);
