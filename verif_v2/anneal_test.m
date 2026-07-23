%ANNEAL_TEST  Module V2.1: critical-annealing law  beta*(n) -> H^- .
%
%   Theory (critical annealing): for a memory-r source the log-loss
%   optimal temperature approaches the phase boundary from below,
%   beta*(n) = H - Theta(1/log n), with excess loss Theta(1/log n).
%   This test measures beta*(n) and excess(n) over a dyadic n-grid.

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;
[~, ~, H] = gen_markov(10, A, delta);
betas = [0.40 0.50 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.10 1.20 1.40];
seeds = [1 2 3];
ns    = [1000 2000 4000 8000 16000];

fprintf('[V2.1] Critical annealing:  delta = %.2f, H = %.4f nats\n', delta, H);
fprintf('   %7s | %8s %10s | %10s %10s\n', ...
        'n', 'beta*', 'H-beta*', 'loss(b*)', 'excess');
for i = 1:numel(ns)
  n  = ns(i);
  t0 = n - min(400, floor(n/5)) + 1;
  lossacc = zeros(1, numel(betas));
  for s = seeds
    x = gen_markov(n, A, delta, 5000*i + s);
    lossacc = lossacc + beta_sweep(x, betas, A, nu, t0);
  end
  loss = lossacc / numel(seeds);
  [lmin, imin] = min(loss);
  if imin > 1 && imin < numel(betas)
    b3 = betas(imin-1:imin+1);  l3 = loss(imin-1:imin+1);
    c  = polyfit(b3, l3, 2);
    bstar = -c(2)/(2*c(1));
    if bstar < b3(1) || bstar > b3(3), bstar = betas(imin); end
    lmin = polyval(c, bstar);
  else
    bstar = betas(imin);
  end
  fprintf('   %7d | %8.3f %10.3f | %10.4f %10.4f\n', ...
          n, bstar, H - bstar, lmin, lmin - H);
end
fprintf('   [theory: H-beta* and excess both decrease ~ 1/log n]\n');
