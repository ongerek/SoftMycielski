%ORDER2_TEST  Module V1.5: mechanism test for the phase-transition claim.
%
%   For the order-2 source of GEN_MARKOV2 (order-1 statistics exactly
%   uniform), the smoothing phase beta < H mixes in useless low-order
%   estimates; the claim predicts beta* shifts markedly upward relative
%   to the order-1 source at the SAME entropy rate.

clear all; more off;
A  = 4;
nu = 1e-3;
betas = [0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 5 6 8 10 12 Inf];
nfin  = sum(isfinite(betas));
seeds = [1 2];
n = 2500;  t0 = 1250;
deltas = [0.20 0.35];

fprintf('[V1.5] Order-2 source (order-1 stats uniform): beta* vs H\n');
fprintf('   %6s %8s | %8s %10s | %10s %10s | %12s\n', ...
        'delta', 'H', 'beta*', 'beta*/H', 'loss(b*)', 'loss(inf)', 'loss(0.25)');
for d = 1:numel(deltas)
  [~, H] = gen_markov2(10, A, deltas(d));
  lossacc = zeros(1, numel(betas));
  for s = seeds
    x = gen_markov2(n, A, deltas(d), 2000*d + s);
    loss = beta_sweep(x, betas, A, nu, t0);
    lossacc = lossacc + loss;
  end
  loss = lossacc / numel(seeds);
  [lmin, imin] = min(loss(1:nfin));
  if imin > 1 && imin < nfin
    b3 = betas(imin-1:imin+1);  l3 = loss(imin-1:imin+1);
    c  = polyfit(b3, l3, 2);
    bstar = -c(2)/(2*c(1));
    if bstar < b3(1) || bstar > b3(3), bstar = betas(imin); end
  else
    bstar = betas(imin);
  end
  fprintf('   %6.2f %8.4f | %8.3f %10.3f | %10.4f %10.4f | %12.4f\n', ...
          deltas(d), H, bstar, bstar/H, lmin, loss(end), loss(1));
  fprintf('   full loss curve (finite betas):\n');
  for b = 1:nfin
    fprintf('      beta = %5.2f   loss = %.4f\n', betas(b), loss(b));
  end
end
fprintf('   (compare order-1 at same H: beta* = 0.605 (d=0.20), 0.731 (d=0.35))\n');
