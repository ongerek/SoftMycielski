%RUNME  Module V1: numerical verification for the soft-Mycielski project.
%
%   V1.1  Abel-summation identity (Proposition 1)      -> machine precision
%   V1.2  Zero-temperature limit (Theorem 1)           -> e^{-beta} rate
%   V1.3  Temperature-entropy law (informs Theorem 2)  -> beta* vs H table
%   V1.4  Match-length law L_t ~ log(t)/H (Wyner-Ziv / Ornstein-Weiss)
%   V1.5  Order-2 mechanism test for the phase-transition claim
%
%   Toolbox-free, GNU Octave compatible.  Runtime: a few minutes.


% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
clear all; more off;
A  = 4;
nu = 1e-3;

fprintf('==================================================================\n');
fprintf(' Module V1 -- soft-Mycielski verification  (A = %d, nu = %.0e)\n', A, nu);
fprintf('==================================================================\n\n');

% ---------- V1.1 Abel identity ------------------------------------------
fprintf('[V1.1] Abel-summation identity, 3 chains, beta = 1.5\n');
maxerr = 0;
for s = 1:3
  x = gen_markov(2000, A, 0.15, 100 + s);
  relerr = abel_check(x, 2000, 1.5, A);
  fprintf('   chain %d: max relative error = %.3e\n', s, relerr);
  maxerr = max(maxerr, relerr);
end
if maxerr < 1e-10
  fprintf('   PASS (identity holds to machine precision)\n\n');
else
  fprintf('   FAIL\n\n');
end

% ---------- V1.2 zero-temperature limit ---------------------------------
fprintf('[V1.2] Zero-temperature limit  ||p_beta - p_inf||_1  vs  e^{-beta}\n');
x    = gen_markov(2000, A, 0.15, 7);
ell  = matchlens(x, 2000);
cont = x(1:2000);
pinf = soft_myc_dist(ell, cont, Inf, A, nu);
for beta = [5 10 15 20 25 30]
  pb = soft_myc_dist(ell, cont, beta, A, nu);
  fprintf('   beta = %4.1f :  L1 distance = %.3e   (e^{-beta} = %.3e)\n', ...
          beta, sum(abs(pb - pinf)), exp(-beta));
end
fprintf('\n');

% ---------- V1.3 temperature-entropy law --------------------------------
fprintf('[V1.3] Temperature-entropy law: beta* vs H  (n = 2500, eval half)\n');
deltas = [0.02 0.05 0.10 0.20 0.35 0.50];
betas  = [0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 5 6 8 10 12 Inf];
nfin   = sum(isfinite(betas));
seeds  = [1 2];
n  = 2500;  t0 = 1250;

fprintf('   %6s %8s | %8s %10s | %10s %10s %10s\n', ...
        'delta', 'H', 'beta*', 'beta*/H', 'loss(b*)', 'loss(inf)', 'loss(0.25)');
res = zeros(numel(deltas), 6);
Lall = [];  Hall = [];
for d = 1:numel(deltas)
  [~, ~, H] = gen_markov(10, A, deltas(d));
  lossacc = zeros(1, numel(betas));
  for s = seeds
    x = gen_markov(n, A, deltas(d), 1000*d + s);
    [loss, Lrec] = beta_sweep(x, betas, A, nu, t0);
    lossacc = lossacc + loss;
    Lall = [Lall, Lrec(2, end-99:end)];          %#ok<AGROW>
    Hall = [Hall, H * ones(1, 100) ./ log(Lrec(1, end-99:end))]; %#ok<AGROW>
  end
  loss = lossacc / numel(seeds);
  [lmin, imin] = min(loss(1:nfin));
  % quadratic refinement of beta* on the three points around the minimum
  if imin > 1 && imin < nfin
    b3 = betas(imin-1:imin+1);  l3 = loss(imin-1:imin+1);
    c  = polyfit(b3, l3, 2);
    bstar = -c(2) / (2*c(1));
    if bstar < b3(1) || bstar > b3(3), bstar = betas(imin); end
  else
    bstar = betas(imin);
  end
  res(d, :) = [deltas(d), H, bstar, bstar/H, lmin, loss(end)];
  fprintf('   %6.2f %8.4f | %8.3f %10.3f | %10.4f %10.4f %10.4f\n', ...
          deltas(d), H, bstar, bstar/H, lmin, loss(end), loss(1));
end
fprintf('\n   (loss in nats/symbol; loss(inf) = hard Mycielski,');
fprintf(' loss(0.25) = near-marginal)\n\n');

% ---------- V1.4 match-length law ---------------------------------------
fprintf('[V1.4] Match-length law:  mean of  L_t * H / log(t)  over late t\n');
ratio = Lall .* Hall;                 % L_t * H / log(t)
fprintf('   mean ratio = %.3f   (std %.3f, N = %d)   [theory -> 1]\n\n', ...
        mean(ratio), std(ratio), numel(ratio));

fprintf('==================================================================\n');
fprintf(' Summary table (delta, H, beta*, beta*/H, loss*, loss_hard):\n');
disp(res);
save('-text', 'results_v1.txt', 'res');

% ---------- V1.5 order-2 mechanism test ---------------------------------
order2_test;

fprintf(' Done.\n');
