%V3_ANNEAL_SELF  Modules V3.1 and V3.2: average-case adaptivity.
%
%   V3.1  Self-critical schedule beta_t = log t/(L_t+1) vs the fixed
%         near-oracle temperatures beta*(n) of Module V2.1, on the
%         SAME sequences and windows.  Prediction of the profile
%         model: both excesses are Theta(1/log n) with the self-tuned
%         constant strictly smaller, i.e. an O(1) relative gain.
%
%   V3.2  Losses binned by the realized match length L_t.  Prediction
%         (supercriticality mechanism): a fixed temperature is LOCALLY
%         SUPERCRITICAL exactly on long-match steps (local entropy
%         log t/(L_t+1) below beta), so the self-tuned schedule should
%         win most where L_t is large.

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;
[~, ~, H] = gen_markov(10, A, delta);

% ---------------- V3.1 ---------------------------------------------------
ns     = [1000 2000 4000 8000 16000];
bstars = [0.714 0.729 0.747 0.760 0.776];      % measured in Module V2.1
seeds  = [1 2 3];
fprintf('[V3.1] Self-critical vs oracle over n  (delta=%.2f, H=%.4f)\n', ...
        delta, H);
fprintf('   %7s | %10s %10s | %10s %10s | %7s\n', 'n', ...
        'exc(orac)', 'exc(self)', 'loss(orac)', 'loss(self)', 'ratio');
for i = 1:numel(ns)
  n  = ns(i);
  t0 = n - min(400, floor(n/5)) + 1;
  Lo = 0;  Ls = 0;  cnt = 0;
  for s = seeds
    x = gen_markov(n, A, delta, 5000*i + s);    % same seeds as V2.1
    for t = t0:(n-1)
      ell  = matchlens(x, t);
      cont = x(1:t);
      Lt   = max(ell);
      po = soft_myc_dist(ell, cont, bstars(i), A, nu);
      ps = soft_myc_dist(ell, cont, log(t)/(Lt+1), A, nu);
      Lo = Lo - log(po(x(t+1)));
      Ls = Ls - log(ps(x(t+1)));
      cnt = cnt + 1;
    end
  end
  Lo = Lo/cnt;  Ls = Ls/cnt;
  fprintf('   %7d | %10.4f %10.4f | %10.4f %10.4f | %7.2f\n', ...
          n, Lo - H, Ls - H, Lo, Ls, (Lo - H)/(Ls - H));
end
fprintf('\n');

% ---------------- V3.2 ---------------------------------------------------
fprintf('[V3.2] Loss binned by realized L_t  (n = 4000, beta_fix = 0.731)\n');
n = 4000;  t0 = 2001;  bfix = 0.731;
edges = [0 4 5 6 7 8 9 Inf];    % bins: <=4, 5, 6, 7, 8, 9, >=10
nb = numel(edges) - 1;
Bo = zeros(1, nb);  Bs = zeros(1, nb);  Bc = zeros(1, nb);
for s = seeds
  x = gen_markov(n, A, delta, 9000 + s);
  for t = t0:(n-1)
    ell  = matchlens(x, t);
    cont = x(1:t);
    Lt   = max(ell);
    b    = min(nb, sum(Lt >= edges(1:end-1)) );
    po = soft_myc_dist(ell, cont, bfix, A, nu);
    ps = soft_myc_dist(ell, cont, log(t)/(Lt+1), A, nu);
    Bo(b) = Bo(b) - log(po(x(t+1)));
    Bs(b) = Bs(b) - log(ps(x(t+1)));
    Bc(b) = Bc(b) + 1;
  end
end
lab = {'<=4','5','6','7','8','9','>=10'};
fprintf('   %6s %8s | %10s %10s | %10s\n', ...
        'L_t', 'count', 'fixed', 'self', 'gain');
for b = 1:nb
  if Bc(b) > 0
    fprintf('   %6s %8d | %10.4f %10.4f | %+10.4f\n', ...
            lab{b}, Bc(b), Bo(b)/Bc(b), Bs(b)/Bc(b), (Bo(b)-Bs(b))/Bc(b));
  end
end
fprintf('   [gain > 0 means self-tuned wins in that bin]\n');
