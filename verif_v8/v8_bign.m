%V8_BIGN  Module V8.2: real-system normalized excesses at large n,
%         to set against the idealized-model constants of V8.1.
%
%   c2(n) = log(n) * excess(self);  c1cr(n) = sqrt(log n) * excess(H);
%   c1(n) = sqrt(log n) * excess(oracle over local sweep).
%   Model comparison at matched L ~ log(n)/H.

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;  seeds = [1 2 3];
[~, PT, H] = gen_markov(10, A, delta);
ns = [100000 250000];
sweep = [0.79 0.81 0.83];

% Rao-Blackwellized estimator (mandatory at large n: the excess
% ~ c/log n is smaller than the composition noise of raw windows):
% per-step conditional KL against the known source row.
for i = 1:numel(ns)
  n = ns(i);  t0 = n - 400 + 1;
  Lb = zeros(1, numel(sweep));  Lh = 0;  Ls = 0;  cnt = 0;
  for s = seeds
    x = gen_markov(n, A, delta, 99000 + 10*i + s);
    for t = t0:(n-1)
      ell  = matchlens(x, t);
      cont = x(1:t);
      Lt   = max(ell);
      Pt   = PT(x(t), :);
      for j = 1:numel(sweep)
        p = soft_myc_dist(ell, cont, sweep(j), A, nu);
        Lb(j) = Lb(j) + sum(Pt .* log(Pt ./ p(:)'));
      end
      p  = soft_myc_dist(ell, cont, H, A, nu);
      Lh = Lh + sum(Pt .* log(Pt ./ p(:)'));
      p  = soft_myc_dist(ell, cont, log(t)/(Lt+1), A, nu);
      Ls = Ls + sum(Pt .* log(Pt ./ p(:)'));
      cnt = cnt + 1;
    end
  end
  Lb = Lb/cnt;  Lh = Lh/cnt;  Ls = Ls/cnt;
  [eo, jo] = min(Lb);
  fprintf('[V8.2] n = %d (%d steps, log n = %.2f, L ~ %.1f)\n', ...
          n, cnt, log(n), log(n)/H);
  fprintf('   exc: self %.4f  H %.4f  oracle %.4f (beta* ~ %.2f)\n', ...
          Ls, Lh, eo, sweep(jo));
  fprintf('   c2(n) = %.3f   c1cr(n) = %.3f   c1(n) = %.3f\n', ...
          Ls*log(n), Lh*sqrt(log(n)), eo*sqrt(log(n)));
end
fprintf('   [model limits: c2 = 0.16(3), c1cr = 0.222(3), c1 = 0.10(1)]\n');
