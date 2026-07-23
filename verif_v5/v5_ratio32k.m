%V5_RATIO32K  Module V5.4: adaptivity ratio at n = 32000.
%
%   Separation theory: fixed-beta excess ~ 1/sqrt(log n) (head atom,
%   free-walk linear small-max law) vs self-critical ~ 1/log n
%   (bridge, quadratic law), so the excess ratio should keep growing
%   ~ sqrt(log n) up to logs.  V3.1 measured 1.58 -> 2.40 over
%   n = 1000 -> 16000; this module adds n = 32000.

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;  n = 32000;  seeds = [1 2];
[~, ~, H] = gen_markov(10, A, delta);
t0 = n - 400 + 1;
bgrid = [0.76 0.78 0.80 0.82];      % local sweep around extrapolated beta*

Lb = zeros(1, numel(bgrid));  Ls = 0;  cnt = 0;
for s = seeds
  x = gen_markov(n, A, delta, 5000*6 + s);   % continue V2.1 seed pattern
  for t = t0:(n-1)
    ell  = matchlens(x, t);
    cont = x(1:t);
    Lt   = max(ell);
    for j = 1:numel(bgrid)
      p = soft_myc_dist(ell, cont, bgrid(j), A, nu);
      Lb(j) = Lb(j) - log(p(x(t+1)));
    end
    p  = soft_myc_dist(ell, cont, log(t)/(Lt+1), A, nu);
    Ls = Ls - log(p(x(t+1)));
    cnt = cnt + 1;
  end
end
Lb = Lb/cnt;  Ls = Ls/cnt;
[lo, jo] = min(Lb);
fprintf('[V5.4] n = 32000, delta = %.2f, H = %.4f  (%d eval steps)\n', ...
        delta, H, cnt);
fprintf('   sweep: ');  fprintf('%.2f->%.4f  ', [bgrid; Lb]);
fprintf('\n   oracle beta* ~ %.2f: excess %.4f;  self: excess %.4f;  ratio %.2f\n', ...
        bgrid(jo), lo - H, Ls - H, (lo - H)/(Ls - H));
fprintf('   [V3.1 ratios: 1.58 1.92 1.63 1.79 2.40 for n = 1k..16k]\n');
