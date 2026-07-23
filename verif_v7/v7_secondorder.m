%V7_SECONDORDER  Module V7.1: the second-order misallocation
%                principle behind the polylog-free rates.
%
%   Conditionally on the past, E[q(X)/P(X)] = 1 for ANY
%   past-measurable q, so mixing weight w of a wrong predictor into
%   the true conditional costs
%       excess(w) = E KL( P || (1-w)P + w q )  ~  (w^2/2) chi2t,
%   chi2t = sum_a (P(a)-q(a))^2 / P(a)   (note: /P, not /q),
%   QUADRATIC in w uniformly on [0, 3/4] -- not linear.
%
%   Testbed: P = true conditional of the circulant chain, q = the
%   running empirical marginal.  Exact predictions (delta = 0.35,
%   q -> uniform):  excess/w^2 = 0.345, 0.340, 0.333, 0.324, 0.325
%   at w = 0.05, 0.1, 0.2, 0.4, 0.6  (chi2t/2 = 0.352 at w -> 0),
%   while excess/w varies by a factor ~6 over the same grid.

clear all; more off;
A = 4;  delta = 0.35;  n = 4000;  seeds = [1 2 3];
[~, P, H] = gen_markov(10, A, delta);
ws = [0.05 0.10 0.20 0.40 0.60];
t0 = 2001;

% Rao-Blackwellized: the lemma concerns the CONDITIONAL KL, so
% compute it exactly per step (the raw Monte-Carlo estimator is
% ill-conditioned at small w precisely BECAUSE the linear term
% cancels: composition noise in the realized transition mix then
% dominates the w^2 signal).
acc = zeros(1, numel(ws));  cnt = 0;
for s = seeds
  x = gen_markov(n, A, delta, 81000 + s);
  cnts = zeros(1, A);
  for t = 1:(n-1)
    cnts(x(t)) = cnts(x(t)) + 1;
    if t < t0, continue; end
    q0 = cnts / t;
    Pt = P(x(t), :);
    for j = 1:numel(ws)
      w = ws(j);
      acc(j) = acc(j) + sum(Pt .* log(Pt ./ ((1-w)*Pt + w*q0)));
    end
    cnt = cnt + 1;
  end
end
exc = acc/cnt;
fprintf('[V7.1] Second-order misallocation (delta=%.2f, chi2t/2 = 0.352)\n', delta);
fprintf('   %6s | %8s | %9s | %9s\n', 'w', 'excess', 'exc/w^2', 'exc/w');
pred = [0.345 0.340 0.333 0.324 0.325];
for j = 1:numel(ws)
  e = exc(j);
  fprintf('   %6.2f | %8.5f | %9.3f | %9.3f\n', ws(j), e, e/ws(j)^2, e/ws(j));
end
fprintf('   [predicted exc/w^2: ');  fprintf('%.3f ', pred);
fprintf(' -- flat = quadratic law]\n');
