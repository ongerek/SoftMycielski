%V9_HMM  Module V9.1: beyond-Markov testbed (hidden Markov source).
%
%   Source: 2 hidden states, sticky (p_stay = 0.95); emissions
%   e1 = [.7 .1 .1 .1], e2 = [.1 .3 .3 .3].  The observable process
%   has infinite Markov order; the exact predictive P(a|x_1^t) is
%   computable by the forward filter, giving (i) the entropy rate to
%   MC precision and (ii) a Rao-Blackwellized per-step excess
%   (conditional KL) for every expert.
%
%   Experts (all streaming, O(K) or O(A^2) per symbol):
%     E1  M78-R self-critical (the paper's streaming realization)
%     E2  bigram KT (add-1/2)
%   Combiners: Bayes mixture (Theorem: <= min + log K total),
%   fixed share (alpha = 1e-3), sigmoid gate eq:gate with
%   (alpha,gamma,c) = (1, 1/H, 0).

clear all; more off;
A = 4;  nu = 1e-3;  Kcap = 24;  seeds = [1 2 3];
T2 = [0.95 0.05; 0.05 0.95];
E  = [0.7 0.1 0.1 0.1; 0.1 0.3 0.3 0.3];

% ---- entropy rate by long forward run -----------------------------------
rand('state', 12345);
nH = 300000;  s = 1;  pi_f = [0.5 0.5];  accH = 0;
for t = 1:nH
  s = 1 + (rand() > T2(s,1));
  a = 1 + sum(rand() > cumsum(E(s,1:A-1)));
  pr = pi_f * T2;  pv = pr * E;               % predictive over symbols
  if t > 1000, accH = accH - log(pv(a)); end
  po = (pr .* E(:,a)');  pi_f = po / sum(po);
end
Hh = accH / (nH - 1000);
fprintf('[V9.1] HMM testbed: H = %.4f (forward MC, n = %d)\n', Hh, nH);

n = 30000;  t0w = 500;  te = 5000;            % mixture warmup; eval start
alfs = 1e-3;
labs = {'M78R-self', 'bigramKT', 'Bayes-mix', 'fix-share', 'sig-gate'};
X2 = zeros(1,5);  Lc2 = zeros(1,2);  cnt = 0;  Lmix = 0;
for sd = seeds
  rand('state', 55000 + sd);
  % generate
  x = zeros(1,n);  s = 1;
  for t = 1:n
    s = 1 + (rand() > T2(s,1));
    x(t) = 1 + sum(rand() > cumsum(E(s,1:A-1)));
  end
  % expert/combiner state
  child = zeros(n+2, A);  cntm = zeros(n+2, A);  nn = 1;  pcur = 1;
  nodes = zeros(1, Kcap+1);  nodes(1) = 1;
  Cb = zeros(A, A);
  Lcum = zeros(1,2);  v = ones(1,2)/2;
  pi_f = [0.5 0.5];
  for t = 1:(n-1)
    a = x(t);
    % consume x(t): M78R updates
    live = find(nodes > 0);
    for j = 1:numel(live), cntm(nodes(live(j)), a) = cntm(nodes(live(j)), a) + 1; end
    c = child(pcur, a);
    if c > 0, pcur = c; else, nn = nn+1; child(pcur,a) = nn; pcur = 1; end
    for k = min(Kcap,t):-1:1
      if nodes(k) > 0, nodes(k+1) = child(nodes(k), a); else, nodes(k+1) = 0; end
    end
    if t > 1, Cb(x(t-1), a) = Cb(x(t-1), a) + 1; end
    % forward filter through x(t)
    pr = pi_f * T2;  po = pr .* E(:,a)';  pi_f = po / sum(po);
    if t < t0w, continue; end
    % ---- predict x(t+1) ----
    live = find(nodes > 0);  d = live(end) - 1;
    b = log(t) / (d + 1);  fac = 1 - exp(-b);
    p1 = exp(-b*d) * cntm(1,:);
    for j = 2:numel(live)
      k = live(j) - 1;
      p1 = p1 + fac * exp(b*(k-d)) * cntm(nodes(live(j)), :);
    end
    ssum = sum(p1);  if ssum > 0, p1 = p1/ssum; else, p1 = ones(1,A)/A; end
    p1 = (1-nu)*p1 + nu/A;
    p2 = (Cb(a,:) + 0.5) / (sum(Cb(a,:)) + A/2);
    P2 = [p1; p2];
    w  = exp(-(Lcum - min(Lcum)));  w = w/sum(w);
    pm = w * P2;
    pf = v * P2;
    lam = 1 / (1 + exp(-(d - log(t)/Hh)));
    pg = lam*p1 + (1-lam)*p2;
    % truth predictive
    pr = pi_f * T2;  pv = pr * E;
    if t >= te
      X2(1) = X2(1) + sum(pv .* log(pv ./ p1));
      X2(2) = X2(2) + sum(pv .* log(pv ./ p2));
      X2(3) = X2(3) + sum(pv .* log(pv ./ pm));
      X2(4) = X2(4) + sum(pv .* log(pv ./ pf));
      X2(5) = X2(5) + sum(pv .* log(pv ./ pg));
      cnt = cnt + 1;
      Lmix = Lmix - log(pm(x(t+1)));
    end
    % realized-loss updates for combiners
    l = -log(P2(:, x(t+1)))';
    Lcum = Lcum + l;
    if t >= te, Lc2 = Lc2 + l; end
    v = v .* exp(-l);  v = v / sum(v);
    v = (1-alfs)*v + alfs/2;
  end
end
X2 = X2 / cnt;
fprintf('   RB excess (nats):');
for i = 1:5, fprintf('  %s %.4f', labs{i}, X2(i)); end
fprintf('\n   oracle check: L_mix - min_i L_i = %.3f  (bound log 2 = 0.693,\n', ...
        Lmix - min(Lc2));
fprintf('    per-step %.6f over %d eval steps)\n', (Lmix - min(Lc2))/cnt, cnt);
