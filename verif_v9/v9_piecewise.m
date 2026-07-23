%V9_PIECEWISE  Module V9.2: piecewise-stationary testbed.
%
%   Source: circulant chain with delta alternating 0.10 / 0.50 every
%   Tseg = 2500 symbols (n = 20000: 8 segments, 7 switches).  The
%   per-step truth P_{delta(t)}(.|x_t) is known, so the excess is the
%   exact conditional KL (Rao-Blackwellized).
%
%   Pre-registered: global experts (M78R, bigram) hit the
%   regime-averaging floor  E_reg KL(P_delta || P_avg) = 0.102 nats
%   (occurrence-averaged successor prob 0.7); windowed bigram
%   (W = 1000) ~ 0.03-0.05; Bayes mixture tracks the best; fixed
%   share helps only if the expert ranking flips between regimes.

clear all; more off;
A = 4;  nu = 1e-3;  Kcap = 24;  seeds = [1 2 3];
ds = [0.10 0.50];  Tseg = 2500;  n = 20000;  W = 1000;
Pd = zeros(2, A, A);
for r = 1:2
  P = (ds(r)/(A-1)) * ones(A);
  for i = 1:A, P(i, mod(i,A)+1) = 1 - ds(r); end
  Pd(r,:,:) = P;
end
t0w = 300;  te = 2501;  alfs = 1e-3;
labs = {'M78R-self', 'bigram-gl', 'bigram-W', 'Bayes-mix', 'fix-share', 'sig-gate'};
Href = 0.5*(0.4349 + 1.2425);
X = zeros(2, 6);  cntr = zeros(1,2);  Lc3 = zeros(1,3);  Lmix = 0;  Lfs = 0;
for sd = seeds
  rand('state', 66000 + sd);
  x = zeros(1,n);  x(1) = ceil(rand()*A);
  for t = 2:n
    r = 1 + (mod(floor((t-1)/Tseg), 2) == 1);
    row = squeeze(Pd(r, x(t-1), :))';
    x(t) = 1 + sum(rand() > cumsum(row(1:A-1)));
  end
  child = zeros(n+2, A);  cntm = zeros(n+2, A);  nn = 1;  pcur = 1;
  nodes = zeros(1, Kcap+1);  nodes(1) = 1;
  Cb = zeros(A, A);  Cw = zeros(A, A);
  Lcum = zeros(1,3);  v = ones(1,3)/3;
  for t = 1:(n-1)
    a = x(t);
    live = find(nodes > 0);
    for j = 1:numel(live), cntm(nodes(live(j)), a) = cntm(nodes(live(j)), a) + 1; end
    c = child(pcur, a);
    if c > 0, pcur = c; else, nn = nn+1; child(pcur,a) = nn; pcur = 1; end
    for k = min(Kcap,t):-1:1
      if nodes(k) > 0, nodes(k+1) = child(nodes(k), a); else, nodes(k+1) = 0; end
    end
    if t > 1
      Cb(x(t-1), a) = Cb(x(t-1), a) + 1;
      Cw(x(t-1), a) = Cw(x(t-1), a) + 1;
      if t > W + 1, Cw(x(t-W-1), x(t-W)) = Cw(x(t-W-1), x(t-W)) - 1; end
    end
    if t < t0w, continue; end
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
    p3 = (Cw(a,:) + 0.5) / (sum(Cw(a,:)) + A/2);
    P3 = [p1; p2; p3];
    w  = exp(-(Lcum - min(Lcum)));  w = w/sum(w);
    pm = w * P3;
    pf = v * P3;
    lam = 1 / (1 + exp(-(d - log(t)/Href)));
    pg = lam*p1 + (1-lam)*p3;
    r = 1 + (mod(floor(t/Tseg), 2) == 1);      % regime of step t+1
    pv = squeeze(Pd(r, a, :))';
    if t >= te
      ex = [sum(pv.*log(pv./p1)) sum(pv.*log(pv./p2)) sum(pv.*log(pv./p3)) ...
            sum(pv.*log(pv./pm)) sum(pv.*log(pv./pf)) sum(pv.*log(pv./pg))];
      X(r,:) = X(r,:) + ex;  cntr(r) = cntr(r) + 1;
    end
    % full-horizon cumulatives (the theorem's statement) from t0w on
    Lmix = Lmix - log(pm(x(t+1)));
    Lfs  = Lfs  - log(pf(x(t+1)));
    l = -log(P3(:, x(t+1)))';
    Lcum = Lcum + l;
    Lc3 = Lc3 + l;
    v = v .* exp(-l);  v = v/sum(v);
    v = (1-alfs)*v + alfs/3;
  end
end
fprintf('[V9.2] Piecewise testbed (delta 0.10/0.50, Tseg = %d)\n', Tseg);
fprintf('   %-10s | %9s %9s | %9s\n', 'expert', 'reg d=.10', 'reg d=.50', 'overall');
for i = 1:6
  fprintf('   %-10s | %9.4f %9.4f | %9.4f\n', labs{i}, X(1,i)/cntr(1), ...
          X(2,i)/cntr(2), (X(1,i)+X(2,i))/sum(cntr));
end
fprintf('   [floor prediction for global experts: 0.116/.087 -> 0.102 avg]\n');
fprintf('   oracle (full horizon, per seed): L_mix - min_i L_i = %.3f (bound 3*log3 = 3.30)\n', ...
        Lmix - min(Lc3));
fprintf('   tracking gain: (L_mix - L_fs)/steps = %.4f nats/step\n', ...
        (Lmix - Lfs)/sum(cntr));
