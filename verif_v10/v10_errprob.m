%V10_ERRPROB  Module V10.1: error-probability (0-1 loss) companion
%             to the log-loss tables.
%
%   For a predictor with distribution p(.|x_t), the plug-in decision
%   is argmax_a p(a); its error probability, measured against the
%   TRUE conditional P(.|x_t) of the testbed, is
%
%       err = E_t [ 1 - P(argmax_a p(a) | x_t) ],
%
%   and the Bayes floor is  err* = E_t [ 1 - max_a P(a|x_t) ].  We
%   report the EXCESS error  (err - err*)  in the same spirit as the
%   log-loss tables report excess over H.  Computed Rao-Blackwellized
%   (against P, not the realized symbol) so the estimator has no
%   next-symbol sampling noise; ties in argmax broken uniformly.
%
%   Covers the three settings of Tables III/IV/VI:
%     (a) annealing sweep beta*(n), delta = 0.35;
%     (b) fixed 0.5 / fixed 1.0 / oracle / self-critical, three sources;
%     (c) oracle vs self across n (for the paired plot).
%   Saves per-n records to errprob_data.mat for V10_PLOTS.

clear all; more off;
A = 4;  nu = 1e-3;

function [errx, err, bayes] = err_excess(x, P, betas, A, nu, t0, selfcol)
  % returns excess error (err-bayes) per beta column, plus raw err, bayes
  n = numel(x);  nb = numel(betas);
  err = zeros(1, nb + selfcol);  bayes = 0;  cnt = 0;
  for t = t0:(n-1)
    ell  = matchlens(x, t);
    cont = x(1:t);
    Pt   = P(x(t), :);
    [pm, ~] = max(Pt);
    bayes = bayes + (1 - pm);
    for b = 1:nb
      p = soft_myc_dist(ell, cont, betas(b), A, nu);
      err(b) = err(b) + errprob_one(p, Pt);
    end
    if selfcol
      L = max(ell);
      p = soft_myc_dist(ell, cont, log(t)/(L+1), A, nu);
      err(nb+1) = err(nb+1) + errprob_one(p, Pt);
    end
    cnt = cnt + 1;
  end
  err = err / cnt;  bayes = bayes / cnt;
  errx = err - bayes;
end

function e = errprob_one(p, Pt)
  % Rao-Blackwellized 0-1 loss of the plug-in decision argmax p
  m = max(p);
  win = find(p >= m - 1e-12);          % tie set
  % expected error averaging ties uniformly:  1 - mean_{a in win} P(a)
  e = 1 - mean(Pt(win));
end

% ---------------- (a) annealing sweep, delta = 0.35 ----------------------
fprintf('[V10.1a] Error-probability along the annealing sweep (delta=0.35)\n');
delta = 0.35;  [~, P, H] = gen_markov(10, A, delta);
ns    = [1000 2000 4000 8000 16000];
bstar = [0.714 0.729 0.747 0.760 0.776];
fprintf('  %6s | %8s %8s %9s | %9s\n', 'n', 'err(b*)', 'bayes', 'excess', 'exc(logloss)');
excll = [0.1073 0.0856 0.0904 0.0975 0.0855];  % from tab:selfn oracle col (ref)
Erec = zeros(numel(ns), 3);
for i = 1:numel(ns)
  n = ns(i);  t0 = n - min(4000, floor(n/2)) + 1;
  acc = 0;  accb = 0;  cnt = 0;
  for s = 1:3
    x = gen_markov(n, A, delta, 100000 + 10*i + s);
    [ex, er, by] = err_excess(x, P, bstar(i), A, nu, t0, 0);
    acc = acc + er;  accb = accb + by;  cnt = cnt + 1;
  end
  er = acc/cnt;  by = accb/cnt;
  Erec(i,:) = [n, er - by, by];
  fprintf('  %6d | %8.4f %8.4f %9.4f | %9.4f\n', n, er, by, er-by, excll(i));
end

% ---------------- (b) self-tune comparison, three sources ----------------
fprintf('\n[V10.1b] Error-probability: fixed / oracle / self-critical\n');
srcs = {{1, 0.10}, {1, 0.35}, {2, 0.20}};   % {order, delta}
bstar_b = [0.42 0.73 0.60];                  % approx oracle beta* per source
fprintf('  %10s | %7s %7s %7s %7s | %7s\n', 'source', 'fix.5', 'fix1', ...
        'oracle', 'self', 'bayes');
for is = 1:3
  ordr = srcs{is}{1};  delta = srcs{is}{2};
  if ordr == 1
    [~, P, H] = gen_markov(10, A, delta);  gen = @(nn,sd) gen_markov(nn,A,delta,sd);
    P2 = []; use2 = false;
  else
    [~, H] = gen_markov2(10, A, delta); gen = @(nn,sd) gen_markov2(nn,A,delta,sd);
    % build P2(i,j,a) = Pr(x=a | prev2=i, prev1=j) from the successor rule
    P2 = zeros(A, A, A);
    for i = 1:A
      for j = 1:A
        s = mod(i + j, A) + 1;
        P2(i, j, :) = delta/(A-1);
        P2(i, j, s) = 1 - delta;
      end
    end
    use2 = true;
  end
  n = 4000;  t0 = 2001;
  betas = [0.5 1.0 bstar_b(is)];
  acc = zeros(1,4); accb = 0; cnt = 0;
  for s = 1:3
    x = gen(n, 110000 + 100*is + s);
    if use2
      [ex, er, by] = err_excess2(x, P2, betas, A, nu, t0);
    else
      [ex, er, by] = err_excess(x, P, betas, A, nu, t0, 1);
    end
    acc = acc + er; accb = accb + by; cnt = cnt + 1;
  end
  er = acc/cnt; by = accb/cnt;
  nm = sprintf('o%d d=%.2f', ordr, delta);
  fprintf('  %10s | %7.4f %7.4f %7.4f %7.4f | %7.4f\n', nm, ...
          er(1), er(2), er(3), er(4), by);
end

% ---------------- (c) oracle vs self across n (for plot) -----------------
fprintf('\n[V10.1c] Oracle vs self error excess across n (delta=0.35)\n');
Crec = zeros(numel(ns), 3);
fprintf('  %6s | %9s %9s %6s\n', 'n', 'exc(orac)', 'exc(self)', 'ratio');
for i = 1:numel(ns)
  n = ns(i);  t0 = n - min(4000, floor(n/2)) + 1;
  ao = 0; as_ = 0; cnt = 0;
  for s = 1:3
    x = gen_markov(n, A, delta, 120000 + 10*i + s);
    [ex, er, by] = err_excess(x, P, bstar(i), A, nu, t0, 1);
    ao = ao + (er(1) - by);  as_ = as_ + (er(2) - by);  cnt = cnt + 1;
  end
  ao = ao/cnt; as_ = as_/cnt;
  Crec(i,:) = [n, ao, as_];
  fprintf('  %6d | %9.4f %9.4f %6.2f\n', n, ao, as_, ao/max(as_,1e-6));
end

save('-mat', 'errprob_data.mat', 'Erec', 'Crec');
fprintf('\n  [saved errprob_data.mat for V10_PLOTS]\n');
