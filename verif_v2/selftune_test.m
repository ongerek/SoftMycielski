%SELFTUNE_TEST  Module V2.2: self-critical temperature schedules.
%
%   The Ornstein-Weiss law  L_t/log t -> 1/H  lets the predictor set
%   its own temperature from its own match lengths:
%
%     schedule (a):  beta_t = log(t) / (L_t + 1)
%                    ~ H - H^2/log t          (marginal: eps*log t -> H^2)
%     schedule (b):  beta_t = log(t) / (L_t + sqrt(L_t + 1))
%                    ~ H - H^{3/2}/sqrt(log t) (safe:  eps*log t -> inf)
%
%   Compared against fixed betas and the oracle beta* from Module V1.

clear all; more off;
A = 4;  nu = 1e-3;  n = 4000;  t0 = 2001;  seeds = [1 2 3];

srcs = {  % {name, generator handle, delta, oracle beta* (Module V1)}
  {'order-1, d=0.10', @(nn,ss) gen_markov (nn, A, 0.10, ss), 0.10, 0.536},
  {'order-1, d=0.35', @(nn,ss) gen_markov (nn, A, 0.35, ss), 0.35, 0.731},
  {'order-2, d=0.20', @(nn,ss) gen_markov2(nn, A, 0.20, ss), 0.20, 0.731}
};

fprintf('[V2.2] Self-critical schedules  (n = %d, eval second half)\n', n);
fprintf('   %-18s %8s | %8s %8s %8s | %8s %8s\n', 'source', 'H', ...
        'b=0.5', 'b=1.0', 'oracle', 'self-a', 'self-b');
for is = 1:numel(srcs)
  name  = srcs{is}{1};  gen = srcs{is}{2};  bo = srcs{is}{4};
  if is <= 2, [~, ~, H] = gen_markov(10, A, srcs{is}{3});
  else,       [~, H]    = gen_markov2(10, A, srcs{is}{3});
  end
  L5 = zeros(1, 5);           % [fixed .5, fixed 1.0, oracle, self-a, self-b]
  for s = seeds
    x = gen(n, 6000*is + s);
    for t = t0:(n-1)
      ell  = matchlens(x, t);
      cont = x(1:t);
      Lt   = max(ell);
      ba   = log(t) / (Lt + 1);
      bb   = log(t) / (Lt + sqrt(Lt + 1));
      bs   = [0.5, 1.0, bo, ba, bb];
      for j = 1:5
        p = soft_myc_dist(ell, cont, bs(j), A, nu);
        L5(j) = L5(j) - log(p(x(t+1)));
      end
    end
  end
  L5 = L5 / (numel(seeds) * (n - t0));
  fprintf('   %-18s %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f\n', ...
          name, H, L5(1), L5(2), L5(3), L5(4), L5(5));
end
fprintf('   [self-tuned schedules should track the oracle without tuning]\n');
