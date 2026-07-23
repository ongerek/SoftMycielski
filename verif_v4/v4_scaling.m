%V4_SCALING  Modules V4.1 and V4.2: sharp-rate mechanism tests.
%
%   V4.1  beta*(n) and excess(n) at delta = 0.10 and 0.50, extending
%         Table VIII beyond the single entropy of Module V2.1.
%         Reported scaling variables:  y = (H - beta*) log(n) / H  and
%         excess * log(n).  Discriminates the H- vs H^2-scaling of the
%         optimal temperature offset and tests the sign of approach
%         (delta = 0.10 starts SUPERcritical at small n).
%
%   V4.2  Distribution of the realized order-0 weight W_0 at fixed
%         beta = H - eps.  Excursion mechanism prediction: W_0 carries
%         a Theta(eps)-probability atom at Theta(1), so
%         E[W_0^2]/E[W_0] = Theta(1)  and  P[W_0 > 0.3] = Theta(eps);
%         the idealized geometric profile instead gives
%         W_0 \approx 1-e^{-eps} deterministically, i.e. ratio ~ eps
%         and P[W_0 > 0.3] = 0.

clear all; more off;
A = 4;  nu = 1e-3;  seeds = [1 2 3];

% ---------------- V4.1 ---------------------------------------------------
dset = {0.10, 0.50};
grids = { [0.30 0.36 0.42 0.46 0.50 0.54 0.58 0.62 0.68 0.75 0.85 1.00], ...
          [0.50 0.60 0.68 0.74 0.80 0.86 0.92 0.98 1.05 1.15 1.30 1.50] };
ns = [1000 2000 4000 8000 16000];
fprintf('[V4.1] Scaling of the optimal temperature across entropies\n');
for id = 1:2
  delta = dset{id};  betas = grids{id};
  [~, ~, H] = gen_markov(10, A, delta);
  fprintf('  delta = %.2f, H = %.4f\n', delta, H);
  fprintf('   %7s | %7s %8s | %8s | %7s %8s\n', 'n', 'beta*', ...
          'H-beta*', 'excess', 'y', 'exc*logn');
  for i = 1:numel(ns)
    n  = ns(i);
    t0 = n - min(400, floor(n/5)) + 1;
    lossacc = zeros(1, numel(betas));
    for s = seeds
      x = gen_markov(n, A, delta, 20000*id + 5000*i + s);
      lossacc = lossacc + beta_sweep(x, betas, A, nu, t0);
    end
    loss = lossacc / numel(seeds);
    [lmin, imin] = min(loss);
    if imin > 1 && imin < numel(betas)
      b3 = betas(imin-1:imin+1);  l3 = loss(imin-1:imin+1);
      c  = polyfit(b3, l3, 2);
      bstar = -c(2)/(2*c(1));
      if bstar < b3(1) || bstar > b3(3), bstar = betas(imin); end
      lmin = polyval(c, bstar);
    else
      bstar = betas(imin);
    end
    fprintf('   %7d | %7.3f %8.3f | %8.4f | %7.2f %8.3f\n', n, bstar, ...
            H-bstar, lmin-H, (H-bstar)*log(n)/H, (lmin-H)*log(n));
  end
end
fprintf('  [delta=0.35 reference (V2.1): y = 2.13, 2.23, 2.29, 2.37, 2.40;\n');
fprintf('   exc*logn = 1.01, 0.82, 1.00, 0.99, 0.83]\n\n');

% ---------------- V4.2 ---------------------------------------------------
fprintf('[V4.2] Head-weight atom test  (delta = 0.35, n = 4000)\n');
delta = 0.35;  [~, ~, H] = gen_markov(10, A, delta);
n = 4000;  t0 = 2001;
epss = [0.05 0.10 0.20 0.30];
fprintf('   %6s %7s | %9s %9s %9s | %9s %9s\n', 'eps', 'beta', ...
        'E[W0]', 'E[W0^2]', 'ratio', 'P[W0>.3]', 'ideal');
for ie = 1:numel(epss)
  eps = epss(ie);  b = H - eps;
  m1 = 0;  m2 = 0;  patom = 0;  cnt = 0;
  for s = seeds
    x = gen_markov(n, A, delta, 31000 + s);
    for t = t0:(n-1)
      ell = matchlens(x, t);
      L   = max(ell);
      h   = histc(ell, 0:L);              % h(k+1) = #{m: ell = k}
      w   = h .* exp(b * ((0:L) - L));    % overflow-safe masses
      W0  = w(1) / sum(w);
      m1 = m1 + W0;  m2 = m2 + W0^2;
      patom = patom + (W0 > 0.3);
      cnt = cnt + 1;
    end
  end
  m1 = m1/cnt;  m2 = m2/cnt;  patom = patom/cnt;
  fprintf('   %6.2f %7.3f | %9.4f %9.4f %9.3f | %9.4f %9.4f\n', ...
          eps, b, m1, m2, m2/m1, patom, 1-exp(-eps));
end
fprintf(['   [mechanism: ratio ~ const, P[W0>.3] ~ eps;' ...
         '  idealized: ratio ~ eps, P = 0]\n']);
