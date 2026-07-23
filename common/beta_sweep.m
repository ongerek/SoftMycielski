function [loss, Lrec] = beta_sweep(x, betas, A, nu, t0)
%BETA_SWEEP  Average sequential log-loss of soft-Mycielski over a beta grid.
%
%   [loss, Lrec] = BETA_SWEEP(x, betas, A, nu, t0)
%
%   Predicts x(t+1) from x(1..t) for t = t0 .. n-1 and every beta in
%   `betas`, reusing the (beta-independent) match lengths at each t.
%
%   loss : 1 x numel(betas), mean of -log p_beta(x(t+1)) in NATS
%   Lrec : [t ; L_t] record of maximal match lengths (for the
%          Wyner-Ziv / Ornstein-Weiss check  L_t ~ log(t)/H )

n  = numel(x);
nb = numel(betas);
loss  = zeros(1, nb);
cnt   = 0;
Lrec  = zeros(2, n - t0);

for t = t0:(n-1)
  ell  = matchlens(x, t);
  cont = x(1:t);
  cnt  = cnt + 1;
  Lrec(:, cnt) = [t; max(ell)];
  truth = x(t+1);
  for b = 1:nb
    p = soft_myc_dist(ell, cont, betas(b), A, nu);
    loss(b) = loss(b) - log(p(truth));
  end
end
loss = loss / cnt;
Lrec = Lrec(:, 1:cnt);
end
