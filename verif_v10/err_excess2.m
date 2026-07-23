function [errx, err, bayes] = err_excess2(x, P2, betas, A, nu, t0)
%ERR_EXCESS2  Error-probability excess for an order-2 source.
%   P2 is A x A x A with P2(i,j,a) = Pr(x=a | prev2=i, prev1=j).
%   Column layout: [betas..., self-critical].  Rao-Blackwellized.

n = numel(x);  nb = numel(betas);
err = zeros(1, nb + 1);  bayes = 0;  cnt = 0;
for t = t0:(n-1)
  ell  = matchlens(x, t);
  cont = x(1:t);
  Pt   = squeeze(P2(x(t-1), x(t), :))';
  bayes = bayes + (1 - max(Pt));
  for b = 1:nb
    p = soft_myc_dist(ell, cont, betas(b), A, nu);
    m = max(p);  win = find(p >= m - 1e-12);
    err(b) = err(b) + (1 - mean(Pt(win)));
  end
  L = max(ell);
  p = soft_myc_dist(ell, cont, log(t)/(L+1), A, nu);
  m = max(p);  win = find(p >= m - 1e-12);
  err(nb+1) = err(nb+1) + (1 - mean(Pt(win)));
  cnt = cnt + 1;
end
err = err / cnt;  bayes = bayes / cnt;
errx = err - bayes;
end
