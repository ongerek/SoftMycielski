function p = soft_myc_dist(ell, cont, beta, A, nu)
%SOFT_MYC_DIST  Soft-Mycielski next-symbol distribution.
%
%   p = SOFT_MYC_DIST(ell, cont, beta, A, nu)
%
%   ell  : 1 x t match lengths ell(m,t), m = 0..t-1   (from MATCHLENS)
%   cont : 1 x t continuations, cont(m+1) = x(m+1)
%   beta : inverse temperature (beta = Inf gives the hard Mycielski
%          predictor with uniform tie-averaging)
%   A    : alphabet size
%   nu   : additive smoothing weight; p <- (1-nu)*p + nu/A
%
%       p_beta(a) propto sum_m 1[x(m+1)=a] * exp(beta*ell(m,t))
%
%   Weights are computed as exp(beta*(ell - max(ell))) for overflow safety.

if isinf(beta)
  L = max(ell);
  w = double(ell == L);
else
  w = exp(beta * (ell - max(ell)));
end

p = zeros(1, A);
for a = 1:A
  p(a) = sum(w(cont == a));
end
p = p / sum(p);
p = (1 - nu) * p + nu / A;
end
