function relerr = abel_check(x, t, beta, A)
%ABEL_CHECK  Verify the Abel-summation identity behind soft-Mycielski.
%
%   For each symbol a, compare
%     LHS(a) = sum_{m<t} 1[x(m+1)=a] * exp(beta*ell(m,t))
%   with
%     RHS(a) = N_0(a) + (1 - exp(-beta)) * sum_{k>=1} exp(beta*k) * N_k(a)
%   where N_k(a) = #{ m < t : ell(m,t) >= k, x(m+1) = a }.
%
%   Returns the maximum relative error over a; should be O(machine eps).

ell  = matchlens(x, t);
cont = x(1:t);
Lmax = max(ell);

LHS = zeros(1, A);
RHS = zeros(1, A);
for a = 1:A
  sel    = (cont == a);
  LHS(a) = sum(exp(beta * ell(sel)));
  Nk     = zeros(1, Lmax + 1);          % Nk(k+1) = N_k(a), k = 0..Lmax
  for k = 0:Lmax
    Nk(k+1) = sum(sel & (ell >= k));
  end
  RHS(a) = Nk(1) + (1 - exp(-beta)) * sum(exp(beta * (1:Lmax)) .* Nk(2:end));
end
relerr = max(abs(LHS - RHS) ./ max(abs(LHS), realmin));
end
