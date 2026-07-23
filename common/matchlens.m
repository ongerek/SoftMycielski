function ell = matchlens(x, t)
%MATCHLENS  Suffix match lengths ell(m,t) for all candidate positions m.
%
%   ell = MATCHLENS(x, t) returns a 1 x t vector where ell(m+1) is the
%   length of the longest common suffix between x(1..m) and x(1..t):
%
%       ell(m,t) = max{ L >= 0 : x(m-L+1..m) == x(t-L+1..t) },  0 <= m <= t-1
%
%   Position m=0 (empty history, ell=0) is included so that every symbol
%   x(m+1), m = 0..t-1, is a candidate continuation.  Vectorized level-by-
%   level: survivors of level k satisfy ell >= k; total cost O(t * L_t).

ell   = zeros(1, t);          % ell(1) is m=0, stays 0
alive = 1:(t-1);              % candidate m values (m >= 1)
k = 0;
while ~isempty(alive)
  k = k + 1;
  alive = alive(alive >= k);              % need x(m-k+1) to exist
  if isempty(alive), break; end
  keep  = x(alive - k + 1) == x(t - k + 1);
  alive = alive(keep);
  ell(alive + 1) = k;                     % survivors have ell >= k
end
end
