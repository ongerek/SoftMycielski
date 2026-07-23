function [x, Hrate] = gen_markov2(n, A, delta, seed)
%GEN_MARKOV2  Order-2 Markov chain with UNINFORMATIVE order-1 statistics.
%
%   Successor rule: s(i,j) = mod(i+j, A) + 1.  Given the pair
%   (x(t-1), x(t)) = (i,j), the next symbol is s(i,j) with probability
%   (1-delta), and each other symbol with probability delta/(A-1).
%
%   Because i -> s(i,j) is a bijection for fixed j, the pair chain is
%   doubly stochastic: the stationary law over pairs is uniform, hence
%     * entropy rate:  H = -(1-delta)log(1-delta) - delta*log(delta/(A-1))
%       (identical closed form to GEN_MARKOV), and
%     * the order-1 conditional P(x(t+1) | x(t)) is EXACTLY uniform:
%       single-symbol contexts carry zero predictive information.
%
%   This isolates the effect of source memory on the optimal temperature.

if nargin >= 4
  rand('state', seed);
end

Hrate = -((1-delta)*log(1-delta) + delta*log(delta/(A-1)));

x = zeros(1, n);
x(1) = ceil(rand()*A);
x(2) = ceil(rand()*A);
others = @(k) [1:(k-1), (k+1):A];
for t = 3:n
  i = x(t-2);  j = x(t-1);
  s = mod(i + j, A) + 1;
  if rand() < 1 - delta
    x(t) = s;
  else
    o = others(s);
    x(t) = o(ceil(rand()*(A-1)));
  end
end
end
