function [x, P, Hrate] = gen_markov(n, A, delta, seed)
%GEN_MARKOV  Order-1 Markov chain over {1,...,A} with closed-form entropy rate.
%
%   From state i, the chain moves to succ(i) = mod(i,A)+1 with probability
%   (1-delta), and to each of the remaining A-1 states with probability
%   delta/(A-1).  The transition matrix is a doubly stochastic circulant,
%   hence the stationary distribution is uniform and the entropy rate is
%
%       H = -(1-delta)*log(1-delta) - delta*log(delta/(A-1))   [nats/symbol]
%
%   Inputs : n (length), A (alphabet size), delta in (0,(A-1)/A], seed (opt)
%   Outputs: x (1 x n sequence), P (A x A transitions), Hrate (nats/symbol)

if nargin >= 4
  rand('state', seed);
end

P = (delta/(A-1)) * ones(A, A);
for i = 1:A
  j = mod(i, A) + 1;
  P(i, j) = 1 - delta;
end

Hrate = -((1-delta)*log(1-delta) + delta*log(delta/(A-1)));

x = zeros(1, n);
x(1) = ceil(rand() * A);
cP = cumsum(P, 2);
for t = 2:n
  u = rand();
  x(t) = find(cP(x(t-1), :) >= u, 1);
end
end
