function [loss, nphrases, tper] = m78_predict(x, A, beta, nu, t0, mode)
%M78_PREDICT  Streaming soft Mycielski-78 predictor on an LZ78 trie.
%
%   [loss, nphrases, tper] = M78_PREDICT(x, A, beta, nu, t0, mode)
%
%   Maintains an incrementally grown LZ78 dictionary (trie).  The node
%   path root -> current spells the suffix consumed since the last
%   phrase boundary; each node stores continuation counts.  Prediction
%   is the geometric order-mixture of the Abel identity restricted to
%   the trie path, computed in O(depth) per symbol:
%
%     p(a) propto C_0(a) + (1 - e^{-beta}) sum_{k=1..d} e^{beta k} C_k(a)
%
%   mode = 'lz'  : classic LZ78 update, current node only, O(1)/symbol
%   mode = 'path': update all path nodes, O(depth)/symbol
%
%   Outputs: mean log-loss (nats) over t = t0..n-1 predicting x(t),
%   number of dictionary phrases, and mean wall time per symbol [s].

n     = numel(x);
child = zeros(n + 2, A);
cnt   = zeros(n + 2, A);
nn    = 1;                       % node count; node 1 = root
cur   = 1;  d = 0;
path  = zeros(1, 512);  path(1) = 1;
loss  = 0;  npred = 0;
fac   = 1 - exp(-beta);

tstart = tic();
for t = 1:n
  % ---- predict x(t) from the current path -----------------------------
  if t > t0
    p = exp(-beta * d) * cnt(path(1), :);          % order 0, rescaled
    for k = 1:d
      p = p + fac * exp(beta * (k - d)) * cnt(path(k+1), :);
    end
    s = sum(p);
    if s > 0, p = p / s; else, p = ones(1, A) / A; end
    p = (1 - nu) * p + nu / A;
    loss  = loss - log(p(x(t)));
    npred = npred + 1;
  end
  % ---- consume x(t): update counts, descend or spawn ------------------
  a = x(t);
  if strcmp(mode, 'path')
    for k = 0:d
      cnt(path(k+1), a) = cnt(path(k+1), a) + 1;
    end
  else
    cnt(cur, a) = cnt(cur, a) + 1;
  end
  c = child(cur, a);
  if c > 0
    cur = c;  d = d + 1;  path(d+1) = cur;
  else
    nn = nn + 1;  child(cur, a) = nn;              % new phrase
    cur = 1;  d = 0;                               % context reset
  end
end
tper     = toc(tstart) / n;
loss     = loss / max(npred, 1);
nphrases = nn - 1;
end
