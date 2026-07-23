function [loss, nphrases, tper] = m78r_predict(x, A, beta, nu, t0, Kcap, selftune)
%M78R_PREDICT  Context-restart soft Mycielski-78 (per-order suffix tracking).
%
%   The dictionary GROWS by the standard LZ78 parse (phrase count and
%   O(1) amortized growth unchanged), but prediction contexts are the
%   TRUE suffixes of every order k = 0..Kcap, tracked incrementally by
%   the identity
%
%       node_k(t+1) = child( node_{k-1}(t), x(t) ),   node_0 = root,
%
%   i.e., the length-k suffix at time t+1 extends the length-(k-1)
%   suffix at time t by the new symbol.  Cost O(Kcap) per symbol.
%   Counts at node_k are updated whenever the length-k suffix is
%   representable, so cnt(w,a) = #{times w was the length-|w| suffix
%   and a followed}: the exact N_k(a) restricted to trie contexts.
%   No context reset ever occurs.
%
%   selftune = 0 : fixed beta;  1 : beta_t = log(t)/(d_t+1), d_t the
%   deepest representable order (truncated at Kcap).

n     = numel(x);
child = zeros(n + 2, A);
cnt   = zeros(n + 2, A);
nn    = 1;                         % root = 1
pcur  = 1;                         % LZ78 parse pointer (growth only)
nodes = zeros(1, Kcap + 1);        % nodes(k+1) = length-k suffix node (0 if absent)
nodes(1) = 1;
loss  = 0;  npred = 0;

tstart = tic();
for t = 1:n
  live = find(nodes > 0);          % existing orders + 1
  d    = live(end) - 1;            % deepest representable order
  % ---- predict x(t) ----------------------------------------------------
  if t > t0
    if selftune, b = log(t) / (d + 1); else, b = beta; end
    fac = 1 - exp(-b);
    p = exp(-b * d) * cnt(1, :);                    % order 0
    for j = 2:numel(live)
      k = live(j) - 1;
      p = p + fac * exp(b * (k - d)) * cnt(nodes(live(j)), :);
    end
    s = sum(p);
    if s > 0, p = p / s; else, p = ones(1, A) / A; end
    p = (1 - nu) * p + nu / A;
    loss  = loss - log(p(x(t)));
    npred = npred + 1;
  end
  % ---- consume x(t) ----------------------------------------------------
  a = x(t);
  for j = 1:numel(live)
    cnt(nodes(live(j)), a) = cnt(nodes(live(j)), a) + 1;
  end
  % dictionary growth via LZ78 parse
  c = child(pcur, a);
  if c > 0
    pcur = c;
  else
    nn = nn + 1;  child(pcur, a) = nn;  pcur = 1;
  end
  % incremental suffix-node shift (descending k avoids overwrite)
  for k = min(Kcap, t):-1:1
    if nodes(k) > 0
      nodes(k+1) = child(nodes(k), a);
    else
      nodes(k+1) = 0;
    end
  end
end
tper     = toc(tstart) / n;
loss     = loss / max(npred, 1);
nphrases = nn - 1;
end
