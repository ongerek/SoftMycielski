%V11_FIGDATA  Module V11.1: verified data for the explanatory figures.
%
%   The tutorial figures must illustrate TRUE numbers.  This script
%   computes, for the running example x = "abracadabra":
%     (i)   the match lengths ell(m,t) at t = 11 and the continuations,
%     (ii)  the soft Mycielski distribution at several beta, including
%           the beta -> 0 limit (order-0 marginal) and beta -> Inf
%           (hard Mycielski),
%     (iii) the LZ78 phrase parse and the resulting trie, plus the
%           per-order suffix nodes of the M78-R recursion (eq:shift),
%     (iv)  a realized fluctuation path D_k with its chord and bridge
%           B_k = D_k - (k/L) D_L, exported as TikZ coordinates.

clear all; more off;
sym = 'abracadabra';
map = 'abcdr';                 % alphabet, sorted
A   = numel(map);
x   = zeros(1, numel(sym));
for i = 1:numel(sym), x(i) = find(map == sym(i)); end
t   = numel(x);

fprintf('[V11.1] example x = %s   (t = %d, A = %d, alphabet %s)\n\n', ...
        sym, t, A, map);

% ---------- (i) match lengths at t ---------------------------------------
fprintf('  m  x(m)  ell(m,t)  cont x(m+1)\n');
ell = zeros(1, t-1);  cont = zeros(1, t-1);
for m = 1:(t-1)
  L = 0;
  while (m - L >= 1) && (t - L >= 1) && (x(m-L) == x(t-L))
    L = L + 1;
  end
  ell(m) = L;  cont(m) = x(m+1);
  fprintf('  %2d   %s      %d         %s\n', m, map(x(m)), L, map(cont(m)));
end
Lmax = max(ell);
fprintf('\n  longest match L_t = %d at m = %s -> hard prediction "%s"\n\n', ...
        Lmax, mat2str(find(ell == Lmax)), map(cont(find(ell == Lmax, 1))));

% ---------- (ii) soft distributions --------------------------------------
betas = [0 0.5 1 2];
fprintf('  beta  |');  for a = 1:A, fprintf('   p(%s)', map(a)); end
fprintf('\n');
for b = betas
  w = exp(b * ell);
  p = zeros(1, A);
  for a = 1:A, p(a) = sum(w(cont == a)); end
  p = p / sum(p);
  fprintf('  %4.1f  |', b);  fprintf('  %.3f', p);  fprintf('\n');
end
w = double(ell == Lmax);  p = zeros(1, A);
for a = 1:A, p(a) = sum(w(cont == a)); end
p = p / sum(p);
fprintf('   Inf  |');  fprintf('  %.3f', p);  fprintf('   (hard Mycielski)\n');
% order-0 check: beta = 0 must equal the empirical marginal of x(2..t)
mar = zeros(1, A);
for a = 1:A, mar(a) = sum(x(2:t) == a); end
mar = mar / sum(mar);
fprintf('   emp. marginal of x(2..t): ');  fprintf(' %.3f', mar);
fprintf('   <- beta=0 row must match\n\n');

% ---------- (iii) LZ78 parse and trie ------------------------------------
fprintf('  LZ78 parse:  ');
child = zeros(64, A);  nn = 1;  pcur = 1;  phrases = {};  cur = '';
par = zeros(1,64);  lab = zeros(1,64);
for i = 1:t
  a = x(i);  cur = [cur map(a)];
  c = child(pcur, a);
  if c > 0
    pcur = c;
  else
    nn = nn + 1;  child(pcur, a) = nn;  par(nn) = pcur;  lab(nn) = a;
    phrases{end+1} = cur;  cur = '';  pcur = 1;
  end
end
if ~isempty(cur), phrases{end+1} = [cur ' (incomplete)']; end
for i = 1:numel(phrases), fprintf('%s|', phrases{i}); end
fprintf('\n  trie nodes (node: parent -%s-> ):\n', 'sym');
for i = 2:nn
  fprintf('    node %d: parent %d, edge %s\n', i, par(i), map(lab(i)));
end

% M78-R per-order suffix nodes after consuming the whole string:
% node_k tracks the node reached by the length-k suffix of x_1..t.
fprintf('  M78-R suffix nodes after t = %d:\n', t);
for k = 1:4
  s = x(t-k+1:t);  nd = 1;  ok = true;
  for j = 1:k
    c = child(nd, s(j));
    if c == 0, ok = false; break; end
    nd = c;
  end
  if ok
    fprintf('    k = %d  suffix "%s" -> node %d\n', k, map(s), nd);
  else
    fprintf('    k = %d  suffix "%s" -> (absent)\n', k, map(s));
  end
end

% ---------- (iv) bridge path for the schematic ---------------------------
A2 = 4;  delta = 0.35;
H  = -((1-delta)*log(1-delta) + delta*log(delta/(A2-1)));
xiv = [H+log(1-delta), H+log(delta/(A2-1))];
rand('state', 20260723);
L = 14;
xi = xiv(1 + (rand(1,L) > (1-delta)));
D  = [0 cumsum(xi)];
B  = D - (0:L)/L * D(end);
fprintf('\n  bridge schematic (L = %d, H = %.4f):\n', L, H);
fprintf('    D coords: ');
for k = 0:L, fprintf('(%d,%.3f) ', k, D(k+1)); end
fprintf('\n    B coords: ');
for k = 0:L, fprintf('(%d,%.3f) ', k, B(k+1)); end
fprintf('\n    chord end D_L = %.3f;  max D = %.3f at k = %d;  max B = %.3f at k = %d\n', ...
        D(end), max(D), find(D == max(D)) - 1, max(B), find(B == max(B)) - 1);
