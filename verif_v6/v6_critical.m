%V6_CRITICAL  Modules V6.1-V6.3: matching upper bounds and the
%             rotation (Vervaat) identity.
%
%   V6.1  Excess of the UNTUNED critical fixed temperature beta = H
%         against the oracle beta*(n).  Theorem (fixed upper bound):
%         both are Theta~(1/sqrt(log n)), so the ratio must stay
%         bounded (predicted ~1.5-2.2, no growth).
%
%   V6.2  Rotation identity under the self-critical tilt:
%         P[max B <= s]  =  E[ N_s/(L+1) ],   N_s = #{k: B_k >= M-s},
%         exact for exchangeable increments; expect agreement within
%         ~25% (cycle-boundary, gamma_k and R_k corrections).
%         (No identity is expected under the fixed tilt: free walk.)
%
%   V6.3  End-window argmax probability P[argmax >= L-1] (w = 2):
%         arcsine ~ sqrt(w/L) under beta = H  vs  uniform ~ w/L
%         under the self tilt; the ratio should grow ~ sqrt(log n).

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;  seeds = [1 2 3];
[~, ~, H] = gen_markov(10, A, delta);

% ---------------- V6.1 ---------------------------------------------------
ns    = [1000 2000 4000 8000 16000];
bstar = [0.714 0.729 0.747 0.760 0.776];
fprintf('[V6.1] Untuned critical beta = H vs oracle beta*(n)\n');
fprintf('   %6s | %8s %8s | %6s\n', 'n', 'exc(H)', 'exc(b*)', 'ratio');
for i = 1:numel(ns)
  n = ns(i);  t0 = n - min(400, floor(n/5)) + 1;
  acc = zeros(1,2);
  for s = seeds
    x = gen_markov(n, A, delta, 61000 + 10*i + s);
    acc = acc + beta_sweep(x, [H bstar(i)], A, nu, t0);
  end
  acc = acc/numel(seeds) - H;
  fprintf('   %6d | %8.4f %8.4f | %6.2f\n', n, acc(1), acc(2), acc(1)/acc(2));
end
fprintf('   [prediction: ratio bounded ~1.5-2.2, not growing]\n\n');

% ---------------- V6.2 + V6.3 --------------------------------------------
ns2   = [1000 4000 16000];
svals = [0.5 1.0 2.0];
fprintf('[V6.2/V6.3] Rotation identity and end-window argmax\n');
for i = 1:numel(ns2)
  n = ns2(i);  t0 = floor(n/2) + 1;
  Pm = zeros(1,3);  Nn = zeros(1,3);          % self: P[M<=s], E[N_s/(L+1)]
  Ef = 0;  Es = 0;  cnt = 0;                  % end-window argmax, both tilts
  for s = seeds
    x = gen_markov(n, A, delta, 71000 + 10*i + s);
    for t = t0:(n-1)
      ell = matchlens(x, t);
      L   = max(ell);
      if L < 3, continue; end
      h  = histc(ell, 0:L);
      bs = [H, log(t)/(L+1)];
      for j = 1:2
        w = h .* exp(bs(j) * ((0:L) - L));
        B = log(w / w(1));
        [M, km] = max(B);
        if j == 1
          Ef = Ef + (km >= L);                % km is 1-based: orders km-1
        else
          Es = Es + (km >= L);
          Pm = Pm + (M <= svals);
          Nn = Nn + sum(B(:) >= M - svals, 1) / (L+1);
        end
      end
      cnt = cnt + 1;
    end
  end
  fprintf('  n = %5d:  identity (self): ', n);
  fprintf('s=%.1f: %.3f vs %.3f   ', [svals; Pm/cnt; Nn/cnt]);
  fprintf('\n              end-window: fixed(H) %.3f  self %.3f  ratio %.2f\n', ...
          Ef/cnt, Es/cnt, Ef/max(Es,1));
end
fprintf(['   [predictions: identity columns agree within ~25%%;\n' ...
         '    end-window fixed ~ 1/sqrt(logn) > self ~ 1/logn, ratio growing]\n']);
