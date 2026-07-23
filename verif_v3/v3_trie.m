%V3_TRIE  Module V3.3: context-restart trie vs plain M78 vs exact.
%
%   Same setup as Module V2.3 (delta = 0.35, n = 8000, beta = 0.7,
%   seeds 7001-7002) so the numbers are directly comparable.

clear all; more off;
A = 4;  nu = 1e-3;  delta = 0.35;  n = 8000;  t0 = 4000;
beta = 0.7;  seeds = [1 2];  Kcap = 16;
[~, ~, H] = gen_markov(10, A, delta);

fprintf('[V3.3] Context-restart trie  (delta=%.2f, H=%.4f, beta=%.2f)\n', ...
        delta, H, beta);
Lex=0; Tex=0; Llz=0; Tlz=0; Lr=0; Tr=0; Lrs=0; Trs=0; NP=0;
for s = seeds
  x = gen_markov(n, A, delta, 7000 + s);
  % exact
  tstart = tic();  le = 0;
  for t = t0:(n-1)
    ell = matchlens(x, t);
    p   = soft_myc_dist(ell, x(1:t), beta, A, nu);
    le  = le - log(p(x(t+1)));
  end
  Tex = Tex + toc(tstart)/(n-t0);  Lex = Lex + le/(n-t0);
  % plain M78 (v2)
  [ll, np, tp] = m78_predict(x, A, beta, nu, t0, 'lz');
  Llz = Llz + ll;  Tlz = Tlz + tp;  NP = NP + np;
  % context-restart, fixed beta
  [lr, ~, tp] = m78r_predict(x, A, beta, nu, t0, Kcap, 0);
  Lr = Lr + lr;  Tr = Tr + tp;
  % context-restart, self-tuned beta
  [ls, ~, tp] = m78r_predict(x, A, beta, nu, t0, Kcap, 1);
  Lrs = Lrs + ls;  Trs = Trs + tp;
end
ns = numel(seeds);
fprintf('   %-28s %10s %14s\n', 'predictor', 'loss', 'time/symbol');
fprintf('   %-28s %10.4f %12.2e s\n', 'exact soft Mycielski',       Lex/ns, Tex/ns);
fprintf('   %-28s %10.4f %12.2e s\n', 'M78 plain (v2, lz)',         Llz/ns, Tlz/ns);
fprintf('   %-28s %10.4f %12.2e s\n', 'M78-R restart, fixed beta',  Lr/ns,  Tr/ns);
fprintf('   %-28s %10.4f %12.2e s\n', 'M78-R restart, self-tuned',  Lrs/ns, Trs/ns);
fprintf('   phrases c(n) = %d;  H = %.4f\n', round(NP/ns), H);
fprintf('   gap closure: plain gap %.4f -> restart gap %.4f (%.0f%%)\n', ...
        Llz/ns - Lex/ns, Lr/ns - Lex/ns, ...
        100*(1 - (Lr/ns - Lex/ns)/(Llz/ns - Lex/ns)));
