%V8_TILTSWEEP  Module V8.1b: tuned-fixed constant c1 via a tilt
%              sweep in the idealized model, plus refined c2 and
%              kappa_h statistics at large L.
%
%   Fixed beta = H - eps, weights propto e^{D_k - eps k}; sweep the
%   scaled tilt z = eps*sqrt(L)/sv over a grid covering the
%   fluctuation window (real-system optimum sits near z ~ sqrt(H)
%   = 1.02).  c1(L) = sqrt(L) * min_z E[KL(z)].

clear all; more off;
A = 4;  delta = 0.35;
H  = -((1-delta)*log(1-delta) + delta*log(delta/(A-1)));
sv  = sqrt((1-delta)*log(1-delta)^2 + delta*log(delta/(A-1))^2 - H^2);
xiv = [H+log(1-delta), H+log(delta/(A-1))];
P   = [1-delta, repmat(delta/(A-1),1,A-1)];
piv = ones(1,A)/A;
Ncap = 300;
zs   = [0 0.3 0.6 0.9 1.2 1.6 2.2];
Ls   = [100 200 400 800];
reps = [24000 16000 10000 7000];

fprintf('[V8.1b] Tilt sweep, refined constants (fresh MC block)\n');
fprintf('  %5s | %7s %8s | %8s | z-grid sqrt(L)*E[KL]:\n', ...
        'L', 'kap_h', 'c2(L)', 'c1(L)');
for il = 1:numel(Ls)
  L = Ls(il);  R = reps(il);
  rand('state', 95000 + L);
  accz = zeros(1, numel(zs));  acc2 = 0;  acck = 0;
  for rep = 1:R
    xi = xiv(1 + (rand(1,L) > (1-delta)));
    D  = [0 cumsum(xi)];
    B  = D - (0:L)/L * D(end);
    lN = max((L-(0:L))*(H - D(end)/L) + B, 0);
    invN = exp(-min(lN, 700));
    lowm = (lN <= log(Ncap));  lowm(1) = false;
    mmax = Ncap;
    st   = 1 + sum(repmat(rand(1,mmax)',1,A-1) > ...
                   repmat(cumsum(P(1:A-1)),mmax,1), 2)';
    cum  = zeros(A, mmax+1);
    for a = 1:A, cum(a,:) = [0 cumsum(st == a)]; end
    qs = zeros(sum(lowm), A);  ii = 0;
    for k = find(lowm)
      m = max(1, round(exp(lN(k))));  ii = ii + 1;
      qs(ii,:) = cum(:, m+1)' / m;
    end
    for tt = 0:numel(zs)
      if tt == 0, u = B; else, u = D - zs(tt)*sv/sqrt(L)*(0:L); end
      u = u - max(u);  w = exp(u);  w = w / sum(w);
      qt = w(1)*piv + (sum(w(2:end)) - sum(w(lowm)))*P;
      ii = 0;
      for k = find(lowm)
        ii = ii + 1;  qt = qt + w(k)*qs(ii,:);
      end
      kl = sum(P .* log(P ./ qt));
      wk = w(2:end);  iv = invN(2:end);  lm = lowm(2:end);
      tailall = fliplr(cumsum(fliplr(wk))) - wk;
      qall = sum(wk.^2 .* iv) + 2*sum(wk .* iv .* tailall);
      wl = wk .* lm;
      taill = fliplr(cumsum(fliplr(wl))) - wl;
      qlow = sum(wl.^2 .* iv) + 2*sum(wl .* iv .* taill);
      kl = kl + (A-1)/2 * (qall - qlow);
      if tt == 0
        acc2 = acc2 + kl;  acck = acck + exp(2*(B(1) - max(B)));
      else
        accz(tt) = accz(tt) + kl;
      end
    end
  end
  cz = sqrt(L) * accz / R;
  fprintf('  %5d | %7.3f %8.3f | %8.3f | ', L, L*acck/R, L*acc2/R, min(cz));
  fprintf('%.3f ', cz);  fprintf('\n');
end
fprintf('  [z grid: ');  fprintf('%.1f ', zs);  fprintf(']\n');
