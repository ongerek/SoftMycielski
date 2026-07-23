%V8_CONSTANTS  Module V8.1: the conjectured constants, computed by
%              Monte Carlo in the idealized model at large L.
%
%   The constants c1, c2 of Conjecture 1 are DEFINED as limits of
%   the idealized channel model (Remark on constants).  The model
%   needs only the surprisal walk, not the sequence, so L = 800
%   (equivalent sequence length e^{HL} ~ 10^358) is directly
%   simulable.  For the circulant testbed the walk increments are
%   exactly i.i.d., so the model has no cycle-aggregation error.
%
%   Model per realization (delta = 0.35, r = 1):
%     D_k   = centered surprisal walk,  B_k = D_k - (k/L) D_L
%     N_k   = max(1, exp((L-k)(H - D_L/L) + B_k))   [N_0 = t, N_L ~ 1]
%     w     propto e^{B_k} (self tilt)  or  e^{D_k} (fixed beta = H)
%     q_0   = pi (uniform);  q_k = nested multinomial empiricals
%             (sampled for N_k <= Ncap, analytic quadratic beyond)
%     cost  = KL(P || sum_k w_k q_k), exact per realization
%   Outputs:  kappa_h(L) = L E[e^{-2 max B}]  (predicted -> 1/sv^2
%   = 1.490);  c2(L) = L E[KL_self];  c1crit(L) = sqrt(L) E[KL_H].

clear all; more off;
A = 4;  delta = 0.35;
H  = -((1-delta)*log(1-delta) + delta*log(delta/(A-1)));
sv2 = (1-delta)*log(1-delta)^2 + delta*log(delta/(A-1))^2 - H^2;
xiv = [H+log(1-delta), H+log(delta/(A-1))];    % increment values
P   = [1-delta, repmat(delta/(A-1),1,A-1)];    % true conditional
piv = ones(1,A)/A;                             % order-0 prediction
Ncap = 300;
Ls   = [12 25 50 100 200 400 800];
reps = [20000 20000 20000 12000 8000 6000 4000];

fprintf('[V8.1] Idealized-model constants (delta=%.2f, 1/sv^2=%.3f)\n', ...
        delta, 1/sv2);
fprintf('  %5s | %7s | %8s %8s | %8s | %8s\n', 'L', 'kap_h', ...
        'c2(L)', 'head%', 'c1cr(L)', 'reps');
for il = 1:numel(Ls)
  L = Ls(il);  R = reps(il);
  rand('state', 90000 + L);
  acc2 = 0;  acc1 = 0;  acck = 0;  acchead = 0;
  for rep = 1:R
    xi = xiv(1 + (rand(1,L) > (1-delta)));
    D  = [0 cumsum(xi)];
    B  = D - (0:L)/L * D(end);
    lN = max((L-(0:L))*(H - D(end)/L) + B, 0);   % log N_k, floored
    invN = exp(-min(lN, 700));
    lowm = (lN <= log(Ncap));  lowm(1) = false;  % k>=1, small counts
    % nested continuation stream (shared prefixes across orders)
    mmax = Ncap;
    st   = 1 + sum(repmat(rand(1,mmax)',1,A-1) > ...
                   repmat(cumsum(P(1:A-1)),mmax,1), 2)';
    cum  = zeros(A, mmax+1);
    for a = 1:A, cum(a,:) = [0 cumsum(st == a)]; end
    for tilt = 1:2
      if tilt == 1, u = B - max(B); else, u = D - max(D); end
      w = exp(u);  w = w / sum(w);
      qt = w(1) * piv;
      qhi = sum(w(2:end)) - sum(w(lowm));        % high-count weight
      qt = qt + qhi * P;
      for k = find(lowm)
        m  = max(1, round(exp(lN(k))));
        qt = qt + w(k) * (cum(:, m+1)' / m);
      end
      kl = sum(P .* log(P ./ qt));
      % analytic quadratic for pairs involving a high-count order:
      % (A-1)/2 * [ quad(all k>=1) - quad(low only) ],
      % quad(S) = sum_{k in S} w_k^2/N_k + 2 sum_{k<k', both in S} w_k w_k' / N_k
      wk = w(2:end);  iv = invN(2:end);  lm = lowm(2:end);
      tailall = fliplr(cumsum(fliplr(wk))) - wk;
      qall = sum(wk.^2 .* iv) + 2*sum(wk .* iv .* tailall);
      wl = wk .* lm;
      taill = fliplr(cumsum(fliplr(wl))) - wl;
      qlow = sum(wl.^2 .* iv) + 2*sum(wl .* iv .* taill);
      kl = kl + (A-1)/2 * (qall - qlow);
      if tilt == 1
        acc2 = acc2 + kl;
        acck = acck + exp(2*u(1) - 2*0);          % e^{-2 max B} = e^{2 u(1)}... u(1)=B_0-max=-max
        acchead = acchead + w(1)^2;
      else
        acc1 = acc1 + kl;
      end
    end
  end
  c2 = L*acc2/R;  c1 = sqrt(L)*acc1/R;  kap = L*acck/R;
  hs = L*(0.7033/2)*acchead/R / c2 * 100;
  fprintf('  %5d | %7.3f | %8.3f %7.1f%% | %8.3f | %8d\n', ...
          L, kap, c2, hs, c1, R);
end
fprintf(['  [predictions: kap_h -> 1.490; c2, c1cr converge with',...
         ' ~1/sqrt(L) drift]\n']);
