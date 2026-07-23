%V10_PLOTS_PNG  Module V10.3: Octave-rendered PNG versions of the
%               three manuscript figures.
%
%   Companion to V10_MAKETIKZ (which emits TikZ bodies compiled by the
%   manuscript).  This script renders the SAME data with Octave's own
%   graphics, writing standalone PNG files for inspection or manual
%   insertion via \includegraphics.  The LaTeX source is not touched.
%
%   Requires a graphics toolkit (gnuplot or qt) and the FreeFont/DejaVu
%   fonts; on a bare container:
%       apt-get install -y gnuplot-nox fonts-freefont-otf
%
%   Outputs (in ../figures/png/):  fig_approach.png
%                         fig_separation.png
%                         fig_errprob.png

clear all; more off;

tk = available_graphics_toolkits();
if isempty(tk)
  error(['no graphics toolkit available: install gnuplot-nox (and ' ...
         'fonts-freefont-otf) or use v10_maketikz.m instead']);
end
graphics_toolkit(tk{1});
fprintf('[V10.3] rendering with the %s toolkit\n', tk{1});

outdir = '../figures/png';
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---------------- data (identical to the TikZ figures) ------------------
ns      = [1000 2000 4000 8000 16000];
logn    = log(ns);
H035    = 1.0320;
bstar   = [0.714 0.729 0.747 0.760 0.776];
exc_or  = [0.1464 0.1074 0.1201 0.1098 0.0855];
exc_slf = [0.0926 0.0558 0.0739 0.0612 0.0356];

if exist('errprob_data.mat', 'file')
  load('errprob_data.mat');            % Erec, Crec
else
  % measured values from Module V10.1 (fallback if the .mat is absent)
  Crec = [ 1000 0.0189 0.0089
           2000 0.0176 0.0055
           4000 0.0184 0.0052
           8000 0.0208 0.0064
          16000 0.0198 0.0046];
end

BLUE = [0.10 0.30 0.70];
RED  = [0.75 0.15 0.12];
GRN  = [0.15 0.45 0.18];
VIO  = [0.45 0.20 0.55];
DPI  = '-r200';

set(0, 'defaultaxesfontsize', 10);
set(0, 'defaulttextfontsize', 10);
set(0, 'defaultlinelinewidth', 1.6);

% ================= Figure 1: approach to criticality ====================
f1 = figure('visible', 'off', 'position', [0 0 640 520]);

subplot(2,1,1);
plot(logn, bstar, 'o-', 'color', BLUE, 'markerfacecolor', BLUE, ...
     'markersize', 5); hold on;
plot([logn(1)-0.2 logn(end)+0.2], H035*[1 1], 'k--', 'linewidth', 1.2);
grid on;
xlabel('log n'); ylabel('beta*(n)');
xlim([logn(1)-0.2 logn(end)+0.2]); ylim([0.68 1.08]);
legend({'measured', 'H = 1.0320'}, 'location', 'northwest');
legend boxoff;
title('(a) optimal inverse temperature rises toward H from below');

subplot(2,1,2);
plot(1./logn, exc_or, 's-', 'color', RED, 'markerfacecolor', RED, ...
     'markersize', 5); hold on;
c = exc_or(end)*logn(end);
plot(1./logn, c./logn, 'k:', 'linewidth', 1.3);
grid on;
xlabel('1 / log n'); ylabel('excess (nats)');
legend({'oracle', 'proportional to 1/log n'}, 'location', 'northwest');
legend boxoff;
title('(b) oracle excess log-loss against the 1/log n law');

print(f1, '-dpng', DPI, fullfile(outdir, 'fig_approach.png'));
close(f1);

% ================= Figure 2: two-sided separation =======================
f2 = figure('visible', 'off', 'position', [0 0 640 480]);

c1 = exc_or(1)*sqrt(logn(1));
c2 = exc_slf(1)*logn(1);
loglog(logn, exc_or,  's-', 'color', RED,  'markerfacecolor', RED, ...
       'markersize', 6); hold on;
loglog(logn, exc_slf, 'o-', 'color', BLUE, 'markerfacecolor', BLUE, ...
       'markersize', 6);
loglog(logn, c1./sqrt(logn), 'k--', 'linewidth', 1.2);
loglog(logn, c2./logn,       'k:',  'linewidth', 1.4);
grid on;
xlabel('log n'); ylabel('excess log-loss (nats)');
xlim([logn(1)*0.96 logn(end)*1.04]);
legend({'fixed beta*(n)', 'self-critical', ...
        'c / sqrt(log n)', 'c / log n'}, 'location', 'southwest');
legend boxoff;
title('two-sided separation: fixed ~ 1/sqrt(log n), self ~ 1/log n');

print(f2, '-dpng', DPI, fullfile(outdir, 'fig_separation.png'));
close(f2);

% ================= Figure 3: error probability ==========================
f3 = figure('visible', 'off', 'position', [0 0 640 520]);

subplot(2,1,1);
semilogx(Crec(:,1), Crec(:,2), 's-', 'color', RED, ...
         'markerfacecolor', RED, 'markersize', 5); hold on;
semilogx(Crec(:,1), Crec(:,3), 'o-', 'color', BLUE, ...
         'markerfacecolor', BLUE, 'markersize', 5);
grid on;
xlabel('n'); ylabel('excess error prob.');
ylim([0 0.025]);
legend({'fixed beta*(n)', 'self-critical'}, 'location', 'northeast');
legend boxoff;
title('(a) 0-1 loss: the self-critical advantage widens with n');

subplot(2,1,2);
r01 = Crec(:,2) ./ max(Crec(:,3), 1e-9);
rll = exc_or(:) ./ exc_slf(:);
semilogx(Crec(:,1), r01, 'd-', 'color', GRN, 'markerfacecolor', GRN, ...
         'markersize', 5); hold on;
semilogx(ns, rll, '^-', 'color', VIO, 'markerfacecolor', VIO, ...
         'markersize', 5);
grid on;
xlabel('n'); ylabel('adaptivity ratio');
ylim([1 5]);
legend({'0-1 loss', 'log-loss'}, 'location', 'northwest');
legend boxoff;
title('(b) adaptivity gain is larger for decisions than for codelength');

print(f3, '-dpng', DPI, fullfile(outdir, 'fig_errprob.png'));
close(f3);

fprintf('  wrote %s/fig_approach.png, fig_separation.png, fig_errprob.png\n', outdir);
fprintf('  0-1 adaptivity ratio: ');  fprintf('%.2f ', r01);  fprintf('\n');
fprintf('  log-loss ratio:       ');  fprintf('%.2f ', rll);  fprintf('\n');
