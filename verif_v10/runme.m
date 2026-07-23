%RUNME  Entry point for verification package v10.
%  V10.1: error-probability (0-1 loss) columns with excess over the
%         Bayes floor (writes errprob_data.mat).
%  V10.2: emit self-contained TikZ/pgfplots figures into ../figures/tikz/.
%  V10.3: Octave-rendered PNG versions of the same three figures
%         (needs a graphics toolkit; writes ../figures/png/).

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v10_errprob;
v10_maketikz;
if ~isempty(available_graphics_toolkits()), v10_plots_png; end
