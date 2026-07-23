%RUNME  Entry point for verification package v9.
%  V9.1: hidden-Markov testbed (streaming experts, Bayes/fixed-share/
%        sigmoid gates, forward-filter Rao-Blackwellized excess);
%  V9.2: piecewise-stationary testbed (regime-averaging floor,
%        windowed expert, fixed-share tracking).

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v9_hmm;
v9_piecewise;
