%RUNME  Entry point for verification package v6.
%  V6.1: untuned critical beta = H vs oracle;  V6.2: rotation
%  identity check;  V6.3: arcsine vs uniform argmax end-window laws.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v6_critical;
