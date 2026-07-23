%RUNME  Entry point for verification package v8.
%  V8.1:  idealized-model constants MC (base grid);
%  V8.1b: tilt sweep for tuned c1 + refined block;
%  V8.1c: precision block at L = 800, 1600;
%  V8.2:  Rao-Blackwellized real-system anchors at n = 1e5, 2.5e5.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v8_constants;
v8_tiltsweep;
v8_precision;
v8_bign;
