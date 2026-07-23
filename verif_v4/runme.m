%RUNME  Entry point for verification package v4.
%  V4.1: cross-entropy scaling of beta*(n); V4.2: head-weight atom test.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v4_scaling;
