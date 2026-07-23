%RUNME  Entry point for verification package v7.
%  V7.1: exact conditional-KL test of the second-order
%  misallocation principle (quadratic law to three decimals).

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v7_secondorder;
