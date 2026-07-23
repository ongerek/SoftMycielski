%RUNME  Entry point for verification package v2.
%  V2.1: critical-annealing law beta*(n) -> H from below;
%  V2.2: self-critical temperature schedules;
%  V2.3/V2.4: Mycielski-78 streaming predictor and order-weight
%             profiles across the phase boundary.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
anneal_test;
selftune_test;
m78_and_profile_test;
