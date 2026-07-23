%RUNME  Entry point for verification package v5.
%  V5.1-V5.3: bridge representation, endpoint pinning, pinned-atom
%  and small-max law diagnostics;  V5.4: adaptivity ratio at n=32000.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v5_bridge;
v5_ratio32k;
