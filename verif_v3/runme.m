%RUNME  Entry point for verification package v3.
%  V3.1/V3.2: adaptivity across n and binned mechanism; V3.3: M78-R trie.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v3_anneal_self;
v3_trie;
