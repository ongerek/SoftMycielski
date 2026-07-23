%RUNME  Entry point for verification package v11.
%  V11.1: verified data behind the explanatory figures --
%    match lengths and soft distributions on x = "abracadabra",
%    the LZ78 parse/trie and M78-R suffix nodes, and the realized
%    fluctuation path with its chord and bridge.
%  See ../figures/checkfig.py for the companion geometry checker,
%  which reports text overlaps and extents in a compiled figure.

% Make the shared routines in ../common visible, whichever
% directory this script is launched from.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
v11_figdata;
