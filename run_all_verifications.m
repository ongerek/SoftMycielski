function run_all_verifications(which)
%RUN_ALL_VERIFICATIONS  Reproduce the numerical claims of the paper.
%
%   run_all_verifications        runs all eleven verification packages
%   run_all_verifications(1:4)   runs a chosen subset
%   run_all_verifications(7)     runs a single package
%
%   Every package is a directory verif_v1 ... verif_v11 containing a
%   runme.m entry point; shared routines live in ./common and are put
%   on the path automatically.  Results are printed to the console:
%   the packages are verification scripts, not a build, so the output
%   IS the deliverable.  Compare it against the tables of the paper
%   (the mapping is in README.md).
%
%   RUNTIME.  The full sweep is dominated by three packages -- v8
%   (Monte-Carlo constants at large horizon), v9 (beyond-Markov
%   testbeds) and v2/v3 (temperature sweeps).  Expect on the order of
%   one to two hours single-threaded; the individual packages are
%   independent and can be run separately.
%
%   REQUIREMENTS.  GNU Octave 8 or later, or MATLAB.  No toolboxes.
%   Package v10 additionally emits figures; its PNG module is skipped
%   automatically when no graphics toolkit is present.
%
%   Tested with GNU Octave 8.4 on Ubuntu 24.04.

  root = fileparts(mfilename('fullpath'));
  addpath(fullfile(root, 'common'));

  if nargin < 1 || isempty(which)
    which = 1:11;
  end

  titles = { ...
    'Abel identity, zero-temperature limit, beta sweeps', ...
    'Critical annealing, self-tuning, M78 and profiles', ...
    'Adaptivity across n, M78-R context restart', ...
    'Cross-entropy scaling and the head-weight atom', ...
    'Bridge representation and pinned-atom diagnostics', ...
    'Critical fixed temperature, rotation identity, argmax laws', ...
    'Second-order misallocation law', ...
    'Model constants by Monte Carlo, large-n anchors', ...
    'Hidden-Markov and piecewise-stationary hybrids', ...
    'Error probability (0-1 loss) and figure generation', ...
    'Worked example behind the explanatory figures'};

  fprintf('\n');
  fprintf('==============================================================\n');
  fprintf('  Soft Mycielski prediction -- numerical verification\n');
  fprintf('  %d package(s) selected; shared routines on the path.\n', ...
          numel(which));
  fprintf('==============================================================\n');

  status = cell(numel(which), 1);
  secs   = zeros(numel(which), 1);
  here   = pwd();

  for ii = 1:numel(which)
    k   = which(ii);
    pkg = sprintf('verif_v%d', k);
    dst = fullfile(root, pkg);
    fprintf('\n--------------------------------------------------------------\n');
    fprintf('  [%2d/%2d]  %s : %s\n', ii, numel(which), pkg, titles{k});
    fprintf('--------------------------------------------------------------\n');
    if ~exist(dst, 'dir')
      fprintf('  ** directory %s not found, skipping\n', pkg);
      status{ii} = 'MISSING';  continue;
    end
    t0 = tic();
    try
      run_one(dst);            % isolated scope: the packages clear all
      status{ii} = 'ok';
    catch err
      status{ii} = 'FAILED';
      fprintf('  ** %s failed: %s\n', pkg, err.message);
    end
    cd(here);
    secs(ii) = toc(t0);
    fprintf('  [%s: %s in %.1f s]\n', pkg, status{ii}, secs(ii));
  end

  fprintf('\n==============================================================\n');
  fprintf('  summary\n');
  fprintf('==============================================================\n');
  for ii = 1:numel(which)
    fprintf('   verif_v%-3d %-8s %8.1f s   %s\n', which(ii), status{ii}, ...
            secs(ii), titles{which(ii)});
  end
  nbad = sum(~strcmp(status, 'ok'));
  fprintf('   ----------------------------------------------------------\n');
  fprintf('   %d of %d package(s) completed, total %.1f s\n', ...
          numel(which) - nbad, numel(which), sum(secs));
  if nbad > 0
    fprintf('   %d package(s) did not complete; see the messages above.\n', ...
            nbad);
  end
  fprintf('\n');
end

% -------------------------------------------------------------------
function run_one(dst)
%RUN_ONE  Execute one package's entry point in a private workspace.
%   The verification scripts begin with "clear all" for reproducibility.
%   Calling them from here means that statement empties this function's
%   workspace rather than the caller's loop state.
  cd(dst);
  runme;
end
