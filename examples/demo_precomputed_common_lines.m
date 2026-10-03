function [estimated_rots, voted_angles] = demo_precomputed_common_lines(se_code_dir, method)
%DEMO_PRECOMPUTED_COMMON_LINES Run Soft-Voting on a small saved test fixture.
%   [R,A] = DEMO_PRECOMPUTED_COMMON_LINES() uses the linear variant and
%   looks for ../3DReconstruction_SE/Code beside the Soft-Voting repo.
%   [R,A] = DEMO_PRECOMPUTED_COMMON_LINES(SE_CODE_DIR,METHOD) accepts an
%   explicit Lu spherical-embedding Code directory and one of:
%   'linear', 'sigmoid', or 'tanh'.

if nargin < 2 || isempty(method)
    method = 'linear';
end
if ~any(strcmp(method, {'linear', 'sigmoid', 'tanh'}))
    error('Method must be linear, sigmoid, or tanh.');
end

example_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(example_dir);

if nargin < 1 || isempty(se_code_dir)
    candidate = fullfile(fileparts(repo_root), '3DReconstruction_SE', 'Code');
    if isfolder(candidate)
        se_code_dir = candidate;
    else
        existing_fn = which('computeMaxhandSts_ang');
        if ~isempty(existing_fn)
            se_code_dir = fileparts(existing_fn);
        else
            error(['Lu dependency not found. Clone 3DReconstruction_SE beside ', ...
                'this repository or pass the absolute path to its Code folder.']);
        end
    end
end

required = {'computeMaxhandSts_ang.m', 'ComputeAngleCorr.m', ...
    'elliptic_embed_unitShpere.m', 'sphricalDist.m', ...
    'MatchPerpedicularAngles.m', 'FindGammaAngles.m', ...
    'calcuRotationMatrix.m', 'NormRMSError.m'};
for i = 1:numel(required)
    if exist(fullfile(se_code_dir, required{i}), 'file') ~= 2
        error('Missing Lu dependency: %s', fullfile(se_code_dir, required{i}));
    end
end

old_path = path;
restore_path = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(se_code_dir);
addpath(fullfile(repo_root, 'src'));

fixture = load(fullfile(example_dir, 'data', 'demo_clcorr_24.mat'), ...
    'clstack', 'corrstack');
clstack = fixture.clstack;
corrstack = fixture.corrstack;
K = size(clstack, 1);
if ~isequal(size(clstack), [K K]) || ~isequal(size(corrstack), [K K])
    error('The demo fixture must contain matching square matrices.');
end

% These are illustrative demonstration settings, not manuscript tuning.
p = 10;
sigma = 0.1; % Retained in the historical function signature.
a = 10;     % Used by sigmoid/tanh; ignored by linear.

timer = tic;
[estimated_rots, voted_angles] = ...
    ComputeAnglesFromProjs_SE_wvote_nonlinear( ...
    clstack, corrstack, p, method, sigma, a);
elapsed = toc(timer);

if ~isequal(size(estimated_rots), [3 3 K]) || ...
        ~isequal(size(voted_angles), [K K]) || ...
        ~all(isfinite(estimated_rots(:))) || ...
        ~all(isfinite(voted_angles(:)))
    error('Soft-Voting returned an invalid output shape or nonfinite value.');
end

max_orthogonality_error = 0;
max_absolute_determinant_error = 0;
for k = 1:K
    R = estimated_rots(:,:,k);
    max_orthogonality_error = max(max_orthogonality_error, ...
        norm(R' * R - eye(3), 'fro'));
    max_absolute_determinant_error = max( ...
        max_absolute_determinant_error, abs(abs(det(R)) - 1));
end
if max_orthogonality_error > 1e-6 || ...
        max_absolute_determinant_error > 1e-6
    error('Estimated transforms failed the orthogonality check.');
end

fprintf(['PASS: %s, K=%d, elapsed=%.3f s, max orthogonality ', ...
    'error=%.3g, handedness sign=%+d\n'], method, K, elapsed, ...
    max_orthogonality_error, sign(det(estimated_rots(:,:,1))));
end
