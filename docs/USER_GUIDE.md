# Soft-Voting MATLAB user guide

## What the package does

Soft-Voting starts from pairwise common-line estimates and their correlation scores. It uses a uniform-vote initialization, updates the pairwise reliability scores by continuous-confidence voting, and passes the final scores to the two spherical embeddings used for orientation recovery. The method estimates orientations; it does not perform 2D classification, particle picking, map alignment, or FSC calculation.

## Input and output contract

Call:

```matlab
[PredRotations, votedAngle] = ComputeAnglesFromProjs_SE_wvote_nonlinear( ...
    clstack, corrstack, p, method, sigma, a);
```

- `clstack`: real `K x K` matrix. Entry `(i,j)` is the common-line angle on image `i` for image pair `(i,j)`, measured in degrees. The matrix is directed and is generally not symmetric.
- `corrstack`: real `K x K` matrix of pairwise common-line correlation scores. It should be populated for both `(i,j)` and `(j,i)`; diagonal entries are ignored.
- `p`: percentile used for reliability thresholding and embedding edge filtering; the example uses `10`.
- `method`: the manuscript variants are `'linear'`, `'sigmoid'`, and `'tanh'`.
- `a`: steepness of the sigmoid or tanh map. It has no effect for `'linear'`.
- `sigma`: required by the historical function signature. It has no effect for the three manuscript variants. The source also contains older experimental `gaussian` and `power` branches, which are not part of this example.
- `PredRotations`: `3 x 3 x K` orthogonal orientation transforms. A global reflected handedness branch may be selected, so determinants may be `-1` rather than `+1`.
- `votedAngle`: `K x K` voted pairwise dihedral angles in radians.

The `p`, `a`, and `sigma` values in the small example are illustrative. To reproduce a specific manuscript experiment, use its reported parameter settings and input matrices; the 24-image fixture is not a substitute for the 1,000-image synthetic dataset.

## Dependencies and installation

1. Install MATLAB with the Statistics and Machine Learning Toolbox. The example was tested on MATLAB R2018b.
2. Download [Lu et al.'s 3DReconstruction_SE code](https://github.com/yluATlzu/3DReconstruction_SE). Keep its `Code` folder intact.
3. Place `3DReconstruction_SE` beside this repository, or pass the absolute `Code` path to `demo_precomputed_common_lines`.
4. Download [ASPIRE-MATLAB 0.14](https://github.com/PrincetonUniversity/aspire) only if you need to generate common lines from projection images or reconstruct volumes. The precomputed-matrix demonstration does not call ASPIRE.

The MATLAB example adds the Soft-Voting `src` folder and the supplied Lu `Code` folder to the path for the duration of the call, then restores the previous path.

## Run the small example

From the repository root in MATLAB:

```matlab
addpath('examples')
[R, A] = demo_precomputed_common_lines([], 'linear');
size(R)          % [3 3 24]
size(A)          % [24 24]
```

Change `'linear'` to `'sigmoid'` or `'tanh'` to exercise the other variants. The example verifies finite output, orthogonality, and absolute determinant one. It does not evaluate angular error, because there is no global reference alignment in this quick-start test.

## Starting from projection images

For a full experiment, preprocess the images, estimate their common lines with ASPIRE-MATLAB, then call Soft-Voting. The following shows the interface used by the existing experiment scripts; adapt parameters to your data:

```matlab
% ASPIRE-MATLAB 0.14 must already be initialized on the MATLAB path.
projs = ReadMRC('your_class_averages.mrcs');
n = size(projs, 1);
n_theta = 360;
npf = cryo_pft(projs, n, n_theta, 'single');
[clstack, corrstack] = cryo_clmatrix(npf, -1, 1, 5, 1);
corrstack = corrstack + corrstack';

[R, votedAngle] = ComputeAnglesFromProjs_SE_wvote_nonlinear( ...
    clstack, corrstack, 10, 'linear', 0.1, 10);
```

The above is an integration example, not the small automated test. Common-line estimation and 3D reconstruction can take substantially longer than the 24-image demo. The common-line angular resolution, shifts, image centering, and pixel size must be chosen consistently with the dataset and reconstruction workflow.

## Data provenance

The real-data class-average ZIP archives are on the [Soft-Voting data-v1 release](https://github.com/Yikai-Xu628/Soft-Voting/releases/tag/data-v1). These are class averages, not the original particles. The 531 EMPIAR-10028 averages and 390 EMPIAR-10328 averages follow the published preprocessing of Lu et al. (2022) and were also used by Wang et al. (2024). Please cite the original studies when reusing them. The original particles remain available through EMPIAR.

## Troubleshooting

- **Undefined `computeMaxhandSts_ang`, `ComputeAngleCorr`, or `elliptic_embed_unitShpere`:** download the Lu dependency and pass its `Code` directory to the demo, or add that directory to the MATLAB path.
- **Undefined `prctile`:** install/enable the Statistics and Machine Learning Toolbox.
- **`ReadMRC` or `cryo_clmatrix` missing:** initialize ASPIRE-MATLAB 0.14. These functions are not needed for the included precomputed-matrix example.
- **Estimated transforms have determinant `-1`:** the method may select the globally mirrored branch. This is not a failed orthogonality check. Resolve global handedness consistently before comparison with reference rotations or maps.
- **Slow execution on large datasets:** voting evaluates projection triples and the embedding performs iterative optimization. First verify installation with the 24-image fixture; a 1,000-image run is not expected to finish as quickly.
