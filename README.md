# Soft-Voting for cryo-EM common-line assessment

This repository contains the MATLAB implementation of Soft-Voting, a continuous-confidence voting method for assessing common lines before spherical-embedding orientation estimation. The three variants evaluated in the manuscript are `linear`, `sigmoid`, and `tanh`.

The release includes a small, runnable example based on precomputed common-line matrices. It does not require downloading raw cryo-EM particles or running a 3D reconstruction.

## Repository contents

| Path | Purpose |
| --- | --- |
| `src/ComputeAnglesFromProjs_SE_wvote_nonlinear.m` | Soft-Voting plus spherical-embedding orientation estimation |
| `src/computeMaxh_nonlinear.m` | Iterative continuous-confidence voting |
| `examples/demo_precomputed_common_lines.m` | Minimum working example for all three variants |
| `examples/data/demo_clcorr_24.mat` | Small precomputed common-line test fixture (24 images) |
| `docs/USER_GUIDE.md` | Input formats, parameters, external dependencies, and troubleshooting |

## Requirements

- MATLAB R2018b was used to test the example. The code also requires the Statistics and Machine Learning Toolbox for `prctile`.
- The [spherical-embedding code of Lu et al.](https://github.com/yluATlzu/3DReconstruction_SE) is an external dependency. Its `Code` directory provides `computeMaxhandSts_ang`, `ComputeAngleCorr`, `elliptic_embed_unitShpere`, `sphricalDist`, `MatchPerpedicularAngles`, `FindGammaAngles`, `calcuRotationMatrix`, and `NormRMSError`. These third-party files are not copied into this repository; retain their original notices and license.
- [ASPIRE-MATLAB 0.14](https://github.com/PrincetonUniversity/aspire) is needed only when starting from projection images (common-line detection or FIRM reconstruction). It is **not** needed for the included precomputed-matrix example.

## Install and run the minimum working example

Clone this repository and the Lu et al. dependency into sibling folders:

```text
workspace/
  Soft-Voting/
  3DReconstruction_SE/
```

Then run in MATLAB:

```matlab
cd('path/to/workspace/Soft-Voting')
addpath('examples')
[R, voted_angles] = demo_precomputed_common_lines();
```

The default is the `linear` variant. To run another variant:

```matlab
[R_sigmoid, voted_sigmoid] = demo_precomputed_common_lines([], 'sigmoid');
[R_tanh, voted_tanh] = demo_precomputed_common_lines([], 'tanh');
```

If the dependency is not in a sibling folder, provide the absolute path to its `Code` directory:

```matlab
[R, voted_angles] = demo_precomputed_common_lines( ...
    'D:\path\to\3DReconstruction_SE\Code', 'linear');
```

The example prints `PASS` with its runtime, image count, and an orthogonality check. It returns estimated orientation transforms of size `3 x 3 x 24` and voted dihedral angles of size `24 x 24`; it does not write output files. The supplied fixture is a 24-image subset of precomputed synthetic common-line data from the EMD-33026/SNR=1 experiment. Its purpose is a fast functionality check, **not** reproduction of the manuscript's 1,000-image accuracy or runtime tables.

## Core function

```matlab
[R, voted_angles] = ComputeAnglesFromProjs_SE_wvote_nonlinear( ...
    clstack, corrstack, p, method, sigma, a);
```

`clstack` and `corrstack` are `K x K` common-line angle and correlation matrices. `clstack(i,j)` is the directed common-line angle on image `i`, in degrees; it need not equal `clstack(j,i)`. `corrstack` holds the corresponding pairwise correlation scores. `p` is a percentile threshold in `[0,100]`. Choose `method` from `linear`, `sigmoid`, or `tanh`; `a` controls the sigmoid/tanh steepness and is ignored for `linear`. `sigma` is retained in the function signature but is not used by these three manuscript variants. See the [user guide](docs/USER_GUIDE.md) for the complete contract and an ASPIRE input example.

The implementation uses weighted voting to update common-line reliability, then uses the final histogram peak heights as weights in the two spherical embeddings. Because cryo-EM common-line reconstruction has a global handedness ambiguity, the returned `3 x 3 x K` orthogonal transforms may use a reflected global branch. Do not compare them directly with ground-truth rotations without a consistent global alignment/handedness convention.

## Experimental data

The class-average data used for the real-data experiments are distributed as ZIP assets in the [data-v1 release](https://github.com/Yikai-Xu628/Soft-Voting/releases/tag/data-v1): EMPIAR-10028 (531 class averages) and EMPIAR-10328 (390 class averages). The original particle datasets are [EMPIAR-10028](https://www.ebi.ac.uk/empiar/EMPIAR-10028/) and [EMPIAR-10328](https://www.ebi.ac.uk/empiar/EMPIAR-10328/). The class-average preparation was described by [Lu et al. (2022)](https://doi.org/10.1038/s42003-022-03255-6) and the same sets were used by [Wang et al. (2024)](https://doi.org/10.1109/TCBB.2024.3476619).

The EMPIAR-10328 class averages were downsampled to `128 x 128` pixels for the reported reconstructions; users should check the image dimensions in the downloaded archive before running a reconstruction. The small fixture in `examples/data` is independent of those real-data archives.

## Scope and attribution

This repository distributes the Soft-Voting extension and a runnable example. The spherical-embedding routines are maintained by Lu et al. and must be obtained from their repository. The included example begins with precomputed common-line matrices; generating those matrices from images and reconstructing volumes additionally require the ASPIRE-MATLAB workflow. See `docs/USER_GUIDE.md` for details.
