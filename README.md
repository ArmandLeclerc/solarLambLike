# Lamb-like waves in the Sun
This repository holds the data and scripts used to support the article written by A. Le Saux, A. Leclerc, G. Laibe, P. Delplace, A. Venaille.

profile_1D_sun_r15140_Z0p02.data.FGONG is a 1D model of the interior of the Sun.

solar_linear_modes.ipynb is a Jupyter notebook which computes the linear modes using Dedalus (v2).

rotation_kernels.ipynb is a Jupyter notebook which computes the rotation kernel and the frequency splitting of the solar mixed mode. The 
linear modes are already computed and loaded from the file mixedMode_acousticMode_forKernels.npy.

The output of the MUSIC simulation used to support this article is available in release v1.0.0 of this repository. A Python script
named plot_PS_MUSIC.py is provided to open and read these outputs (power spectra). They show the displacements of the f-modes,
which are compared to the ones calculated by the linear theory and read from fmodes_simu_0.15-0.9_ells_4-13_r15140_Z0p02.npy.
