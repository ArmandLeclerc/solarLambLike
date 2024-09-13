# Lamb-like waves in the Sun
This repository holds the data and scripts used to support the article written by A. Le Saux, A. Leclerc, G. Laibe, P. Delplace, A. Venaille.

profile_1D_sun_r15140.data.FGONG is a 1D model of the interior of the Sun.

solar_linear_modes.ipynb is a Jupyter notebook which computes the linear modes using Dedalus.

rotation_kernels.ipynb is a Jupyter notebook which computes the rotation kernel and the frequency splitting of the solar mixed mode. The 
linear modes are already computed and loaded from the file mixedMode_acousticMode_forKernels.npy.
