### Plot the power spectrum of the MUSIC simulation output

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py

### ----------------- as a function of radius ---------------###
filename = "./RvsFreq_solarLL_HSE_cont2dtrec2e2_vel2_0_1e4D_r1-4.h5"

def get_idx(val, vec):
    return (np.abs(np.array(vec)-val)).argmin()


# Load data from MUSIC simulation
f1 = h5py.File(filename,"r")
rads = np.array(f1.get('radius'))
ells = np.array(f1.get('ell_eff'))
freqs = np.array(f1.get('frequency'))
spectrum = np.array(f1.get('spectrum'))

# Load eigenfunction of f modes, computed with Dedalus
file_dedalus = "./fmodes_simu_0.15-0.9_ells_4-13_r15140_Z0p02.npy"
r, xir_4, xih_4, xir_6, xih_6, xir_9, xih_9, xir_11, xih_11, xir_13, xih_13 = np.load(file_dedalus)

# Radius to normalise amplitude to 1
rnorm = 0.66
Rstar = 7e10
idxrnorm = get_idx(rnorm,rads/Rstar)
idxrnormDed = get_idx(rnorm,r)

# Selects frequencies to plot
idx_f1 = get_idx(344.5e-6, freqs) #ell 4.327 HSE mixed mode
idx_f2 = get_idx(422.4e-6, freqs) #ell 6.723 HSE
idx_f3 = get_idx(497.3e-6, freqs) #ell 9.12 HSE
idx_f4 = get_idx(581.7e-6, freqs) #ell 11.518

# Average over 3 freqyency bin
delta_f = 1

# Selects index in effective angular degrees array to plot
idxell1 = 3
idxell2 = 5
idxell3 = 8
idxell4 = 9

# Plots the displacement and eigenfunctions
plt.figure(figsize=(5.8,3.9))
plt.plot(rads/Rstar,np.mean(np.sqrt(spectrum[idx_f1-delta_f:idx_f1+(delta_f+1),:,idxell1]),axis=0)/(np.mean(np.sqrt(spectrum[idx_f1-delta_f:idx_f1+(delta_f+1),idxrnorm,idxell1]),axis=0)),
            color='deepskyblue',label=r'($\omega$, $\tilde{{\ell}}$) = ({:.1f},{:.1f})'.format(freqs[idx_f1]*1e6, ells[idxell1]))
plt.plot(rads/Rstar,np.mean(np.sqrt(spectrum[idx_f2-delta_f:idx_f2+(delta_f+1),:,idxell2]),axis=0)/(np.mean(np.sqrt(spectrum[idx_f2-delta_f:idx_f2+(delta_f+1),idxrnorm,idxell2]),axis=0)),
            color='gold',label=r'($\omega$, $\tilde{{\ell}}$) = ({:.1f},{:.1f})'.format(freqs[idx_f2]*1e6, ells[idxell2]))
plt.plot(rads/Rstar,np.mean(np.sqrt(spectrum[idx_f3-delta_f:idx_f3+(delta_f+1),:,idxell3]),axis=0)/(np.mean(np.sqrt(spectrum[idx_f3-delta_f:idx_f3+(delta_f+1),idxrnorm,idxell3]),axis=0)),
            color='coral',label=r'($\omega$, $\tilde{{\ell}}$) = ({:.1f},{:.1f})'.format(freqs[idx_f3]*1e6, ells[idxell3]))
plt.plot(rads/Rstar,np.mean(np.sqrt(spectrum[idx_f4-delta_f:idx_f4+(delta_f+1),:,idxell4]),axis=0)/(np.mean(np.sqrt(spectrum[idx_f4-delta_f:idx_f4+(delta_f+1),idxrnorm,idxell4]),axis=0)),
            color='red',label=r'($\omega$, $\tilde{{\ell}}$) = ({:.1f},{:.1f})'.format(freqs[idx_f4]*1e6, ells[idxell4]))
plt.plot(r, np.abs(xih_4/(xih_4[idxrnormDed])),':', color='navy',)
plt.plot(r, np.abs(xih_6/(xih_6[idxrnormDed])),':', color='darkgreen')
plt.plot(r, np.abs(xih_9/(xih_9[idxrnormDed])),':', color='brown')
plt.plot(r, np.abs(xih_11/(xih_11[idxrnormDed])),':', color='black')
plt.ylim(-0.05,2.9)
plt.xlim(0.145,0.9)
plt.xlabel(r'$r/R_{\rm star}$')
plt.ylabel(r'$\xi_{\rm h}$/$\xi_{\rm h, 0}$')
plt.legend()

plt.show()


### ----------------- as a function of angular degree and frequency ---------------###
filename = "./FreqVSell_solarLL_HSE_cont2dtrec2e2_vel2_0_2e4D_r0p6avg5.h5"

#load data
f1 = h5py.File(filename,"r")
rads = np.array(f1.get('radius'))
ells = np.array(f1.get('ell_eff'))
freqs = np.array(f1.get('frequency'))
spectrum = np.array(f1.get('spectrum'))

# Values of effective angular degree comes as doublet
# we thus sum these two close values for clarity
ell_eff = []
ell_eff.append(ells[0])
for ll in range(1,len(ells)-1,2):
    ell_eff.append((ells[ll]+ells[ll+1])/2)

new_spectrum = np.ones((len(freqs),len(ell_eff)))
print(np.shape(new_spectrum))
new_spectrum[:,0] = spectrum[:,0]
k = 1
for ll in range(1,len(ells)-1,2):
    new_spectrum[:,k] = (spectrum[:,ll] + spectrum[:,ll])
    k += 1

# Plots the power spectrum
plt.figure(figsize=(5.8,3.9))
plt.pcolormesh(
        np.array(ell_eff),
        freqs*1e6,
        np.array(new_spectrum),
        cmap="inferno",
        norm=mpl.colors.LogNorm(),
        shading="nearest",
        vmin=np.quantile(new_spectrum, 0.15),
        vmax=np.quantile(new_spectrum, .99),
    )
plt.xlim(0.99,40)
plt.ylim(0,1200)
plt.colorbar().set_label(r'$P[\hat{\rm v}_r^2] ~[\rm cm^2.s^{-2}]$')
plt.ylabel(r'$\omega$ ($\rm \mu$Hz)')
plt.xlabel(r'$\tilde{\ell}$')
plt.show()


### ----------------- as a function of radius and frequency ---------------###
filename = "./RvsFreq_solarLL_HSE_cont2dtrec2e2_vel2_0_1e4D_r1-4.h5"

#load data
f1 = h5py.File(filename,"r")
rads = np.array(f1.get('radius'))
ells = np.array(f1.get('ell_eff'))
freqs = np.array(f1.get('frequency'))
spectrum = np.array(f1.get('spectrum'))

# plots the power spectra of the first 20 values of ell effective
for ell_idx in range(20):
    plt.figure(figsize=(6.8,4))
    plt.title(r'$\ell_{{\rm eff}}$ = {:.3f}'.format(ells[ell_idx]))
    plt.pcolormesh(
        freqs*1e6,
        rads,
        np.abs(np.real(spectrum[:,:,ell_idx].T)),
        cmap="inferno",
        norm=mpl.colors.LogNorm(),
        vmin=np.quantile(spectrum[:,:,ell_idx].T , 0.05),
        vmax=np.quantile(spectrum[:,:,ell_idx].T , 0.99),
    )

    plt.colorbar().set_label(r'$P[\hat{\rm v}_i^2] ~[\rm cm^2.s^{-2}]$')
    plt.xlabel(r'Frequency ($\mu$Hz)')
    plt.ylabel(r'r/$R_{\rm star}$')

plt.show()