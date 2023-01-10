import numpy as np
import matplotlib.pyplot as pt 
from scipy.interpolate import interp1d


linear = "/home/astro/variu/phd/voids/chengscodes/CosmoGAME/EisensteinHu_Pnw.txt"
k_lin_nw, pk_linw_nw = np.loadtxt(linear, usecols=(0, 1), unpack=True)
fint_lin_nw = interp1d(np.log(k_lin_nw), np.log(pk_linw_nw), kind="cubic")
f_lin_nw = lambda k: np.exp(fint_lin_nw(np.log(k)))


k0, pk0 = np.loadtxt('/scratch/variu/clustering/patchy_cmass_subset/box1/redshift_recon/pkvoids/avg_495_recon.pspec', usecols=(0, 1), unpack=True)
k1, pk1 = np.loadtxt('/scratch/variu/clustering/patchy_cmass_subset/box1/real/pkvoids_g/avg_500.pspec', usecols=(0, 1), unpack=True)
k2, pk2 = np.loadtxt("/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.pspec", usecols=(0, 1), unpack=True)
k3, pk3 = np.loadtxt("/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_voids/stitched_G2048_50_G512_2000.txt", usecols=(0, 1), unpack=True)

range_2 = k3 <= np.max(k0)
k2t, pk2t = k3[range_2], pk3[range_2]

range_ = slice(0, k2.size, 5)
k2, pk2 = k2[range_], pk2[range_]

range_ = slice(0, k3.size, 5)
k3, pk3 = k3[range_], pk3[range_]

fig, ax = pt.subplots(2, 1, figsize=(5, 4), sharex=True, gridspec_kw={"hspace":0, "height_ratios":[3, 1]})

ax[0].plot(k0, pk0 / f_lin_nw(k0), label="post-recon", color="grey")
ax[0].plot(k1, 1.2 * pk1 / f_lin_nw(k1), label="pre-recon", color="k")
ax[0].plot(k2, 1.45 * pk2 / f_lin_nw(k2), label="CG", color="red", ls="", marker="o", ms=1.5)
ax[0].plot(k3, 1.3 * pk3 / f_lin_nw(k3), label="SK", color="blue", ls="", marker="x", ms=1.5)

ax[0].set_xlim([0, 0.6])
ax[0].set_ylim([0, 17])

ax[1].plot(k2t, pk0  / pk2t)

ax[0].legend()
fig.savefig("/home/astro/variu/test.png")
# pt.show()   