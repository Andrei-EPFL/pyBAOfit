#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from getdist import plots, MCSamples, loadMCSamples
import sys, os

if len(sys.argv) != 2 and len(sys.argv) != 3:
  print('Usage: {} FILE_ROOT [C_MAX]'.format(sys.argv[0]), file=sys.stderr)
  sys.exit(1)

# Setting parameter names and ranges
npar = 3
names = ['alpha', 'B', 'Snl']
labels = [r'$\alpha$', r'$B$', r'$\Sigma_{\rm nl}$']
lowbound = ['N'] * npar
upbound = ['N'] * npar
if len(sys.argv) == 3:
  upbound[npar - 1] = sys.argv[2]

# I/O files
fileroot = sys.argv[1]
path, name = os.path.split(fileroot)
if path == '':
  fileroot = './' + fileroot
chains = fileroot + '.txt'
fparam = fileroot + '.paramnames'
frange = fileroot + '.ranges'
ofile = fileroot + 'getdist.pdf'
if not os.path.isfile(chains):
  print('Error: cannot access {}'.format(chains), file=sys.stderr)
  sys.exit(1)

np.savetxt(fparam, np.transpose([names, labels]), fmt='%s')
np.savetxt(frange, np.transpose([names, lowbound, upbound]), fmt='%s')

plt.rcParams['text.usetex'] = True

# Load sample from FILE_ROOT.txt
sample = loadMCSamples(fileroot, \
    settings={'fine_bins_2D':1024,'fine_bins':8192})

g = plots.getSubplotPlotter()
g.settings.lab_fontsize = 16
g.triangle_plot(sample, filled='True', \
    line_args={'lw':2,'color':'#006FED'})

print('Statistics from getdist ...')
stats = sample.getMargeStats()
best = np.zeros(npar)
lower = np.zeros(npar)
upper = np.zeros(npar)
mean = np.zeros(npar)
sigma = np.zeros(npar)
for i in range(npar):
  par = stats.parWithName(names[i])
  best[i] = par.bestfit_sample
  mean[i] = par.mean
  sigma[i] = par.err
  lower[i] = par.limits[0].lower
  upper[i] = par.limits[0].upper
  print('{0:s}: {1:.5f} + {2:.6f} - {3:.6f}, or {4:.5f} +- {5:.6f}'.format( \
        names[i], best[i], upper[i]-best[i], best[i]-lower[i], mean[i], \
        sigma[i]))

#  ax = g.subplots[i,i]
#  ax.axvline(best[i], color='gray', ls='--', lw=0.5)
#  ax.axvline(lower[i], color='silver', ls=':', lw=0.5)
#  ax.axvline(upper[i], color='silver', ls=':', lw=0.5)

print('My statistics ...')
d = np.loadtxt(chains)
for i in range(npar):
  j = i + 2
  b = list(zip(d[:,0], d[:,j]))
  b.sort(key=lambda x: x[1])
  b = np.array(b)
  b[:,0] = b[:,0].cumsum()
  sig1 = 0.5 + 0.6826 / 2.
  bi = lambda x: np.interp(x, b[:,0], b[:,1], left=b[0,1], right=b[-1,1])

  low1 = bi(1 - sig1)
  high1 = bi(sig1)
  median = bi(0.5)
  print('{0:s}: {1:.5f} + {2:.6f} - {3:.6f}'.format( \
        names[i], median, high1-median, median-low1))

#  ax = g.subplots[i,i]
#  ax.axvline(median, color='blue', ls='--', lw=1)
#  ax.axvline(low1, color='red', ls=':', lw=1)
#  ax.axvline(high1, color='red', ls=':', lw=1)

# Reset axis ticks for the c parameter
ax = g.subplots[2,2]
xmin, xmax = ax.get_xlim()
dx = xmax - xmin
dx_mag = 10**(int(np.log10(dx)))
if dx / dx_mag > 5:
  dx_mag *= 2
elif dx / dx_mag < 2:
  dx_mag /= 2.5
elif dx / dx_mag < 3:
  dx_mag /= 2
xmin_new = int(xmin / dx_mag) * dx_mag
xmax_new = int(xmax / dx_mag) * dx_mag
nx = int(np.round((xmax_new - xmin_new) / dx_mag + 1))
if (xmax - xmax_new) / dx < 0.05:
  xmax_new -= dx_mag
  nx -= 1
if (xmin_new - xmin) / dx < 0.05:
  xmin_new += dx_mag
  nx -= 1
ax.set_xticks(np.round(np.linspace(xmin_new, xmax_new, nx), 0))
for ax in g.subplots[2,:2]:
  ax.set_yticks(np.round(np.linspace(xmin_new, xmax_new, nx), 0))

for ax in g.subplots[:,0]:
  ax.axvline(1, color='gray', ls='--')

#g.export(ofile)
g.fig.savefig(ofile, bbox_inches='tight', transparent=True)

