#!/usr/bin/env python3
import numpy as np
from getdist import loadMCSamples
import sys, os

if len(sys.argv) != 2:
  print('Usage: {} FILE_ROOT'.format(sys.argv[0]), file=sys.stderr)
  sys.exit(1)

# Setting parameter names and ranges
npar = 4
names = ['alpha', 'B', 'Snl', 'c']
labels = [r'$\alpha$', r'$B$', r'$\Sigma_{\rm nl}$', r'$c$']
lowbound = ['N'] * npar
upbound = ['N'] * npar

# I/O files
fileroot = sys.argv[1]
path, name = os.path.split(fileroot)
if path == '':
  fileroot = './' + fileroot
chains = fileroot + '.txt'
fstats = fileroot + 'stats.dat'
fparam = fileroot + '.paramnames'
frange = fileroot + '.ranges'
ofile = fileroot + 'mystats.txt'
if not os.path.isfile(chains):
  print('Error: cannot access {}'.format(chains), file=sys.stderr)
  sys.exit(1)

np.savetxt(fparam, np.transpose([names, labels]), fmt='%s')
np.savetxt(frange, np.transpose([names, lowbound, upbound]), fmt='%s')

# Load sample from FILE_ROOT.txt
sample = loadMCSamples(fileroot, \
    settings={'fine_bins_2D':1024,'fine_bins':8192})

stats = sample.getMargeStats()
par = stats.parWithName(names[0])
lower = par.limits[0].lower
upper = par.limits[0].upper
sigma = (upper - lower) * 0.5
#best = par.bestfit_sample

# Get the maximum likelihood value
d = np.loadtxt(chains, usecols=(0,2))
sel = (d[:,1] > lower) & (d[:,1] < upper)
d = d[sel,:]
d = d[d[:,1].argsort()]

ngrid = 100
agrid = np.linspace(lower, upper, ngrid+1)
a = (agrid[:-1] + agrid[1:]) * 0.5
L = np.zeros(ngrid)

j = 0
for i in range(len(d)):
  while d[i,1] >= agrid[j+1]:
    j += 1
  L[j] += d[i,0]
best = a[np.argmax(L)]

with open(fstats, "r") as f:
  f.readline()
  line = f.readline()
  evi = float(line.split()[5])

chi2_arr = np.loadtxt(chains, usecols=(1), unpack=True)
min_chi2 = np.min(chi2_arr)

with open(ofile, "w") as f:
  f.write('{0:.5f} {1:.6f} {2:.5f} {3:.5f}'.format(best, sigma, evi, min_chi2))

