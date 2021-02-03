#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from getdist import plots, MCSamples, loadMCSamples
import sys, os


def my_sample(fileroot, npar):
  # Setting parameter names and ranges
  if(npar==4):
    names = ['alpha', 'B', 'Snl', 'c']
    labels = [r'$\alpha$', r'$B$', r'$\Sigma_{\rm nl}$', r'$c$']
  elif(npar==3):
    names = ['alpha', 'B', 'Snl']
    labels = [r'$\alpha$', r'$B$', r'$\Sigma_{\rm nl}$']
  else:
    sys.exit()

  lowbound = ['N'] * npar
  upbound = ['N'] * npar
  if len(sys.argv) == 2:
    upbound[npar - 1] = sys.argv[1]

  # I/O files
  path, name = os.path.split(fileroot)
  if path == '':
    fileroot = './' + fileroot
  chains = fileroot + '.txt'
  fparam = fileroot + '.paramnames'
  frange = fileroot + '.ranges'
  if not os.path.isfile(chains):
    print('Error: cannot access {}'.format(chains), file=sys.stderr)
    sys.exit(1)

  np.savetxt(fparam, np.transpose([names, labels]), fmt='%s')
  np.savetxt(frange, np.transpose([names, lowbound, upbound]), fmt='%s')

  plt.rcParams['text.usetex'] = True

  # Load sample from FILE_ROOT.txt
  sample = loadMCSamples(fileroot, \
      settings={'fine_bins_2D':1024,'fine_bins':8192})
  
  return sample


def plot_corner(samples, ofile):

  g = plots.getSubplotPlotter()
  g.settings.lab_fontsize = 16
  #g.triangle_plot(samples, filled='True', legend_labels=['Parabola', 'Galaxy', 'Void', 'VoidGL', 'VoidGL2'], \
  #    line_args=[{'lw':2, 'color':'green'}, {'lw':2, 'color':'blue'}, {'lw':2, 'color':'red'}, {'lw':2, 'color':'orange'}, {'lw':2, 'color':'magenta'}], contour_colors=['green', 'blue', 'red', 'orange', 'magenta'])
  g.triangle_plot(samples, filled='True', legend_labels=['Parabola', 'G1024CIC', 'G1024CICB'], \
      line_args=[{'lw':0.7, 'color':'orange'}, {'lw':0.7, 'color':'red'}, {'lw':0.7, 'color':'blue'}], contour_colors=['orange','red', 'blue'])
  #sample_par, sample_gal, sample_voiGL, sample_voiCIC, sample_voiCIC512
  # Reset axis ticks for the c parameter
  ax = g.subplots[3,3]
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
  for ax in g.subplots[3,:3]:
    ax.set_yticks(np.round(np.linspace(xmin_new, xmax_new, nx), 0))

  for ax in g.subplots[:,0]:
    ax.axvline(1, color='gray', ls='--')

  #g.export(ofile)
  g.fig.savefig(ofile, bbox_inches='tight', transparent=True)


def main():
  inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/avg_fit/"
  #fileroot='BAOfit_voids_Box_z0.465600_16R.2pcf_'
  fileroot='BAOfit_avg_500_16R.2pcf_'
  
  range_list = ["50_150", "50_170"]
  for range_ in range_list:

    ofile = inpath + fileroot +range_+'_getdist_60_150.pdf'
    sample_par = my_sample(inpath + "/parab_"+range_+"_fast/" + fileroot, 4)
    sample_voiCIC = my_sample(inpath + "/G1024CIC_"+range_+"_2000/" + fileroot, 3)
    sample_voiCICB = my_sample(inpath + "/G1024CIC_"+range_+"_B2000/" + fileroot, 3)
    plot_corner([sample_par, sample_voiCIC, sample_voiCICB], ofile)

if __name__== '__main__':
  if len(sys.argv) != 1 and len(sys.argv) != 2:
    print('Usage: python {} [C_MAX]'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)
  main()
