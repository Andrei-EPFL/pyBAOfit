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
  g.triangle_plot(samples, filled='True', legend_labels=['100', '500', '2000 filtered'], \
      line_args=[{'lw':0.7, 'color':'green'}, {'lw':0.7, 'color':'magenta'}, {'lw':0.7, 'color':'orange'}], contour_colors=['green','magenta', 'orange'])
#sample_par, sample_gal, sample_voiGL, sample_voiCIC, sample_voiCIC512
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


def main():
  inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/"
  fileroot='BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_'
  #fileroot='BAOfit_avg_16R.2pcf_'
  
  ofile = '/home/astro/variu/' + fileroot +'getdist_60_150.pdf'
  
  
  sample_gal = my_sample(inpath + "/void_G1024CIC_60_150_mocks_100temp/" + fileroot, 3)
  sample_par = my_sample(inpath + "/void_G1024CIC_60_150_mocks_500temp/" + fileroot, 3)
  sample_voiGL = my_sample(inpath + "/G1024CIC_60_150_m_2000_filtemp/" + fileroot, 3)
  #sample_voiCIC512 = my_sample(inpath + "/void_G512CIC_60_150/" + fileroot, 3)
  
  # sample_gal = my_sample(inpath + "/galaxy_60_150_fast/" + fileroot, 3)
  # sample_par = my_sample(inpath + "/parab_60_150/" + fileroot, 4)
  # sample_voi = my_sample(inpath + "/voidtemplate_conv_60_150/" + fileroot, 3)
  # sample_voi_GL = my_sample(inpath + "/voidtemplate_60_150_G2048L512_1000bins/" + fileroot, 3)
  # sample_voi_GL2 = my_sample(inpath + "/voidtemplate_60_150_G2048L512_1000bins_2/" + fileroot, 3)
  
  #plot_corner([sample_par, sample_gal, sample_voi, sample_voi_GL, sample_voi_GL2], ofile)
  #plot_corner([sample_par, sample_voi, sample_voi_GL2], ofile)
  plot_corner([sample_gal, sample_par,  sample_voiGL], ofile)

if __name__== '__main__':
  if len(sys.argv) != 1 and len(sys.argv) != 2:
    print('Usage: python {} [C_MAX]'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)
  main()
