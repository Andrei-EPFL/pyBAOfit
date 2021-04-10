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


def plot_corner(samples, ofile, colors, legends, line_args):

  g = plots.getSubplotPlotter()
  g.settings.lab_fontsize = 16
  #g.settings.alpha_filled_add = 0.5
  #g.settings.alpha_factor_contour_lines = 0.5
  
  g.triangle_plot(samples, filled='True', legend_labels=legends, \
      line_args=line_args, contour_colors=colors)

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
  inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit//"
  #fileroot='BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_'
  fileroot = "BAOfit_avg_16R.2pcf_" #'BAOfit_BigMD_vv.2pcf_' #'BAOfit_avg_16R_vhxcf.xcf_' #'BAOfit_BigMD_gv.xcf_' #
  
  ofile1 = '/home/astro/variu/' + fileroot +'_fast2.0_R16_getdist_60_150_myst.pdf'
  

  sample_1 = my_sample(inpath + "/G1024CIC_60_150_2000F/" + fileroot, 3)
  sample_2 = my_sample(inpath + "/stitched_subvolumes/" + fileroot, 3)
  sample_3 = my_sample(inpath + "/lin_sav_gol_71_5_stitched_subvolumes/" + fileroot, 3)
  sample_4 = my_sample(inpath + "/stitched_G2048_50_G512_2000/" + fileroot, 3)
  sample_5 = my_sample(inpath + "/sm_stitched_G2048_50_G512_2000/" + fileroot, 3)
  #sample_4 = my_sample(inpath + "/FFTlog/powspecFT_R-scaled2.37-int_flat_0.95/" + fileroot, 3)
  #sample_ = my_sample("../output/BAOfit_daniels2pcf.2pcf_", 3)
  plot_corner([sample_1, sample_2, sample_3, sample_4, sample_5], ofile1, ["green", "magenta", "blue", "orange", "cyan"], ["avg2000", "stitch subvolume", "stitch subvolume sm", "stitch 2000 50", "stitch 2000 50 sm"], [{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}, {"lw":0.7, "color":"blue","ls":"--"}, {"lw":0.7, "color":"orange","ls":"--"}, {"lw":0.7, "color":"cyan","ls":"--"}])
  #plot_corner([sample_], "../output/BAOfit.posterior.pdf", ["green"], ["Test"], [{"lw":0.7, "color":"green","ls":"--"}])

  # sample_list = []
  # legend_list = []
  # line_args = []

  # fracs = ["10","20", "50", "100", "500","1000", "1500"]
  
  # cm = plt.get_cmap('gist_rainbow')
  # NUM_COLORS = len(fracs) + 1
  # color_list=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]

  # for i, frac in enumerate(fracs):
  #   sample_ = my_sample(inpath + "//powspecFT_R-16-int_flat_1.00-N{}/".format(frac) + fileroot, 3)
  #   sample_list.append(sample_)
  #   legend_list.append("powspecFT-R-16-int-flat-1.00-N{}".format(frac))
  #   line_args.append({"lw":0.7, "color":color_list[i],"ls":"--"})
  
  # sample_ = my_sample(inpath + "/mysmoothed_temp_R16_fa/" + fileroot, 3)
  # sample_list.append(sample_)
  # legend_list.append("mysmoothed-temp-R16-fa")
  # line_args.append({"lw":0.7, "color":color_list[7],"ls":"--"})
  

  # plot_corner(sample_list, ofile1, color_list, legend_list, line_args)
  
if __name__== '__main__':
  if len(sys.argv) != 1 and len(sys.argv) != 2:
    print('Usage: python {} [C_MAX]'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)
  main()
