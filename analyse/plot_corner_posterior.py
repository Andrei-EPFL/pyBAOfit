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
  # g.settings.line_styles = "tab20"
  g.triangle_plot(samples, filled='True', legend_labels=legends, \
      line_args=line_args, contour_colors=colors)#colors='tab20')#, contour_colors=colors)

  # Reset axis ticks for the c parameter
  ax = g.subplots[2,2]
  # ax = g.subplots[3,3]

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
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit//"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/avg_fit/"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/range_test/"
  inpath = "//scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/range_test/"
  #fileroot='BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_'
  # fileroot = 'BAOfit_avg_16R.2pcf_' # 'BAOfit_avg_16R_vhxcf.xcf_' #'BAOfit_BigMD_vv.2pcf_' #"BAOfit_avg_16R.2pcf_" # 'BAOfit_BigMD_gv.xcf_' #
  # fileroot = "BAOfit_avg_16R_vhxcf.xcf_"
  fileroot = "BAOfit_avg_1000_16R.2pcf_"
  ofile1 = inpath + fileroot +'_temp2.pdf'
  sample_1 = my_sample(inpath + "/parab_60_150/" + fileroot, 4)
  sample_2 = my_sample(inpath + "/stitched_16R_G2048_50_G512_2000_CG_LC_60_150/" + fileroot, 3)
  sample_3 = my_sample(inpath + "/avg_16R_2000_SK_60_150/" + fileroot, 3)
  
  plot_corner([sample_1, sample_2, sample_3], ofile1, ["green", "magenta", "blue"], ["parab", "CGLC", "SKB"], [{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}, {"lw":0.7, "color":"blue","ls":"--"}])
  

  exit()
  fileroot = "BAOfit_CATALPTCICz0.466G960S1014964334.dat.2pcf_"
  ofile1 = inpath + fileroot +'_temp.pdf'
  
  
  sample_1, sigmanl1 = my_sample(inpath + "/stitched_16R_G2048_50_G512_2000_CG_B/" + fileroot, 3)
  sample_2, sigmanl2 = my_sample(inpath + "/stitched_16R_G2048_50_G512_2000_CG_LC/" + fileroot, 3)
  sample_3, sigmanl3 = my_sample(inpath + "/avg_16R_2000_SK_60_150/" + fileroot, 3)
  sample_4, sigmanl4 = my_sample(inpath + "/parab_60_150_fixc/" + fileroot, 4)
  sample_5, sigmanl5 = my_sample(inpath + "/parab_60_150_gauss/" + fileroot, 4)
  
  

  # plot_corner([sample_4, sample_7], ofile1, ["green", "magenta"], ["18R 60 150", "18R 60 170"], [{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}])
  # plot_corner([sample_1, sample_2, sample_3, sample_4], ofile1, ["green", "magenta", "blue", "orange"], ["parab init", "parab large c", "parablargec2", "parablargec3"], [{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}, {"lw":0.7, "color":"blue","ls":"--"}, {"lw":0.7, "color":"orange","ls":"--"}])
  plot_corner([sample_5, sample_4, sample_1, sample_2, sample_3], [sigmanl5, sigmanl4, sigmanl1, sigmanl2, sigmanl3], ofile1, ["orange", "red", "green", "magenta", "blue"], ["parabg", "fixc", "CGB", "CGLC", "SKB"], [{"lw":0.7, "color":"orange","ls":"--"}, {"lw":0.7, "color":"red","ls":"--"}, {"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}, {"lw":0.7, "color":"blue","ls":"--"}])
  # plot_corner([sample_1, sample_2, sample_3, sample_4, sample_5, sample_6, sample_7, sample_8, sample_9, sample_10, sample_11, sample_12], ofile1, ["green", "magenta", "blue", "green", "magenta", "blue", "green", "magenta", "blue", "red", "orange", "cyan"], ["18RF0.075 40 150", "16RF0.3 40 150", "16RF0.075 40 150", "18RF0.075 60 150", "16RF0.3 60 150", "16RF0.075 60 150", "18RF0.075 60 170", "16RF0.3 60 170", "16RF0.075 60 170", "18RF0.075 40 170", "16RF0.3 40 170", "16RF0.075 40 170"], [{"lw":0.7, "ls":"-"}, {"lw":0.7, "ls":"-"}, {"lw":0.7, "ls":"-"},{"lw":0.7, "ls":"--"}, {"lw":0.7, "ls":"--"}, {"lw":0.7, "ls":"--"},{"lw":0.7, "ls":"-"}, {"lw":0.7, "ls":"-"}, {"lw":0.7, "ls":"-"}, {"lw":0.7, "ls":"--"}, {"lw":0.7, "ls":"--"}, {"lw":0.7, "ls":"--"}])
  # plot_corner([sample_1, sample_2, sample_3, sample_4, sample_5, sample_6, sample_7, sample_8, sample_9, sample_10, sample_11, sample_12], ofile1, ["green", "magenta", "blue", "green", "magenta", "blue", "green", "magenta", "blue", "red", "orange", "cyan"], ["18RF0.075 40 150", "16RF0.3 40 150", "16RF0.075 40 150", "18RF0.075 60 150", "16RF0.3 60 150", "16RF0.075 60 150", "18RF0.075 60 170", "16RF0.3 60 170", "16RF0.075 60 170", "18RF0.075 40 170", "16RF0.3 40 170", "16RF0.075 40 170"], [{"lw":0.7, "color":"green","ls":"-"}, {"lw":0.7, "color":"magenta","ls":"-"}, {"lw":0.7, "color":"blue","ls":"-"},{"lw":0.7, "color":"green","ls":":"}, {"lw":0.7, "color":"magenta","ls":":"}, {"lw":0.7, "color":"blue","ls":":"},{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}, {"lw":0.7, "color":"blue","ls":"--"}, {"lw":0.7, "color":"red","ls":"--"}, {"lw":0.7, "color":"orange","ls":"--"}, {"lw":0.7, "color":"cyan","ls":"--"}])
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
