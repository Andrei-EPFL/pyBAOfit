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
    names = ['alpha', 'B_V', 'Snl_VV', 'c_VV']
    labels = ["\\alpha", "B_{\\mathrm{V}}", "\\Sigma_{\\mathrm{nl}, VV}", "c_{VV}"]
  elif(npar==3):
    names = ['alpha', 'B_V', 'Snl_VV']
    labels = ["\\alpha", "B", "\\Sigma_{\\mathrm{nl}}"]
  elif(npar==6):
    labels = ["\\alpha", "B_{\\mathrm{V}}", "\\Sigma_{\\mathrm{nl}, VV}", "a_{0, {\\mathrm{VV}}}", "a_{1, {\\mathrm{VV}}}", "a_{2, {\\mathrm{VV}}}"]
    names = ["alpha", "B_V", "Snl_VV","a0_VV", "a1_VV", "a2_VV"] 
  elif(npar==7):
    labels = ["\\alpha", "B_{\\mathrm{V}}", "\\Sigma_{\\mathrm{nl}, VV}", "c_{VV}", "a_{0, {\\mathrm{VV}}}", "a_{1, {\\mathrm{VV}}}", "a_{2, {\\mathrm{VV}}}"]
    names = ["alpha", "B_V", "Snl_VV", "c_VV","a0_VV", "a1_VV", "a2_VV"] 
  else:
    print("choose npar correctly")
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

  # plt.rcParams['text.usetex'] = True

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
  # ax = g.subplots[2,2]
  # # ax = g.subplots[3,3]

  # xmin, xmax = ax.get_xlim()
  # dx = xmax - xmin
  # dx_mag = 10**(int(np.log10(dx)))
  # if dx / dx_mag > 5:
  #   dx_mag *= 2
  # elif dx / dx_mag < 2:
  #   dx_mag /= 2.5
  # elif dx / dx_mag < 3:
  #   dx_mag /= 2
  # xmin_new = int(xmin / dx_mag) * dx_mag
  # xmax_new = int(xmax / dx_mag) * dx_mag
  # nx = int(np.round((xmax_new - xmin_new) / dx_mag + 1))
  # if (xmax - xmax_new) / dx < 0.05:
  #   xmax_new -= dx_mag
  #   nx -= 1
  # if (xmin_new - xmin) / dx < 0.05:
  #   xmin_new += dx_mag
  #   nx -= 1
  # ax.set_xticks(np.round(np.linspace(xmin_new, xmax_new, nx), 0))
  # for ax in g.subplots[2,:2]:
  #   ax.set_yticks(np.round(np.linspace(xmin_new, xmax_new, nx), 0))

  # for ax in g.subplots[:,0]:
  #   ax.axvline(1, color='gray', ls='--')

  #g.export(ofile)
  g.fig.savefig(ofile, bbox_inches='tight', transparent=False)


def main():
  f = mpl.font_manager.findSystemFonts(fontpaths="/home/astro/variu//fonts/")
  # print(f)
  for font_file in f:
      mpl.font_manager.fontManager.addfont(font_file)

  plt.rcParams.update({'font.family': "serif"})
  plt.rcParams.update({'font.serif': "Times New Roman"})
  mpl.rcParams.update({'font.size': 16})

  mpl.rcParams['mathtext.fontset'] = 'cm'
  mpl.rcParams['mathtext.rm'] = 'serif'

  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit//"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/avg_fit/"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/range_test/"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/sickle_realisations_fit/"
  # inpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/sickle_realisations_fit/avg/vhxcf/"
  #fileroot='BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_'
  # fileroot = 'BAOfit_avg_16R.2pcf_' # 'BAOfit_avg_16R_vhxcf.xcf_' #'BAOfit_BigMD_vv.2pcf_' #"BAOfit_avg_16R.2pcf_" # 'BAOfit_BigMD_gv.xcf_' #
  # fileroot = "BAOfit_avg_16R_vhxcf.xcf_"
  # fileroot = "BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_"
  # fileroot = "BAOfit_avg_16R_vhxcf.xcf_"
  # fileroot = "BAOfit_CATALPTCICz0.466G960S1014964334.dat.2pcf_"
  
  # ofile1 = './' + fileroot +'_temp2.pdf'
  # sample_1 = my_sample(inpath0 + fileroot, 4)
  # sample_2 = my_sample(inpath1 + fileroot, 4)
  # sample_3 = my_sample(inpath + "/avg_16R_2000_SK_60_150/" + fileroot, 3)
  
  # plot_corner([sample_1, sample_2], ofile1, ["green", "magenta"], ["parab0", "parab1"], [{"lw":0.7, "color":"green","ls":"--"}, {"lw":0.7, "color":"magenta","ls":"--"}])
  

  # exit()
  
  # ofile = "/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/BAOfit_avg_sc_parab_nuisance.png"
  # models = ["/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/parab/BAOfit_avg_sc_", "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/range_test/scaled_covariance/parab_60_150/BAOfit_avg_16R.2pcf_"]
  
  
  # ofile = "/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/BAOfit_CATALPTCICz0.466G960S1057271749.VOID_parab_nuisance.png"
  # models = ["/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/parab/BAOfit_CATALPTCICz0.466G960S1057271749.VOID_", "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/parab_60_150_fixc/BAOfit_CATALPTCICz0.466G960S1057271749.VOID.dat.2pcf_"]
  
  # models = ["/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/cg/BAOfit_CATALPTCICz0.466G960S1057271749.VOID_", "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG/BAOfit_CATALPTCICz0.466G960S1057271749.VOID.dat.2pcf_"]
  # ofile = "/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/BAOfit_CATALPTCICz0.466G960S1057271749.VOID_CG_nuisance.png"
  
  models = ["/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/cg/BAOfit_avg_sc_", "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/range_test/vv2pcf/scaled_covariance/stitched_16R_G2048_50_G512_2000_CG_60_150/BAOfit_avg_16R.2pcf_"]
  ofile = "/home/astro/variu/phd/voids/chengscodes/BAOflit_new/tmp/BAOfit_avg_sc_CG_nuisance.png"
  
  
  sample_list = []
  line_args = []
  legend_list = ["C", "py"]
  color_list = ["red", "blue"]
  # nparam = [6, 4]
  nparam = [6, 3]

  # # cm = plt.get_cmap('gist_rainbow')
  # # NUM_COLORS = len(models)
  # # color_list = [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)]
  
  for i, model in enumerate(models):
    sample_ = my_sample(model, nparam[i])
    sample_list.append(sample_)
    line_args.append({"lw":0.7, "color":color_list[i],"ls":"--"})
  
  
  plot_corner(sample_list, ofile, color_list, legend_list, line_args)
  
if __name__== '__main__':
  if len(sys.argv) != 1 and len(sys.argv) != 2:
    print('Usage: python {} [C_MAX]'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)
  main()
