#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as pt
import numpy as np
import sys
import glob
import matplotlib.backends.backend_pdf

        
def plot_avg(inpath, respath, ax, partname="partname", type_="type_", color="color", xleft=60, xright=150, DOF=-1, data=False):
    
    if data:
        s_data, cf_data, std_data = np.loadtxt(inpath + "/"+partname, usecols=(0, 1, 2), unpack=True)
        ax.errorbar(s_data, s_data*s_data*cf_data, yerr=s_data*s_data*std_data, color='k',ls='-', label="Data")
        range_ = np.where(np.logical_and(s_data <=175, s_data >=35))
        ymax = 1.2 * np.max(cf_data[range_] * s_data[range_]**2)
        ymin = 1.2 * np.min(cf_data[range_] * s_data[range_]**2)
        ax.set_ylim([ymin,ymax])

    
    sb, cf_best = np.loadtxt(respath + type_ + "/BAOfit_"+partname + "_best.dat", usecols=(0, 1), unpack=True)
    alpha_m, chi2, evi = np.loadtxt(respath + type_ + "/BAOfit_"+partname + "_mystats.txt", usecols=(0, 1, 2), unpack=True)
    ax.plot(sb, sb * sb * cf_best, color=color, label=r'$\chi^2$=' + str(chi2) + '; ' + str(DOF) + ' DOF\nln(ev)=' + str(evi) + "\n")
     
    
    ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
    ax.set_ylabel(r"s$^2\xi$(s) [h$^{-2}$ Mpc$^2$]")

    #ax.set_ylim([-100,50])
    ax.set_xlim(left=35, right=175)
    
    ax.axvline(x=xleft, color='grey', linestyle='--')
    ax.axvline(x=xright, color='grey', linestyle='--')
    
    ax.legend(bbox_to_anchor=[0.8, 0.74], loc='center', framealpha=0.5)
    #ax.legend(loc='upper right', framealpha=1)


def main():
    pt.rcParams['font.size'] = 14

    inpath_avg500 = "/scratch/variu/clustering/patchy_cmass_subset/lightcone_box1/real/vv2pcf/16R/xi/"
    bestfitpath_avg500 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/avg_fit/"
        
    
    outpath = "/home/astro/variu/"
     
    fig, ax = pt.subplots(2,3, figsize=(30,20))
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][0], partname="avg_500_16R.2pcf", type_="/G1024CIC_40_150_2000/", color="red", xleft=40, xright=150, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][0], partname="avg_500_16R.2pcf", type_="/G1024CIC_40_150_B2000/", color="blue", xleft=40, xright=150)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][0], partname="avg_500_16R.2pcf", type_="/parab_40_150_fast/", color="orange", xleft=40, xright=150)
    
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][1], partname="avg_500_16R.2pcf", type_="/G1024CIC_50_150_2000/", color="red", xleft=50, xright=150, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][1], partname="avg_500_16R.2pcf", type_="/G1024CIC_50_150_B2000/", color="blue", xleft=50, xright=150)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][1], partname="avg_500_16R.2pcf", type_="/parab_50_150_fast/", color="orange", xleft=50, xright=150)

    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][2], partname="avg_500_16R.2pcf", type_="/G1024CIC_60_150_2000/", color="red", xleft=60, xright=150, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][2], partname="avg_500_16R.2pcf", type_="/G1024CIC_60_150_B2000/", color="blue", xleft=60, xright=150)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[0][2], partname="avg_500_16R.2pcf", type_="/parab_60_150_fast/", color="orange", xleft=60, xright=150)
    
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][0], partname="avg_500_16R.2pcf", type_="/G1024CIC_40_170_2000/", color="red", xleft=40, xright=170, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][0], partname="avg_500_16R.2pcf", type_="/G1024CIC_40_170_B2000/", color="blue", xleft=40, xright=170)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][0], partname="avg_500_16R.2pcf", type_="/parab_40_170_fast/", color="orange", xleft=40, xright=170)
    
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][1], partname="avg_500_16R.2pcf", type_="/G1024CIC_50_170_2000/", color="red", xleft=50, xright=170, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][1], partname="avg_500_16R.2pcf", type_="/G1024CIC_50_170_B2000/", color="blue", xleft=50, xright=170)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][1], partname="avg_500_16R.2pcf", type_="/parab_50_170_fast/", color="orange", xleft=50, xright=170)

    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][2], partname="avg_500_16R.2pcf", type_="/G1024CIC_60_170_2000/", color="red", xleft=60, xright=170, data=True)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][2], partname="avg_500_16R.2pcf", type_="/G1024CIC_60_170_B2000/", color="blue", xleft=60, xright=170)
    plot_avg(inpath_avg500, bestfitpath_avg500, ax[1][2], partname="avg_500_16R.2pcf", type_="/parab_60_170_fast/", color="orange", xleft=60, xright=170)
    fig.savefig(bestfitpath_avg500 + "/best_fit_avg_500_16R.2pcf.png")
    
    #pt.show()



if __name__== '__main__':
    main()
