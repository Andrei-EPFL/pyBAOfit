#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as pt
import numpy as np
import sys
import glob
import matplotlib.backends.backend_pdf

        
def plot_avg(inpath, respath, ax, partname="partname", type_="type_", color="color", xleft=60, xright=150, DOF=-1, data=False, col2=1, label="label"):
    
    if data:
        s_data, cf_data, std_data = np.loadtxt(inpath + "/"+partname, usecols=(0, 1, 2), unpack=True)
        ax.errorbar(s_data, s_data*s_data*cf_data, yerr=s_data*s_data*std_data, color='k',ls='-', label="Data")

        #s_data, cf_data = np.loadtxt(inpath + "/"+partname, usecols=(0, col2), unpack=True)
        #ax.plot(s_data, s_data*s_data*cf_data, color='k',ls='-', label="Data")
        range_ = np.where(np.logical_and(s_data <=175, s_data >=35))
        ymax = 1.2 * np.max(cf_data[range_] * s_data[range_]**2)
        ymin = 1.2 * np.min(cf_data[range_] * s_data[range_]**2)
        ax.set_ylim([ymin,ymax])

    
    sb, cf_best = np.loadtxt(respath + type_ + "/BAOfit_"+partname + "_best.dat", usecols=(0, 1), unpack=True)
    #alpha_m, sigma, evi, chi2 = np.loadtxt(respath + type_ + "/BAOfit_"+partname + "_mystats.txt", usecols=(0, 1, 2, 3), unpack=True)
    ax.plot(sb, sb * sb * cf_best, color=color, label=label)#+ '\nln(ev)=' + '{:.4g}'.format(evi) +"\n"+ r'$\chi^2$=' + '{:.4g}'.format(chi2/DOF) +r"; $\nu=$" +'{}'.format(DOF))# + '; ' + '{:.4g}'.format(alpha_m) + ' alpha\nln(ev)=' + '{:.4g}'.format(evi) + "\n" + label)
     
    
    ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
    ax.set_ylabel(r"s$^2\xi$(s) [h$^{-2}$ Mpc$^2$]")

    #ax.set_ylim([-100,50])
    ax.set_xlim(left=35, right=175)
    
    ax.axvline(x=xleft, color='grey', linestyle='--')
    ax.axvline(x=xright, color='grey', linestyle='--')
    
    ax.legend(bbox_to_anchor=[0.8, 0.74], loc='center', framealpha=0.5)
    #ax.legend(loc='upper right', framealpha=1)


def test_range():

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


def chengsplots():
    bestfitpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/"
    inpath = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/"
    opath = "/home/astro/variu/chengsplots/"

    fig, ax = pt.subplots(figsize=(10,10))
    plot_avg(inpath+ "vv2pcf/16R/", bestfitpath+"vv2pcf/avg_fit/", ax, partname="avg_16R.2pcf", type_="/parab_60_150_fast2F/", color="green", xleft=60, xright=150, data=True, label="Parabolic", DOF=11)
    plot_avg(inpath+ "vv2pcf/16R/", bestfitpath+"vv2pcf/avg_fit/", ax, partname="avg_16R.2pcf", type_="/galaxy_60_150/", color="magenta", xleft=60, xright=150, data=False, label="Galaxy", DOF=12)
    plot_avg(inpath+ "vv2pcf/16R/", bestfitpath+"vv2pcf/avg_fit/", ax, partname="avg_16R.2pcf", type_="/G1024CIC_60_150_2000F/", color="blue", xleft=60, xright=150, data=False, label="Template", DOF=12)
    fig.tight_layout()
    fig.savefig(opath + "best_fit_vv2pcf_avg_16R.2pcf.png")

    fig, ax = pt.subplots(figsize=(10,10))
    plot_avg(inpath+ "vhxcf/16R/", bestfitpath+"vhxcf/avg_fit/", ax, partname="avg_16R_vhxcf.xcf", type_="/parab_60_150_fast2F/", color="green", xleft=60, xright=150, data=True, label="Parabolic", DOF=11)
    plot_avg(inpath+ "vhxcf/16R/", bestfitpath+"vhxcf/avg_fit/", ax, partname="avg_16R_vhxcf.xcf", type_="/galaxy_60_150/", color="magenta", xleft=60, xright=150, data=False, label="Galaxy", DOF=12)
    plot_avg(inpath+ "vhxcf/16R/", bestfitpath+"vhxcf/avg_fit/", ax, partname="avg_16R_vhxcf.xcf", type_="/G1024CIC_60_150_2001F/", color="blue", xleft=60, xright=150, data=False, label="Template", DOF=12)
    fig.tight_layout()
    fig.savefig(opath + "best_fit_vhxcf_avg_16R_vhxcf.xcf.png")

    # inpath = "/hpcstorage/variu/BigMD/clustering/"
    # fig, ax = pt.subplots(figsize=(10,10))
    # plot_avg(inpath, bestfitpath+"vv2pcf/bigmd_fit/", ax, partname="BigMD_vv.2pcf", type_="/parab_60_150_fast2F/", color="green", xleft=60, xright=150, data=True, col2=3, label="Parabolic", DOF=11)
    # plot_avg(inpath, bestfitpath+"vv2pcf/bigmd_fit/", ax, partname="BigMD_vv.2pcf", type_="/galaxy_60_150/", color="magenta", xleft=60, xright=150, data=False, col2=3, label="Galaxy", DOF=12)
    # plot_avg(inpath, bestfitpath+"vv2pcf/bigmd_fit/", ax, partname="BigMD_vv.2pcf", type_="/G1024CIC_60_150_2000F/", color="blue", xleft=60, xright=150, data=False, col2=3, label="Template", DOF=12)
    # fig.tight_layout()
    # fig.savefig(opath + "best_fit_vv2pcf_BigMD_vv.2pcf.png")

    # fig, ax = pt.subplots(figsize=(10,10))
    # plot_avg(inpath, bestfitpath+"vhxcf/bigmd_fit/", ax, partname="BigMD_gv.xcf", type_="/parab_60_150_fast2F/", color="green", xleft=60, xright=150, data=True, col2=3, label="Parabolic", DOF=11)
    # plot_avg(inpath, bestfitpath+"vhxcf/bigmd_fit/", ax, partname="BigMD_gv.xcf", type_="/galaxy_60_150/", color="magenta", xleft=60, xright=150, data=False, col2=3, label="Galaxy", DOF=12)
    # plot_avg(inpath, bestfitpath+"vhxcf/bigmd_fit/", ax, partname="BigMD_gv.xcf", type_="/G1024CIC_60_150_2001F/", color="blue", xleft=60, xright=150, data=False, col2=3, label="Template", DOF=12)
    # fig.tight_layout()
    # fig.savefig(opath + "best_fit_vhxcf_BigMD_gv.xcf.png")

def main():
    pt.rcParams['font.size'] = 18
    sdata, xidata = np.loadtxt("/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/avg_16R.2pcf", usecols=(0,1), unpack=True)
    sbest1, xibest1 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit//G1024CIC_60_150_2000F/BAOfit_avg_16R.2pcf_best.dat", usecols=(0,1),unpack=True)
    sbest2, xibest2 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/stitched_subvolumes/BAOfit_avg_16R.2pcf_best.dat", usecols=(0,1),unpack=True)
    sbest3, xibest3 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/lin_sav_gol_71_5_stitched_subvolumes/BAOfit_avg_16R.2pcf_best.dat", usecols=(0,1),unpack=True)
    sbest4, xibest4 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/stitched_G2048_50_G512_2000/BAOfit_avg_16R.2pcf_best.dat", usecols=(0,1),unpack=True)
    sbest5, xibest5 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/sm_stitched_G2048_50_G512_2000/BAOfit_avg_16R.2pcf_best.dat", usecols=(0,1),unpack=True)
    
    pt.plot(sdata, xidata*sdata*sdata, color="k", marker=".", ls='-')
    pt.plot(sbest1, xibest1*sbest1*sbest1, color="green", label=f"avg2000: median={1.00017}; avg={1.00008}; sigma={0.006809}", ls="--")
    pt.plot(sbest2, xibest2*sbest2*sbest2, color="magenta", label=f"stitch subvolume: median={1.00032}; avg={1.00055}; sigma={0.007015}", ls="--")
    pt.plot(sbest3, xibest3*sbest3*sbest3, color="blue", label=f"stitch subvolume smooth: median={1.00024}; avg={1.00034}; sigma={0.007070}", ls="--")
    pt.plot(sbest4, xibest4*sbest4*sbest4, color="orange", label=f"stitch G2048; G512: median={1.00026}; avg={1.00025}; sigma={0.007110}", ls="--")
    pt.plot(sbest5, xibest5*sbest5*sbest5, color="cyan", label=f"stitch G2048; G512 smooth: median={1.00011}; avg={1.00032}; sigma={0.007005}", ls="--")
    pt.xlim([55, 155])
    pt.ylim([-10, 26])
    pt.axvline(60, color="grey", ls="--")
    pt.axvline(150, color="grey", ls="--")
    pt.legend()
    pt.show()
    exit()
    chengsplots()
    exit()

    inpath_avg500 = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/"
    bestfitpath_avg500 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/"
    partname = "avg_16R.2pcf"
    outpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/cases_best_fit_avg_500_16R.2pcf"
    

    ### powspecFT_R-16-int_flat_1.00-N{}
    fig, ax = pt.subplots(figsize=(20,20))
    fracs4 = ["10","20", "50", "100", "500","1000", "1500"]
    cm = pt.get_cmap('gist_rainbow')
    NUM_COLORS = len(fracs4) + 1
    ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    for i, frac in enumerate(fracs4):
        if i == 0:
            plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspecFT_R-16-int_flat_1.00-N{}/".format(frac), color="red", xleft=60, xright=150, data=True)
        else:
            plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspecFT_R-16-int_flat_1.00-N{}/".format(frac), color="red", xleft=60, xright=150, data=False)

    plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/mysmoothed_temp_R16_fa/", color="red", xleft=60, xright=150, data=False)

    fig.savefig(outpath + "_4.png")
    

    # fig, ax = pt.subplots(figsize=(20,20))
    
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/fast_2.0/powspec_direct_R16/", color="green", xleft=60, xright=150, data=True)
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/fast_2.0/powspecFT_R-scaled2.37-int_flat_{}/".format("1.00-R16"), color="violet", xleft=60, xright=150, data=False)
    
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/fast_2.0/mysmoothed_temp_R16/", color="blue", xleft=60, xright=150, data=False)

    # #plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/fast_2.0/powspec_R-scaled2.37-int_flat_{}/".format("0.95"), color="green", xleft=60, xright=150, data=True)
    # #plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/fast_2.0/powspecFT_R-scaled2.37-int_flat_{}/".format("0.95"), color="violet", xleft=60, xright=150, data=False)
    # fig.savefig(outpath+"_fast_2.0_my.png")
    
    
    # fig, ax = pt.subplots(figsize=(20,20))
    
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/FFTlog/powspec_direct_R16/", color="blue", xleft=60, xright=150, data=True)
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/FFTlog/powspecFT_R-scaled2.37-int_flat_{}/".format("1.00-R16"), color="orange", xleft=60, xright=150, data=False)

    # #plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/FFTlog/powspec_R-scaled2.37-int_flat_{}/".format("0.95"), color="blue", xleft=60, xright=150, data=False)
    # #plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/FFTlog/powspecFT_R-scaled2.37-int_flat_{}/".format("0.95"), color="orange", xleft=60, xright=150, data=False)
    # plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/FFTlog/mysmoothed_temp_R16/", color="blue", xleft=60, xright=150, data=False)

    # fig.savefig(outpath+"_FFTlog_my.png")
    
    # exit()

    # ### powspec_R-dr{}-int
    # fracs = ["0.75", "0.80", "0.85", "0.90", "0.95", "1.00", "1.05", "1.10", "1.15"]
    # cm = pt.get_cmap('gist_rainbow')
    # NUM_COLORS = len(fracs)
    # ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    # for i, frac in enumerate(fracs):
    #     if i == 0:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspec_R-dr{}-int/".format(frac), color="red", xleft=60, xright=150, data=True)
    #     else:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspec_R-dr{}-int/".format(frac), color="red", xleft=60, xright=150, data=False)

    # fig.savefig(outpath+".png")
    ###

    # ### powspec_R-scaled2.37-int_flat_{}
    # fig, ax = pt.subplots(figsize=(20,20))
    # #fracs2 = ["0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40","0.45","0.50","0.55", "0.60", "0.65", "0.70","0.75", "0.80", "0.85", "0.90", "0.95"]
    # fracs2 = ["0.50","0.55", "0.60", "0.65", "0.70","0.75", "0.80", "0.85", "0.90", "0.95"]
    # cm = pt.get_cmap('gist_rainbow')
    # NUM_COLORS = len(fracs2)
    # ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    # for i, frac in enumerate(fracs2):
    #     if i == 0:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspec_R-scaled2.37-int_flat_{}/".format(frac), color="red", xleft=60, xright=150, data=True)
    #     else:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspec_R-scaled2.37-int_flat_{}/".format(frac), color="red", xleft=60, xright=150, data=False)

    # fig.savefig(outpath + "_2.png")
    # ###

    # ### powspecFT_R-scaled2.37-int_flat_{}
    # fig, ax = pt.subplots(figsize=(20,20))
    # #fracs3 = ["0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40","0.45","0.50","0.55", "0.60", "0.65", "0.70","0.75", "0.80", "0.85", "0.90", "0.95", "1.00", "1.00-R16"]
    # fracs3 = ["0.50","0.55", "0.60", "0.65", "0.70","0.75", "0.80", "0.85", "0.90", "0.95","1.00-R16"]
    # cm = pt.get_cmap('gist_rainbow')
    # NUM_COLORS = len(fracs3)
    # ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    # for i, frac in enumerate(fracs3):
    #     if i == 0:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspecFT_R-scaled2.37-int_flat_{}/".format(frac), color="red", xleft=60, xright=150, data=True)
    #     else:
    #         plot_avg(inpath_avg500, bestfitpath_avg500, ax, partname=partname, type_="/powspecFT_R-scaled2.37-int_flat_{}/".format(frac), color="red", xleft=60, xright=150, data=False)

    # fig.savefig(outpath + "_3.png")
    # ###



if __name__== '__main__':
    main()
