import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pt
import numpy as np
import sys
import glob
import os
from matplotlib.backends.backend_pdf import PdfPages

        
def myplot(bestpath="best", datapath="test", type_="Test", outpath="test"):
    files = glob.glob(datapath + "/CA*2pcf")
    print(len(files))
    
    pp = PdfPages(outpath)
    std = np.loadtxt(datapath + "/avg_1000_16R.2pcf", usecols=(2), unpack=True)

    for i, file_ in enumerate(files):
        print(i)
        sd, cf_data = np.loadtxt(file_, usecols=(0, 1), unpack=True)
        filename = os.path.basename(file_)
    
        sb, cf_best = np.loadtxt(bestpath + f"/{type_}/BAOfit_"+filename + "_best_nw.dat", usecols=(0, 1), unpack=True)
        sigma, alpha, bias, sigmanl = np.loadtxt(bestpath + f"/{type_}/BAOfit_"+filename + "_mystats.txt", usecols=(1, 5, 6, 7), unpack=True)
        
        sb1, cf_best1 = np.loadtxt(bestpath + f"/stitched_16R_G2048_50_G512_2000_CG_B/BAOfit_"+filename + "_best_nw.dat", usecols=(0, 1), unpack=True)
        sigma1, alpha1, bias1, sigmanl1 = np.loadtxt(bestpath + f"/stitched_16R_G2048_50_G512_2000_CG_B/BAOfit_"+filename + "_mystats.txt", usecols=(1, 5, 6, 7), unpack=True)

        sb2, cf_best2 = np.loadtxt(bestpath + f"/avg_16R_2000_SK_60_150/BAOfit_"+filename + "_best_nw.dat", usecols=(0, 1), unpack=True)
        sigma2, alpha2, bias2, sigmanl2 = np.loadtxt(bestpath + f"/avg_16R_2000_SK_60_150/BAOfit_"+filename + "_mystats.txt", usecols=(1, 5, 6, 7), unpack=True)

        
        # if np.abs(2*(sigma - sigma1)/(sigma + sigma1)) > 0.1:
        if sigmanl < 9:# or sigmanl1 < 9 or sigmanl2 > 16:
            fig, ax = pt.subplots()
            print(filename)
            ax.set_title(filename)
            ax.errorbar(sd, sd * sd * cf_data, yerr=sd * sd * std, color="k", ls="", marker="o", label="data")
            ax.plot(sb, sb * sb * cf_best, color="red", ls="--", label=f"best CG_LC {alpha} {bias} {sigmanl}")
            
            ax.plot(sb1, sb1 * sb1 * cf_best1, color="blue", ls="--", label=f"best CG_B {alpha1} {bias1} {sigmanl1}")
            ax.plot(sb2, sb2 * sb2 * cf_best2, color="orange", ls="--", label=f"best SK_B {alpha2} {bias2} {sigmanl2}")
            ax.legend(bbox_to_anchor=[0.8, 0.74], loc='center', framealpha=0.)

            ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
            ax.set_ylabel(r"s$^2\xi$(s) [h$^{-2}$ Mpc$^2$]")

            ax.set_ylim([-25, 50])
            ax.set_xlim(left=35, right=175)
            ax.axvline(60, color="grey", ls="--")
            ax.axvline(150, color="grey", ls="--")

            pp.savefig(fig)
            pt.close()
        # pt.clf()

    pp.close()
    
def main():
    # outpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/"
    # datapath = "/scratch/variu/clustering/patchy_cmass_subset/lightcone_box1/real/vv2pcf/16R_1000cf/"
    # bestpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/"
    # myplot(datapath=datapath, bestpath=bestpath, type_ = "/stitched_16R_G2048_50_G512_2000_CG_LC/", outpath=outpath + "stitched_16R_G2048_50_G512_2000_CG_LC_best.pdf")
    # myplot(datapath=datapath, bestpath=bestpath, type_ = "/stitched_16R_G2048_50_G512_2000_CG_B/", outpath=outpath + "stitched_16R_G2048_50_G512_2000_CG_B_best.pdf")
    # myplot(datapath=datapath, bestpath=bestpath, type_ = "/parab_60_150_fixc/", outpath=outpath + "parab_60_150_fixc.pdf")

    # bestpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/"
    # outpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/"
    # datapath = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/"
    # myplot(datapath=datapath, bestpath=bestpath, type_ = "/stitched_16R_G2048_50_G512_2000_CG/", outpath=outpath + "stitched_16R_G2048_50_G512_2000_CG_B_best_B.pdf")


    bestpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/fixSigmaNL/"
    outpath = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/fixSigmaNL/"
    datapath = "/scratch/variu/clustering/patchy_cmass_subset/lightcone_box1/real/vv2pcf/16R_1000cf/"    
    myplot(datapath=datapath, bestpath=bestpath, type_ = "/stitched_16R_G2048_50_G512_2000_CG_LC/", outpath=outpath + "more.pdf")


if __name__== '__main__':
    main()
