import numpy as np
import matplotlib.pyplot as pt
import matplotlib as mpl


def ax_style(ax, title_):
    ax.set_xlabel(r"s$_\mathrm{max}$")
    ax.set_ylabel(r"s$_\mathrm{min}$" + ", " + title_)

    smin = np.array([50, 70])
    smax = np.array([130, 150, 170])

    ax.set_xticks([0., 2., 4.])
    ax.set_yticks([1., 3.])
    ax.set_xticklabels(smax)
    ax.set_yticklabels(smin)

def obtain_matrix(fig, ax, smin, smax, inpath="inpath", type_="type", index=None, title_="title", filestyle="BAOfit_avg_16R.2pcf_mystats.txt"):
    alpham = np.zeros((len(smin), len(smax)))
    sigma = np.zeros((len(smin), len(smax)))

    for j in range(len(smax)):
        for i in range(len(smin)):
            file_m = inpath + type_.format(smin[i], smax[j]) + "/" + filestyle
            sigma[i][j], alpham[i][j] = np.loadtxt(file_m, usecols=(1, 4), unpack=True)
      

    # h0 = ax[0].imshow(alpham-1, vmin=-0.005, vmax=0.005, origin="lower", cmap='seismic')
    # h0 = ax[0].imshow(alpham-1, vmin=-0.015, vmax=0.015, origin="lower", cmap='seismic')
    h1 = ax.imshow((alpham-1) / (sigma), vmin=-1, vmax=1, origin="lower", cmap='seismic')
    
    
    
    # ax_style(ax[0], smin, smax)
    ax_style(ax, title_)
   
    # ax[0].text(5.5, 1.5, title_, rotation='vertical')
    # ax.text(5.7, 1.5, title_, rotation='vertical')
    
    return sigma, alpham, h1

def main():
    pt.rcParams.update({'font.family': "serif"})
    pt.rcParams.update({'font.serif': "Times New Roman"})
    # pt.rcParams['mathtext.fontset'] = 'Times New Roman'
    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['mathtext.rm'] = 'serif'
    
    smin = np.array([40, 50, 60, 70, 80])
    smax = np.array([130, 140, 150, 160, 170, 180])

    path_g_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/range_test/vv2pcf/unscaled_covariance/"
    path_g_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/range_test/vhxcf/unscaled_covariance/"

    path_p_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/range_test/unscaled_covariance/"
    path_p_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vhxcf/range_test/unscaled_covariance/"
    
    path_cg_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/avg_fit/"
    path_cg_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/avg_fit/"

    path_sk_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/avg_fit/range_test/"
    path_sk_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/avg_fit/range_test/"
    
    # path_4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf/avg_fit/range_test/"
    # path_5 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_CG/avg_fit/"

    path_lc = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/range_test/"


    ###
    fig, ax = pt.subplots(4, 2, figsize=(6, 6), sharex=True, sharey="col", gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})

    _, _, h1 = obtain_matrix(None, ax[0][0], smin, smax, path_g_vv, type_="galaxy_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="A, GAL")
    obtain_matrix(None, ax[1][0], smin, smax, path_p_vv, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="A, PAR")
    obtain_matrix(None, ax[2][0], smin, smax, path_sk_vv, type_="G1024CIC_{}_{}_m_2000F", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="A, SK")
    obtain_matrix(None, ax[3][0], smin, smax, path_cg_vv, type_="stitched_16R_G2048_50_G512_2000_CG_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="A, CG")

    obtain_matrix(None, ax[0][1], smin, smax, path_g_vh, type_="galaxy_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="X, GAL")
    obtain_matrix(None, ax[1][1], smin, smax, path_p_vh, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="X, PAR")
    obtain_matrix(None, ax[2][1], smin, smax, path_sk_vh, type_="avg_2001_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="X, SK")
    obtain_matrix(None, ax[3][1], smin, smax, path_cg_vh, type_="stitched_16R_G2048_50_G512_2000_CG_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="X, CG")

    fig.colorbar(h1, ax=ax[:,1], orientation='vertical')
    fig.savefig("/home/astro/variu/void_project_figs/AX_B_CG_SK_GAL_PAR_fitting_interval.pdf", bbox_inches='tight')

    ###
    fig, ax = pt.subplots(4, 2, figsize=(6, 6), sharex=True, sharey="col", gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})

    _, _, h1 = obtain_matrix(None, ax[0][0], smin, smax, path_cg_vv, type_="stitched_16R_G2048_50_G512_2000_CG_imp_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_=r"A, CG$_\mathrm{def}$")
    obtain_matrix(None, ax[1][0], smin, smax, path_cg_vv, type_="stitched_16R_G2048_50_G512_2000_CG_80_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_=r"A, CG$_\mathrm{80}$")
    obtain_matrix(None, ax[2][0], smin, smax, path_cg_vv, type_="stitched_16R_G2048_50_G512_2000_CG_120_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_=r"A, CG$_\mathrm{120}$")
    obtain_matrix(None, ax[3][0], smin, smax, path_sk_vv, type_="G1024F0.3CIC_{}_{}_m_2000F", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_=r"A, SK$_\mathrm{def}$")

    obtain_matrix(None, ax[0][1], smin, smax, path_cg_vh, type_="stitched_16R_G2048_50_G512_2000_CG_imp_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_=r"X, CG$_\mathrm{def}$")
    obtain_matrix(None, ax[1][1], smin, smax, path_cg_vh, type_="stitched_16R_G2048_50_G512_2000_CG_80_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_=r"X, CG$_\mathrm{80}$")
    obtain_matrix(None, ax[2][1], smin, smax, path_cg_vh, type_="stitched_16R_G2048_50_G512_2000_CG_120_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_=r"X, CG$_\mathrm{120}$")
    obtain_matrix(None, ax[3][1], smin, smax, path_sk_vh, type_="G1024CICF0.3_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_=r"X, SK$_\mathrm{def}$")

    fig.colorbar(h1, ax=ax[:, 1], orientation='vertical')
    fig.savefig("/home/astro/variu/void_project_figs/AX_B_CGimp_GC80_GC120_SKimp_fitting_interval.pdf", bbox_inches='tight')

    ###
    fig, ax = pt.subplots(5, 1, figsize=(3, 4.5), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.00, "bottom":0.11})

    _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_lc, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="PAR")
    obtain_matrix(None, ax[1], smin, smax, path_lc, type_="avg_16R_2000_SK_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="SK")
    obtain_matrix(None, ax[2], smin, smax, path_lc, type_="stitched_16R_G2048_50_G512_2000_CG_LC_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="CG")
    
    obtain_matrix(None, ax[3], smin, smax, path_lc, type_="stitched_G2048_50_G512_2000_SK_B_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_=r"SK$_\mathrm{B}$")
    obtain_matrix(None, ax[4], smin, smax, path_lc, type_="stitched_16R_G2048_50_G512_2000_CG_B_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_=r"CG$_\mathrm{B}$")

    # fig.colorbar(h0, ax=ax[:,0], orientation='vertical')
    fig.colorbar(h1, ax=ax[:], orientation='vertical')
    fig.savefig("/home/astro/variu/void_project_figs/A_LC_CG_SK_CGB_SKB_par_fitting_interval.pdf", bbox_inches='tight')

    exit()

    ###
    # fig, ax = pt.subplots(4, 1, figsize=(3, 6), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})

    # _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_1, type_="galaxy_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="GAL")
    # obtain_matrix(None, ax[1], smin, smax, path_1, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="PAR")
    # obtain_matrix(None, ax[2], smin, smax, path_1, type_="G1024CIC_{}_{}_m_2000F", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="SK")
    # obtain_matrix(None, ax[3], smin, smax, path_0, type_="stitched_16R_G2048_50_G512_2000_CG_{}_{}", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="CG")

    # fig.colorbar(h1, ax=ax[:], orientation='vertical')
    # fig.savefig("/home/astro/variu/void_project_figs/A_B_CG_SK_GAL_PAR_avg_alpha.pdf")

    ####
    # fig, ax = pt.subplots(4, 1, figsize=(3, 6), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})

    # _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_0, type_="stitched_16R_G2048_50_G512_2000_CG_imp_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="CGimp")
    # obtain_matrix(None, ax[1], smin, smax, path_0, type_="stitched_16R_G2048_50_G512_2000_CG_80_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="CG80")
    # obtain_matrix(None, ax[2], smin, smax, path_0, type_="stitched_16R_G2048_50_G512_2000_CG_120_{}_{}", index=0, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="CG120")
    # obtain_matrix(None, ax[3], smin, smax, path_1, type_="G1024F0.3CIC_{}_{}_m_2000F", index=1, filestyle="BAOfit_avg_16R.2pcf_mystats.txt", title_="SKimp")

    # # fig.colorbar(h0, ax=ax[:,0], orientation='vertical')
    # fig.colorbar(h1, ax=ax[:], orientation='vertical')
    # fig.savefig("/home/astro/variu/void_project_figs/A_B_CGimp_GC80_GC120_SKimp_avg_alpha.pdf")

    
    ###
    # fig, ax = pt.subplots(4, 1, figsize=(3, 6), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})

    # _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_2, type_="galaxy_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="GAL")
    # obtain_matrix(None, ax[1], smin, smax, path_2, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="PAR")
    # obtain_matrix(None, ax[2], smin, smax, path_2, type_="avg_2001_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="SK")
    # obtain_matrix(None, ax[3], smin, smax, path_3, type_="stitched_16R_G2048_50_G512_2000_CG_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="CG")

    # # fig.colorbar(h0, ax=ax[:,0], orientation='vertical')
    # fig.colorbar(h1, ax=ax[:], orientation='vertical')
    # fig.savefig("/home/astro/variu/void_project_figs/X_B_CG_SK_GAL_PAR_avg_alpha.pdf")


    ###
    # fig, ax = pt.subplots(4, 1, figsize=(3, 6), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.1, "bottom":0.08, "right":0.85})
    # _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_3, type_="stitched_16R_G2048_50_G512_2000_CG_imp_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="CGimp")
    # obtain_matrix(None, ax[1], smin, smax, path_3, type_="stitched_16R_G2048_50_G512_2000_CG_80_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="CG80")
    # obtain_matrix(None, ax[2], smin, smax, path_3, type_="stitched_16R_G2048_50_G512_2000_CG_120_{}_{}", index=0, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="CG120")
    # obtain_matrix(None, ax[3], smin, smax, path_2, type_="G1024CICF0.3_{}_{}", index=1, filestyle="BAOfit_avg_16R_vhxcf.xcf_mystats.txt", title_="SKimp")

    # # fig.colorbar(h0, ax=ax[:,0], orientation='vertical')
    # fig.colorbar(h1, ax=ax[:], orientation='vertical')
    # fig.savefig("/home/astro/variu/void_project_figs/X_B_CGimp_GC80_GC120_SKimp_avg_alpha.pdf")

    ###
    fig, ax = pt.subplots(3, 1, figsize=(3, 4.5), sharex=True, sharey=True, gridspec_kw={"hspace":0, "top":0.96, "left":0.00, "bottom":0.11})

    # _, _, h0, h1 = obtain_matrix(None, [ax[0][0], ax[0][1]], smin, smax, path_2, type_="galaxy_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="GAL")
    _, _, h1 = obtain_matrix(None, ax[0], smin, smax, path_5, type_="parab_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="PAR")
    obtain_matrix(None, ax[1], smin, smax, path_4, type_="stitched_G2048_50_G512_2000_SK_LC_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="SK")
    obtain_matrix(None, ax[2], smin, smax, path_0, type_="stitched_16R_G2048_50_G512_2000_CG_LC_{}_{}", index=1, filestyle="BAOfit_avg_1000_16R.2pcf_mystats.txt", title_="CG")

    # fig.colorbar(h0, ax=ax[:,0], orientation='vertical')
    fig.colorbar(h1, ax=ax[:], orientation='vertical')
    fig.savefig("/home/astro/variu/void_project_figs/A_LC_CG_SK_GAL_PAR_avg_alpha.pdf")

    
    
    pt.show()


if __name__== '__main__':
      main()
