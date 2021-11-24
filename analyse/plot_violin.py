#import matplotlib as mpl
#mpl.use('Agg')
import os, re, sys
import glob
import numpy as np
import matplotlib.pyplot as pt
from getdist import plots, MCSamples, loadMCSamples
from matplotlib.backends.backend_pdf import PdfPages


def plot_data_mystats(inpath, endcffile=".dat.2pcf"):
    files = glob.glob(inpath + "*mystats*")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)

    print(len(files))
    alpha = np.zeros(len(files))
    alpham = np.zeros(len(files))
    sigma = np.zeros(len(files))
    for i, file_ in enumerate(files):
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + endcffile + "_mystats.txt"

        alpha[i], sigma[i], alpham[i]= np.loadtxt(file_m, usecols=(0, 1,4), unpack=True)
 

    return alpha, sigma, alpham

def plot_all_comb(ax, alpha1=None, alpha2=None, alpha3=None, label1="label1", label2="label2", label3="label3", mode=0):
    ax[0].plot(alpha3, alpha1, marker=".", ls="")
    ax[0].plot(np.linspace(np.min(alpha3), np.max(alpha3), 10), np.linspace(np.min(alpha3), np.max(alpha3), 10), color="grey", ls="--")
    ax[0].set_ylabel(label1, fontsize=16)
    ax[0].set_xlabel(label3, fontsize=16)
    
    ax[1].plot(alpha3, alpha2, marker=".", ls="")
    ax[1].plot(np.linspace(np.min(alpha3), np.max(alpha3), 10), np.linspace(np.min(alpha3), np.max(alpha3), 10), color="grey", ls="--")
    ax[1].set_ylabel(label2, fontsize=16)
    ax[1].set_xlabel(label3, fontsize=16)

    if mode == 0:
        # ax[1].set_xlim([0.79, 1.21])
        # ax[1].set_ylim([0.79, 1.21])
        # ax[0].set_ylim([0.79, 1.21])
        # ax[0].set_xlim([0.79, 1.21])

        ax[1].set_xlim([0.93, 1.07])
        ax[1].set_ylim([0.93, 1.07])
        ax[0].set_ylim([0.93, 1.07])
        ax[0].set_xlim([0.93, 1.07])

def plot_alphabest_alpha(ax, alphabest=None, alpha=None, labelbest="label1", label="label2"):
    ax.plot(alphabest, alpha, marker=".", ls="")
    ax.plot(np.linspace(0.8, 1.2, 10), np.linspace(0.8, 1.2, 10), color="grey", ls="--")
    ax.set_xlabel(labelbest, fontsize=16)
    ax.set_ylabel(label, fontsize=16)
    ax.set_ylim([0.93, 1.07])
    ax.set_xlim([0.93, 1.07])
    
def plot_three_cases(inpath1, inpath2, inpath3, dict_labels, endcffile=".dat.2pcf"):
    fig0, ax0 = pt.subplots(figsize=(4,3))
    fig1, ax1 = pt.subplots(figsize=(4,3))
    fig2, ax2 = pt.subplots(figsize=(4,3))

    # fig3, ax3 = pt.subplots(2, 2, figsize=(20,20), sharex=True, gridspec_kw={"hspace":0})
    # fig4, ax4 = pt.subplots(1, 2, figsize=(20,10), sharex=True, sharey=True)
    # fig5, ax5 = pt.subplots(1, 3, figsize=(30,10))
    # fig6, ax6 = pt.subplots(2, 2, figsize=(20,20))
    
    alpha_b1, chi2_b1, evi_1 = np.loadtxt(inpath1+"/DAT_alpha_chi2_500.dat",usecols=(0,1,2), unpack=True)
    alpha_a1, sigma_1, alpha_m1 = plot_data_mystats(inpath1, endcffile=endcffile)

    alpha_b2, chi2_b2, evi_2 = np.loadtxt(inpath2+"/DAT_alpha_chi2_500.dat",usecols=(0,1,2), unpack=True)
    alpha_a2, sigma_2, alpha_m2 = plot_data_mystats(inpath2, endcffile=endcffile)

    alpha_b3, chi2_b3, evi_3 = np.loadtxt(inpath3+"/DAT_alpha_chi2_500.dat",usecols=(0,1,2), unpack=True)
    alpha_a3, sigma_3, alpha_m3 = plot_data_mystats(inpath3, endcffile=endcffile)

    ax0.axhline(1, color="grey", ls=":")
    parts_a = ax0.violinplot([alpha_a1, alpha_a2, alpha_a3], [1, 2, 3], vert=True, showmeans=False, showmedians=False, showextrema=False)
    
    for b in parts_a['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further right than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
        b.set_color('r')

    parts_m = ax0.violinplot([alpha_m1, alpha_m2, alpha_m3], [1, 2, 3], vert=True, showmeans=True, showmedians=True, showextrema=True)
    
    for b in parts_m['bodies']:
        # get the center
        m = np.mean(b.get_paths()[0].vertices[:, 0])
        # modify the paths to not go further left than the center
        b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
        b.set_color('b')

    parts_m['cmeans'].set_color("b")
    parts_m['cmedians'].set_color("b")
    
    ax0.legend([parts_a['bodies'][0],parts_m['bodies'][0]],['avg', 'med'])

    ax0.set_xticks([1,2,3])
    ax0.set_xticklabels([dict_labels["label1"], dict_labels["label2"], dict_labels["label3"]])
    ax0.set_ylabel(r"$\alpha$")
    # plot_all_comb([ax3[0][0],ax3[1][0]], alpha1=alpha_b1, alpha2=alpha_b2, alpha3=alpha_b3, label1=dict_labels['alpha1best'], label2=dict_labels['alpha2best'], label3=dict_labels['alpha3best'])
    # plot_all_comb([ax3[0][1],ax3[1][1]], alpha1=alpha_a1, alpha2=alpha_a2, alpha3=alpha_a3, label1=dict_labels['alpha1'], label2=dict_labels['alpha2'], label3=dict_labels['alpha3'])
    # plot_all_comb(ax4, alpha1=sigma_1, alpha2=sigma_2, alpha3=sigma_3, label1=dict_labels['sigma1'], label2=dict_labels['sigma2'], label3=dict_labels['sigma3'], mode=1)

    # plot_alphabest_alpha(ax5[0], alphabest=alpha_b1, alpha=alpha_a1, labelbest=dict_labels['alpha1best'], label=dict_labels['alpha1'])
    # plot_alphabest_alpha(ax5[1], alphabest=alpha_b2, alpha=alpha_a2, labelbest=dict_labels['alpha2best'], label=dict_labels['alpha2'])
    # plot_alphabest_alpha(ax5[2], alphabest=alpha_b3, alpha=alpha_a3, labelbest=dict_labels['alpha3best'], label=dict_labels['alpha3'])


    fig0.tight_layout()
    fig0.savefig(dict_labels['filename']+"temp1.pdf")
    fig1.tight_layout()
    fig1.savefig(dict_labels['filename']+"temp2.pdf")
    fig2.tight_layout()
    fig2.savefig(dict_labels['filename']+"temp3.pdf")
    # fig3.tight_layout()
    # fig3.savefig(dict_labels['filename']+"temp4.png")
    # fig4.tight_layout()
    # fig4.savefig(dict_labels['filename']+"temp5.png")
    # fig5.tight_layout()
    # fig5.savefig(dict_labels['filename']+"temp6.png")
    # fig6.tight_layout()
    # fig6.savefig(dict_labels['filename']+"temp7.png")
    

def main():
    ### Test the effect of the template
    # inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_100_temp_exte/"
    # inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_500_temp_exte/"
    # inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_short/"
    # inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_exte_filt/"
    ###
    # plot_all_bestfit()
    # exit()

    inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/galaxy_60_150_m/"
    inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/parab_60_150_m_fast/"
    inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000/"

    inpath20 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000E/"
    inpath21 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000F/"
    inpath22 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_500/"
    inpath23 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_500F/"
    inpath24 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000_old/"
    
    inpath25 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/galaxy_60_150_m/"
    inpath26 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/parab_60_150_m_fast2F/"
    inpath27 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024CIC_60_150_m_2001F/"
    inpath37 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024F0.3CIC_60_150_m_2000F/"
    inpath38 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024CIC18R_60_150_m_2000F/"

    inpath30 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/stitched_subvolumes/"
    inpath31 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/lin_sav_gol_71_5_stitched_subvolumes/"

    inpath32 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/stitched_G2048_50_G512_2000/"
    inpath33 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/sm_stitched_G2048_50_G512_2000/"
    
    inpath34 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000FS/"
    inpath35 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024F0.3CIC_60_150_m_2000F/"
    inpath36 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC18R_60_150_m_2000F/"

    inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/galaxy_60_150_m/"
    inpath5 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/parab_60_150_m_fast/"
    inpath6 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/parab_60_150_m_fast2/"
    inpath7 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000_old/"
    inpath8 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000/"
    inpath9 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_B2000/"
    
    inpath10 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000D/"
    inpath11 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000E/"
    inpath12 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000F/"
    
    inpath13 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_500/"
    inpath14 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_500F/"

    #plot_all_bestfit()
    # obtain_chi2_alpha_evi(inpath35, endcffile=".VOID.dat.2pcf")
    # obtain_chi2_alpha_evi(inpath36, endcffile=".VOID.dat.2pcf")
    # obtain_chi2_alpha_evi(inpath37, endcffile=".dat.dr.2pcf")
    # obtain_chi2_alpha_evi(inpath38, endcffile=".dat.dr.2pcf")
    #obtain_chi2_alpha_evi(inpath34, endcffile=".VOID.dat.2pcf")
    #exit()

    dict_labels = {"label1":"F0.075 18R", "label3":"F0.075 16R", "label2":"FAC0.3 16R",
                    "alpha1best":r"$\alpha_\mathrm{best, 18R}$", "alpha3best":r"$\alpha_\mathrm{best, ST}$", "alpha2best":r"$\alpha_\mathrm{best, F03}$",
                    "alpha1":r"$\alpha_\mathrm{18R}$", "alpha3":r"$\alpha_\mathrm{ST}$", "alpha2":r"$\alpha_\mathrm{F03}$",
                    "sigma1":r"$\sigma_\mathrm{18R}$", "sigma3":r"$\sigma_\mathrm{ST}$", "sigma2":r"$\sigma_\mathrm{F03}$",
                    "filename":"/home/astro/variu/temp/X_B_18R_ST_F03_"}
    plot_three_cases(inpath38, inpath37, inpath37, dict_labels, endcffile=".dat.dr.2pcf")


    # dict_labels = {"label1":"G1024FAC0.075 18R", "label3":"G1024FAC0.075 st_2000_50 fast", "label2":"G1024FAC0.3 fast",
    #                 "alpha1best":r"$\alpha_\mathrm{best, 18R}$", "alpha3best":r"$\alpha_\mathrm{best, ST}$", "alpha2best":r"$\alpha_\mathrm{best, F03}$",
    #                 "alpha1":r"$\alpha_\mathrm{18R}$", "alpha3":r"$\alpha_\mathrm{ST}$", "alpha2":r"$\alpha_\mathrm{F03}$",
    #                 "sigma1":r"$\sigma_\mathrm{18R}$", "sigma3":r"$\sigma_\mathrm{ST}$", "sigma2":r"$\sigma_\mathrm{F03}$",
    #                 "filename":"/home/astro/variu/temp/B_18R_ST_F03_"}
    # plot_three_cases(inpath36, inpath35, inpath32, dict_labels, endcffile=".VOID.dat.2pcf")

    # dict_labels = {"label1":"Original fast", "label2":"stitch_2000_50 fast", "label3":"Shifted 0.0055 fast",
    #                 "alpha1best":r"$\alpha_\mathrm{best, OF}$", "alpha2best":r"$\alpha_\mathrm{best, ST}$", "alpha3best":r"$\alpha_\mathrm{best, SF}$",
    #                 "alpha1":r"$\alpha_\mathrm{OF}$", "alpha2":r"$\alpha_\mathrm{ST}$", "alpha3":r"$\alpha_\mathrm{SF}$",
    #                 "sigma1":r"$\sigma_\mathrm{OF}$", "sigma2":r"$\sigma_\mathrm{ST}$", "sigma3":r"$\sigma_\mathrm{SF}$",
    #                 "filename":"/home/astro/variu/temp/B_OF_ST_SF_"}
    # plot_three_cases(inpath21, inpath32, inpath34, dict_labels, endcffile=".VOID.dat.2pcf")


    # dict_labels = {"label1":"stitch subvolume", "label2":"sm stitch subvolume", "label3":"avg2000",
    #                 "alpha1best":r"$\alpha_\mathrm{best, ST}$", "alpha2best":r"$\alpha_\mathrm{best, STSM}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{ST}$", "alpha2":r"$\alpha_\mathrm{STSM}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{ST}$", "sigma2":r"$\sigma_\mathrm{STSM}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_STSUB_STSUBSM_TEM_"}
    # plot_three_cases(inpath30, inpath31, inpath21, dict_labels, endcffile=".VOID.dat.2pcf")

    
    # dict_labels = {"label1":"stitch_2000_50", "label2":"sm stitch_2000_50", "label3":"avg2000",
    #                 "alpha1best":r"$\alpha_\mathrm{best, ST}$", "alpha2best":r"$\alpha_\mathrm{best, STSM}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{ST}$", "alpha2":r"$\alpha_\mathrm{STSM}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{ST}$", "sigma2":r"$\sigma_\mathrm{STSM}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_ST200050_ST200050SM_TEM_"}
    # plot_three_cases(inpath32, inpath33, inpath21, dict_labels, endcffile=".VOID.dat.2pcf")


    # dict_labels = {"label1":"Galaxy: 12 DOF", "label2":"Parabola: 11 DOF", "label3":"Template: 12 DOF",
    #                 "alpha1best":r"$\alpha_\mathrm{best, GAL}$", "alpha2best":r"$\alpha_\mathrm{best, PAR}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{GAL}$", "alpha2":r"$\alpha_\mathrm{PAR}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{GAL}$", "sigma2":r"$\sigma_\mathrm{PAR}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/chengsplots/B_GAL_PAR_TEM_"}
    # plot_three_cases(inpath1, inpath2, inpath21, dict_labels, endcffile=".VOID.dat.2pcf")

    # dict_labels = {"label1":"Galaxy: 12 DOF", "label2":"Parabola: 11 DOF", "label3":"Template: 12 DOF",
    #                 "alpha1best":r"$\alpha_\mathrm{best, GAL}$", "alpha2best":r"$\alpha_\mathrm{best, PAR}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{GAL}$", "alpha2":r"$\alpha_\mathrm{PAR}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{GAL}$", "sigma2":r"$\sigma_\mathrm{PAR}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_GAL_PAR_TEM_"}
    # plot_three_cases(inpath1, inpath2, inpath3, dict_labels, endcffile=".VOID.dat.2pcf")


    # dict_labels = {"label1":"Galaxy: 12 DOF", "label2":"Parabola 2: 11 DOF", "label3":"Template: 12 DOF",
    #                 "alpha1best":r"$\alpha_\mathrm{best, GAL}$", "alpha2best":r"$\alpha_\mathrm{best, PAR2}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{GAL}$", "alpha2":r"$\alpha_\mathrm{PAR2}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{GAL}$", "sigma2":r"$\sigma_\mathrm{PAR2}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/LC_GAL_PAR2_TEM_"}
    # plot_three_cases(inpath4, inpath6, inpath8, dict_labels)

    # dict_labels = {"label1":"Template old: 12 DOF", "label2":"Template Box: 12 DOF", "label3":"Template new: 12 DOF",
    #                 "alpha1best":r"$\alpha_\mathrm{best, OLD}$", "alpha2best":r"$\alpha_\mathrm{best, BOX}$", "alpha3best":r"$\alpha_\mathrm{best, NEW}$",
    #                 "alpha1":r"$\alpha_\mathrm{OLD}$", "alpha2":r"$\alpha_\mathrm{BOX}$", "alpha3":r"$\alpha_\mathrm{NEW}$",
    #                 "sigma1":r"$\sigma_\mathrm{OLD}$", "sigma2":r"$\sigma_\mathrm{BOX}$", "sigma3":r"$\sigma_\mathrm{NEW}$",
    #                 "filename":"/home/astro/variu/temp/LC_OLD_BOX_NEW_"}
    # plot_three_cases(inpath7, inpath9, inpath8, dict_labels)


    # dict_labels = {"label1":"Galaxy: 12 DOF", "label2":"Parabola: 11 DOF", "label3":"Parabola 2: 11 DOF",
    #                 "alpha1best":r"$\alpha_\mathrm{best, GAL}$", "alpha2best":r"$\alpha_\mathrm{best, PAR}$", "alpha3best":r"$\alpha_\mathrm{best, PAR2}$",
    #                 "alpha1":r"$\alpha_\mathrm{GAL}$", "alpha2":r"$\alpha_\mathrm{PAR}$", "alpha3":r"$\alpha_\mathrm{PAR2}$",
    #                 "sigma1":r"$\sigma_\mathrm{GAL}$", "sigma2":r"$\sigma_\mathrm{PAR}$", "sigma3":r"$\sigma_\mathrm{PAR2}$",
    #                 "filename":"/home/astro/variu/temp/LC_GAL_PAR_PAR2_"}
    # plot_three_cases(inpath4, inpath5, inpath6, dict_labels)


    # dict_labels = {"label1":"Temp fast 1.5", "label2":"Temp fast 2.0", "label3":"Temp FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEMe}$", "alpha2best":r"$\alpha_\mathrm{best, TEMf}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEMe}$", "alpha2":r"$\alpha_\mathrm{TEMf}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEMe}$", "sigma2":r"$\sigma_\mathrm{TEMf}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/LC_TEMe_TEMf_TEM_"}
    # plot_three_cases(inpath11, inpath12, inpath8, dict_labels)
    

    # dict_labels = {"label1":"Temp fast 1.5", "label2":"Temp fast 1.0", "label3":"Temp FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEMe}$", "alpha2best":r"$\alpha_\mathrm{best, TEMd}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEMe}$", "alpha2":r"$\alpha_\mathrm{TEMd}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEMe}$", "sigma2":r"$\sigma_\mathrm{TEMd}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/LC_TEMe_TEMd_TEM_"}
    # plot_three_cases(inpath11, inpath10, inpath8, dict_labels)
    

    # dict_labels = {"label1":"Temp 500 FFTlinlog", "label2":"Temp 500 fast 2.0", "label3":"Temp 2000 FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEM500FFT}$", "alpha2best":r"$\alpha_\mathrm{best, TEM500F}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEM500FFT}$", "alpha2":r"$\alpha_\mathrm{TEM500F}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEM500FFT}$", "sigma2":r"$\sigma_\mathrm{TEM500F}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/LC_TEM500FFT_TEM500F_TEM_"}
    # plot_three_cases(inpath13, inpath14, inpath8, dict_labels)


    # dict_labels = {"label1":"Temp 500 FFTlinlog", "label2":"Temp 500 fast 2.0", "label3":"Temp 2000 fast 2.0",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEM500FFT}$", "alpha2best":r"$\alpha_\mathrm{best, TEM500F}$", "alpha3best":r"$\alpha_\mathrm{best, TEMF}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEM500FFT}$", "alpha2":r"$\alpha_\mathrm{TEM500F}$", "alpha3":r"$\alpha_\mathrm{TEMF}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEM500FFT}$", "sigma2":r"$\sigma_\mathrm{TEM500F}$", "sigma3":r"$\sigma_\mathrm{TEMF}$",
    #                 "filename":"/home/astro/variu/temp/LC_TEM500FFT_TEM500F_TEMF_"}
    # plot_three_cases(inpath13, inpath14, inpath12, dict_labels)
    

    # dict_labels = {"label1":"Temp 500 FFTlinlog", "label2":"Temp 500 fast 2.0", "label3":"Temp 2000 fast 2.0",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEM500FFT}$", "alpha2best":r"$\alpha_\mathrm{best, TEM500F}$", "alpha3best":r"$\alpha_\mathrm{best, TEMF}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEM500FFT}$", "alpha2":r"$\alpha_\mathrm{TEM500F}$", "alpha3":r"$\alpha_\mathrm{TEMF}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEM500FFT}$", "sigma2":r"$\sigma_\mathrm{TEM500F}$", "sigma3":r"$\sigma_\mathrm{TEMF}$",
    #                 "filename":"/home/astro/variu/temp/B_TEM500FFT_TEM500F_TEMF_"}
    # plot_three_cases(inpath22, inpath23, inpath21, dict_labels, endcffile=".VOID.dat.2pcf")


    # dict_labels = {"label1":"Temp 2000 fast 1.5", "label2":"Temp 2000 fast 2.0", "label3":"Temp 2000 FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEME}$", "alpha2best":r"$\alpha_\mathrm{best, TEMF}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEME}$", "alpha2":r"$\alpha_\mathrm{TEMF}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEME}$", "sigma2":r"$\sigma_\mathrm{TEMF}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_TEME_TEMF_TEM_"}
    # plot_three_cases(inpath20, inpath21, inpath3, dict_labels, endcffile=".VOID.dat.2pcf")

    

    # dict_labels = {"label1":"Temp 2000 fast 2.0", "label2":"Temp 2000 FFTlinlog old", "label3":"Temp 2000 FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEMF}$", "alpha2best":r"$\alpha_\mathrm{best, TEMold}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEMF}$", "alpha2":r"$\alpha_\mathrm{TEMold}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEMF}$", "sigma2":r"$\sigma_\mathrm{TEMold}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_TEMF_TEMold_TEM_"}
    # plot_three_cases(inpath21, inpath24, inpath3, dict_labels, endcffile=".VOID.dat.2pcf")


    # dict_labels = {"label1":"Temp 500 fast 2.0", "label2":"Temp 2000 FFTlinlog old", "label3":"Temp 2000 FFTlinlog",
    #                 "alpha1best":r"$\alpha_\mathrm{best, TEM500F}$", "alpha2best":r"$\alpha_\mathrm{best, TEMold}$", "alpha3best":r"$\alpha_\mathrm{best, TEM}$",
    #                 "alpha1":r"$\alpha_\mathrm{TEM500F}$", "alpha2":r"$\alpha_\mathrm{TEMold}$", "alpha3":r"$\alpha_\mathrm{TEM}$",
    #                 "sigma1":r"$\sigma_\mathrm{TEM500F}$", "sigma2":r"$\sigma_\mathrm{TEMold}$", "sigma3":r"$\sigma_\mathrm{TEM}$",
    #                 "filename":"/home/astro/variu/temp/B_TEM500F_TEMold_TEM_"}
    # plot_three_cases(inpath23, inpath24, inpath3, dict_labels, endcffile=".VOID.dat.2pcf")


    exit()

    
    #compute_covariance([inpath5, inpath, inpath6])
    #exit()
    
    fig, ax = pt.subplots(1, 2, figsize=(20,10))
    #plot_chi2_alpha(inpath, "", ax, color="red", label="G1024CIC")
    #plot_chi2_alpha(inpath5, "", ax, color="green", label="G1024IC_100")
    plot_chi2_alpha(inpath7, "", ax, color="red", label="G1024CIC_2001")
    
    #plot_chi2_alpha(inpath2, "_voidGL", ax, color="blue", label="G2048L512")
    #plot_chi2_alpha(inpath3, "_galaxy", ax, color="orange", label="galaxy")

    fig.savefig("/home/astro/variu/vhxcf_temporar4.pdf")

if __name__== '__main__':
    main()