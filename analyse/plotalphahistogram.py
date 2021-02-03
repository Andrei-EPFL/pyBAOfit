#import matplotlib as mpl
#mpl.use('Agg')
import os, re, sys
import glob
import numpy as np
import matplotlib.pyplot as pt
from getdist import plots, MCSamples, loadMCSamples
from matplotlib.backends.backend_pdf import PdfPages

def plot_errorbar(path, leglabel, ax, col, line, mode="cf", bias=1, plotBOOL=True):
    '''
    Plots the powerspectrum with error bar for all 100 halo/void files!
    '''
    files = glob.glob(path)
    print("For: {}; There are {}".format(leglabel, len(files)))
    x_aux, y_aux = np.loadtxt(files[0], comments="#", usecols=(0,1), unpack=True)
    y_arr = np.zeros(( len(files), len(y_aux) ))
    
    for i, filename in enumerate(files):
        try:
            x, y_tmp = np.loadtxt(filename, comments="#", usecols=(0,1), unpack=True)
            #ax.plot(k, pk_tmp, color='grey')
            y_arr[i] = y_tmp
        except:
            continue
    y_avg = np.mean(y_arr, 0)
    y_std = np.std(y_arr, 0)
    
    if(mode=="pk"):
        if plotBOOL:
            ax.errorbar(x, y_avg, yerr=y_std, label=leglabel, linestyle=line, color=col)
            ax.set_yscale('log')
            ax.set_xscale('log')
    elif(mode=="cf"):
        if plotBOOL:
            ax.errorbar(x, x*x*y_avg, yerr=x*x*y_std, label=leglabel, linestyle=line, color=col)
    else:
        print("Error: you have to choose between cf or pk as modes")
    return x, y_avg, y_std

def plot_chi2_alpha(inpath, endname, ax, color="red", label="label", plot=True):
    alpha, chi2 = np.loadtxt(inpath+"/DAT_alpha_chi2_500"+endname+".dat",usecols=(0,1), unpack=True)
    print(np.min(alpha), np.max(alpha))
    
    range_ = alpha>0.1
    alpha = alpha[range_]
    chi2 = chi2[range_]
    bins = np.linspace(0.75, 1.25, 51)
    print(bins)
    if plot == True:
        ax[0].plot(alpha, chi2, color=color, marker='.', ls='', label=label+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)))
        ax[0].set_xlabel(r"$\alpha_\mathrm{best}$")
        ax[0].set_ylabel(r"$\chi^2_\mathrm{best}$")
        ax[0].legend()
        
        #bins = np.linspace(0.8, 1.2, 161)
        ax[1].hist(alpha, bins=bins, alpha=1, histtype="step", color=color, label=label)#+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)))#, fill=False, edgecolor=color)
        ax[1].axvline(np.mean(alpha), color=color, ls=":")
        ax[1].axvline(1, color="grey", ls="--")
        
        #ax[1].set_xlim([0.96, 1.04])
        ax[1].set_xlabel(r"$\alpha_\mathrm{best}$")
        ax[1].set_ylabel(r"counts")
        ax[1].legend(loc="upper left")

    return alpha, chi2

def plot_chi2_alpha_evi(inpath, endname, ax, color="red", label="label", plot=True, n=1):
    alpha, chi2, evi = np.loadtxt(inpath+"/DAT_alpha_chi2_500"+endname+".dat",usecols=(0,1,2), unpack=True)
    chi2 = chi2/n

    range_ = alpha>0.1
    chi2 = chi2[range_]
    evi = evi[range_]
    alpha = alpha[range_]

    if plot == True:
        ax[0].plot(alpha, chi2, color=color, marker='.', ls='', label=label)
        ax[0].axvline(np.mean(alpha), color=color, ls="--")
        ax[0].set_xlabel(r"$\alpha_\mathrm{best}$")
        ax[0].set_ylabel(r"$\chi^2_\mathrm{best}/\nu$")
        ax[0].legend()
    
        ax[1].plot(evi, chi2, color=color, marker='.', ls='', label=label)
        ax[1].axvline(np.mean(evi), color=color, ls="--")
        ax[1].axhline(np.mean(chi2), color=color, ls="--")
        ax[1].set_xlabel(r"ln(ev)")
        ax[1].set_ylabel(r"$\chi^2_\mathrm{best}/\nu$")
        ax[1].legend()

    return alpha, chi2

def plot_data_mystats(inpath, ax, label="label", color="red"):
    files = glob.glob(inpath + "*mystats*")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)

    print(len(files))
    alpha = np.zeros(len(files))
    sigma = np.zeros(len(files))
    for i, file_ in enumerate(files):
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".dat.2pcf_mystats.txt"

        alpha[i], sigma[i], _ = np.loadtxt(file_m, usecols=(0,1,2), unpack=True)
        #if(sigma[i]>0.015):
        #    print(file_)
        #    plot_best(file_, inpath, label)
    
    ax[0].plot(alpha, sigma, ls="", marker=".", label=label, color=color)
    #ax[0].set_ylim([0.0045, 0.013])
    ax[0].set_xlabel(r"$\alpha$")
    ax[0].set_ylabel(r"$\sigma_{\alpha}$")
    ax[0].legend()
    
    print(np.min(sigma), np.max(sigma))
    bins = np.linspace(0.007, 0.197, 96)
    print(bins)
    ax[1].hist(sigma, alpha=1, histtype="step", bins=bins, label=label+": mean= {}; std={}".format(np.mean(sigma), np.std(sigma)), color=color)#, bins=bins)
    ax[1].set_xlabel(r"$\sigma_{\alpha}$")
    ax[1].set_ylabel("counts")
    ax[1].legend()

    bins = np.linspace(0.75, 1.25, 51)
    ax[2].hist(alpha, alpha=1,histtype="step", bins=bins, label=label+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)), color=color)
    ax[2].set_xlabel(r"$\alpha$")
    ax[2].set_ylabel("counts")
    ax[2].legend()

    return alpha, sigma

def obtain_chi2_alpha_evi(inpath):
    files = glob.glob(inpath + "*_.txt")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)

    print(len(files))
    chi2_arr = np.zeros(len(files))
    alpha_arr = np.zeros(len(files))
    evi_arr = np.zeros(len(files))
    
    for i, file_ in enumerate(files):
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".dat.2pcf_.txt"

        if os.path.isfile(file_m):
            chi2_tmp, alpha_tmp = np.loadtxt(file_m, usecols=(1, 2), unpack=True)
            
            file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".dat.2pcf_mystats.txt"
            if os.path.isfile(file_m):
            
                evi_tmp = np.loadtxt(file_m, usecols=(2), unpack=True)
                
                min_pos = chi2_tmp==np.min(chi2_tmp)
                
                chi2_arr[i] = chi2_tmp[min_pos]
                alpha_arr[i] = alpha_tmp[min_pos]
                evi_arr[i] = evi_tmp
            else:
                print(file_m)
        else:
            print(file_m)
    np.savetxt(inpath + "/DAT_alpha_chi2_500.dat", np.array([alpha_arr, chi2_arr, evi_arr]).T)

def provide_bestfit_alpha(path_fit, filename_cf):
    best_file = path_fit + "/BAOfit_" + filename_cf+"_best.dat"
    chi2_file = path_fit + "/BAOfit_" + filename_cf + "_.txt"
    mystats_f = path_fit + "/BAOfit_" + filename_cf + "_mystats.txt"

    sb, cfb = np.loadtxt(best_file, usecols=(0, 1), unpack=True)

    chi2_arr, alpha_arr = np.loadtxt(chi2_file, usecols=(1, 2), unpack=True)
    min_pos = chi2_arr==np.min(chi2_arr)

    alpha, sigma, evi = np.loadtxt(mystats_f, usecols=(0,1,2), unpack=True)

    return sb, cfb, chi2_arr[min_pos], alpha_arr[min_pos], alpha, sigma, evi

def plot_1file_chi2alpha(infile, npar):
    chi2, alpha, B, Snl = np.loadtxt(infile, usecols=(1, 2, 3, 4), unpack=True)
    
    alphalim1 = np.logical_and(alpha<1.08, alpha>1.05)
    alphalim2 = np.logical_and(alpha<1.05, alpha>1.027)
    alphalim3 = np.logical_and(alpha<1.027, alpha>0.99)

    pos1 = chi2 == np.min(chi2[alphalim1])
    pos2 = chi2 == np.min(chi2[alphalim2])
    pos3 = chi2 == np.min(chi2[alphalim3])

    print(np.min(chi2))
    print(np.min(chi2[alphalim1]), alpha[pos1])
    print(np.min(chi2[alphalim2]), alpha[pos2])
    print(np.min(chi2[alphalim3]), alpha[pos3])
    
    pt.plot(alpha[alphalim1], np.exp(-chi2[alphalim1]), ls="", marker=".")
    pt.axvline(alpha[pos1], color="grey")
    pt.plot(alpha[alphalim2], np.exp(-chi2[alphalim2]), ls="", marker=".")
    pt.axvline(alpha[pos2], color="grey")
    pt.plot(alpha[alphalim3], np.exp(-chi2[alphalim3]), ls="", marker=".")
    pt.axvline(alpha[pos3], color="grey")
    #pt.ylim([0.4, 3.6])
    pt.xlim([0.98, 1.11])
    pt.savefig("./temporaryalpha.pdf")
       
def plot_best(filepath_cf, filename_cf, pdf, std):
    path_fit = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/"    
    s_data, cf_data = np.loadtxt(filepath_cf + filename_cf, usecols=(0,1), unpack=True)

    sbG, cfbG, chi2G, alpha_bG, alpha_G, sigma_G, evi_G = provide_bestfit_alpha(path_fit + "/galaxy_60_150_m/", filename_cf)
    sbP, cfbP, chi2P, alpha_bP, alpha_P, sigma_P, evi_P = provide_bestfit_alpha(path_fit + "/parab_60_150_m_fast/", filename_cf)
    sbT, cfbT, chi2T, alpha_bT, alpha_T, sigma_T, evi_T = provide_bestfit_alpha(path_fit + "/G1024CIC_60_150_m_2000/", filename_cf)
    if alpha_bT > 1.05:
        fig, ax = pt.subplots(figsize=(10, 10))
        
        ax.errorbar(s_data, cf_data*s_data**2, yerr=s_data*s_data*std, label="data: "+ filename_cf, color="k")
        
        ax.plot(sbG, sbG*sbG*cfbG, label="Galaxy: " + r"$\chi2$=%.4g; $\alpha_\mathrm{best}$=%.4g; $\alpha$=%.4g; $\sigma$=%.4g; ln(ev)=%.4g" % (chi2G, alpha_bG, alpha_G, sigma_G, evi_G), ls="--", color="magenta")
        ax.plot(sbP, sbP*sbP*cfbP, label="Parabola: " + r"$\chi2$=%.4g; $\alpha_\mathrm{best}$=%.4g; $\alpha$=%.4g; $\sigma$=%.4g; ln(ev)=%.4g" % (chi2P, alpha_bP, alpha_P, sigma_P, evi_P), ls="--", color="green")
        ax.plot(sbT, sbT*sbT*cfbT, label="Template: " + r"$\chi2$=%.4g; $\alpha_\mathrm{best}$=%.4g; $\alpha$=%.4g; $\sigma$=%.4g; ln(ev)=%.4g" % (chi2T, alpha_bT, alpha_T, sigma_T, evi_T), ls="--", color="blue")
        
        ax.legend()

        range_ = np.where(np.logical_and(s_data <=165, s_data >=45))
        ymax = 1.1 * np.max(cf_data[range_] * s_data[range_]**2)
        ymin = 1.1 *np.min(cf_data[range_] * s_data[range_]**2)
        
        ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
        ax.set_ylabel(r"s$^2\xi$(s)")

        ax.set_ylim([ymin,ymax])
        ax.set_xlim(left=45, right=165)
        
        ax.axvline(x=60, color='grey', linestyle='--')
        ax.axvline(x=150, color='grey', linestyle='--')
        
        pdf.savefig(fig)
        pt.close(fig)

def compute_covariance(inpath_list):
    alpha_list = []
    sigma_list = []
    for inpath in inpath_list:
        files = glob.glob(inpath + "*mystats*")
        print(len(files))
        alpha = np.zeros(len(files))
        sigma = np.zeros(len(files))
        for i, file_ in enumerate(files):
            alpha[i], sigma[i], _ = np.loadtxt(file_, usecols=(0,1,2), unpack=True)
        
        alpha_list.append(alpha)
        sigma_list.append(sigma)
    
    alpha_cov = np.corrcoef(np.array(alpha_list))
    print(alpha_cov.shape)
    pt.imshow(alpha_cov)
    pt.colorbar()
    pt.savefig("./output/alpha_cov.pdf")
 
def plot_tests_effect_template(inpath1, inpath2, inpath3, inpath4):
    ### Fig 1
    fig, ax = pt.subplots(2, 1, figsize=(30, 10))
    alpha_100, chi2_100 = plot_chi2_alpha(inpath1, "", ax, color="red", label="100")
    alpha_500, chi2_500 = plot_chi2_alpha(inpath2, "", ax, color="blue", label="500")
    alpha_2000, chi2_2000 = plot_chi2_alpha(inpath3, "", ax, color="orange", label="2000")
    alpha_2000f, chi2_2000f = plot_chi2_alpha(inpath4, "", ax, color="green", label="2000filt")
    fig.savefig("/home/astro/variu/temp.pdf")
    
    ### Fig 2
    fig, ax = pt.subplots(2, 2, figsize=(30, 10))
    ax[0][0].plot(alpha_2000, alpha_100, color="red", marker=".", ls='', label=np.corrcoef(alpha_2000, alpha_100)[0, 1])
    ax[0][0].plot(alpha_2000, alpha_500, color="blue", marker=".", ls='', label=np.corrcoef(alpha_2000, alpha_500)[0, 1])
    ax[0][0].plot(alpha_2000, alpha_2000, color="orange", marker=".", ls='', label=np.corrcoef(alpha_2000, alpha_2000)[0, 1])
    ax[0][0].plot(alpha_2000, alpha_2000f, color="green", marker=".", ls='', label=np.corrcoef(alpha_2000, alpha_2000f)[0, 1])
    
    ax[1][0].plot(alpha_2000, 100*(alpha_2000-alpha_100)/alpha_2000, color="red", marker=".", ls='')
    ax[1][0].plot(alpha_2000, 100*(alpha_2000-alpha_500)/alpha_2000, color="blue", marker=".", ls='')
    ax[1][0].plot(alpha_2000, 100*(alpha_2000-alpha_2000)/alpha_2000, color="orange", marker=".", ls='')
    ax[1][0].plot(alpha_2000, 100*(alpha_2000-alpha_2000f)/alpha_2000, color="green", marker=".", ls='')
    ax[1][0].set_xlabel("alpha_2000")
    ax[1][0].set_ylabel("%")

    ax[0][1].plot(chi2_2000, chi2_100, color="red", marker=".", ls='', label=np.corrcoef(chi2_2000, chi2_100)[0, 1])
    ax[0][1].plot(chi2_2000, chi2_500, color="blue", marker=".", ls='', label=np.corrcoef(chi2_2000, chi2_500)[0, 1])
    ax[0][1].plot(chi2_2000, chi2_2000, color="orange", marker=".", ls='', label=np.corrcoef(chi2_2000, chi2_2000)[0, 1])
    ax[0][1].plot(chi2_2000, chi2_2000f, color="green", marker=".", ls='', label=np.corrcoef(chi2_2000, chi2_2000f)[0, 1])
    
    ax[1][1].plot(chi2_2000, 100*(chi2_2000-chi2_100)/chi2_2000, color="red", marker=".", ls='')
    ax[1][1].plot(chi2_2000, 100*(chi2_2000-chi2_500)/chi2_2000, color="blue", marker=".", ls='')
    ax[1][1].plot(chi2_2000, 100*(chi2_2000-chi2_2000)/chi2_2000, color="orange", marker=".", ls='')
    ax[1][1].plot(chi2_2000, 100*(chi2_2000-chi2_2000f)/chi2_2000, color="green", marker=".", ls='')
    ax[1][1].set_xlabel("chi2_2000")
    ax[1][1].set_ylabel("%")
    
    ax[0][0].legend()
    ax[0][1].legend()
    fig.savefig("/home/astro/variu/temp_2.pdf")
    
    ### Fig 3
    fig, ax = pt.subplots(1, 3, figsize=(30,10))
    alpham_100, sigma_100 = plot_data_mystats(inpath1, ax, label="G1024CIC_100", color="red")
    alpham_500, sigma_500 = plot_data_mystats(inpath2, ax, label="G1024CIC_500", color="blue")
    alpham_2000, sigma_2000 = plot_data_mystats(inpath3, ax, label="G1024CIC_2000", color="orange")
    alpham_2000f, sigma_2000f = plot_data_mystats(inpath4, ax, label="G1024CIC_2000filt", color="green")
    fig.savefig("/home/astro/variu/temp_3.pdf")

    ### Fig 4
    fig, ax = pt.subplots(2, 2, figsize=(30,10), gridspec_kw={"hspace":0}, sharex="col")
    ax[0][0].plot(alpham_2000, alpham_100, color="red", marker=".", ls='', label=np.corrcoef(alpham_2000, alpham_100)[0, 1])
    ax[0][0].plot(alpham_2000, alpham_500, color="blue", marker=".", ls='', label=np.corrcoef(alpham_2000, alpham_500)[0, 1])
    ax[0][0].plot(alpham_2000, alpham_2000, color="orange", marker=".", ls='', label=np.corrcoef(alpham_2000, alpham_2000)[0, 1])
    ax[0][0].plot(alpham_2000, alpham_2000f, color="green", marker=".", ls='', label=np.corrcoef(alpham_2000, alpham_2000f)[0, 1])
    ax[0][0].set_xlabel("alpha_2000")
    
    ax[1][0].plot(alpham_2000, 100*(alpham_2000-alpham_100)/alpham_2000, color="red", marker=".", ls='')
    ax[1][0].plot(alpham_2000, 100*(alpham_2000-alpham_500)/alpham_2000, color="blue", marker=".", ls='')
    ax[1][0].plot(alpham_2000, 100*(alpham_2000-alpham_2000)/alpham_2000, color="orange", marker=".", ls='')
    ax[1][0].plot(alpham_2000, 100*(alpham_2000-alpham_2000f)/alpham_2000, color="green", marker=".", ls='')
    ax[1][0].set_xlabel("alpha_2000")
    ax[1][0].set_ylabel("%")

    ax[0][1].plot(sigma_2000, sigma_100, color="red", marker=".", ls='', label=np.corrcoef(sigma_2000, sigma_100)[0, 1])
    ax[0][1].plot(sigma_2000, sigma_500, color="blue", marker=".", ls='', label=np.corrcoef(sigma_2000, sigma_500)[0, 1])
    ax[0][1].plot(sigma_2000, sigma_2000, color="orange", marker=".", ls='', label=np.corrcoef(sigma_2000, sigma_2000)[0, 1])
    ax[0][1].plot(sigma_2000, sigma_2000f, color="green", marker=".", ls='', label=np.corrcoef(sigma_2000, sigma_2000f)[0, 1])
    
    ax[1][1].plot(sigma_2000, 100*(sigma_2000-sigma_100)/sigma_2000, color="red", marker=".", ls='')
    ax[1][1].plot(sigma_2000, 100*(sigma_2000-sigma_500)/sigma_2000, color="blue", marker=".", ls='')
    ax[1][1].plot(sigma_2000, 100*(sigma_2000-sigma_2000)/sigma_2000, color="orange", marker=".", ls='')
    ax[1][1].plot(sigma_2000, 100*(sigma_2000-sigma_2000f)/sigma_2000, color="green", marker=".", ls='')
    ax[1][1].set_xlabel("sigma_2000")
    ax[1][1].set_ylabel("%")

    fig.savefig("/home/astro/variu/temp_4.pdf")

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
        ax[1].set_xlim([0.79, 1.21])
        ax[1].set_ylim([0.79, 1.21])
        ax[0].set_ylim([0.79, 1.21])
        ax[0].set_xlim([0.79, 1.21])

def plot_alphabest_alpha(ax, alphabest=None, alpha=None, labelbest="label1", label="label2"):
    ax.plot(alphabest, alpha, marker=".", ls="")
    ax.plot(np.linspace(0.8, 1.2, 10), np.linspace(0.8, 1.2, 10), color="grey", ls="--")
    ax.set_xlabel(labelbest, fontsize=16)
    ax.set_ylabel(label, fontsize=16)
    
def plot_three_cases(inpath1, inpath2, inpath3):
    ### Fig 1
    fig0, ax0 = pt.subplots()
    fig1, ax1 = pt.subplots()
    fig2, ax2 = pt.subplots(1, 3, figsize=(30,10))
    fig3, ax3 = pt.subplots(2, 2, figsize=(20,20), sharex=True, gridspec_kw={"hspace":0})
    fig4, ax4 = pt.subplots(1, 2, figsize=(20,10), sharex=True, sharey=True)
    fig5, ax5 = pt.subplots(1, 3, figsize=(30,10))
    fig6, ax6 = pt.subplots()
    fig7, ax7 = pt.subplots()
    #alpha_500, chi2_500 = plot_chi2_alpha_evi(inpath1, "", [ax0, ax1], color="magenta", label="Galaxy: 12 DOF", n=12)
    #alpha_2000, chi2_2000 = plot_chi2_alpha_evi(inpath2, "", [ax0, ax1], color="green", label="Parabola: 11 DOF", n=11)
    #alpha_100, chi2_100 = plot_chi2_alpha_evi(inpath3, "", [ax0, ax1], color="blue", label="Template: 12 DOF", n=12)

    alpha_bG, _ = plot_chi2_alpha(inpath1, "", [ax0, ax1], color="magenta", label="Galaxy: 12 DOF")
    plot_chi2_alpha_evi(inpath1, "", [ax6, ax7], color="magenta", label="Galaxy: 12 DOF")
    alpha_aG, sigma_G = plot_data_mystats(inpath1, ax2, label="Galaxy", color="magenta")

    alpha_bP, _ = plot_chi2_alpha(inpath2, "", [ax0, ax1], color="green", label="Parabola: 11 DOF")
    plot_chi2_alpha_evi(inpath2, "", [ax6, ax7], color="green", label="Parabola: 11 DOF")
    alpha_aP, sigma_P = plot_data_mystats(inpath2, ax2, label="Parabola", color="green")
    
    alpha_bT, _ = plot_chi2_alpha(inpath3, "", [ax0, ax1], color="blue", label="Template: 12 DOF")
    plot_chi2_alpha_evi(inpath3, "", [ax6, ax7], color="blue", label="Template: 12 DOF")
    alpha_aT, sigma_T = plot_data_mystats(inpath3, ax2, label="Template", color="blue")

    plot_all_comb([ax3[0][0],ax3[1][0]], alpha1=alpha_bG, alpha2=alpha_bP, alpha3=alpha_bT, label1=r"$\alpha_\mathrm{best, GAL}$", label2=r"$\alpha_\mathrm{best, PAR}$", label3=r"$\alpha_\mathrm{best, TEM}$")
    plot_all_comb([ax3[0][1],ax3[1][1]], alpha1=alpha_aG, alpha2=alpha_aP, alpha3=alpha_aT, label1=r"$\alpha_\mathrm{GAL}$", label2=r"$\alpha_\mathrm{PAR}$", label3=r"$\alpha_\mathrm{TEM}$")
    plot_all_comb(ax4, alpha1=sigma_G, alpha2=sigma_P, alpha3=sigma_T, label1=r"$\sigma_\mathrm{GAL}$", label2=r"$\sigma_\mathrm{PAR}$", label3=r"$\sigma_\mathrm{TEM}$", mode=1)

    plot_alphabest_alpha(ax5[0], alphabest=alpha_bG, alpha=alpha_aG, labelbest=r"$\alpha_\mathrm{best, GAL}$", label=r"$\alpha_\mathrm{GAL}$")
    plot_alphabest_alpha(ax5[1], alphabest=alpha_bP, alpha=alpha_aP, labelbest=r"$\alpha_\mathrm{best, PAR}$", label=r"$\alpha_\mathrm{PAR}$")
    plot_alphabest_alpha(ax5[2], alphabest=alpha_bT, alpha=alpha_aT, labelbest=r"$\alpha_\mathrm{best, TEM}$", label=r"$\alpha_\mathrm{TEM}$")

    fig0.tight_layout()
    fig0.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp4.png")
    fig1.tight_layout()
    fig1.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp5.png")
    fig2.tight_layout()
    fig2.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp6.png")
    fig3.tight_layout()
    fig3.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp7.png")
    fig4.tight_layout()
    fig4.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp8.png")
    fig5.tight_layout()
    fig5.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp9.png")
    fig6.tight_layout()
    fig6.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp10.png")
    fig7.tight_layout()
    fig7.savefig("/home/astro/variu/GAL_PAR2_TEMn_temp11.png")

def plot_three_cases_2(inpath1, inpath2, inpath3):
    ### Fig 1
    fig0, ax0 = pt.subplots()
    fig1, ax1 = pt.subplots()
    fig2, ax2 = pt.subplots(1, 3, figsize=(30,10))
    fig3, ax3 = pt.subplots(2, 2, figsize=(20,20), sharex=True, gridspec_kw={"hspace":0})
    fig4, ax4 = pt.subplots(1, 2, figsize=(20,10), sharex=True, sharey=True)
    fig5, ax5 = pt.subplots(1, 3, figsize=(30,10))
    fig6, ax6 = pt.subplots()
    fig7, ax7 = pt.subplots()
    
    alpha_bG, _ = plot_chi2_alpha(inpath1, "", [ax0, ax1], color="magenta", label="Template old: 12 DOF")
    plot_chi2_alpha_evi(inpath1, "", [ax6, ax7], color="magenta", label="Template old: 12 DOF")
    alpha_aG, sigma_G = plot_data_mystats(inpath1, ax2, label="Template old", color="magenta")

    alpha_bP, _ = plot_chi2_alpha(inpath2, "", [ax0, ax1], color="green", label="Template Box: 12 DOF")
    plot_chi2_alpha_evi(inpath2, "", [ax6, ax7], color="green", label="Template Box: 12 DOF")
    alpha_aP, sigma_P = plot_data_mystats(inpath2, ax2, label="Template Box", color="green")
    
    alpha_bT, _ = plot_chi2_alpha(inpath3, "", [ax0, ax1], color="blue", label="Template new: 12 DOF")
    plot_chi2_alpha_evi(inpath3, "", [ax6, ax7], color="blue", label="Template new: 12 DOF")
    alpha_aT, sigma_T = plot_data_mystats(inpath3, ax2, label="Template new", color="blue")

    plot_all_comb([ax3[0][0],ax3[1][0]], alpha1=alpha_bG, alpha2=alpha_bP, alpha3=alpha_bT, label1=r"$\alpha_\mathrm{best, OLD}$", label2=r"$\alpha_\mathrm{best, BOX}$", label3=r"$\alpha_\mathrm{best, NEW}$")
    plot_all_comb([ax3[0][1],ax3[1][1]], alpha1=alpha_aG, alpha2=alpha_aP, alpha3=alpha_aT, label1=r"$\alpha_\mathrm{OLD}$", label2=r"$\alpha_\mathrm{BOX}$", label3=r"$\alpha_\mathrm{NEW}$")
    plot_all_comb(ax4, alpha1=sigma_G, alpha2=sigma_P, alpha3=sigma_T, label1=r"$\sigma_\mathrm{OLD}$", label2=r"$\sigma_\mathrm{BOX}$", label3=r"$\sigma_\mathrm{NEW}$", mode=1)

    plot_alphabest_alpha(ax5[0], alphabest=alpha_bG, alpha=alpha_aG, labelbest=r"$\alpha_\mathrm{best, OLD}$", label=r"$\alpha_\mathrm{OLD}$")
    plot_alphabest_alpha(ax5[1], alphabest=alpha_bP, alpha=alpha_aP, labelbest=r"$\alpha_\mathrm{best, BOX}$", label=r"$\alpha_\mathrm{BOX}$")
    plot_alphabest_alpha(ax5[2], alphabest=alpha_bT, alpha=alpha_aT, labelbest=r"$\alpha_\mathrm{best, NEW}$", label=r"$\alpha_\mathrm{NEW}$")

    fig0.tight_layout()
    fig0.savefig("/home/astro/variu/OLD_BOX_NEW_temp4.png")
    fig1.tight_layout()
    fig1.savefig("/home/astro/variu/OLD_BOX_NEW_temp5.png")
    fig2.tight_layout()
    fig2.savefig("/home/astro/variu/OLD_BOX_NEW_temp6.png")
    fig3.tight_layout()
    fig3.savefig("/home/astro/variu/OLD_BOX_NEW_temp7.png")
    fig4.tight_layout()
    fig4.savefig("/home/astro/variu/OLD_BOX_NEW_temp8.png")
    fig5.tight_layout()
    fig5.savefig("/home/astro/variu/OLD_BOX_NEW_temp9.png")
    fig6.tight_layout()
    fig6.savefig("/home/astro/variu/OLD_BOX_NEW_temp10.png")
    fig7.tight_layout()
    fig7.savefig("/home/astro/variu/OLD_BOX_NEW_temp11.png")

def plot_three_cases_3(inpath1, inpath2, inpath3):
    ### Fig 1
    fig0, ax0 = pt.subplots()
    fig1, ax1 = pt.subplots()
    fig2, ax2 = pt.subplots(1, 3, figsize=(30,10))
    fig3, ax3 = pt.subplots(2, 2, figsize=(20,20), sharex=True, gridspec_kw={"hspace":0})
    fig4, ax4 = pt.subplots(1, 2, figsize=(20,10), sharex=True, sharey=True)
    fig5, ax5 = pt.subplots(1, 3, figsize=(30,10))
    fig6, ax6 = pt.subplots()
    fig7, ax7 = pt.subplots()
    
    alpha_bG, _ = plot_chi2_alpha(inpath1, "", [ax0, ax1], color="magenta", label="Galaxy: 12 DOF")
    plot_chi2_alpha_evi(inpath1, "", [ax6, ax7], color="magenta", label="Galaxy: 12 DOF")
    alpha_aG, sigma_G = plot_data_mystats(inpath1, ax2, label="Galaxy", color="magenta")

    alpha_bP, _ = plot_chi2_alpha(inpath2, "", [ax0, ax1], color="green", label="Parabola: 11 DOF")
    plot_chi2_alpha_evi(inpath2, "", [ax6, ax7], color="green", label="Parabola: 11 DOF")
    alpha_aP, sigma_P = plot_data_mystats(inpath2, ax2, label="Parabola", color="green")
    
    alpha_bT, _ = plot_chi2_alpha(inpath3, "", [ax0, ax1], color="blue", label="Parabola 2: 11 DOF")
    plot_chi2_alpha_evi(inpath3, "", [ax6, ax7], color="blue", label="Parabola 2: 11 DOF")
    alpha_aT, sigma_T = plot_data_mystats(inpath3, ax2, label="Parabola 2", color="blue")

    plot_all_comb([ax3[0][0],ax3[1][0]], alpha1=alpha_bG, alpha2=alpha_bP, alpha3=alpha_bT, label1=r"$\alpha_\mathrm{best, GAL}$", label2=r"$\alpha_\mathrm{best, PAR}$", label3=r"$\alpha_\mathrm{best, PAR2}$")
    plot_all_comb([ax3[0][1],ax3[1][1]], alpha1=alpha_aG, alpha2=alpha_aP, alpha3=alpha_aT, label1=r"$\alpha_\mathrm{GAL}$", label2=r"$\alpha_\mathrm{PAR}$", label3=r"$\alpha_\mathrm{PAR2}$")
    plot_all_comb(ax4, alpha1=sigma_G, alpha2=sigma_P, alpha3=sigma_T, label1=r"$\sigma_\mathrm{GAL}$", label2=r"$\sigma_\mathrm{PAR}$", label3=r"$\sigma_\mathrm{PAR2}$", mode=1)

    plot_alphabest_alpha(ax5[0], alphabest=alpha_bG, alpha=alpha_aG, labelbest=r"$\alpha_\mathrm{best, GAL}$", label=r"$\alpha_\mathrm{GAL}$")
    plot_alphabest_alpha(ax5[1], alphabest=alpha_bP, alpha=alpha_aP, labelbest=r"$\alpha_\mathrm{best, PAR}$", label=r"$\alpha_\mathrm{PAR}$")
    plot_alphabest_alpha(ax5[2], alphabest=alpha_bT, alpha=alpha_aT, labelbest=r"$\alpha_\mathrm{best, PAR2}$", label=r"$\alpha_\mathrm{PAR2}$")

    fig0.tight_layout()
    fig0.savefig("/home/astro/variu/GAL_PAR_PAR2_temp4.png")
    fig1.tight_layout()
    fig1.savefig("/home/astro/variu/GAL_PAR_PAR2_temp5.png")
    fig2.tight_layout()
    fig2.savefig("/home/astro/variu/GAL_PAR_PAR2_temp6.png")
    fig3.tight_layout()
    fig3.savefig("/home/astro/variu/GAL_PAR_PAR2_temp7.png")
    fig4.tight_layout()
    fig4.savefig("/home/astro/variu/GAL_PAR_PAR2_temp8.png")
    fig5.tight_layout()
    fig5.savefig("/home/astro/variu/GAL_PAR_PAR2_temp9.png")
    fig6.tight_layout()
    fig6.savefig("/home/astro/variu/GAL_PAR_PAR2_temp10.png")
    fig7.tight_layout()
    fig7.savefig("/home/astro/variu/GAL_PAR_PAR2_temp11.png")



def plot_all_bestfit():
    pp = PdfPages('/home/astro/variu/bestfit_alpha_1.05.pdf')
    filepath_cf = "/scratch/variu/clustering/patchy_cmass_subset/lightcone_box1/real/vv2pcf/16R/xi/"
    files = glob.glob(filepath_cf + "/CA*2pcf")
    s, avg, std = plot_errorbar(filepath_cf + "/CA*2pcf", "label", None, "grey", "-", plotBOOL=False)

    for i, file_ in enumerate(files):
        print(i)
        filename_cf = os.path.basename(file_)
        #try:
        plot_best(filepath_cf, filename_cf, pp, std)
        #except:
        #    print(file_)
        #    continue
        
    pp.close()

def main():
    ### Test the effect of the template
    # inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_100_temp_exte/"
    # inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_500_temp_exte/"
    # inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_short/"
    # inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_exte_filt/"
    # plot_tests_effect_template(inpath1, inpath2, inpath3, inpath4)
    # exit()
    ###

    inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/galaxy_60_150_m/"
    inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/parab_60_150_m_FFT/"
    inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000/"
    
    inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/galaxy_60_150_m/"
    inpath5 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/parab_60_150_m_fast/"
    inpath6 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/parab_60_150_m_fast2/"
    inpath7 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000_old/"
    inpath8 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_2000/"
    inpath9 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1/real/vv2pcf/G1024CIC_60_150_m_B2000/"
    
    #plot_all_bestfit()

    #obtain_chi2_alpha_evi(inpath7)
    #plot_three_cases(inpath1, inpath2, inpath3)
    #plot_three_cases(inpath4, inpath6, inpath8)
    #plot_three_cases_2(inpath7, inpath9, inpath8)
    plot_three_cases_3(inpath4, inpath5, inpath6)
    

    exit()

    #plot_1file_chi2alpha(inpath2 + "/BAOfit_CATALPTCICz0.466G960S1711891480.VOID.dat.2pcf_", 3)    
    #plot_1file_chi2alpha(inpath + "/BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_.txt", 3)
    
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