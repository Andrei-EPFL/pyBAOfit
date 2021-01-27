#import matplotlib as mpl
#mpl.use('Agg')import os, re, sys
import glob
import numpy as np
import matplotlib.pyplot as pt
from getdist import plots, MCSamples, loadMCSamples

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
    if plot == True:
        ax[0].plot(alpha, chi2, color=color, marker='.', ls='', label=label+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)))
        ax[0].set_xlabel(r"$\alpha$")
        ax[0].set_ylabel(r"$\chi^2$")
        ax[0].legend()
        
        #bins = np.linspace(0.8, 1.2, 161)
        ax[1].hist(alpha, alpha=0.4, facecolor=color, label=label+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)))#, fill=False, edgecolor=color)
        #ax[1].set_xlim([0.96, 1.04])
        ax[1].set_xlabel(r"$\alpha$")
        ax[1].set_ylabel(r"counts")
        ax[1].legend()

    return alpha, chi2

def plot_chi2_alpha_2(inpath, endname, ax, color="red", label="label", plot=True, n=1):
    alpha, chi2, evi = np.loadtxt(inpath+"/DAT_alpha_chi2_500"+endname+".dat",usecols=(0,1,2), unpack=True)
    chi2 = chi2/n
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

def obtain_chi2alpha(inpath):
    files = glob.glob(inpath + "*_.txt")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)

    print(len(files))
    chi2_arr = np.zeros(len(files))
    alpha_arr = np.zeros(len(files))
    evi_arr = np.zeros(len(files))
    
    for i, file_ in enumerate(files):
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".VOID.dat.2pcf_.txt"
        chi2_tmp, alpha_tmp = np.loadtxt(file_m, usecols=(1, 2), unpack=True)
        
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".VOID.dat.2pcf_mystats.txt"
        evi_tmp = np.loadtxt(file_m, usecols=(2), unpack=True)
        
        min_pos = chi2_tmp==np.min(chi2_tmp)
        
        chi2_arr[i] = chi2_tmp[min_pos]
        alpha_arr[i] = alpha_tmp[min_pos]
        evi_arr[i] = evi_tmp

    np.savetxt(inpath + "/DAT_alpha_chi2_500.dat", np.array([alpha_arr, chi2_arr, evi_arr]).T)

    

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
       

def plot_best(file_mystats, inpath_best, part):
    name = os.path.basename(file_mystats)[7:-12]
    mode = re.split('/',inpath_best)[-2]
    print(mode)
    s, cf = np.loadtxt("/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/" + name,usecols=(0,1), unpack=True)
    sb, cfb = np.loadtxt(inpath_best + "/BAOfit_" + name+"_best.dat",usecols=(0,1), unpack=True)
    sb1, cfb1 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_500temp/BAOfit_" + name+"_best.dat",usecols=(0,1), unpack=True)
    sb2, cfb2 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G512CIC_60_150_mocks/BAOfit_" + name+"_best.dat",usecols=(0,1), unpack=True)
    sb3, cfb3 = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_GL_60_150_mocks/BAOfit_" + name+"_best.dat",usecols=(0,1), unpack=True)
    
    fig, ax = pt.subplots()
    s, avg, std = plot_errorbar("/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/CA*2pcf", "label", ax, "grey", "-", plotBOOL=False)
    
    ax.errorbar(s, cf*s*s, yerr=s*s*std, label="data: "+ name, color="k")
    ax.plot(sb, sb*sb*cfb, label="Bestfit with high sigma:" + mode, ls="--")
    #ax.plot(sb1, sb1*sb1*cfb1, label="Bestfit: void_G1024CIC_60_150", ls="--")
    #ax.plot(sb2, sb2*sb2*cfb2, label="Bestfit: void_G512CIC_60_150", ls="--")
    #ax.plot(sb3, sb3*sb3*cfb3, label="Bestfit: void_GL_60_150", ls="--")
    
    ax.legend()

    range_ = np.where(np.logical_and(s <=165, s >=45))
    ymax = 1.1 * np.max(cf[range_] * s[range_]**2)
    ymin = 1.1 *np.min(cf[range_] * s[range_]**2)
    
    ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
    ax.set_ylabel(r"s$^2\xi$(s)")

    ax.set_ylim([ymin,ymax])
    ax.set_xlim(left=45, right=165)
    
    ax.axvline(x=60, color='grey', linestyle='--')
    ax.axvline(x=150, color='grey', linestyle='--')
    
    fig.savefig("./output/temporary"+part+".pdf")

def obtain_data(inpath, ax, label="label", color="red"):
    files = glob.glob(inpath + "*mystats*")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)

    print(len(files))
    alpha = np.zeros(len(files))
    sigma = np.zeros(len(files))
    for i, file_ in enumerate(files):
        file_m = inpath + "/BAOfit_CATALPTCICz0.466G960S" + str(int(randnum[i])) + ".VOID.dat.2pcf_mystats.txt"

        alpha[i], sigma[i], _ = np.loadtxt(file_m, usecols=(0,1,2), unpack=True)
        #if(sigma[i]>0.015):
        #    print(file_)
        #    plot_best(file_, inpath, label)
    
    ax[0].plot(alpha, sigma, ls="", marker=".", label=label, color=color)
    #ax[0].set_ylim([0.0045, 0.013])
    ax[0].set_xlabel("alpha")
    ax[0].set_ylabel("sigma")
    ax[0].legend()
    
    #bins = np.linspace(0.0045, 0.013, 50)
    ax[1].hist(sigma, alpha=0.4, label=label+": mean= {}; std={}".format(np.mean(sigma), np.std(sigma)), color=color)#, bins=bins)
    ax[1].set_xlabel("sigma")
    ax[1].set_ylabel("counts")
    ax[1].legend()

    ax[2].hist(alpha, alpha=0.4, label=label+": mean= {}; std={}".format(np.mean(alpha), np.std(alpha)), color=color)
    ax[2].set_xlabel("alpha")
    ax[2].set_ylabel("counts")
    ax[2].legend()

    return alpha, sigma

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


def plot_four_cases(inpath1, inpath2, inpath3, inpath4):
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
    alpham_100, sigma_100 = obtain_data(inpath1, ax, label="G1024CIC_100", color="red")
    alpham_500, sigma_500 = obtain_data(inpath2, ax, label="G1024CIC_500", color="blue")
    alpham_2000, sigma_2000 = obtain_data(inpath3, ax, label="G1024CIC_2000", color="orange")
    alpham_2000f, sigma_2000f = obtain_data(inpath4, ax, label="G1024CIC_2000filt", color="green")
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

def plot_three_cases(inpath1, inpath2, inpath3):
    ### Fig 1
    fig0, ax0 = pt.subplots()
    fig1, ax1 = pt.subplots()
    alpha_500, chi2_500 = plot_chi2_alpha_2(inpath2, "", [ax0, ax1], color="magenta", label="Galaxy: 12 DOF", n=12)
    alpha_2000, chi2_2000 = plot_chi2_alpha_2(inpath3, "", [ax0, ax1], color="green", label="Parabola: 11 DOF", n=11)
    alpha_100, chi2_100 = plot_chi2_alpha_2(inpath1, "", [ax0, ax1], color="blue", label="Template: 12 DOF", n=12)
    fig0.tight_layout()
    fig0.savefig("/home/astro/variu/temp.png")
    fig1.tight_layout()
    fig1.savefig("/home/astro/variu/temp1.png")

def main():
    #inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_100_temp_exte/"
    #inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_500_temp_exte/"
    #inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_short/"
    #inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks_2000_temp_exte_filt/"
    
    #inpath5 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/parab_60_150_mocks/"
    #inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/galaxy_60_150_fast_mocks/"
    #inpath4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G512CIC_60_150_mocks/"
    #inpath7 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/16R/void_G1024CIC_60_150_mocks_2001_temp_short/"
    
    inpath1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000/"
    inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/galaxy_60_150_m/"
    inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/parab_60_150_m_FFT/"
    

    #obtain_chi2alpha(inpath7)
    #obtain_chi2alpha(inpath1)
    #obtain_chi2alpha(inpath2)
    #obtain_chi2alpha(inpath3)
    #exit()
    #plot_four_cases(inpath1, inpath2, inpath3, inpath4)
    plot_three_cases(inpath1, inpath2, inpath3)
    exit()

    #plot_1file_chi2alpha(inpath2 + "/BAOfit_CATALPTCICz0.466G960S1711891480.VOID.dat.2pcf_", 3)    
    #plot_1file_chi2alpha(inpath + "/BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_.txt", 3)
    
    #compute_covariance([inpath5, inpath, inpath6])
    #exit()

    fig, ax = pt.subplots(1,3, figsize=(30,10))
    
    obtain_data(inpath2, ax, label="Galaxy", color="magenta")
    obtain_data(inpath3, ax, label="Parabola", color="green")
    obtain_data(inpath1, ax, label="Template", color="blue")
    
    fig.savefig("/home/astro/variu/temp_parab.pdf")
    exit()
    
    fig, ax = pt.subplots(1, 2, figsize=(20,10))
    #plot_chi2_alpha(inpath, "", ax, color="red", label="G1024CIC")
    #plot_chi2_alpha(inpath5, "", ax, color="green", label="G1024IC_100")
    plot_chi2_alpha(inpath7, "", ax, color="red", label="G1024CIC_2001")
    
    #plot_chi2_alpha(inpath2, "_voidGL", ax, color="blue", label="G2048L512")
    #plot_chi2_alpha(inpath3, "_galaxy", ax, color="orange", label="galaxy")

    fig.savefig("/home/astro/variu/vhxcf_temporar4.pdf")

if __name__== '__main__':
    main()