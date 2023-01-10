#import matplotlib as mpl
#mpl.use('Agg')
import os, re, sys
import glob
import numpy as np
import matplotlib.pyplot as pt
from getdist import plots, MCSamples, loadMCSamples
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import matplotlib.font_manager as fm
import matplotlib as mpl

def obtain_chi2_alpha_evi(inpath, endcffile=".dat.2pcf", begcffile=""):
    files = glob.glob(inpath + "*_.txt")
    # randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_1000.txt", usecols=(0), unpack=True)

    print(len(files))
    chi2_arr = np.zeros(len(files))
    alpha_arr = np.zeros(len(files))
    evi_arr = np.zeros(len(files))
    
    for i, file_ in enumerate(files):
        file_m = inpath + f"/BAOfit_{begcffile}CATALPTCICz0.466G960S" + str(int(randnum[i])) + endcffile + "_.txt"

        if os.path.isfile(file_m):
            chi2_tmp, alpha_tmp = np.loadtxt(file_m, usecols=(1, 2), unpack=True)
            
            file_m = inpath + f"/BAOfit_{begcffile}CATALPTCICz0.466G960S" + str(int(randnum[i])) + endcffile + "_mystats.txt"
            if os.path.isfile(file_m):
            
                evi_tmp = np.loadtxt(file_m, usecols=(2), unpack=True)
                
                min_pos = np.argmin(chi2_tmp)
                
                chi2_arr[i] = chi2_tmp[min_pos]
                alpha_arr[i] = alpha_tmp[min_pos]
                evi_arr[i] = evi_tmp
            else:
                print(file_m)
        else:
            print(file_m)
    np.savetxt(inpath + "/DAT_alpha_chi2_500.dat", np.array([alpha_arr, chi2_arr, evi_arr]).T)


def read_chi2_alpha_evi(inpath, endname, n=1):
    alpha, chi2, evi = np.loadtxt(inpath+"/DAT_alpha_chi2_500"+endname+".dat",usecols=(0,1,2), unpack=True)

    chi2 = chi2/n
    range_ = alpha>0.1

    chi2 = chi2[range_]
    evi = evi[range_]
    alpha = alpha[range_]

    print("Lines in the the DAT_alpha_chi2_500 file: ", len(alpha))

    return alpha, chi2, evi


def read_data_mystats(inpath, endcffile=".dat.2pcf", begcffile=""):
    files = glob.glob(inpath + "*mystats*")
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_1000.txt", usecols=(0), unpack=True)

    print(f"There are {len(files)} files in the folder")
    alpha = np.zeros(len(files))
    alpham = np.zeros(len(files))
    sigma = np.zeros(len(files))
    sigmanl = np.zeros(len(files))
    for i, file_ in enumerate(files):
        files = glob.glob(inpath + f"/BAOfit_*CATALPTCICz0.466G960S" + str(int(randnum[i])) + endcffile + "_mystats.txt")
        if len(files) == 1:
            file_m = files[0]
        else:
            print(len(files))
            print("there are more files with the same seed")
            sys.exit()
        # alpha[i], sigma[i], alpham[i], sigmanl[i] = np.loadtxt(file_m, usecols=(0, 1, 4, 7), unpack=True)
        alpha[i], sigma[i], alpham[i] = np.loadtxt(file_m, usecols=(0, 1, 4), unpack=True)
    return alpha, sigma, alpham, sigmanl


def determine_xmin_xmax(fraction, xmin, xmax):
    xmin_tmp, xmax_tmp = xmin, xmax
    if xmin > np.min(fraction):
        xmin_tmp = np.min(fraction)
    if xmax < np.max(fraction):
        xmax_tmp = np.max(fraction)
    
    return xmin_tmp, xmax_tmp


def ax_style_my_corner_UL(ax, xmin, xmax, keys, outliers=False):
    outlier = 0
    if outliers:
        outlier = 1

    for row in range(len(keys)):
        for col in range(len(keys)):
            tcol = col + outlier
            ax[row][tcol].set_xlim([0.99 * xmin, 1.01 * xmax])
            if row != col:
                ax[row][tcol].set_ylim([0.99 * xmin, 1.01 * xmax])


def ax_style_my_corner_L(ax, xmin, xmax, keys, outliers=False):
    outlier = 0
    if outliers:
        outlier = 1

    for row in range(len(keys)):
        for col in range(row + 1):
            tcol = col + outlier
            ax[row][tcol].set_xlim([0.99 * xmin, 1.01 * xmax])
            if row != col:
                ax[row][tcol].set_ylim([0.99 * xmin, 1.01 * xmax])


def plot_my_corner_lower(ax, data_dict, dict_labels, keys, outliers=False):
    xmin = 100
    xmax = 0
    min_alpha = 0.7
    max_alpha = 1.3
    outlier = 0
    if outliers:
        min_alpha = 0.9
        max_alpha = 1.1
        outlier = 1

    for i, key in enumerate(keys):
        range_ = (data_dict[key]["alphamed"] > min_alpha) & (data_dict[key]["alphamed"] < max_alpha)
        print("There are  values above 1.1 or below 0.9", data_dict[key]["alphamed"][~range_])
        xmin, xmax = determine_xmin_xmax(data_dict[key]["alphamed"][range_], xmin, xmax)
        
    for i, key in enumerate(keys):
        ti = i + outlier
        # ax[i][ti].text(1, 1, dict_labels[keys[i]]["label"] + " " + dict_labels[keys[i]]["label"])  # this is for debugging
        ax[i][ti].hist(data_dict[key]["alphamed"], bins=dict_labels["alphamed"]["bins"], density=True, histtype="step", color="blue")
        ax[i][ti].axvline(1, color="grey", ls="--")

    for row, key in enumerate(keys):
        for col in range(row):
            tcol = col + outlier
            # ax[row][tcol].text(0.98, 1, dict_labels[keys[col]]["label"] + " " + dict_labels[keys[row]]["label"])  # this is for debugging
            ax[row][tcol].scatter(data_dict[keys[col]]["alphamed"], data_dict[keys[row]]["alphamed"], s=1, color="blue")
            ax[row][tcol].plot([xmin, xmax], [xmin, xmax], color="k", ls="--")
    
    tcol = outlier
    
    for row in range(len(keys) - 1):
        ax[row + 1][tcol].set_yticks(dict_labels["alphamed"]["ticks"])
        ax[row + 1][tcol].set_yticklabels(dict_labels["alphamed"]["ticks"])
   
    for col in range(len(keys)):
        tcol = col + outlier
        ax[-1][tcol].set_xticks(dict_labels["alphamed"]["ticks"])
        ax[-1][tcol].set_xticklabels(dict_labels["alphamed"]["ticks"], rotation=30)

    for i, key in enumerate(keys):
        ti = i + outlier
        ax[i][0].set_ylabel("${0}$".format(dict_labels[keys[i]]["alphamed"]))
        ax[-1][ti].set_xlabel("${0}$".format(dict_labels[keys[i]]["alphamed"]))
    
    ax[0][0].set_ylabel("")

    if outliers:
        for row in range(len(keys) - 1):
            ax[row + 1][0].set_xlim([0.79,0.84])
            ax[row + 1][0].scatter(data_dict[keys[0]]["alphamed"], data_dict[keys[row + 1]]["alphamed"], s=1, color="blue")
        
            ax[-1][0].set_xticks([0.8])
            ax[-1][0].set_xticklabels([0.8], rotation=30)
            print(ax[-1][0].get_xticks())
    
    return xmin, xmax


def plot_my_corner_upper(ax, data_dict, dict_labels, keys, outliers=False):
    xmin = 100
    xmax = 0
    min_alpha = 0.7
    max_alpha = 1.3
    outlier = 0
    if outliers:
        min_alpha = 0.9
        max_alpha = 1.1
        outlier = 1

    for i, key in enumerate(keys):
        range_ = (data_dict[key]["alphamed"] > min_alpha) & (data_dict[key]["alphamed"] < max_alpha)
        print("There are  values above 1.1 or below 0.9", data_dict[key]["alphamed"][~range_])
        xmin, xmax = determine_xmin_xmax(data_dict[key]["alphamed"][range_], xmin, xmax)
        
    
    for i, key in enumerate(keys):
        ti = i + outlier
        # ax[i][ti].text(1, 1, dict_labels[keys[i]]["label"] + " " + dict_labels[keys[i]]["label"])  # this is for debugging
        ax[i][ti].hist(data_dict[key]["alphamed"], bins=dict_labels["alphamed"]["bins"], density=True, histtype="step", color="red")
        ax[i][ti].axvline(1, color="grey", ls="--")

    for row, key in enumerate(keys):
        for col in range(row):
            trow = row + outlier
            # ax[col][trow].text(0.98, 1, "x: "+ dict_labels[keys[row]]["label"] + "; y: " + dict_labels[keys[col]]["label"])  # this is for debugging
            ax[col][trow].scatter(data_dict[keys[row]]["alphamed"], data_dict[keys[col]]["alphamed"], s=1, color="red")
            ax[col][trow].plot([xmin, xmax], [xmin, xmax], color="k", ls="--")

    tcol = outlier    
    for row in range(len(keys) - 1):
        ax[row][-1].set_yticks(dict_labels["alphamed"]["ticks"])
        ax[row][-1].set_yticklabels(dict_labels["alphamed"]["ticks"])
        ax[row][-1].tick_params(axis="y", direction="out", left=False, right=True, labelleft=False, labelright=True)
        ax[row][-1].set_ylabel("${0}$".format(dict_labels[keys[row]]["alphamed"]))
        ax[row][-1].yaxis.set_label_position("right")

    for i in range(len(keys)):
        ti = i + outlier
        ax[0][ti].set_xticks(dict_labels["alphamed"]["ticks"])
        ax[0][ti].set_xticklabels(dict_labels["alphamed"]["ticks"], rotation=30)
        ax[0][ti].tick_params(axis="x", direction="out", bottom=False, top=True, labelbottom=False, labeltop=True)
        ax[0][ti].set_xlabel("${0}$".format(dict_labels[keys[i]]["alphamed"]))
        ax[0][ti].xaxis.set_label_position("top")
       
    return xmin, xmax


def plot_evi_tension(fig, ax, data_dict, dict_labels, keys, evi_chi2="evi"):

    xmin = 100
    xmax = 0
    
    # remove all ticks and tick labels
    for row in range(len(keys)):
        for col in range(len(keys)):
            ax[col][row].set_yticklabels([])
            ax[col][row].set_yticks([])
            ax[col][row].set_xticklabels([])
            ax[col][row].set_xticks([])

    # add ticks and tick labels only for the first row and last rows
    for col in range(len(keys) - 1):
        ax[-1][col].set_xticks(dict_labels["tensionparams"]["ticks"])
        ax[-1][col].set_xticklabels(dict_labels["tensionparams"]["ticks"], rotation=30)
        ax[0][col + 1].set_xticks(dict_labels["bayesfactor"]["ticks"])
        ax[0][col + 1].set_xticklabels(dict_labels["bayesfactor"]["ticks"], rotation=30)
        
    # the diagonal plots
    x_n = np.linspace(-4., 4., 1001)
    for i in range(len(keys)):
        range_1 = data_dict[keys[i]]["sigma_nl"] < 9
        fraction_ = (data_dict[keys[i]]["alphamed"] - 1) / data_dict[keys[i]]["sigma"]
        range_ = (fraction_ > dict_labels["pullone"]["maxval"]) | (fraction_ < dict_labels["pullone"]["minval"])
        counter = len(fraction_[range_])
        print(f"# cases when sigma_nl < 9 = {len(fraction_[range_1])}")
        ax[i][i].hist(fraction_[range_1], color="green", density=True, histtype="step", bins=dict_labels["pullone"]["bins"])#, label=str(counter))
        # ax[i][i].hist(fraction_, color="orange", ls="--", density=True, histtype="step", bins=dict_labels["pullone"]["bins"])#, label=str(counter))
        ax[i][i].plot(x_n, gauss(x_n, 0, 1), ls="--", color="k")

        ax[i][i].tick_params(axis="x", direction="in", bottom=True, top=True, labelbottom=False, labeltop=True)
        ax[i][i].set_xticks(dict_labels["pullone"]["ticks"])

        # ax[i][i].legend(loc="upper right")
    
    ax[-1][-1].tick_params(axis="x", direction="in", bottom=True, top=True, labelbottom=True, labeltop=False)
    ax[0][0].set_xticklabels(dict_labels["pullone"]["ticks"], rotation=30)
    ax[-1][-1].set_xticklabels(dict_labels["pullone"]["ticks"], rotation=30)
        
    xmin_l = 100
    xmax_l = -100
    
    xmin_u = 100
    xmax_u = -100
    
    for row, key in enumerate(keys):
        for col in range(row):
            # lower triangular plot
            # ax[row][col].text(-0.5, 0.1, dict_labels[keys[col]]["label"] + "-" + dict_labels[keys[row]]["label"])  # this is for debugging
            upper_val = (data_dict[keys[col]]["alphamed"] - data_dict[keys[row]]["alphamed"])
            lower_val = np.sqrt(data_dict[keys[col]]["sigma"]**2 + data_dict[keys[row]]["sigma"]**2)
            fraction = upper_val / lower_val
            range_ = (fraction < dict_labels["tensionparams"]["maxval"]) & (fraction > dict_labels["tensionparams"]["minval"])
            
            xmin_l, xmax_l = determine_xmin_xmax(fraction[range_], xmin_l, xmax_l)

            ax[row][col].hist(fraction, bins=dict_labels["tensionparams"]["bins"], color="red", density=True, histtype="step", label=str(len(fraction[~range_])))
            # ax[row][col].legend()
            ax[row][col].axvline(0, color="k", ls="--")

            # upper triangular plot
            # ax[col][row].text(-0.5, 0.1, dict_labels[keys[col]]["label"] + "-" + dict_labels[keys[row]]["label"])  # this is for debugging
            delta_evi = data_dict[keys[col]][evi_chi2] - data_dict[keys[row]][evi_chi2]
            range_ = (delta_evi < dict_labels["bayesfactor"]["maxval"]) & (delta_evi > dict_labels["bayesfactor"]["minval"])

            xmin_u, xmax_u = determine_xmin_xmax(delta_evi[range_], xmin_u, xmax_u)

            ax[col][row].hist(delta_evi, bins=dict_labels["bayesfactor"]["bins"], color="blue", density=True, histtype="step", label=str(len(delta_evi[~range_])) )
            # ax[col][row].legend()
            ax[col][row].tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)
            ax[col][row].tick_params(axis="y", left=False, right=False, labelleft=False, labelright=True)
            ax[col][row].axvline(0, color="k", ls="--")


    for row, key in enumerate(keys):
        for col in range(row):
            ax[row][col].set_xlim([0.995 * xmin_l, xmax_l * 1.005])
            ax[col][row].set_xlim([0.995 * xmin_u, xmax_u * 1.005])

    for i, key in enumerate(keys):
        # lower triangular plot
        ax[i][0].set_ylabel((dict_labels[keys[i]]["label"]))
        ax[-1][i].set_xlabel((dict_labels[keys[i]]["label"]))
    
        # upper triangular plot
        ax[i][-1].set_ylabel((dict_labels[keys[i]]["label"]))
        ax[i][-1].yaxis.set_label_position("right")

        ax[0][i].set_xlabel((dict_labels[keys[i]]["label"]))
        ax[0][i].xaxis.set_label_position("top")

    ax[0][0].set_ylabel("")
    ax[0][0].set_xlabel("")
    ax[-1][-1].set_ylabel("")
    ax[-1][-1].set_xlabel("")
    
    line_l = Line2D([0], [0], color="red", lw=1)
    line_u = Line2D([0], [0], color="blue", lw=1)
    line_g = Line2D([0], [0], color="green", lw=1)
        
    # ax00 = ax[0][0].twinx()
    # ax00.set_yticklabels([])
    # ax00.set_yticks([])
    # ax[0][0].legend([line_g, line_l, line_u], [r"$\frac{\alpha_x-1}{\sigma_x}$", r"$\frac{\alpha_x-\alpha_y}{\sqrt{\sigma_x^2 + \sigma_y^2}}$", r"$\ln\left(\frac{\mathcal{Z}_y}{\mathcal{Z}_x}\right)$"], loc="upper right", bbox_to_anchor=(0, 0.5), fontsize="large", frameon=False)
    ax[0][0].legend([line_g, line_l, line_u], [r"$\frac{\alpha_x-1}{\sigma_x}$", r"$\frac{\alpha_x-\alpha_y}{\sqrt{\sigma_x^2 + \sigma_y^2}}$", r"$\ln\left(\frac{\mathcal{Z}_y}{\mathcal{Z}_x}\right)$"], bbox_to_anchor=(0, 2.5), loc='upper left', ncol=3, fontsize="x-large", frameon=False)


def plot_deltaalpha_deltasigma(ax, data_dict, dict_labels, keys):
    
    bins_1 = np.linspace(-0.05, 0.05, 501)
    bins_2 = np.linspace(-2, 2, 501)
    bins_3 = np.linspace(-2, 2, 201)
    xmin_1, xmax_1 = 100, -100
    xmin_2, xmax_2 = 100, -100
    xmin_3, xmax_3 = 100, -100
    legend=""
    
    x_n = np.linspace(-4., 4., 1001)
    
    for row in range(len(keys) - 1):
        ax[row][0].set_ylabel(dict_labels[keys[row + 1]]["label"])

        # Delta alpha / alpha
        upper_val = (data_dict[keys[row + 1]]["alphamed"] - data_dict[keys[0]]["alphamed"])
        fraction = upper_val / data_dict[keys[0]]["alphamed"]
        range_ = (fraction < dict_labels["deltaalphaoveralpha"]["maxval"]) & (fraction > dict_labels["deltaalphaoveralpha"]["minval"])

        ax[row][0].hist(fraction[range_], bins=dict_labels["deltaalphaoveralpha"]["bins"], color=dict_labels["color"], density=True, histtype="step", label=legend+" "+str(len(fraction[~range_])))
        ax[row][0].axvline(0, color="k", ls="--")
        xmin_1, xmax_1 = determine_xmin_xmax(fraction[range_], xmin_1, xmax_1)
        # ax[row][0].legend()
        
        # Delta sigma / sigma
        upper_val = (data_dict[keys[row + 1]]["sigma"] - data_dict[keys[0]]["sigma"])
        fraction = upper_val / data_dict[keys[0]]["sigma"]
        range_ = (fraction < dict_labels["deltasigmaoversigma"]["maxval"]) & (fraction > dict_labels["deltasigmaoversigma"]["minval"])
        
        ax[row][1].hist(fraction[range_], bins=dict_labels["deltasigmaoversigma"]["bins"], color=dict_labels["color"], density=True, histtype="step", label=str(len(fraction[~range_])))
        ax[row][1].axvline(0, color="k", ls="--")
        xmin_2, xmax_2 = determine_xmin_xmax(fraction[range_], xmin_2, xmax_2)
        # ax[row][1].legend()
        
        # Delta alpha / sqrt(sigma1**2 + sig2 **2) tension
        upper_val = (data_dict[keys[row + 1]]["alphamed"] - data_dict[keys[0]]["alphamed"])
        lower_val = np.sqrt(data_dict[keys[0]]["sigma"]**2 + data_dict[keys[row + 1]]["sigma"]**2)
        fraction = upper_val / lower_val
        range_ = (fraction < dict_labels["deltaalphaoveravgsigma"]["maxval"]) & (fraction > dict_labels["deltaalphaoveravgsigma"]["minval"])
        
        ax[row][2].hist(fraction[range_], bins=dict_labels["deltaalphaoveravgsigma"]["bins"], color=dict_labels["color"], density=True, histtype="step", label=str(len(fraction[~range_])))
        ax[row][2].axvline(0, color="k", ls="--")
        xmin_3, xmax_3 = determine_xmin_xmax(fraction[range_], xmin_3, xmax_3)
        # ax[row][2].legend()
    
        # (alpha - 1) / sigma
        upper_val = (data_dict[keys[row + 1]]["alphamed"] - 1)
        lower_val = data_dict[keys[row + 1]]["sigma"]
        fraction = upper_val / lower_val
        range_ = (fraction < dict_labels["alphamoneoversigma"]["maxval"]) & (fraction > dict_labels["alphamoneoversigma"]["minval"])
        
        ax[row][3].hist(fraction[range_], bins=dict_labels["alphamoneoversigma"]["bins"], color=dict_labels["color"], density=True, histtype="step", label=str(len(fraction[~range_])))
        ax[row][3].plot(x_n, gauss(x_n, 0, 1), ls="--", color="k")
        # xmin_3, xmax_3 = determine_xmin_xmax(fraction[range_], xmin_3, xmax_3)
        # ax[row][3].legend()
    


    ax[-1][0].set_xticks(dict_labels["deltaalphaoveralpha"]["ticks"])
    ax[-1][0].set_xticklabels(dict_labels["deltaalphaoveralpha"]["ticks"], rotation=30)
    
    ax[-1][1].set_xticks(dict_labels["deltasigmaoversigma"]["ticks"])
    ax[-1][1].set_xticklabels(dict_labels["deltasigmaoversigma"]["ticks"], rotation=30)
    
    ax[-1][2].set_xticks(dict_labels["deltaalphaoveravgsigma"]["ticks"])
    ax[-1][2].set_xticklabels(dict_labels["deltaalphaoveravgsigma"]["ticks"], rotation=30)
   
    ax[-1][3].set_xticks(dict_labels["alphamoneoversigma"]["ticks"])
    ax[-1][3].set_xticklabels(dict_labels["alphamoneoversigma"]["ticks"], rotation=30)

    ax[-1][0].set_xlabel("$\\frac{{\\alpha_y}}{{{lower}}} - 1$".format(lower=str(dict_labels[keys[0]]["alphamed"])), fontsize=14)
    ax[-1][1].set_xlabel("$\\frac{{\\sigma_y}}{{{lower}}} - 1$".format(lower=dict_labels[keys[0]]["sigma"]), fontsize=14)
    ax[-1][2].set_xlabel("$\\frac{{\\alpha_y-{upper} }}{{\\sqrt{{\\sigma_y^2 + {lower}^2}} }}$".format(upper=dict_labels[keys[0]]["alphamed"], lower=dict_labels[keys[0]]["sigma"]), fontsize = 14)
    ax[-1][3].set_xlabel("$\\frac{\\alpha_y-1}{\\sigma_y}$", fontsize = 14)

    ax[0][0].set_xlim([1.01 * xmin_1, xmax_1 * 1.01])
    ax[0][1].set_xlim([1.01 * xmin_2, xmax_2 * 1.01])
    ax[0][2].set_xlim([1.01 * xmin_3, xmax_3 * 1.01])
    ax[0][3].set_xlim([dict_labels["alphamoneoversigma"]["minval"] * 1.02, dict_labels["alphamoneoversigma"]["maxval"] * 1.02])
    
    return [xmin_1, xmin_2, xmin_3], [xmax_1, xmax_2, xmax_3]
    
    
def plot_comparison(dict_labels, keys, outliers=False):
    from mycorner_lib import my_corner, my_corner_outliers, my_box

    if outliers:
        size = len(keys)
        fig1, ax1 = my_corner_outliers(size=size)
    else:
        fig1, ax1 = my_corner(size=len(keys))
    fig2, ax2 = my_box(rows=len(keys), cols=len(keys))
    # fig4, ax4 = my_box(rows=(len(keys) - 1), cols=4, sharex="col", sharey="col", figsize=(5, len(keys)*4.4/4.))
    fig4, ax4 = my_box(rows=(len(keys) - 1), cols=4, sharex="col", sharey="col", figsize=(5, 4))

    data_dict = {}

    for i, key in enumerate(keys):
        dict_ = dict_labels[key]
        
        alpha_b, chi2, evi = read_chi2_alpha_evi(dict_["filename"], "", n=1)
        alpha_a, sigma_, alpha_m, sigma_nl = read_data_mystats(dict_["filename"], endcffile=dict_labels["endcffile"], begcffile=dict_labels["begcffile"])

        temp_dict = {"alpha":alpha_a, "alphabest":alpha_b, "alphamed":alpha_m, "sigma":sigma_, "chi2":chi2, "evi":evi, "sigma_nl":sigma_nl}
        data_dict[key] = temp_dict

    xmin_L, xmax_L = plot_deltaalpha_deltasigma(ax4, data_dict, dict_labels, keys)
    
    print(xmin_L, xmax_L)
    ax4[0][0].set_xlim([0.99 * xmin_L[0], xmax_L[0] * 1.01])
    ax4[0][1].set_xlim([0.99 * xmin_L[1], xmax_L[1] * 1.01])
    ax4[0][2].set_xlim([0.99 * xmin_L[2], xmax_L[2] * 1.01])

    xmin, xmax = plot_my_corner_lower(ax1, data_dict, dict_labels, keys, outliers=outliers)
    ax_style_my_corner_L(ax1, xmin, xmax, keys, outliers=outliers)

    plot_evi_tension(fig2, ax2, data_dict, dict_labels, keys)


    # plot_evi_tension(ax3, data_dict, dict_labels, keys, evi_chi2="chi2")

    fig1.tight_layout()
    fig1.savefig(dict_labels['filename']+"_corner1.pdf", bbox_inches='tight')
    
    fig2.tight_layout()
    fig2.savefig(dict_labels['filename']+"_corner2.pdf", bbox_inches='tight')
    
    fig4.tight_layout()
    fig4.savefig(dict_labels['filename']+"_corner4.pdf", bbox_inches='tight')
    

def plot_corner_2outliers(dict_labels, keys):
    from mycorner_lib import my_corner, my_corner_outliers_2sides

    size = len(keys)
    fig, ax = my_corner_outliers_2sides(size=size)


    for row in range(1, len(keys)):
        ax[row][0].set_yticks(dict_labels["alphamed"]["ticks"])
        ax[row][0].set_yticklabels(dict_labels["alphamed"]["ticks"])
   
   
    for col in range(len(keys)):
        if col == 0:
            outlier = 1
        else:
            outlier = 2

        tcol = col + outlier
        ax[-1][tcol].set_xticks(dict_labels["alphamed"]["ticks"])
        ax[-1][tcol].set_xticklabels(dict_labels["alphamed"]["ticks"], rotation=30)

    for i, key in enumerate(keys):
        if i == 0:
            outlier = 1
        else:
            outlier = 2
        ti = i + outlier
        ax[i][0].set_ylabel("${0}$".format(dict_labels[keys[i]]["alphamed"]))
        ax[-1][ti].set_xlabel("${0}$".format(dict_labels[keys[i]]["alphamed"]))
    
    ax[0][0].set_ylabel("")


    data_dict = {}

    for i, key in enumerate(keys):
        dict_ = dict_labels[key]
        
        # alpha_b, chi2, evi = read_chi2_alpha_evi(dict_["filename"], "", n=1)
        alpha_a, sigma_, alpha_m, sigma_nl = read_data_mystats(dict_["filename"], endcffile=dict_labels["endcffile"], begcffile=dict_labels["begcffile"])
        # alpha_a = alpha_m = sigma_ = sigma_nl = np.linspace(0.8, 1.3, 100)
        temp_dict = {"alpha":alpha_a, "alphamed":alpha_m, "sigma":sigma_,  "sigma_nl":sigma_nl}
        data_dict[key] = temp_dict


    xmin = 100
    xmax = 0
    min_alpha = 0.9
    max_alpha = 1.1
    outlier = 1

    for i, key in enumerate(keys):
        range_ = (data_dict[key]["alphamed"] > min_alpha) & (data_dict[key]["alphamed"] < max_alpha)
        print("There are  values above 1.1 or below 0.9", data_dict[key]["alphamed"][~range_])
        xmin, xmax = determine_xmin_xmax(data_dict[key]["alphamed"][range_], xmin, xmax)
    
    print(xmin, xmax)
    # xmin = 0 
    # xmax = 2
    for i, key in enumerate(keys):
        if i == 0:
            outlier = 1
        else:
            outlier = 2
        ti = i + outlier
        # ax[i][ti].text(1, 1, dict_labels[keys[i]]["label"] + " " + dict_labels[keys[i]]["label"])  # this is for debugging
        ax[i][ti].hist(data_dict[key]["alphamed"], bins=dict_labels["alphamed"]["bins"], density=True, histtype="step", color="blue")
        ax[i][ti].axvline(1, color="grey", ls="--")
        ax[i][ti].set_xlim([0.99 * xmin, 1.01 * xmax])


    for row, key in enumerate(keys):
        for col in range(row):
            if col == 0:
                outlier = 1
            else:
                outlier = 2
            
            tcol = col + outlier
            # ax[row][tcol].text(0.9, 1, dict_labels[keys[col]]["label"] + " " + dict_labels[keys[row]]["label"])  # this is for debugging
            ax[row][tcol].scatter(data_dict[keys[col]]["alphamed"], data_dict[keys[row]]["alphamed"], s=1, color="blue")
            ax[row][tcol].plot([xmin, xmax], [xmin, xmax], color="k", ls="--")
            ax[row][tcol].set_xlim([0.99 * xmin, 1.01 * xmax])
            ax[row][tcol].set_ylim([0.99 * xmin, 1.01 * xmax])


    range_ = (data_dict[keys[0]]["alphamed"] >= 0.8) & (data_dict[keys[0]]["alphamed"] <= 0.9)
    print(np.min(data_dict[keys[0]]["alphamed"][range_]), np.max(data_dict[keys[0]]["alphamed"][range_]))

    range_ = (data_dict[keys[0]]["alphamed"] >= 1.1) & (data_dict[keys[0]]["alphamed"] <= 1.2)
    print(np.min(data_dict[keys[0]]["alphamed"][range_]), np.max(data_dict[keys[0]]["alphamed"][range_]))
    
    ax[0][0].hist(data_dict[keys[0]]["alphamed"], bins=dict_labels["alphamed"]["bins"], density=True, histtype="step", color="blue")
    ax[0][2].hist(data_dict[keys[0]]["alphamed"], bins=dict_labels["alphamed"]["bins"], density=True, histtype="step", color="blue")
    ax[0][0].set_xlim([0.79, 0.90])
    ax[0][2].set_xlim([1.14, 1.20])


    for row in range(len(keys) - 1):
        ax[row + 1][0].set_xlim([0.79, 0.90])
        ax[row + 1][2].set_xlim([1.14, 1.20])

        ax[row + 1][0].scatter(data_dict[keys[0]]["alphamed"], data_dict[keys[row + 1]]["alphamed"], s=1, color="blue")
        ax[row + 1][2].scatter(data_dict[keys[0]]["alphamed"], data_dict[keys[row + 1]]["alphamed"], s=1, color="blue")
    
        ax[-1][0].set_xticks([0.85])
        ax[-1][2].set_xticks([1.17])
        ax[-1][0].set_xticklabels([0.85], rotation=30)
        ax[-1][2].set_xticklabels([1.17], rotation=30)
        # print(ax[-1][0].get_xticks())
    
    

    fig.tight_layout()
    fig.savefig(dict_labels['filename']+"_corner0.pdf", bbox_inches='tight')
    

def plot_dampA1_dampA2(dict_labels, dict_labels2, keys="test", outpath="test"):
    fig, ax = pt.subplots(2, 2, figsize=(5, 4), gridspec_kw={"hspace":0, "wspace":0})

    ### remove all ticks and tick labels
    for row in range(len(keys)):
        for col in range(len(keys)):
            ax[col][row].set_yticklabels([])
            ax[col][row].set_yticks([])
            ax[col][row].set_xticklabels([])
            ax[col][row].set_xticks([])

    ### Get the data
    data_dict = {}
    for i, key in enumerate(keys):
        dict_ = dict_labels[key]
        alpha_a, sigma_, alpha_m, _ = read_data_mystats(dict_["filename"], endcffile=dict_labels["endcffile"], begcffile=dict_labels["begcffile"])
        temp_dict = {"alpha_a":alpha_a, "alphamed":alpha_m, "sigma":sigma_}
        data_dict[key] = temp_dict

    data_dict2 = {}
    for i, key in enumerate(keys):
        dict_2 = dict_labels2[key]
        alpha_a, sigma_, alpha_m, _ = read_data_mystats(dict_2["filename"], endcffile=dict_labels2["endcffile"], begcffile=dict_labels2["begcffile"])
        temp_dict2 = {"alpha_a":alpha_a, "alphamed":alpha_m, "sigma":sigma_}
        data_dict2[key] = temp_dict2

    ### the first row plots
    x_n = np.linspace(-4., 4., 1001)
    for i in range(len(keys)):

        fraction_ = (data_dict[keys[i]]["alphamed"] - 1) / data_dict[keys[i]]["sigma"]
        range_ = (fraction_ > dict_labels["pullone"]["maxval"]) | (fraction_ < dict_labels["pullone"]["minval"])
        ax[0][i].hist(fraction_, color=dict_labels["color"], density=True, histtype="step", bins=dict_labels["pullone"]["bins"], label=f"{len(fraction_[range_])}")# {round(np.mean(fraction_[~range_]), 5)}")
        
        
        fraction_2 = (data_dict2[keys[i]]["alphamed"] - 1) / data_dict2[keys[i]]["sigma"]
        range_2 = (fraction_2 > dict_labels2["pullone"]["maxval"]) | (fraction_2 < dict_labels2["pullone"]["minval"])
        ax[0][i].hist(fraction_2, color=dict_labels2["color"], density=True, histtype="step", bins=dict_labels2["pullone"]["bins"], label=f"{len(fraction_2[range_2])}")# {round(np.mean(fraction_2[~range_2]), 5)}")
        

        ax[0][i].plot(x_n, gauss(x_n, 0, 1), ls="--", color="k")
        ax[0][i].set_xticks(dict_labels["pullone"]["ticks"])
        ax[0][i].set_xticklabels(dict_labels["pullone"]["ticks"], rotation=30)
        
        # ax[0][i].legend(loc="upper right")


    # the second row of plots
    upper_val =  data_dict[keys[0]]["sigma"] - data_dict[keys[1]]["sigma"]
    lower_val = (data_dict[keys[0]]["sigma"] + data_dict[keys[1]]["sigma"]) / 2.
    fraction_ = upper_val / lower_val
    range_ = (fraction_ > dict_labels["deltasigmaoversigma"]["maxval"]) | (fraction_ < dict_labels["deltasigmaoversigma"]["minval"])

    upper_val2 =  data_dict2[keys[0]]["sigma"] - data_dict2[keys[1]]["sigma"]
    lower_val2 = (data_dict2[keys[0]]["sigma"] + data_dict2[keys[1]]["sigma"]) / 2.
    fraction_2 = upper_val2 / lower_val2
    range_2 = (fraction_2 > dict_labels2["deltasigmaoversigma"]["maxval"]) | (fraction_2 < dict_labels2["deltasigmaoversigma"]["minval"])
    

    ax[1][1].hist(fraction_  , color=dict_labels["color"], density=True, histtype="step", bins=dict_labels["deltasigmaoversigma"]["bins"],  label=f"{len(fraction_[range_])}")# {round(np.mean(fraction_[~range_]), 5)}")
    ax[1][1].hist(fraction_2, color=dict_labels2["color"], density=True, histtype="step", bins=dict_labels2["deltasigmaoversigma"]["bins"], label=f"{len(fraction_2[range_2])}")# {round(np.mean(fraction_2[~range_2]), 5)}")
    ax[1][1].axvline(0, color="k", ls="--")

    ax[1][1].set_xlim([dict_labels["deltasigmaoversigma"]["minval"], dict_labels["deltasigmaoversigma"]["maxval"]])
    ax[1][1].set_xticks(dict_labels["deltasigmaoversigma"]["ticks"])
    ax[1][1].set_xticklabels(dict_labels["deltasigmaoversigma"]["ticks"], rotation=30)


    # the lower diagonal plots
    upper_val =         data_dict[keys[0]]["alphamed"] - data_dict[keys[1]]["alphamed"]
    lower_val = np.sqrt(data_dict[keys[0]]["sigma"]**2 + data_dict[keys[1]]["sigma"]**2)
    fraction_ = upper_val / lower_val
    range_ = (fraction_ > dict_labels["deltaalphaoveravgsigma"]["maxval"]) | (fraction_ < dict_labels["deltaalphaoveravgsigma"]["minval"])
    
    upper_val2 =         data_dict2[keys[0]]["alphamed"] - data_dict2[keys[1]]["alphamed"]
    lower_val2 = np.sqrt(data_dict2[keys[0]]["sigma"]**2 + data_dict2[keys[1]]["sigma"]**2)
    fraction_2 = upper_val2 / lower_val2
    range_2 = (fraction_2 > dict_labels2["deltaalphaoveravgsigma"]["maxval"]) | (fraction_2 < dict_labels2["deltaalphaoveravgsigma"]["minval"])
    
    ax[1][0].hist(fraction_,  color=dict_labels["color"], density=True, histtype="step",  bins=dict_labels["deltaalphaoveravgsigma"]["bins"],  label=f"{len(fraction_[range_])}")# {round(np.mean(fraction_[~range_]), 5)}")
    ax[1][0].hist(fraction_2, color=dict_labels2["color"], density=True, histtype="step", bins=dict_labels2["deltaalphaoveravgsigma"]["bins"], label=f"{len(fraction_2[range_2])}")# {round(np.mean(fraction_2[~range_2]), 5)}")
    
    ax[1][0].axvline(0, color="k", ls="--")

    ax[1][0].set_xlim([dict_labels["deltaalphaoveravgsigma"]["minval"], dict_labels["deltaalphaoveravgsigma"]["maxval"]])
    ax[1][0].set_xticks(dict_labels["deltaalphaoveravgsigma"]["ticks"])
    ax[1][0].set_xticklabels(dict_labels["deltaalphaoveravgsigma"]["ticks"], rotation=30)
    

    ax[0][0].set_xlabel(r"$(\alpha_1 - 1)/\sigma_1$")
    ax[0][1].set_xlabel(r"$(\alpha_2 - 1)/\sigma_2$")
    ax[0][0].xaxis.set_label_position("top")
    ax[0][1].xaxis.set_label_position("top")

    ax[1][0].set_xlabel(r"$(\alpha_1 - \alpha_2)/(\sqrt{\sigma_1^2 + \sigma_2^2})$")
    ax[1][1].set_xlabel(r"$(\sigma_1 - \sigma_2)/[0.5\times(\sigma_1 + \sigma_2)]$")
    
    # ax[1][0].legend()
    # ax[1][1].legend()

    ax[0][0].tick_params(axis="x", direction="out", bottom=False, top=True, labelbottom=False, labeltop=True)
    ax[0][1].tick_params(axis="x", direction="out", bottom=False, top=True, labelbottom=False, labeltop=True)
    ax[0][1].tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)
    ax[1][1].tick_params(axis="x", direction="out", bottom=True, top=False, labelbottom=True, labeltop=False)

        
    fig.tight_layout()    
    fig.savefig(outpath + "/A_B_DAMP1_DAMP2.pdf", bbox_inches='tight')


    #### second plot
    fig, ax = pt.subplots(2, 2, figsize=(5, 4), gridspec_kw={"hspace":0, "wspace":0})

    ax[0][0].hist( 2 * (data_dict[keys[0]]["alphamed"] - data_dict2[keys[0]]["alphamed"]) / (data_dict[keys[0]]["alphamed"] + data_dict2[keys[0]]["alphamed"]), bins=np.linspace(-0.01, 0.01, 51), color="green", density=True, histtype="step", label="a=1: par-cg")
    ax[0][0].hist( 2 * (data_dict[keys[1]]["alphamed"] - data_dict2[keys[1]]["alphamed"]) / (data_dict[keys[1]]["alphamed"] + data_dict2[keys[1]]["alphamed"]), bins=np.linspace(-0.01, 0.01, 51), color="red", density=True, histtype="step", label="a=2: par-cg")

    ax[0][1].hist(data_dict[keys[0]]["sigma"] - data_dict2[keys[0]]["sigma"], bins=np.linspace(-0.002, 0.002, 51), color="green", density=True, histtype="step")
    ax[0][1].hist(data_dict[keys[1]]["sigma"] - data_dict2[keys[1]]["sigma"], bins=np.linspace(-0.002, 0.002, 51), color="red", density=True, histtype="step")
    
    ax[1][0].hist(data_dict[keys[0]]["alphamed"] - data_dict[keys[1]]["alphamed"], bins=np.linspace(-0.01, 0.01, 101), color=dict_labels[keys[0]]["color"], density=True, histtype="step", label="par")
    ax[1][0].hist(data_dict2[keys[0]]["alphamed"] - data_dict2[keys[1]]["alphamed"], bins=np.linspace(-0.01, 0.01, 101), color=dict_labels2[keys[0]]["color"], density=True, histtype="step", label="cg")
    ax[1][0].set_xlim([-0.007, 0.007])    

    ax[1][1].hist(data_dict[keys[0]]["sigma"] - data_dict[keys[1]]["sigma"], bins=np.linspace(-0.002, 0.002, 101), color=dict_labels[keys[0]]["color"], density=True, histtype="step")
    ax[1][1].hist(data_dict2[keys[0]]["sigma"] - data_dict2[keys[1]]["sigma"], bins=np.linspace(-0.002, 0.002, 101), color=dict_labels2[keys[0]]["color"], density=True, histtype="step")
    ax[1][1].set_xlim([-0.002, 0.002])    
    # ax[1][0].legend()
    
    fig.tight_layout()    
    fig.savefig(outpath + "/A_B_DAMP1_DAMP2_2.pdf", bbox_inches='tight')

def plot_comparison_UL(dict_labels_L, dict_labels_U, keys, figname="test", outliers=False):
    ### declare figures    
    from mycorner_lib import my_box_outliers, my_box
    if outliers:
        fig, ax = my_box_outliers(size=len(keys))
    else:
        fig, ax = my_box(rows=len(keys), cols=len(keys))
    
    # fig1, ax1 = my_box(rows=len(keys) - 1, cols=4, sharex="col", sharey="col", figsize=(5, len(keys)*4.4/4.))
    fig1, ax1 = my_box(rows=len(keys) - 1, cols=4, sharex="col", sharey="col", figsize=(5, 2))
    # fig1, ax1 = my_box(rows=len(keys) - 1, cols=4, sharex="col", sharey="col", figsize=(5, 4))


    ### Obtain data
    data_dict_L = {}
    data_dict_U = {}

    for i, key in enumerate(keys):
        dict_ = dict_labels_L[key]
        alpha_b, chi2, evi = read_chi2_alpha_evi(dict_["filename"], "", n=1)
        alpha_a, sigma_, alpha_m, sigma_nl = read_data_mystats(dict_["filename"], endcffile=dict_labels_L["endcffile"], begcffile=dict_labels_L["begcffile"])
        temp_dict_L = {"alpha":alpha_a, "alphabest":alpha_b, "alphamed":alpha_m, "sigma":sigma_, "chi2":chi2, "evi":evi, "sigma_nl":sigma_nl}
        data_dict_L[key] = temp_dict_L

        dict_ = dict_labels_U[key]
        alpha_b, chi2, evi = read_chi2_alpha_evi(dict_["filename"], "", n=1)
        alpha_a, sigma_, alpha_m, sigma_nl = read_data_mystats(dict_["filename"], endcffile=dict_labels_U["endcffile"], begcffile=dict_labels_U["begcffile"])
        temp_dict_U = {"alpha":alpha_a, "alphabest":alpha_b, "alphamed":alpha_m, "sigma":sigma_, "chi2":chi2, "evi":evi, "sigma_nl":sigma_nl}
        data_dict_U[key] = temp_dict_U

    ### Plot median values of alpha
    ## Lower triangular plot
    xmin_l, xmax_l = plot_my_corner_lower(ax, data_dict_L, dict_labels_L, keys, outliers=outliers)

    ## Upper triangular plot
    xmin_u, xmax_u = plot_my_corner_upper(ax, data_dict_U, dict_labels_U, keys, outliers=outliers)

    ax_style_my_corner_UL(ax, np.min(np.array([xmin_l, xmin_u])), np.max(np.array([xmax_l, xmax_u])), keys, outliers=outliers)

    line_u = Line2D([0], [0], marker='o', color='w', markerfacecolor=dict_labels_U["color"], markersize=5)
    line_l = Line2D([0], [0], marker='o', color='w', markerfacecolor=dict_labels_L["color"], markersize=5)

    outlier = 0
    if outliers:
        outlier = 1
    # ax[0][outlier + 1].legend([line_u], [dict_labels_U['label']],  fontsize="small", frameon=False, loc='upper right')
    # ax[1][outlier].legend([line_l], [dict_labels_L['label']],  fontsize="small", frameon=False, loc='upper right')
    ax[0][-1].legend([line_u, line_l], [dict_labels_U['label'], dict_labels_L['label']],  fontsize="small", bbox_to_anchor=(1, 1), frameon=False, loc='lower left')

    
    ### Plot pull, tension deltasigma and delta alpha
    xmin_L, xmax_L = plot_deltaalpha_deltasigma(ax1, data_dict_L, dict_labels_L, keys)
    xmin_U, xmax_U = plot_deltaalpha_deltasigma(ax1, data_dict_U, dict_labels_U, keys)
    line_l = Line2D([0], [0], marker='o', color='w', markerfacecolor=dict_labels_L["color"], markersize=5)
    line_u = Line2D([0], [0], marker='o', color='w', markerfacecolor=dict_labels_U["color"], markersize=5)
    # ax1[0][0].legend([line_l, line_u], [dict_labels_L["label"], dict_labels_U["label"]], loc="best", bbox_to_anchor=(0.5, 0.5), fontsize="small", frameon=False)

    ax1[0][0].set_xlim([1.01 * np.min([xmin_L[0], xmin_U[0]]), np.max([xmax_L[0], xmax_U[0]]) * 1.01])
    ax1[0][1].set_xlim([1.01 * np.min([xmin_L[1], xmin_U[1]]), np.max([xmax_L[1], xmax_U[1]]) * 1.01])
    ax1[0][2].set_xlim([1.01 * np.min([xmin_L[2], xmin_U[2]]), np.max([xmax_L[2], xmax_U[2]]) * 1.01])
    
    ###
    # fig2, ax2 = my_box(rows=len(keys), cols=len(keys), figsize=None)
    fig2, ax2 = my_box(rows=len(keys), cols=len(keys), figsize=(5, 4))
    plot_evi_tension(fig2, ax2, data_dict_U, dict_labels_U, keys)

    # fig3, ax3 = my_box(rows=len(keys), cols=len(keys), figsize=None)
    fig3, ax3 = my_box(rows=len(keys), cols=len(keys), figsize=(5, 4))
    plot_evi_tension(fig3, ax3, data_dict_L, dict_labels_L, keys)


    ###
    fig.tight_layout()
    fig.savefig(figname+"_corner.pdf", bbox_inches='tight')

    fig1.tight_layout()
    fig1.savefig(figname+"_corner1.pdf", bbox_inches='tight')

    fig2.tight_layout()
    fig2.savefig(figname + "_" + dict_labels_U["label"] + "_corner2.pdf", bbox_inches='tight')

    fig3.tight_layout()
    fig3.savefig(figname + "_" + dict_labels_L["label"] + "_corner3.pdf", bbox_inches='tight')


def gauss(x, mu, sig):
    return (1 / np.sqrt(2 * np.pi * sig**2) * np.exp(- (x - mu)**2 / (2 * sig**2)))


def plot_rec_comparison(dict_labels, keys, outliers=False):

    fig, ax = pt.subplots(figsize=(5, 4))

    data_dict = {}
    bins = np.linspace(-4.5, 4.5, 21)
    print(bins)
    colors = ["red", "blue"]
    for i, key in enumerate(keys):
        dict_ = dict_labels[key]
        
        # alpha_b, chi2, evi = read_chi2_alpha_evi(dict_["filename"], "", n=1)
        alpha_a, sigma_, alpha_m, _ = read_data_mystats(dict_["filename"], endcffile=dict_labels["endcffile"], begcffile=dict_labels["begcffile"])
        print(np.min((alpha_m - 1) / sigma_), np.max((alpha_m - 1) / sigma_))
        ax.hist((alpha_m - 1) / sigma_, histtype="step", color=colors[i], label=dict_["label"], density=True, bins=bins)
        ax.set_xlabel("$\\frac{\\alpha-1}{\\sigma_\\alpha}$", fontsize=14)
        # ax.axvline(0, ls="--", color="k")
    
    x_n = np.linspace(-4.5, 4.5, 1001)
    ax.plot(x_n, gauss(x_n, 0, 1), ls="--", color="k")
    # print(bins[:-1] + 0.09)
    ax.legend()
    fig.tight_layout()
    fig.savefig(dict_labels['filename']+".pdf", bbox_inches='tight')


def get_best_fit_3(filepath):
    s, cf, cf_nw = np.loadtxt(filepath, usecols=(0, 1, 2), unpack=True)
    range_ = (s >= 60) & (s <= 150)
    return s[range_], cf[range_], cf_nw[range_]

def best_fit_threepanels(outpath="test"):
    path = "/scratch/variu/phd_fitOut/patchy_cmass_subset/avg_PATCHY_BOX/"

    fig, ax = pt.subplots(3, 1, figsize=(5, 4), gridspec_kw={"hspace":0, "wspace":0}, sharex=True)

    sP, cf_PA_nw, std_nw = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/avg_PATCHY_BOX/avg_100_patchy_box_void_nw_pre.2pcf", usecols=(0, 1, 2), unpack=True)
    sP, cf_PA,    std    = np.loadtxt("/scratch/variu/phd_fitOut/patchy_cmass_subset/avg_PATCHY_BOX/avg_100_patchy_box_void_pre.2pcf"   , usecols=(0, 1, 2), unpack=True)

    spb, cf_p,  cf_p_nw  = get_best_fit_3(path + "/parab_60_150_fast/BAOfit_avg_100_patchy_box_void_pre.2pcf_best_nw.dat")  
    scb, cf_CG, cf_CG_nw = get_best_fit_3(path + "/stitched_16R_G2048_50_G512_2000_CG/BAOfit_avg_100_patchy_box_void_pre.2pcf_best_nw.dat")
    ssb, cf_SK, cf_SK_nw = get_best_fit_3(path + "/stitched_G2048_50_G512_2000/BAOfit_avg_100_patchy_box_void_pre.2pcf_best_nw.dat")

    ax[0].errorbar(sP, sP * sP * cf_PA, yerr=sP * sP * std, label="PATCHY", color="k", ls="", marker=".")
    ax[0].plot(spb, spb * spb * cf_p, label="PAR$_\\mathrm{{G}}$", color="red", ls="--")
    ax[0].plot(scb, scb * scb * cf_CG, label="CG$_\\mathrm{{B}}$", color="blue", ls=":")
    ax[0].plot(ssb, ssb * ssb * cf_SK, label="SK$_\\mathrm{{B}}$", color="orange", ls=":")

    ax[1].errorbar(sP, sP * sP * (cf_PA - cf_PA_nw), yerr=sP * sP * np.sqrt(std**2 + std_nw**2), label="PATCHY", color="k", ls="", marker=".")
    ax[1].plot(spb, spb * spb * (cf_p - cf_p_nw), label="PAR$_\\mathrm{{G}}$", color="red", ls="--")
    ax[1].plot(scb, scb * scb * (cf_CG - cf_CG_nw), label="CG$_\\mathrm{{B}}$", color="blue", ls=":")
    ax[1].plot(ssb, ssb * ssb * (cf_SK - cf_SK_nw), label="SK$_\\mathrm{{B}}$", color="orange", ls=":")

    ax[2].errorbar(sP, sP * sP * (cf_PA_nw), yerr=sP * sP * std_nw, label="PATCHY", color="k",  ls="", marker=".")
    ax[2].plot(spb, spb * spb * (cf_p_nw), label="PAR$_\\mathrm{{G}}$", color="red", ls="--")
    ax[2].plot(scb, scb * scb * (cf_CG_nw), label="CG$_\\mathrm{{B}}$", color="blue", ls=":")
    ax[2].plot(ssb, ssb * ssb * (cf_SK_nw), label="SK$_\\mathrm{{B}}$", color="orange", ls=":")

    ax[0].set_ylim([-13, 22])
    ax[1].set_ylim([-13, 22])
    ax[2].set_ylim([-13, 22])

    ax[0].set_xlim([50, 160])

    ax[2].set_xlabel(r"s[$h^{-1}$Mpc]")
    
    ax[0].set_ylabel(r"s$^2\xi(s)$  ")
    ax[1].set_ylabel(r"s$^2\xi(s)^\mathrm{BAO}$")
    ax[2].set_ylabel(r"s$^2\xi(s)^\mathrm{nw}$")

    ax[2].legend(ncol=2)
    
    fig.tight_layout()
    fig.savefig(outpath + "A_B_CG_SK_PAR_best_fit.pdf", bbox_inches='tight')

def best_fit_onepanel(outpath="test"):

    data_file = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vhxcf/16R/avg_16R_vhxcf.xcf"
    cgfile    = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/avg_fit/stitched_16R_G2048_50_G512_2000_CG_60_150/BAOfit_avg_16R_vhxcf.xcf_best_nw.dat"
    skfile    = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/avg_fit/range_test/avg_2001_60_150/BAOfit_avg_16R_vhxcf.xcf_best_nw.dat"
    best_fit_onepanel_aux(data_file=data_file, cgfile=cgfile, skfile=skfile, ylim=[-25, 18], xlim=[50, 160], outfile=outpath + "X_B_CG_SK_best_fit.pdf")
    
    data_file = "/scratch/variu/clustering/patchy_cmass_subset/box1/redshift_recon/vv2pcf/avg_16R_500_recon.2pcf"
    cgfile    = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/redshift_recon/vv2pcf_CG/avg_fit/stitched_16R_G2048_50_G512_2000_CG_60_150/BAOfit_avg_16R_500_recon.2pcf_best_nw.dat"
    skfile    = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/redshift_recon/vv2pcf/avg_fit/stitched_G2048_50_G512_2000_60_150/BAOfit_avg_16R_500_recon.2pcf_best_nw.dat"
    best_fit_onepanel_aux(data_file=data_file, cgfile=cgfile, skfile=skfile, ylim=[-17, 35], xlim=[50, 160], outfile=outpath + "A_B_CG_SK_best_fit_recon.pdf")
    

def best_fit_onepanel_aux(data_file="test", cgfile="test", skfile="test", ylim=None, xlim=None, outfile="test"):
    fig, ax = pt.subplots(figsize=(5, 4))
    
    sa, cfa, stda = np.loadtxt(data_file, usecols=(0,1,2), unpack=True)
    sb, cfb, cfb_nw = np.loadtxt(cgfile, usecols=(0,1,2), unpack=True)
    sk, cfk, cfk_nw = np.loadtxt(skfile, usecols=(0,1,2), unpack=True)
    
    ax.errorbar(sa, sa*sa*cfa, yerr=sa*sa*stda, label="PATCHY", color="k", ls="", marker=".")
    
    ax.plot(sb, sb*sb*cfb, label="CG$_\\mathrm{{B}}$", color="blue", ls="--")
    ax.plot(sb, sb*sb*cfb_nw, color="blue", ls=":")
    
    ax.plot(sk, sk*sk*cfk, label="SK$_\\mathrm{{B}}$", color="orange", ls="--")
    ax.plot(sk, sk*sk*cfk_nw, color="orange", ls=":")
    
    ax.legend()
    
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    ax.set_xlabel(r"s[$h^{-1}$Mpc]")    
    ax.set_ylabel(r"s$^2\xi(s)$  ")
    
    fig.tight_layout()
    fig.savefig(outfile, bbox_inches='tight')


def plot_template_aux(ax, f_lin, path="test", label="label", marker=".", markersize=3, markerfacecolor=(1, 0, 0, 0.3), markeredgecolor=(1, 0, 0, 1), A=1.):
    k, pk = np.loadtxt(path, usecols=(0, 1), unpack=True)
    
    range_ = slice(0, k.size, 5)
    k, pk = k[range_], pk[range_]
    ax.plot(k, A * pk / f_lin(k), ls="", label=label, marker=marker, markersize=markersize, markerfacecolor=markerfacecolor, markeredgecolor=markeredgecolor)
    # ax[1].plot(k, A * k * pk, ls="", marker=marker, markersize=markersize, markerfacecolor=markerfacecolor, markeredgecolor=markeredgecolor)
    
def poly2(x, a, b):
    return a + b * x**2

def main_template(outpath="test"):

    dict_1 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.pspec",
            "label":r"CG$_\mathrm{B,vv}$",
            "facecolor":(1, 0, 0, 0.3),
            "edgecolor":(1, 0, 0, 1),
            "A":1.15,
            "marker":"o",
            "markersize":3
            }

    dict_2 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.xpspec",
            "label":r"CG$_\mathrm{B,gv}$",
            "facecolor":(1, 0, 1, 0.3),
            "edgecolor":(1, 0, 1, 1),
            "A":1.,
            "marker":"o",
            "markersize":3
            }

    dict_3 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_voids/stitched_G2048_50_G512_2000.txt",
            "label":r"SK$_\mathrm{B,vv}$",
            "facecolor":(0, 0, 1, 0.3),
            "edgecolor":(0, 0, 1, 1),
            "A":1.05,
            "marker":"s",
            "markersize":3
            }

    dict_4 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_cross/avg_16R_2001.pspec",
            "label":r"SK$_\mathrm{B,gv}$",
            "facecolor":(1, 1, 0, 0.3),
            "edgecolor":(1, 0.7, 0, 1),
            "A":1.,
            "marker":"s",
            "markersize":3
            }
    
    # plot_template(dict_1, dict_2, dict_3, dict_4, outpath=outpath, figurename="correct_template.pdf")


    ######
    dict_1 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_80.pspec",
            "label":r"CG$_\mathrm{B,80,vv}$",
            "facecolor":(0.1, 1, 0.1, 0.3),
            "edgecolor":(0.1, 1, 0.1, 1),
            "A":1.25,
            "marker":"o",
            "markersize":3
            }

    dict_2 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_80.xpspec",
            "label":r"CG$_\mathrm{B,80,gv}$",
            "facecolor":(0.1, 1, 0.1, 0.3),
            "edgecolor":(0.1, 1, 0.1, 1),
            "A":1.,
            "marker":"s",
            "markersize":3
            }

    dict_3 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_120.pspec",
            "label":r"CG$_\mathrm{B,120,vv}$",
            "facecolor":(160/255., 32/255., 240/255., 0.3),
            "edgecolor":(160/255., 32/255., 240/255., 1),
            "A":1.05,
            "marker":"o",
            "markersize":3
            }

    dict_4 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_120.xpspec",
            "label":r"CG$_\mathrm{B,120,gv}$",
            "facecolor":(160/255., 32/255., 240/255., 0.3),
            "edgecolor":(160/255., 32/255., 240/255., 1),
            "A":1.,
            "marker":"s",
            "markersize":3
            }

    # plot_template(dict_1, dict_2, dict_3, dict_4, outpath=outpath, figurename="wrong_80_120_template.pdf")


    ###### 
    dict_1 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_imp.pspec",
            "label":r"CG$_\mathrm{B,def,vv}$",
            "facecolor":(1, 0, 0, 0.3),
            "edgecolor":(1, 0, 0, 1),
            "A":1.,
            "marker":"o",
            "markersize":3
            }

    dict_2 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG_imp.xpspec",
            "label":r"CG$_\mathrm{B,def,gv}$",
            "facecolor":(1, 0, 1, 0.3),
            "edgecolor":(1, 0, 1, 1),
            "A":1.,
            "marker":"o",
            "markersize":3
            }

    dict_3 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024FAC0.3/stitched_G2048_50_G512_2000_FAC0.3_G1024.txt",
            "label":r"SK$_\mathrm{B,def,vv}$",
            "facecolor":(0, 0, 1, 0.3),
            "edgecolor":(0, 0, 1, 1),
            "A":1.3,
            "marker":"s",
            "markersize":3
            }

    dict_4 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024FAC0.3/stitched_G2048_50_G512_2000_FAC0.3_G1024_xpspec.txt",
            "label":r"SK$_\mathrm{B,def,gv}$",
            "facecolor":(1, 1, 0, 0.3),
            "edgecolor":(1, 0.7, 0, 1),
            "A":1.,
            "marker":"s",
            "markersize":3
            }
    
    # plot_template(dict_1, dict_2, dict_3, dict_4, outpath=outpath, figurename="wrong_imp_template.pdf")


    ###### 
    dict_1 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.pspec",
            "label":r"CG$_\mathrm{B,vv}$",
            "facecolor":(1, 0, 0, 0.3),
            "edgecolor":(1, 0, 0, 1),
            "A":0.6,
            "marker":"o",
            "markersize":3
            }

    dict_2 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/lightcone_box1_large_v/stitched_16R_G2048_50_G512_20_CG_large_v.pspec",
            "label":r"CG$_\mathrm{LC,vv}$",
            "facecolor":(0.1, 1, 0.1, 0.3),
            "edgecolor":(0.1, 1, 0.1, 1),
            "A":1.03,
            "marker":"o",
            "markersize":3
            }

    dict_3 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_voids/stitched_G2048_50_G512_2000.txt",
            "label":r"SK$_\mathrm{B,vv}$",
            "facecolor":(0, 0, 1, 0.3),
            "edgecolor":(0, 0, 1, 1),
            "A":0.55,
            "marker":"s",
            "markersize":3
            }

    dict_4 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/lightcone_box1/G1024/16R/large_v/stitched_16R_G2048_50_G512_20_CIC_large_v.pspec",
            "label":r"SK$_\mathrm{LC,vv}$",
            "facecolor":(1, 1, 0, 0.3),
            "edgecolor":(1, 0.7, 0, 1),
            "A":1.03,
            "marker":"s",
            "markersize":3
            }
    
    # plot_template(dict_1, dict_2, dict_3, dict_4, outpath=outpath, figurename="LC_template.pdf", box=False)


    ###### 
    dict_1 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.pspec",
            "label":r"CG$_\mathrm{B,vv}$",
            "facecolor":(1, 0, 0, 0.3),
            "edgecolor":(1, 0, 0, 1),
            "A":1.45,
            "marker":"o",
            "markersize":3
            }

    # dict_2 = {"filename": "/scratch/variu/clustering/linearnw_COSMOGAME/patchy_cmass_like/fix_ampl/box1/avg/stitched_16R_G2048_50_G512_2000_CG.xpspec",
    #         "label":r"CG$_\mathrm{B,gv}$",
    #         "facecolor":(1, 0, 1, 0.3),
    #         "edgecolor":(1, 0, 1, 1),
    #         "A":1.,
    #         "marker":"o",
    #         "markersize":3
    #         }

    dict_3 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_voids/stitched_G2048_50_G512_2000.txt",
            "label":r"SK$_\mathrm{B,vv}$",
            "facecolor":(0, 0, 1, 0.3),
            "edgecolor":(0, 0, 1, 1),
            "A":1.30,
            "marker":"s",
            "markersize":3
            }

    # dict_4 = {"filename": "/scratch/variu/clustering/linearnw_CIC/patchy_cmass_like/fix_ampl/box1/G1024/16R_1000bins/avg_cross/avg_16R_2001.pspec",
    #         "label":r"SK$_\mathrm{B,gv}$",
    #         "facecolor":(1, 1, 0, 0.3),
    #         "edgecolor":(1, 0.7, 0, 1),
    #         "A":1.,
    #         "marker":"s",
    #         "markersize":3
            # }
    
    # plot_template(dict_1, dict_2, dict_3, dict_4, outpath=outpath, figurename="correct_template_post_recon_data.pdf", box=False)
    plot_template(dict_1, None, dict_3, None, outpath=outpath, figurename="correct_template_post_recon_data.pdf", box=False, recon=True)


def plot_template(dict_1, dict_2, dict_3, dict_4, outpath="test", figurename="figname", box=True, recon=False):
    from scipy.interpolate import interp1d

    linear = "/home/astro/variu/phd/voids/chengscodes/CosmoGAME/EisensteinHu_Pnw.txt"
    k_lin_nw, pk_linw_nw = np.loadtxt(linear, usecols=(0, 1), unpack=True)
    fint_lin_nw = interp1d(np.log(k_lin_nw), np.log(pk_linw_nw), kind="cubic")
    f_lin_nw = lambda k: np.exp(fint_lin_nw(np.log(k)))

    fig, ax = pt.subplots(figsize=(5, 4), sharex=True, gridspec_kw={"wspace":0.4})

    if box:
        reference_a = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/pkvoids_g/avg_500.pspec"
    else:
        reference_a = "/scratch/variu/clustering/patchy_cmass_subset/lightcone_box1/real_large_V_pspec/avg_80.pspec"
        if recon:
            reference_a = "/scratch/variu/clustering/patchy_cmass_subset/box1/redshift_recon/pkvoids/avg_495_recon.pspec"
    
    k_ra, pk_ra = np.loadtxt(reference_a, usecols=(0, 1), unpack=True)
    print(np.max(k_ra))

    
    if box:
        reference_x = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/pkcross_g/avg_500.xpspec"
        k_rx, pk_rx = np.loadtxt(reference_x, usecols=(0, 1), unpack=True)
        print(np.max(k_rx))

       
        
    ######################################## Parabolic START
    # from scipy.optimize import curve_fit
    # range_ = k_ra < 0.1

    # for box
    # c = 1064
    # A = 0.8

    # for LC
    # c = 1495.20356669250032
    # A = 0.4

    # ax.axhline(( pk_ra / f_lin_nw(k_ra))[0], color="grey", ls="--", label="DW")

    # ax.plot(k_ra, A * (1 + c * k_ra ** 2), color="green", ls="--", label="PAR")

    # popt, pcov = curve_fit(poly2, k_ra[range_], (pk_ra / f_lin_nw(k_ra))[range_])
    # print(popt)
    # ax.plot(k_ra, poly2(k_ra, *popt), color="grey", ls="--")

    # c = 219
    # A = 1
    # ax.axhline((pk_rx / f_lin_nw(k_rx))[0], color="grey", ls="--", label="DW")

    # ax.plot(k_rx, - A * (1 + c * k_rx ** 2), color="green", ls="--", label="PAR")
    # popt, pcov = curve_fit(poly2, k_rx[range_], (pk_rx / f_lin_nw(k_rx))[range_])
    # print(popt)
    # ax.plot(k_rx, poly2(k_rx, *popt), color="grey", ls="--")
    ######################################## Parabolic END

    
    plot_template_aux(ax, f_lin_nw, path=dict_1["filename"], label=dict_1["label"], marker=dict_1["marker"], markersize=dict_1["markersize"], markerfacecolor=dict_1["facecolor"], markeredgecolor=dict_1["edgecolor"], A=dict_1["A"])
    # plot_template_aux(ax, f_lin_nw, path=dict_2["filename"], label=dict_2["label"], marker=dict_2["marker"], markersize=dict_2["markersize"], markerfacecolor=dict_2["facecolor"], markeredgecolor=dict_2["edgecolor"], A=dict_2["A"])
    plot_template_aux(ax, f_lin_nw, path=dict_3["filename"], label=dict_3["label"], marker=dict_3["marker"], markersize=dict_3["markersize"], markerfacecolor=dict_3["facecolor"], markeredgecolor=dict_3["edgecolor"], A=dict_3["A"])
    # plot_template_aux(ax, f_lin_nw, path=dict_4["filename"], label=dict_4["label"], marker=dict_4["marker"], markersize=dict_4["markersize"], markerfacecolor=dict_4["facecolor"], markeredgecolor=dict_4["edgecolor"], A=dict_4["A"])
  

    if recon:
        ax.plot(k_ra, pk_ra / f_lin_nw(k_ra), label=r"recon mock $P_\mathrm{vv}$", color="k")
    else:
        ax.plot(k_ra, pk_ra / f_lin_nw(k_ra), label=r"mock $P_\mathrm{vv}$", color="k")

    # ax[1].plot(k_ra, k_ra * pk_ra, color="k")

    if box:
        ax.plot(k_rx, pk_rx / f_lin_nw(k_rx), label=r"mock $P_\mathrm{gv}$", color="grey")
        # ax[1].plot(k_rx, k_rx * pk_rx, color="grey")
    
    if box:
        ax.set_xlim([0, 0.60])
        ax.set_ylim([-4.0, 15])
        # ax[1].set_yticks([-1000, 0, 1000, 2000, 3000])
        # ax[1].set_yticklabels([-1, 0, 1, 2, 3])
    else:
        ax.set_xlim([0, 0.45])
        ax.set_ylim([-0.1, 6.0])
        if recon:
            ax.set_xlim([0, 0.60])
            ax.set_ylim([-0.1, 20])
        

        # ax[1].set_yticks([0, 400, 800, 1200, 1600])
        # ax[1].set_yticklabels([0, 4, 8, 12, 16])


    ax.set_xlabel(r"$k[h~\mathrm{Mpc}^{-1}]$")
    ax.set_ylabel(r"$P_\mathrm{mock/template}(k)/P_\mathrm{lin,nw}(k)$")
    
    # ax[1].set_xlabel(r"$k[h~\mathrm{Mpc}^{-1}]$")
    # ax[1].set_ylabel(r"$kP_\mathrm{mock/template}(k)$")
    
    ax.set_yticks([0])
    ax.set_yticklabels([0])
        
    ax.axhline(0, ls=":", color="grey")
    ax.legend(bbox_to_anchor=(0, 1.0), loc='lower left', ncol=3)

    fig.tight_layout()
    fig.savefig(outpath + figurename, bbox_inches='tight')

def main():
    pt.rcParams.update({'font.family': "serif"})
    pt.rcParams.update({'font.serif': "Times New Roman"})

    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['mathtext.rm'] = 'serif'
    # print(pt.rcParams.keys())
    
    outpath = "/home/astro/variu/void_project_figs/"
    # outpath = "/home/astro/variu/temp/"
     
    # main_template(outpath=outpath)
    # exit()
    # best_fit_threepanels(outpath=outpath)
    # best_fit_onepanel(outpath=outpath)
    # exit()
    ## BOX vv2pcf
    # inpath2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/parab_60_150_m_fast/"
    # inpath2c = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/parab_60_150_m_fast_fixc/"
    # inpath3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000/"
    
    # inpath20 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000E/"
    # inpath21 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000F/"
    # inpath22 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_500/"
    # inpath23 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_500F/"
    # inpath24 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000_old/"
    # inpath241 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_100F/"

    # inpath30 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/stitched_subvolumes/"
    # inpath31 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/lin_sav_gol_71_5_stitched_subvolumes/"

    # inpath33 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/sm_stitched_G2048_50_G512_2000/"
    
    # inpath34 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC_60_150_m_2000FS/"
    # inpath36 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024CIC18R_60_150_m_2000F/"

    inpath_gal_vv =        "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/galaxy_60_150_m/"
    inpath_SK_vv =         "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/stitched_G2048_50_G512_2000/"
    inpath_SK_imp_vv =     "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/G1024F0.3CIC_60_150_m_2000F/"
    inpath_parab_vv =      "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/parab_60_150_gauss/"
    inpath_parab_fixc_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vv2pcf/parab_60_150_fixc/"
    
    inpath_CG_vv =     "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG/"
    inpath_CG_imp_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG_imp/"
    inpath_CG_80_vv =  "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG_80/"
    inpath_CG_120_vv = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG_120/"
   
    ## recon BOX vv2pcf
    inpath104 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/redshift_recon/vv2pcf_CG/stitched_G2048_50_G512_2000/"
    inpath105 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/redshift_recon/vv2pcf/stitched_G2048_50_G512_2000/"
 
    ## BOX vhxcf
    # inpath26 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/parab_60_150_m_fast2F/"
    # inpath26_c = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/parab_60_150_m_fast_largec/"
    # inpath26_fixc = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/parab_60_150_m_fast_fixc/"
    # inpath38 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024CIC18R_60_150_m_2000F/"

    inpath_gal_vh =        "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/galaxy_60_150_m/"
    inpath_SK_vh =         "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024CIC_60_150_m_2001F/"
    inpath_SK_imp_vh =     "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf/G1024F0.3CIC_60_150_m_2000F/"
    inpath_parab_vh =      "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vhxcf/parab_60_150_gauss/"
    inpath_parab_fixc_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/parab/vhxcf/parab_60_150_fixc/"
    
    inpath_CG_vh =     "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/stitched_16R_G2048_50_G512_2000_CG/"
    inpath_CG_imp_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/stitched_16R_G2048_50_G512_2000_CG_imp/"
    inpath_CG_80_vh =  "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/stitched_16R_G2048_50_G512_2000_CG_80/"
    inpath_CG_120_vh = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vhxcf_CG/stitched_16R_G2048_50_G512_2000_CG_120/"
   
    ## LC vv2pcf
    # inpath1000 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG/"
    # inpath1001 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG_LC/"
    # inpath1002 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_CG/galaxy_60_150_m/"
    # inpath1003 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_CG/parab_60_150_m/"
    # inpath1004 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf/avg_16R_2000.pspec/"
    # inpath1005 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf/stitched_G2048_50_G512_2000/"


    # obtain_chi2_alpha_evi(inpath_CG_80_vh, endcffile=".dat.dr.2pcf")
    # obtain_chi2_alpha_evi(inpath1002, endcffile=".dat.2pcf")
    # obtain_chi2_alpha_evi(inpath2c, endcffile=".VOID.dat.2pcf")
    # exit()

    
    ### Plot 1: AX_B_CG_SK_GAL_fixc
    # dict_labels_A = {
    #     "first":{"label":"fix c", "alphabest":r"$\alpha_\mathdefault{best, fix c}$", "alphamed":"\\alpha_\\mathdefault{{fix~c}}", "alpha":r"$\alpha_\mathdefault{fix c}$", "sigma":"\\sigma_\\mathdefault{fix c}", "filename":inpath_parab_fixc_vv, "color":"red"},
    #     "second":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, CG}$", "alphamed":"\\alpha_{{\\mathdefault{{CG}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{CG}$", "sigma":"\\sigma_\\mathdefault{{CG}}", "filename":inpath_CG_vv, "color":"green"},
    #     "third":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, SK}$", "alphamed":"\\alpha_{{\\mathdefault{{SK}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{SK}$", "sigma":"\\sigma_\\mathdefault{{SK}}", "filename":inpath_SK_vv,"color":"blue"},
    #     "forth":{"label":"DW", "alphabest":r"$\alpha_\mathdefault{best, GAL}$", "alphamed":"\\alpha_\\mathdefault{{DW}}", "alpha":r"$\alpha_\mathdefault{GAL}$", "sigma":"\\sigma_\\mathdefault{{GAL}}", "filename":inpath_gal_vv, "color":"magenta"},
    #     "filename":outpath + "A_B_CG_SK_GAL_fixc",
    #     "label":"A",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.98, 1.0, 1.02], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-15, "maxval":15, "ticks":[-10, -5, 0, 5, 10], "bins":np.linspace(-15,15, 51)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-2, "maxval":0.8, "ticks":[-1.6, -1.2, -0.9, -0.6, -0.3, 0, 0.3], "bins":np.linspace(-2, 2, 81)},
    #     "deltaalphaoveralpha":{"minval":-0.01, "maxval":0.02, "ticks":[-0.004, 0, 0.004, 0.008, 0.012], "bins":np.linspace(-0.05, 0.05, 201)},
    #     "deltasigmaoversigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3], "bins":np.linspace(-2, 2, 201)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":1.5, "ticks":[-0.4, 0, 0.4, 0.8, 1.2], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }
    # dict_labels_X = {
    #     "first":{"label":"fix c", "alphabest":r"$\alpha_\mathdefault{best, fix c}$", "alphamed":"\\alpha_\\mathdefault{{fix~c}}", "alpha":r"$\alpha_\mathdefault{fix c}$", "sigma":"\\sigma_\\mathdefault{fix c}", "filename":inpath_parab_fixc_vh, "color":"cyan"},
    #     "second":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, CG}$", "alphamed":"\\alpha_{{\\mathdefault{{CG}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{CG}$", "sigma":"\\sigma_\\mathdefault{{CG}}", "filename":inpath_CG_vh, "color":"green"},
    #     "third":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, SK}$", "alphamed":"\\alpha_{{\\mathdefault{{SK}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{SK}$", "sigma":"\\sigma_\\mathdefault{{SK}}", "filename":inpath_SK_vh,"color":"blue"},
    #     "forth":{"label":"DW", "alphabest":r"$\alpha_\mathdefault{best, GAL}$", "alphamed":"\\alpha_\\mathdefault{{DW}}", "alpha":r"$\alpha_\mathdefault{GAL}$", "sigma":"\\sigma_\\mathdefault{{GAL}}", "filename":inpath_gal_vh, "color":"magenta"},
    #     "filename":outpath + "X_B_CG_SK_GAL_fixc",
    #     "label":"X",
    #     "color":"red",
    #     "alphamed":{"ticks":[0.98, 1.0, 1.02], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-15, "maxval":15, "ticks":[-10, -5, 0, 5, 10], "bins":np.linspace(-15,15, 51)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-2, "maxval":0.8, "ticks":[-1.6, -1.2, -0.9, -0.6, -0.3, 0, 0.3], "bins":np.linspace(-2, 2, 81)},
    #     "deltaalphaoveralpha":{"minval":-0.01, "maxval":0.02, "ticks":[-0.004, 0, 0.004, 0.008, 0.012], "bins":np.linspace(-0.05, 0.05, 201)},
    #     "deltasigmaoversigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3], "bins":np.linspace(-2, 2, 201)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":1.5, "ticks":[-0.4, 0, 0.4, 0.8, 1.2], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.dr.2pcf",
    #     "begcffile":""
    #     }
    # plot_comparison_UL(dict_labels_A, dict_labels_X, keys=["first","second","third", "forth"], figname=outpath + "/AX_B_CG_SK_GAL_fixc", outliers=False)
    

    ### Plot 2: AX_B_CG_SK_GAL_PAR_fixc
    # dict_labels_A = {
    #     "first":{"label":"PAR$_\\mathdefault{{G}}$", "alphabest":r"$\alpha_\mathdefault{best, PAR}$", "alphamed":"\\alpha_{{\\mathdefault{{PAR}}_\\mathrm{{G}} }}", "alpha":r"$\alpha_\mathdefault{PAR}$", "sigma":"\\sigma_\\mathdefault{{PAR}}", "filename":inpath_parab_vv, "color":"orange"},
    #     "second":{"label":"fix c", "alphabest":r"$\alpha_\mathdefault{best, fix c}$", "alphamed":"\\alpha_\\mathdefault{{fix~c}}", "alpha":r"$\alpha_\mathdefault{fix c}$", "sigma":"\\sigma_\\mathdefault{{fix c}}", "filename":inpath_parab_fixc_vv, "color":"red"},
    #     "third":{"label":"CG$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, CG}$", "alphamed":"\\alpha_{{\\mathdefault{{CG}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{CG}$", "sigma":"\\sigma_\\mathdefault{{CG}}", "filename":inpath_CG_vv, "color":"green"},
    #     "forth":{"label":"SK$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, SK}$", "alphamed":"\\alpha_{{\\mathdefault{{SK}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{SK}$", "sigma":"\\sigma_\\mathdefault{{SK}}", "filename":inpath_SK_vv,"color":"blue"},
    #     "fifth":{"label":"DW", "alphabest":r"$\alpha_\mathdefault{best, GAL}$", "alphamed":"\\alpha_\\mathdefault{{DW}}", "alpha":r"$\alpha_\mathdefault{GAL}$", "sigma":"\\sigma_\\mathdefault{{GAL}}", "filename":inpath_gal_vv, "color":"magenta"},
    #     "filename":outpath + "A_B_CG_SK_GAL_PAR_fixc",
    #     "label":"A",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.975, 1.0, 1.025], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-15, "maxval":15, "ticks":[-10, -5, 0, 5, 10], "bins":np.linspace(-15,15, 51)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-2, "maxval":0.8, "ticks":[-1.6, -1.2, -0.8, -0.4, 0, 0.4], "bins":np.linspace(-2, 2, 51)},
    #     "deltaalphaoveralpha":{"minval":-0.01, "maxval":0.02, "ticks":[-0.004, 0, 0.004, 0.008], "bins":np.linspace(-0.05, 0.05, 201)},
    #     "deltasigmaoversigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.3, -0.15, -0.1, 0, 0.15], "bins":np.linspace(-2, 2, 201)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":1.5, "ticks":[-0.4, 0, 0.4, 0.8, 1.2], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }
    # dict_labels_X = {
    #     "first":{"label":"PAR$_\\mathdefault{{G}}$", "alphabest":r"$\alpha_\mathdefault{best, PAR}$", "alphamed":"\\alpha_{{\\mathdefault{{PAR}}_\\mathrm{{G}} }}", "alpha":r"$\alpha_\mathdefault{PAR}$", "sigma":"\\sigma_\\mathdefault{{PAR}}", "filename":inpath_parab_vh, "color":"orange"},
    #     "second":{"label":"fix c", "alphabest":r"$\alpha_\mathdefault{best, fix c}$", "alphamed":"\\alpha_\\mathdefault{{fix~c}}", "alpha":r"$\alpha_\mathdefault{fix c}$", "sigma":"\\sigma_\\mathdefault{{fix c}}", "filename":inpath_parab_fixc_vh, "color":"cyan"},
    #     "third":{"label":"CG$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, CG}$", "alphamed":"\\alpha_{{\\mathdefault{{CG}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{CG}$", "sigma":"\\sigma_\\mathdefault{{CG}}", "filename":inpath_CG_vh, "color":"green"},
    #     "forth":{"label":"SK$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathdefault{best, SK}$", "alphamed":"\\alpha_{{\\mathdefault{{SK}}_\\mathrm{{B}}}}", "alpha":r"$\alpha_\mathdefault{SK}$", "sigma":"\\sigma_\\mathdefault{{SK}}", "filename":inpath_SK_vh,"color":"blue"},
    #     "fifth":{"label":"DW", "alphabest":r"$\alpha_\mathdefault{best, GAL}$", "alphamed":"\\alpha_\\mathdefault{{DW}}", "alpha":r"$\alpha_\mathdefault{GAL}$", "sigma":"\\sigma_\\mathdefault{{GAL}}", "filename":inpath_gal_vh, "color":"magenta"},
    #     "filename":outpath + "X_B_CG_SK_GAL_PAR_fixc",
    #     "label":"X",
    #     "color":"red",
    #     "alphamed":{"ticks":[0.98, 1.0, 1.02], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-15, "maxval":15, "ticks":[-10, -5, 0, 5, 10], "bins":np.linspace(-15,15, 51)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-2, "maxval":0.8, "ticks":[-1.6, -1.2, -0.8, -0.4, 0, 0.4], "bins":np.linspace(-2, 2, 51)},
    #     "deltaalphaoveralpha":{"minval":-0.01, "maxval":0.02, "ticks":[-0.004, 0, 0.004, 0.008], "bins":np.linspace(-0.05, 0.05, 201)},
    #     "deltasigmaoversigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.3, -0.15, -0.1, 0, 0.15], "bins":np.linspace(-2, 2, 201)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":1.5, "ticks":[-0.4, 0, 0.4, 0.8, 1.2], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.dr.2pcf",
    #     "begcffile":""
    #     }
    # plot_comparison_UL(dict_labels_A, dict_labels_X, keys=["first","second","third", "forth", "fifth"], figname=outpath + "/AX_B_CG_SK_GAL_PAR_fixc", outliers=True)
    # plot_comparison(dict_labels_A, keys=["first", "second", "third","forth", "fifth"], outliers=True)
    # plot_comparison(dict_labels_X, keys=["first", "second", "third","forth", "fifth"], outliers=False)
    
    
    #### Plot 3: AX_B_CG_CG80_CGimp_CG120
    # dict_labels_A = {
    #     "first":{"label":"CG$_\\mathdefault{{B}}$",                  "alphabest":r"$\alpha_\mathrm{best, CG}$",  "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{B}} }}","filename":inpath_CG_vv, "color":"green"},
    #     "second":{"label":r"CG$_\mathrm{def}$", "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{def}}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{def}} }}", "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{def}}$", "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{def}} }}", "filename":inpath_CG_imp_vv,"color":"blue"},
    #     "third":{"label":r"CG$_\mathrm{80}$",   "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{80}}$",  "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{80}} }}",  "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{80}}$",  "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{80}} }}",  "filename":inpath_CG_80_vv, "color":"magenta"},
    #     "forth":{"label":r"CG$_\mathrm{120}$",  "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{120}}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{120}} }}", "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{120}}$", "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{120}} }}", "filename":inpath_CG_120_vv, "color":"orange"},
    #     "filename":outpath + "A_B_CG_CG80_CGimp_CG120",
    #     "label":"A",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.97, 1.0, 1.03], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 41)},
    #     "tensionparams":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.2, 0, 0.2, 0.4], "bins":np.linspace(-1, 1, 81)},
    #     "deltaalphaoveralpha":{"minval":-0.005, "maxval":0.005, "ticks":[-0.004, -0.002, 0, 0.002, 0.004], "bins":np.linspace(-0.05, 0.05, 401)},
    #     "deltasigmaoversigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.1, -0.05, 0, 0.05, 0.1], "bins":np.linspace(-2, 2, 251)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }

    # dict_labels_X = {
    #     "first":{"label":"CG$_\\mathdefault{{B}}$",                  "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}",                      "alpha":r"$\alpha_\mathrm{CG}$",                "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{B}} }}",                      "filename":inpath_CG_vh, "color":"green"},
    #     "second":{"label":r"CG$_\mathrm{def}$", "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{def}}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{def}} }}", "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{def}}$", "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{def}} }}", "filename":inpath_CG_imp_vh,"color":"blue"},
    #     "third":{"label":r"CG$_\mathrm{80}$",   "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{80}}$",  "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{80}} }}",  "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{80}}$",  "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{80}} }}",  "filename":inpath_CG_80_vh, "color":"magenta"},
    #     "forth":{"label":r"CG$_\mathrm{120}$",  "alphabest":r"$\alpha_{\mathrm{best, CG}_\mathrm{120}}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{120}} }}", "alpha":r"$\alpha_{\mathrm{CG}_\mathrm{120}}$", "sigma":"\\sigma_{{\\mathrm{{CG}}_\\mathrm{{120}} }}", "filename":inpath_CG_120_vh, "color":"orange"},
    #     "filename":outpath + "X_B_CG_CG80_CGimp_CG120",
    #     "label":"X",
    #     "color":"red",
    #     "alphamed":{"ticks":[0.97, 1.0, 1.03], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 41)},
    #     "tensionparams":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.2, 0, 0.2, 0.4], "bins":np.linspace(-1, 1, 81)},
    #     "deltaalphaoveralpha":{"minval":-0.005, "maxval":0.005, "ticks":[-0.004, -0.002, 0, 0.002, 0.004], "bins":np.linspace(-0.05, 0.05, 401)},
    #     "deltasigmaoversigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.1, -0.05, 0, 0.05, 0.1], "bins":np.linspace(-2, 2, 251)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.dr.2pcf",
    #     "begcffile":""
    #     }
    # plot_comparison_UL(dict_labels_A, dict_labels_X, keys=["first","second","third", "forth"], figname=outpath + "/AX_B_CG_CG80_CGimp_CG120")

    
    ### Plot 4: AX_B_SK_SK_imp
    # dict_labels_A = {
    #     "first":{"label":"SK$_\mathrm{{B}}$",                  "alphabest":r"$\alpha_\mathrm{best, SK}$",                "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}}}}",                     "alpha":r"$\alpha_\mathrm{SK}$",                "sigma":"\\sigma_{{\\mathrm{{SK}}_\\mathrm{{B}}}}",                      "filename":inpath_SK_vv, "color":"blue"},
    #     "second":{"label":r"SK$_\mathrm{def}$", "alphabest":r"$\alpha_{\mathrm{best, SK}_\mathrm{def}}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{def}}}}", "alpha":r"$\alpha_{\mathrm{SK}_\mathrm{def}}$", "sigma":"\\sigma_{{\\mathrm{{SK}}_\\mathrm{{def}}}}", "filename":inpath_SK_imp_vv, "color":"red"},
    #     "filename":outpath + "A_B_SK_SKimp",
    #     "label":"A",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.98, 0.99, 1.0, 1.01, 1.02], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-0.8, "maxval":0.8, "ticks":[-0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4], "bins":np.linspace(-1, 1, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.005, "maxval":0.005, "ticks":[-0.004, -0.002, 0], "bins":np.linspace(-0.05, 0.05, 251)},
    #     "deltasigmaoversigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.08, -0.04, 0, 0.04, 0.08], "bins":np.linspace(-2, 2, 251)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }

    # dict_labels_X = {
    #     "first":{"label":"SK$_\mathrm{{B}}$",                  "alphabest":r"$\alpha_\mathrm{best, SK}$",               "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}}}}",                     "alpha":r"$\alpha_\mathrm{SK}$",                "sigma":"\\sigma_{{\\mathrm{{SK}}_\\mathrm{{B}}}}",                     "filename":inpath_SK_vh, "color":"blue"},
    #     "second":{"label":r"SK$_\mathrm{def}$", "alphabest":r"$\alpha_{\mathrm{best, SK}_\mathrm{def}}", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{def}}}}", "alpha":r"$\alpha_{\mathrm{SK}_\mathrm{def}}$", "sigma":"\\sigma_{{\\mathrm{{SK}}_\\mathrm{{def}}}}", "filename":inpath_SK_imp_vh, "color":"red"},
    #     "filename":outpath + "X_B_SK_SKimp",
    #     "label":"X",
    #     "color":"red",
    #     "alphamed":{"ticks":[0.98, 0.99, 1.0, 1.01, 1.02], "bins":np.linspace(0.8, 1.2, 101)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-0.8, "maxval":0.8, "ticks":[-0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4], "bins":np.linspace(-1, 1, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.005, "maxval":0.005, "ticks":[-0.004, -0.002, 0], "bins":np.linspace(-0.05, 0.05, 251)},
    #     "deltasigmaoversigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.08, -0.04, 0, 0.04, 0.08], "bins":np.linspace(-2, 2, 251)},
    #     "deltaalphaoveravgsigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], "bins":np.linspace(-2, 2, 101)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.dr.2pcf",
    #     "begcffile":""
    #     }
    # plot_comparison_UL(dict_labels_A, dict_labels_X, keys=["first","second"], figname=outpath + "/AX_B_SK_SK_imp")
    
    # exit()
    ### Plot 5: A_B_CGrec_SKrec
    # dict_labels = {
    #     "first": {"label":"CG$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_\\mathrm{{CG}}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath104, "color":"green"},
    #     "second":{"label":"SK$_\\mathdefault{{B}}$", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_\\mathrm{{SK}}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}", "filename":inpath105, "color":"blue"},
    #     "filename":outpath + "A_B_CGrec_SKrec",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.95, 1.0, 1.05], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.8, -0.4, 0, 0.4, 0.8], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":"void"
    #     }
    # plot_rec_comparison(dict_labels, keys=["first","second"], outliers=False)

    ###### LC Plots  
    ######

    ### Plot A; First volume
    # inpath_lc = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/"
    # dict_labels = {
    #     "first":{"label":"PAR", "alphabest":r"$\alpha_\mathrm{best, PAR}$", "alphamed":"\\alpha_\\mathrm{{PAR}}", "alpha":r"$\alpha_\mathrm{PAR}$", "sigma":"\\sigma_\\mathrm{{PAR}}", "filename":inpath_lc + "parab_60_150_gauss/", "color":"magenta"},
    #     "second":{"label":"CG", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_\\mathrm{{CG}}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
    #     "third":{"label":"SK", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_\\mathrm{{SK}}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}", "filename":inpath_lc + "avg_16R_2000_SK_60_150/","color":"blue"},
    #     "forth":{"label":"SKB", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_\\mathrm{{SKB}}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc + "stitched_G2048_50_G512_2000_SK_B/", "color":"magenta"},
    #     "fifth":{"label":"CGB", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_\\mathrm{{CGB}}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
    #     "sixth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$", "filename":inpath_lc + "parab_60_150_fixc/", "color":"magenta"},
    #     "filename":outpath + "A_LC_CG_SK_CGB_SKB_par_fixc",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.9, 1.0, 1.1, 1.2], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.8, -0.4, 0, 0.4, 0.8], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.2pcf",
    #     "begcffile":""
    #     }

    # plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth", "sixth"], outliers=False)


    ### Plot B; First volume; fixed Sigma_nl
    # inpath_lc_snl="/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real/vv2pcf_new/fixSigmaNL/"    
    # dict_labels = {
    #     "first":{"label":"PAR", "alphabest":r"$\alpha_\mathrm{best, PAR}$", "alphamed":"\\alpha_\\mathrm{{PAR}}", "alpha":r"$\alpha_\mathrm{PAR}$", "sigma":"\\sigma_\\mathrm{{PAR}}", "filename":inpath_lc_snl + "parab_60_150_gauss/", "color":"magenta"},
    #     "second":{"label":"CG", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_\\mathrm{{CG}}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc_snl + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
    #     "third":{"label":"SK", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_\\mathrm{{SK}}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}", "filename":inpath_lc_snl + "avg_16R_2000_SK_60_150/","color":"blue"},
    #     "forth":{"label":"SKB", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_\\mathrm{{SKB}}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc_snl + "stitched_G2048_50_G512_2000_SK_B/", "color":"magenta"},
    #     "fifth":{"label":"CGB", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_\\mathrm{{CGB}}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc_snl + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
    #     "sixth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$", "filename":inpath_lc_snl + "parab_60_150_fixc/", "color":"magenta"},
    #     "filename":outpath + "A_LC_fixsnl_CG_SK_CGB_SKB_par_fixc",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.9, 1.0, 1.1, 1.2], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.8, -0.4, 0, 0.4, 0.8], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.2pcf",
    #     "begcffile":""
    #     }

    # plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth", "sixth"], outliers=False)


    ### From 27.07.2022
    ### Plot C; Larger volume; fixed Sigma_nl
    # inpath_lc = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100X_2_20X/vv2pcf/"
    # dict_labels = {
    #     "first":{"label":"PAR", "alphabest":r"$\alpha_\mathrm{best, PAR}$", "alphamed":"\\alpha_{{\\mathrm{{PAR}}_\\mathdefault{{G}}}}", "alpha":r"$\alpha_\mathrm{PAR}$", "sigma":"\\sigma_\\mathrm{{PAR}}",           "filename":inpath_lc + "parab_60_150_gauss/", "color":"magenta"},
    #     "second":{"label":"CG$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
    #     "third":{"label":"SK$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}",  "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_LC/","color":"blue"},
    #     "forth":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc + "stitched_G2048_50_G512_2000_SK_B_60_150/", "color":"magenta"},
    #     "fifth":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
    #     "sixth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix~c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$",                       "filename":inpath_lc + "parab_60_150_fixc/", "color":"magenta"},
    #     "filename":outpath + "A_LC_fixsnl_CG_SK_CGB_SKB_par_fixc_largev",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.95, 1.0, 1.05], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.8, -0.4, 0, 0.4, 0.8], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.2pcf",
    #     "begcffile":"void"
    #     }

    # plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth", "sixth"], outliers=False)
    # plot_corner_2outliers(dict_labels, keys=["first", "second", "third", "forth", "fifth", "sixth"])
    
    # ### Plot D; Larger volume; fixed Sigma_nl, no parab with gauss model
    # dict_labels = {
    #     "first":{"label":"CG$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
    #     "second":{"label":"SK$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_LC/","color":"blue"},
    #     "third":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc + "stitched_G2048_50_G512_2000_SK_B_60_150/", "color":"magenta"},
    #     "forth":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
    #     "fifth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix~c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$", "filename":inpath_lc + "parab_60_150_fixc/", "color":"magenta"},
    #     "filename":outpath + "A_LC_fixsnl_CG_SK_CGB_SKB_fixc_largev",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.9, 1.0, 1.1, 1.2], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-2, -1, 0, 1, 2], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.2, 0, 0.2], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.2pcf",
    #     "begcffile":"void"
    #     }

    # plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth"], outliers=False)


    ### From 27.07.2022
    ### Plot E; Larger volume; fixed Sigma_nl
    inpath_lc = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100X_2_20X/vhxcf/"
    
    dict_labels = {
        "first":{"label":"PAR", "alphabest":r"$\alpha_\mathrm{best, PAR}$", "alphamed":"\\alpha_{{\\mathrm{{PAR}}_\\mathdefault{{G}}}}", "alpha":r"$\alpha_\mathrm{PAR}$", "sigma":"\\sigma_\\mathrm{{PAR}}",           "filename":inpath_lc + "parab_60_150_gauss/", "color":"magenta"},
        "second":{"label":"CG$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
        "third":{"label":"SK$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}",  "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_LC/","color":"blue"},
        "forth":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_B/", "color":"magenta"},
        "fifth":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
        "sixth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix~c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$",                       "filename":inpath_lc + "parab_60_150_fixc/", "color":"magenta"},
        "filename":outpath + "X_LC_fixsnl_CG_SK_CGB_SKB_par_fixc_largev",
        "color":"blue",
        "alphamed":{"ticks":[0.95, 1.0, 1.05], "bins":np.linspace(0.8, 1.2, 51)},
        "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-10, 10, 101)},
        "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
        "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.8, -0.4, 0, 0.4, 0.8], "bins":np.linspace(-2, 2, 101)},
        "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
        "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
        "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
        "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
        "endcffile":".dat.xcf",
        "begcffile":"void"
        }

    plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth", "sixth"], outliers=True)
    # plot_corner_2outliers(dict_labels, keys=["first", "second", "third", "forth", "fifth", "sixth"])
    
    ### Plot F; Larger volume; fixed Sigma_nl, no parab with gauss model
    # dict_labels = {
    #     "first":{"label":"CG$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, CG}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{CG}$", "sigma":"\\sigma_\\mathrm{{CG}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_LC/", "color":"green"},
    #     "second":{"label":"SK$_\\mathrm{{LC}}$", "alphabest":r"$\alpha_\mathrm{best, SK}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{LC}} }}", "alpha":r"$\alpha_\mathrm{SK}$", "sigma":"\\sigma_\\mathrm{{SK}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_LC/","color":"blue"},
    #     "third":{"label":"SK$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, SKB}$", "alphamed":"\\alpha_{{\\mathrm{{SK}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{SKB}$", "sigma":"\\sigma_\\mathrm{{SKB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_SK_B/", "color":"magenta"},
    #     "forth":{"label":"CG$_\\mathrm{{B}}$", "alphabest":r"$\alpha_\mathrm{best, CGB}$", "alphamed":"\\alpha_{{\\mathrm{{CG}}_\\mathrm{{B}} }}", "alpha":r"$\alpha_\mathrm{CGB}$", "sigma":"\\sigma_\\mathrm{{CGB}}", "filename":inpath_lc + "stitched_16R_G2048_50_G512_2000_CG_B/", "color":"magenta"},
    #     "fifth":{"label":"fix c", "alphabest":r"$\alpha_\mathrm{best, fix c}$", "alphamed":"\\alpha_\\mathrm{{fix~c}}", "alpha":r"$\alpha_\mathrm{fix c}$", "sigma":"\\sigma_\\mathrm{{fix c}}$", "filename":inpath_lc + "parab_60_150_fixc/", "color":"magenta"},
    #     "filename":outpath + "X_LC_fixsnl_CG_SK_CGB_SKB_fixc_largev",
    #     "color":"blue",
    #     "alphamed":{"ticks":[0.9, 1.0, 1.1, 1.2], "bins":np.linspace(0.8, 1.2, 51)},
    #     "bayesfactor":{"minval":-10, "maxval":10, "ticks":[-2, -1, 0, 1, 2], "bins":np.linspace(-10, 10, 101)},
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-1, "maxval":1, "ticks":[-0.2, 0, 0.2], "bins":np.linspace(-2, 2, 101)},
    #     "deltaalphaoveralpha":{"minval":-0.05, "maxval":0.05, "ticks":[-0.04, -0.02, 0, 0.02, 0.04], "bins":np.linspace(-0.1, 0.1, 201)},
    #     "deltasigmaoversigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "deltaalphaoveravgsigma":{"minval":-1, "maxval":1, "ticks":[-0.6, -0.3, 0, 0.3, 0.6], "bins":np.linspace(-1, 1, 51)},
    #     "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "endcffile":".dat.xcf",
    #     "begcffile":"void"
    #     }

    # plot_comparison(dict_labels, keys=["first", "second", "third","forth", "fifth"], outliers=False)
        

    #####
    # path = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/"

    #### Plot 5: A_B_DAMPA1_DAMPA2
    # dict_labels1 = {
    #     "first": {"label":"a=1", "alphabest":r"$\alpha_\mathrm{best, a1}$", "alphamed":r"$\alpha_\mathrm{{a1}}$", "alpha":r"$\alpha_\mathrm{a1}$", "sigma":r"$\sigma_\mathrm{a1}$", "filename":path + "/damp_a_1/parab_60_150_fixc/", "color":"blue"},
    #     "second":{"label":"a=2", "alphabest":r"$\alpha_\mathrm{best, a2}$", "alphamed":r"$\alpha_\mathrm{{a2}}$", "alpha":r"$\alpha_\mathrm{a2}$", "sigma":r"$\sigma_\mathrm{a2}$", "filename":path + "/parab/vv2pcf/parab_60_150_fixc/", "color":"blue"},
    #     "filename":outpath + "/A_B_DAMP1_DAMP2",
    #     "label":"A",
    #     "color":"blue",
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-0.8, "maxval":0.8, "ticks":[-0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4], "bins":np.linspace(-1, 1, 101)},
    #     "deltaalphaoveravgsigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15], "bins":np.linspace(-2, 2, 201)},
    #     "deltasigmaoversigma":{"minval":-0.1, "maxval":0.1, "ticks":[-0.09, -0.06, -0.03, 0, 0.03, 0.06, 0.09], "bins":np.linspace(-2, 2, 401)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }

    # dict_labels2 = {
    #     "first": {"label":"a=1", "alphabest":r"$\alpha_\mathrm{best, a1}$", "alphamed":r"$\alpha_\mathrm{{a1}}$", "alpha":r"$\alpha_\mathrm{a1}$", "sigma":r"$\sigma_\mathrm{a1}$", "filename":path + "/damp_a_1/stitched_16R_G2048_50_G512_2000_CG/", "color":"orange"},
    #     "second":{"label":"a=2", "alphabest":r"$\alpha_\mathrm{best, a2}$", "alphamed":r"$\alpha_\mathrm{{a2}}$", "alpha":r"$\alpha_\mathrm{a2}$", "sigma":r"$\sigma_\mathrm{a2}$", "filename":path + "/vv2pcf_CG/stitched_16R_G2048_50_G512_2000_CG/", "color":"orange"},
    #     "filename":outpath + "/A_B_DAMP1_DAMP2",
    #     "label":"A",
    #     "color":"orange",
    #     "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 21)},
    #     "tensionparams":{"minval":-0.8, "maxval":0.8, "ticks":[-0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4], "bins":np.linspace(-1, 1, 101)},
    #     "deltaalphaoveravgsigma":{"minval":-0.2, "maxval":0.2, "ticks":[-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15], "bins":np.linspace(-2, 2, 201)},
    #     "deltasigmaoversigma":{"minval":-0.1, "maxval":0.1, "ticks":[-0.09, -0.06, -0.03, 0, 0.03, 0.06, 0.09], "bins":np.linspace(-2, 2, 401)},
    #     "endcffile":".VOID.dat.2pcf",
    #     "begcffile":""
    #     }

    # plot_dampA1_dampA2(dict_labels1, dict_labels2, keys=["first", "second"], outpath=outpath)


if __name__== '__main__':
    # main()
