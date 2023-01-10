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

def obtain_chi2_alpha_evi(inpath, endcffile=".dat.2pcf", begcffile="void_"):
    files = glob.glob(inpath + "*_.txt")
    # randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_500.txt", usecols=(0), unpack=True)
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_30.txt", usecols=(0), unpack=True)

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
    randnum = np.loadtxt("/home/astro/variu/phd/voids/randnum_30.txt", usecols=(0), unpack=True)

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


def main():
    pt.rcParams.update({'font.family': "serif"})
    pt.rcParams.update({'font.serif': "Times New Roman"})

    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['mathtext.rm'] = 'serif'
    # print(pt.rcParams.keys())
    
    outpath = "/home/astro/variu/void_project_figs/"
    # outpath = "/home/astro/variu/temp/"
     
    

    inpath_1 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100_20/vv2pcf/2pcf_r/"
    inpath_2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100_20/vv2pcf/2pcf_init/"
    inpath_3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100_20/vv2pcf/2pcf_t/"
    inpath_4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/lightcone_box1_1000CF/real_large_V_100_20/vv2pcf/2pcf_s/"
    
    
    # obtain_chi2_alpha_evi(inpath_1, endcffile=".dat.2pcf")
    # obtain_chi2_alpha_evi(inpath_2, endcffile=".dat.2pcf")
    # obtain_chi2_alpha_evi(inpath_3, endcffile=".dat.2pcf")
    # obtain_chi2_alpha_evi(inpath_4, endcffile=".dat.2pcf")
    # exit()
    ### Plot 1: AX_B_CG_SK_GAL_fixc
    dict_labels = {
        "first":{"label":"2pcf$_\\mathrm{{r}}$", "alphabest":r"$\alpha_\mathdefault{best, 2pcf r}$",    "alphamed":"\\alpha_{{\\mathdefault{{2pcf}}_\\mathrm{{r}}}}", "alpha":r"$\alpha_\mathdefault{2pcf r}$", "sigma":"\\sigma_\\mathdefault{{2pcf r}}", "filename":inpath_2, "color":"green"},
        "second":{"label":"2pcf init",             "alphabest":r"$\alpha_\mathdefault{best, 2pcf init}$", "alphamed":"\\alpha_\\mathdefault{{2pcf init}}", "alpha":r"$\alpha_\mathdefault{2pcf init}$", "sigma":"\\sigma_\\mathdefault{2pcf init}", "filename":inpath_1, "color":"red"},
        "third":{"label":"2pcf$_\\mathrm{{t}}$",  "alphabest":r"$\alpha_\mathdefault{best, 2pcf t}$",    "alphamed":"\\alpha_{{\\mathdefault{{2pcf}}_\\mathrm{{t}}}}", "alpha":r"$\alpha_\mathdefault{2pcf t}$", "sigma":"\\sigma_\\mathdefault{{2pcf t}}", "filename":inpath_3,"color":"blue"},
        "forth":{"label":"2pcf$_\\mathrm{{s}}$",  "alphabest":r"$\alpha_\mathdefault{best, 2pcf s}$",    "alphamed":"\\alpha_\\mathdefault{{2pcf s}}", "alpha":r"$\alpha_\mathdefault{2pcf s}$", "sigma":"\\sigma_\\mathdefault{{2pcf s}}", "filename":inpath_4, "color":"magenta"},
        "filename":outpath + "A_B_init_r_t_s_",
        "label":"A",
        "color":"blue",
        "alphamed":{"ticks":[0.94, 0.96, 0.98, 1.0, 1.02], "bins":np.linspace(0.8, 1.2, 51)},
        "bayesfactor":{"minval":-15, "maxval":15, "ticks":[-2.0, -1, 0, 1, 2.0], "bins":np.linspace(-15,15, 51)},
        "pullone":{"minval":-4, "maxval":4, "ticks":[-4, -2, 0, 2, 4], "bins":np.linspace(-4, 4, 11)},
        "tensionparams":{"minval":-2, "maxval":0.8, "ticks":[-1.6, -1.2, -0.9, -0.6, -0.3, 0, 0.3], "bins":np.linspace(-2, 2, 81)},
        "deltaalphaoveralpha":{"minval":-0.01, "maxval":0.02, "ticks":[-0.004, 0, 0.004, 0.008, 0.012], "bins":np.linspace(-0.05, 0.05, 101)},
        "deltasigmaoversigma":{"minval":-0.5, "maxval":0.5, "ticks":[-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3], "bins":np.linspace(-2, 2, 101)},
        "deltaalphaoveravgsigma":{"minval":-1, "maxval":1.5, "ticks":[-0.4, 0, 0.4, 0.8, 1.2], "bins":np.linspace(-2, 2, 101)},
        "alphamoneoversigma":{"minval":-4, "maxval":4, "ticks":[-2, 0, 2, 4], "bins":np.linspace(-4, 4, 11)},
        "endcffile":".dat.2pcf",
        "begcffile":"void_"
        }
   
    plot_comparison(dict_labels, keys=["first", "second", "third","forth"], outliers=False)
    
    

if __name__== '__main__':
    main()