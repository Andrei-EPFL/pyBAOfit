import matplotlib.pyplot as pt
import numpy as np


def find_peaks_in_chi2alpha(inpath, partname, chi2_var):
    filename = "/BAOfit_"+partname+"_.txt"
    infile = inpath + filename

    chi2, alpha, B, Snl = np.loadtxt(infile, usecols=(1, 2, 3, 4), unpack=True)
    
    #alphalim1 = np.logical_and(alpha<1.08, alpha>1.05)
    #alphalim2 = np.logical_and(alpha<1.05, alpha>1.027)
    
    #alphalim3 = np.logical_and(alpha<1.027, alpha>0.99)
    
    alphalim1 = np.logical_and(alpha<1.03, alpha>1.015)
    alphalim2 = np.logical_and(alpha<1.015, alpha>0.98)
    
    pos1 = chi2 == np.min(chi2[alphalim1])
    pos2 = chi2 == np.min(chi2[alphalim2])
    #pos3 = chi2 == np.min(chi2[alphalim3])

    print(np.min(chi2))
    print(np.min(chi2[alphalim1]), alpha[pos1], B[pos1], Snl[pos1], chi2_var.chi2_func(alpha[pos1], [B[pos1], Snl[pos1]]))
    print(np.min(chi2[alphalim2]), alpha[pos2], B[pos2], Snl[pos2], chi2_var.chi2_func(alpha[pos2], [B[pos2], Snl[pos2]]))
    #print(np.min(chi2[alphalim3]), alpha[pos3], B[pos3], Snl[pos3], chi2_var.chi2_func(alpha[pos3], [B[pos3], Snl[pos3]]))
    pt.plot(alpha[alphalim1], np.exp(-chi2[alphalim1]), ls="", marker=".")
    pt.axvline(alpha[pos1], color="grey")
    pt.plot(alpha[alphalim2], np.exp(-chi2[alphalim2]), ls="", marker=".")
    pt.axvline(alpha[pos2], color="grey")
    #pt.plot(alpha[alphalim3], np.exp(-chi2[alphalim3]), ls="", marker=".")
    #pt.axvline(alpha[pos3], color="grey")
    #pt.ylim([0.4, 3.6])
    pt.xlim([0.98, 1.11])
    pt.savefig("./output/temporaryalpha.pdf")
    
    np.savetxt("./output/"+filename + "_best11", np.array(chi2_var.best_fit(alpha[pos1], [B[pos1], Snl[pos1]])).T, header="alpha="+str(alpha[pos1]) +"  B=" + str(B[pos1]) + "  Snl=" + str(Snl[pos1]) + "  chi2="+str(np.min(chi2[alphalim1])))
    np.savetxt("./output/"+filename + "_best21", np.array(chi2_var.best_fit(alpha[pos2], [B[pos2], Snl[pos2]])).T, header="alpha="+str(alpha[pos2]) +"  B=" + str(B[pos2]) + "  Snl=" + str(Snl[pos2]) + "  chi2="+str(np.min(chi2[alphalim2])))
    #np.savetxt("./output/"+filename + "_best3", np.array(chi2_var.best_fit(alpha[pos3], [B[pos3], Snl[pos3]])).T, header="alpha="+str(alpha[pos3]) +"  B=" + str(B[pos3]) + "  Snl=" + str(Snl[pos3]) + "  chi2="+str(np.min(chi2[alphalim3])))
    pt.show()

def plot_best_fit(chi2_var, xim_var):
    sb, xib = chi2_var.best_fit_broadband(1.06, [2.5, 23.6])

    
    import matplotlib.pyplot as pt
  
    pt.plot(sb, sb*sb*xib, ls="-", marker="", color="grey")
    
    s, xi = xim_var.xi_model_2([2.50192557, 23.63312126])
    pt.plot(s, s*s*xi, ls="--", marker="", color="red")
    
    s, xi = xim_var.xi_model_2([1.08561133, 12.91535025])
    pt.plot(s, s*s*xi, ls="--", marker="", color="blue")

    pt.ylim([-200, 200])
    pt.savefig("./output/temporarypoly.pdf")