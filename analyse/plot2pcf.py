#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as pt
import numpy as np
import sys
import glob
import matplotlib.backends.backend_pdf

def plot_errorbar(path, leglabel, ax, col="red", line="-", plotBOOL=True, bias=1, type_="cf", col2=1):
       
    files = glob.glob(path)
    print("For: {}; There are {}".format(leglabel, len(files)))
    x_aux, y_aux = np.loadtxt(files[0], comments="#", usecols=(0,col2), unpack=True)
    y_arr = np.zeros(( len(files), len(y_aux) ))
    
    for i, filename in enumerate(files):
        try:
            x, y_tmp = np.loadtxt(filename, comments="#", usecols=(0,col2), unpack=True)
            y_arr[i] = y_tmp
            # if(type_ == "cf"):
            #     ax.errorbar(x, x*x*y_tmp, color="grey")
            # elif(type_ == "pk"):
            #     ax.errorbar(x, y_tmp, color="grey")
            #     ax.set_yscale('log')
            #     ax.set_xscale('log')
        except:
            continue
    y_avg = np.mean(y_arr, 0)
    y_std = np.std(y_arr, 0)

    if(plotBOOL):
        #for i in range(y_arr.shape[1]):
        #    pt.figure(i)
        #    pt.hist(y_arr[:,i], bins=20)
        if(type_ == "pk"):
            ax.errorbar(x, y_avg, yerr=y_std, label=leglabel, linestyle=line, color=col)
            ax.set_yscale('log')
            ax.set_xscale('log')
        elif (type_ == "cf"):
            ax.errorbar(x, x*x*y_avg, yerr=x*x*y_std, label=leglabel, linestyle=line, color=col)
        else:
            print("type_ = pk or type_ = cf")
            exit(1)

    return x, y_avg, y_std, len(files)
        
def plot_3n(inpath, dictionary_list, limits, cf_std, ax, chi2_val=-1, evidence=-1, DOF=-1, mode="avg500", nmocks=1):
    if(mode=='avg500'):
        partname = "avg"
        label = "Average of 500\nPATCHY mocks"
        cf_std = cf_std / np.sqrt(nmocks) 
    elif(mode=="MD"):
        partname = "voids_Box_z0.465600"
        label="N-Body sim."
    else:
        print("The modes are avg500 and MD")
        sys.exit()
    s, cf_data = np.loadtxt(inpath + "/" +limits +"/"+partname+"_"+limits+".2pcf", usecols=(0,1), unpack=True)
    ax.errorbar(s, s*s*cf_data, yerr=s*s*cf_std, color='k',ls='--', label=label)
    
    for dictionary in dictionary_list:
        sb, cf_best = np.loadtxt(dictionary['path'] + "/BAOfit_"+partname+"_" + limits + ".2pcf_best.dat", usecols=(0,1), unpack=True)
        ax.plot(sb, sb * sb * cf_best, color=dictionary['color'], label=dictionary['label']) # \n' + r'$\chi^2$=' + str(chi2_val[0]) + '; ' + str(DOF[0]) + ' DOF\nln(ev)=' + str(evidence[0]) + "\n"))
     
    range_ = np.where(np.logical_and(s <=165, s >=45))
    ymax = 1.2 * np.max(cf_data[range_] * s[range_]**2)
    ymin = 1.2 * np.min(cf_data[range_] * s[range_]**2)
    
    ax.set_xlabel(r"s [h$^{-1}$ Mpc]")
    ax.set_ylabel(r"s$^2\xi$(s) [h$^{-2}$ Mpc$^2$]")

    ax.set_ylim([ymin,ymax])
    #ax.set_ylim([-100,50])
    ax.set_xlim(left=45, right=165)
    
    ax.axvline(x=60, color='grey', linestyle='--')
    ax.axvline(x=150, color='grey', linestyle='--')
    
    ax.legend(bbox_to_anchor=[0.8, 0.74], loc='center', framealpha=1)
    #ax.legend(loc='upper right', framealpha=1)


def main_plot(Rmin_arr, Rmax_arr, dictionary_list, inpath="/path/to/input", inpath_MD="/path/to/input", outpath="/path/to/output",  figname="figname.ext", mode="avg500"):
    outputfile = outpath + figname

    for i, Rmin in enumerate(Rmin_arr):
        Rmax = Rmax_arr[i]
        fig, ax = pt.subplots() 
        
        if(Rmax==1000):
            limits = str(Rmin) + "R"
                  
        else:
            limits = str(Rmin) + "R" + str(Rmax)

        _, _, cf_std, nmocks = plot_errorbar(inpath+"/"+limits + "/CA*.2pcf", limits + " avg 500 files", ax, col="red", line="-", plotBOOL=False )   
        print("The number of mocks is ", nmocks)
        
        #plot_3n(inpath, dictionary_list, limits, cf_std, ax, mode=mode, nmocks=nmocks)
        plot_3n(inpath_MD, dictionary_list, limits, cf_std, ax, mode=mode, nmocks=nmocks)
        
        
        fig.tight_layout()
        fig.savefig(outputfile)  

def temp_plot():
    bestfitpath_avg500_voi3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150_mocks/BAOfit_CATALPTCICz0.466G960S527868848.VOID.dat.2pcf_*best*"
    files = glob.glob(bestfitpath_avg500_voi3)
    sb, cfb = np.loadtxt(files[0], usecols=(0,1), unpack=True)
    s, cf_data = np.loadtxt("/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf//16R/CATALPTCICz0.466G960S527868848.VOID.dat.2pcf", usecols=(0,1), unpack=True)
    pt.plot(s, s*s*cf_data, color="grey")
    pt.plot(sb, sb*sb*cfb, color="red")
    pt.show()    

def main():

    inpath_avg500 = "/scratch/variu/clustering/patchy_cmass_subset/box1/real/overlapping/vv2pcf/"
    inpath_MD="/scratch/variu/clustering/multidark_simulation/box1/real/"
    bestfitpath_avg500_gal = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R//galaxy_60_150_fast/"
    bestfitpath_avg500_par = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/parab_60_150/"
    bestfitpath_avg500_voi = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R//voidtemplate_conv_60_150//"
    bestfitpath_avg500_voi2 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R//voidtemplate_60_150_G2048L512_1000bins/"
    bestfitpath_avg500_voi3 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G1024CIC_60_150/"
    bestfitpath_avg500_voi4 = "/scratch/variu/phd_fitOut/patchy_cmass_subset/box1/real/overlapping/vv2pcf/16R/void_G512CIC_60_150/"
        
    dictionary_gal = { 'path':bestfitpath_avg500_gal, 'color':"magenta", 'label':'Galaxy' }
    dictionary_par = { 'path':bestfitpath_avg500_par, 'color':"green", 'label':'Parabola' }
    dictionary_voi = { 'path':bestfitpath_avg500_voi, 'color':"cyan", 'label':'Void temp. conv' }
    dictionary_voi2 = { 'path':bestfitpath_avg500_voi2, 'color':"orange", 'label':'G2048L512' }
    dictionary_voi3 = { 'path':bestfitpath_avg500_voi3, 'color':"blue", 'label':'G1024CIC'  }
    dictionary_voi4 = { 'path':bestfitpath_avg500_voi4, 'color':"blue", 'label':'G512CIC' }
    
    dictionary_list = [dictionary_gal, dictionary_par, dictionary_voi2, dictionary_voi3]#, dictionary_voi4]
    
    outpath = "/home/astro/variu/"
     
    #Rmin_arr = [4, 6, 8,  10, 12, 14, 16, 18, 20, 22, 16]
    #Rmax_arr = [6, 8, 10, 12, 14, 16, 18, 20, 22, 30, 1000]
    
    Rmin_arr = [16]
    Rmax_arr = [1000]
    pt.rcParams['font.size'] = 14
    main_plot(Rmin_arr, Rmax_arr, dictionary_list, inpath=inpath_avg500, inpath_MD=inpath_MD, outpath=outpath,  figname="bestfitcurve_MD_16R_avg500_60_150.pdf", mode="MD")
    #temp_plot()


if __name__== '__main__':
    main()
