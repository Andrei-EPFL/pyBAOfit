#!/usr/bin/env python3
import sys
import os
import shutil
import json
import argparse
import configparser
import numpy as np
from mydata import XiData
from mycovariance import CovMat
from mymodel import XiModel
from mychi2 import Chi2Class
from mymultinest import MultinestClass
from myplot import find_peaks_in_chi2alpha, plot_best_fit

def obtain_parameters(parser):
    #### Arguments of the code    

    parser.add_argument("--config", type=str, help="ini file holding configuration")
    parser.add_argument("--output_dir", type=str, help="output directory (overrides config file)")
    parser.add_argument("--input_data", type=str, help="input file (file to be fitted) (same)")
    parser.add_argument("--input_mocks", type=str, help="file that contains the paths to all mocks for covariance (same)")
    parser.add_argument("--input_pvoid", type=str, help="void template file (same)")
    parser.add_argument("--model", type=str, help="the required model: pvoid, galaxy, parab (same)")
    parser.add_argument("--method", type=str, help="integration method: FFTlinlog, FFTlog, fast (same)")
    parser.add_argument("--kmin", type=float, help="min k of the P(k) template (same)")
    parser.add_argument("--kmax", type=float, help="max k of the P(k) template (same)")
    parser.add_argument("--num_lnk_bin", type=int, help="number of log spaced k bins for the P(k) template (same)")
    parser.add_argument("--fit_smin", type=float, help="min s range for fit (same)")
    parser.add_argument("--fit_smax", type=float, help="max s range for fit (same)")
    parser.add_argument("--min_s_index", type=int, help="min s index to be used from the data CF (same)")
    parser.add_argument("--damp_a", type=float, help="the damping factor for the fast integration method")

    args = parser.parse_args()

    #### Configuration file
    config_file = args.config # config file
    if config_file is None:
        config_file = '/home/astro/variu/phd/voids/chengscodes/BAOfit/voidnw_FFTlog_myVer/config.ini'
    config_name = os.path.basename(config_file)

    if not os.path.isfile(config_file):
        print("ERROR: The configuration file: " + config_file + " does not exist!")
        sys.exit(1)

    config = configparser.ConfigParser()
    config.read(config_file)
    cosmoparams = {
        "h":float(config['cosmoparams']['h']),
        "Omega_m":float(config['cosmoparams']['Omega_m']),
        "Omega_b":float(config['cosmoparams']['Omega_b']),
        "Tcmb":float(config['cosmoparams']['Tcmb']),
        #"w":float(config['cosmoparams']['w']),
        #"omnuh2":float(config['cosmoparams']['omnuh2']),
        #"omk":float(config['cosmoparams']['omk']),
        #"helium_fraction":float(config['cosmoparams']['helium_fraction']),
        #"massless_neutrinos":float(config['cosmoparams']['massless_neutrinos']),
        #"nu_mass_eigenstates":float(config['cosmoparams']['nu_mass_eigenstates']),
        #"massive_neutrinos":float(config['cosmoparams']['massive_neutrinos']),
        #"nu_mass_fractions":float(config['cosmoparams']['nu_mass_fractions']),
        #"transfer_kmax":float(config['cosmoparams']['transfer_kmax']),
        #"transfer_redshift":float(config['cosmoparams']['transfer_redshift']),
        "ns":float(config['cosmoparams']['ns'])
        #"scalar_amp":float(config['cosmoparams']['scalar_amp'])
    }
    
    #### Output directory
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = config["paths"]["output_dir"]
    
    if not os.path.isdir(output_dir):
        print('ERROR: output_dir "{}" does not exist.'.format(output_dir), file=sys.stderr)
        sys.exit(1)
    
    #### Copy the config file into the output directory
    if not os.path.isfile(output_dir + config_name):
        shutil.copy(config_file, output_dir + config_name)
    
    #### Write the arguments in a file
    if not os.path.isfile(output_dir + 'args.ini'):
        with open(output_dir + 'args.ini', 'w') as f:
            json.dump(args.__dict__, f, indent=2)
    
    #### Input data
    input_data = args.input_data
    if input_data is None:
        input_data = config["paths"]["input_data"]
    
    bname = input_data.split('/')[-1]
    outbase = output_dir + '/BAOfit_' + bname + "_"
    
    #### Define the class variables
    xid_var = XiData(config_file, args)
    xim_var = XiModel(config_file, args, cosmoparams)
    
    if xim_var.npoly >= xid_var.nidx:
        print('ERROR: too many nuisance parameters.', file=sys.stderr)
        sys.exit(1)

    print('STATUS: Read/Compute the covariance matrix.')
    covmat_var = CovMat(config_file, args)
    
    if covmat_var.nmock < xid_var.nidx + 3:
        print('ERROR: the number of mocks is not enough.', file=sys.stderr)
        sys.exit(1)
    if xid_var.ndbin != covmat_var.Rcov.shape[0]:
        print('ERROR: bin size of data and mocks do not match.', file=sys.stderr)
        sys.exit(1)
    
    print("INFO: The number of bins is %i" %xid_var.nidx)
    
    chi2_var = Chi2Class(xim_var, xid_var, covmat_var)
    #find_peaks_in_chi2alpha(output_dir, os.path.basename(input_data), chi2_var)
    #plot_best_fit(chi2_var, xim_var)
    multinest_var = MultinestClass(config_file, outbase, chi2_var)
    #print("test")
    #sys.exit()
    
    multinest_var.run_multinest()
    multinest_var.analyse_multinest()
    
if __name__ == '__main__':
    main()
