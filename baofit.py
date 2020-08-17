#!/usr/bin/env python3
import sys
import os
import configparser
import numpy as np
from mydata import XiData
from mycovariance import CovMat
from mymodel import XiModel
from mychi2 import Chi2Class
from mymultinest import MultinestClass

def main():
    output_dir = sys.argv[1]
    input_data = sys.argv[2]
    input_mocks = sys.argv[3]
    input_pvoid = sys.argv[4]
        
    config_file = '/home/epfl/variu/phd/voids/chengscodes/BAOfit/voidnw_FFTlog_myVer/config.ini'
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
    
    if not os.path.isfile(input_mocks):
        print("ERROR: The file: " + input_mocks + " does not exist!")
        sys.exit(1)

    if not os.path.isdir(output_dir):
        print('ERROR: output_dir "{}" does not exist.'.format(output_dir), file=sys.stderr)
        sys.exit(1)
    
    bname = input_data.split('/')[-1]
    outbase = output_dir + '/BAOfit_' + bname + "_"

    xid_var = XiData(config_file, input_data)
    xim_var = XiModel(config_file, input_pvoid, cosmoparams)
    
    if xim_var.npoly >= xid_var.nidx:
        print('ERROR: too many nuisance parameters.', file=sys.stderr)
        sys.exit(1)

    print('STATUS: Read/Compute the covariance matrix.')
    covmat_var = CovMat(config_file, input_mocks)
    
    if covmat_var.nmock < xid_var.nidx + 3:
        print('ERROR: the number of mocks is not enough.', file=sys.stderr)
        sys.exit(1)
    if xid_var.ndbin != covmat_var.Rcov.shape[0]:
        print('ERROR: bin size of data and mocks do not match.', file=sys.stderr)
        sys.exit(1)
    
    print("INFO: The number of bins is %i" %xid_var.nidx)
    
    chi2_var = Chi2Class(xim_var, xid_var, covmat_var)
    #print(chi2_var.chi2_func(1, [1, 1, 1]))
    
    multinest_var = MultinestClass(config_file, outbase, chi2_var)
    #print(multinest_var.loglike([1,1,1,1], 4, 4))
    #print("test")
    #sys.exit()
    
    multinest_var.run_multinest()
    multinest_var.analyse_multinest()
    
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("USAGE: python " + sys.argv[0] + " outpath avg_file 100_file template_pk_file")
        print("USAGE: If do not need a template file put NONE instead of the template_pk_file")
        sys.exit(2)
    main()
