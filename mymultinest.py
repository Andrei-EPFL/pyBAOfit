import configparser
import numpy as np
import pymultinest as pmn
import sys

class MultinestClass():
    def __init__(self, config_file, outbase, chi2_var):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.live_points = config["multinest"].getint("live_points")
        self.tol = config["multinest"].getfloat("tol")
        self.verbose = config["multinest"].getboolean("verbose")
        
        self.outbase = outbase
        self.chi2_var = chi2_var

        self.parameters = self.chi2_var.list_params
        self.n_params = len(self.parameters)
        print(("INFO: The parameters of the model are ["+', '.join(['%s']*len(self.parameters))+"]") % tuple(self.parameters))
        
        self.prior_params = {}
        for p in self.parameters:
            self.prior_params[p] = tuple(map(float,config.get('priors', p).split(',')))
        print("INFO: The priors of the parameters are " + str(self.prior_params))
        self.c_flat_min = config.getfloat('shapepriors', 'c_flat_min')
        self.c_flat_max = config.getfloat('shapepriors', 'c_flat_max')
        self.c_width = config.getfloat('shapepriors', 'c_width')
                
    def prior(self, cube, ndim, nparams):
        for i in range(ndim):
            val_max = self.prior_params[self.parameters[i]][1]
            val_min = self.prior_params[self.parameters[i]][0]
            cube[i] = cube[i] * (val_max - val_min) + val_min

    def loglike(self, cube, ndim, nparams):
        lnlike = -0.5 * self.chi2_var.chi2_func(cube[0], cube[1:ndim])
        if('c' not in self.parameters or self.c_width == float("inf")):
            return lnlike
        if(cube[3] < self.c_flat_min):
            lnlike = lnlike - 0.5 * ( (cube[3] - self.c_flat_min) / self.c_width )**2
        elif (cube[3] > self.c_flat_max):
            lnlike = lnlike - 0.5 * ( (cube[3] - self.c_flat_max) / self.c_width )**2
        return lnlike

    def run_multinest(self):
        pmn.run(self.loglike, self.prior, self.n_params, outputfiles_basename=self.outbase, resume=True, \
            verbose=self.verbose, n_live_points=self.live_points, evidence_tolerance=self.tol)
        
    def analyse_multinest(self):
        res = pmn.Analyzer(outputfiles_basename=self.outbase, n_params=self.n_params)
        bestpar = res.get_best_fit()['parameters']
        stats = res.get_stats()['marginals']
        median = stats[0]['median']
        print(f"INFO: The median of alpha is {median}")
        sb, bestfit, bestfit_nw, bestfit_broadband = self.chi2_var.best_fit(bestpar[0], bestpar[1:])
        print("INFO: The chi2 for %s is equal to %f " %( self.outbase, self.chi2_var.chi2_func(bestpar[0], bestpar[1:]) ))
        print(("INFO: The best fit parameters are: ["+', '.join(['%f']*len(bestpar))+"]") % tuple(bestpar))
        bestfile = self.outbase + 'best_nw.dat'
        np.savetxt(bestfile, np.transpose([sb, bestfit, bestfit_nw, bestfit_broadband]), fmt='%.8g')
