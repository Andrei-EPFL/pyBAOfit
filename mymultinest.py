import configparser
import numpy as np
import pymultinest as pmn

class MultinestClass():
    def __init__(self, config_file, outbase, chi2_var):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.live_points = config["multinest"].getint("live_points")
        self.tol = config["multinest"].getfloat("tol")
        
        self.outbase = outbase
        self.chi2_var = chi2_var

        self.parameters = ['alpha', 'B', 'Snl']
        self.n_params = len(self.parameters)
        self.amax = config["priors"].getfloat("amax")
        self.amin = config["priors"].getfloat("amin")
        self.bmax = config["priors"].getfloat("bmax")
        self.bmin = config["priors"].getfloat("bmin")
        self.Snlmax = config["priors"].getfloat("Snlmax")
        self.Snlmin = config["priors"].getfloat("Snlmin")
        
    def prior(self, cube, ndim, nparams):
        cube[0] = cube[0] * (self.amax - self.amin) + self.amin
        cube[1] = cube[1] * (self.bmax - self.bmin) + self.bmin
        cube[2] = cube[2] * (self.Snlmax - self.Snlmin) + self.Snlmin

    def loglike(self, cube, ndim, nparams):
        alpha, B, Snl = cube[0], cube[1], cube[2]
        return -0.5 * self.chi2_var.chi2_func([B, Snl], alpha)

    def run_multinest(self):
        pmn.run(self.loglike, self.prior, self.n_params, outputfiles_basename=self.outbase, resume=True, \
            verbose=True, n_live_points=self.live_points, evidence_tolerance=self.tol)
        
    def analyse_multinest(self):
        res = pmn.Analyzer(outputfiles_basename=self.outbase, n_params=self.n_params)
        bestpar = res.get_best_fit()['parameters']
    
        sb, bestfit = self.chi2_var.best_fit(bestpar[1:], bestpar[0])
        print("The chi2 for %s is equal to %f " %( self.outbase, self.chi2_var.chi2_func(bestpar[1:], bestpar[0]) ))
        print("The best fit parameters are:")
        print(bestpar)
        bestfile = self.outbase + 'best.dat'
        np.savetxt(bestfile, np.transpose([sb, bestfit]), fmt='%.8g')
