import os
import sys
import configparser
import numpy as np
import matplotlib.pyplot as pt

class CovMat():
    def __init__(self, config_file, args):
        config = configparser.ConfigParser()
        config.read(config_file)

        self.compute_cov = config["paths"].getboolean("compute_cov")
        self.save_cov = config["paths"].getboolean("save_cov")
        self.cov_file = config["paths"]["cov_file"]
        
        input_mocks = args.input_mocks
        if input_mocks is None:
            input_mocks = config["paths"]["input_mocks"]

        min_s_index = args.min_s_index
        if min_s_index is None:
            min_s_index = config['params'].getint('min_s_index')
        
        self.mock_y_column = config['params'].getint('mock_y_column')
        self.rescale_chi2_by = config["params"].getfloat("rescale_chi2_by")

        self.input_mocks = input_mocks
        if not os.path.isfile(self.input_mocks):
            print("ERROR: The file: " + self.input_mocks + " does not exist!")
            sys.exit(1)

        self.min_s_index = min_s_index
        self.nmock, self.Rcov = self.get_cov()
        self.cov = self.get_icov()
        self.icov = np.linalg.inv(self.cov)

    def get_cov(self):
        '''Read/Compute the pre-processed covariance matrix.
        Return: [Nmock, Rcov], where Rcov is the upper triangular matrix from the
        QR decomposition of the mock matrix.'''
        if self.compute_cov == True:
            # Read the list of 2PCF from mocks
            mocks = []
            with open(self.input_mocks) as f:
                for line in f:
                    fname = line.rstrip('\n')
                    if fname != '':
                        mocks.append(fname)
            Nmock = len(mocks)

            # Read 2PCF of mocks
            ximock = [None] * Nmock
            for i in range(Nmock):
                temp = np.loadtxt(mocks[i], usecols=(self.mock_y_column, ), unpack=True)
                ximock[i] = temp[self.min_s_index: ]

            ximock = np.array(ximock)
            
            # Compute the mock matrix M (C = M^T . M)
            mean = np.mean(ximock, axis=0)
            ximock -= mean
            
            # QR decomposition of M
            Rcov = np.linalg.qr(ximock, mode='r')
            
            if self.save_cov == True:
                np.savetxt(self.cov_file, Rcov, header=str(Nmock))
        else:         # comput_cov = False
            with open(self.cov_file) as f:
                Nmock = int(f.readline())
            Rcov = np.loadtxt(self.cov_file, skiprows=1)

        return [Nmock, Rcov]
        
    def get_icov(self):
        '''Read/Compute the pre-processed covariance matrix.
        Return: [Nmock, Rcov], where Rcov is the upper triangular matrix from the
        QR decomposition of the mock matrix.'''
        if self.compute_cov == True:
            # Read the list of 2PCF from mocks
            mocks = []
            with open(self.input_mocks) as f:
                for line in f:
                    fname = line.rstrip('\n')
                    if fname != '':
                        mocks.append(fname)
            Nmock = len(mocks)

            # Read 2PCF of mocks
            ximock = [None] * Nmock
            for i in range(Nmock):
                temp = np.loadtxt(mocks[i], usecols=(self.mock_y_column, ), unpack=True)
                ximock[i] = temp[self.min_s_index: ]

            ximock = np.array(ximock)
            cov_ = np.cov(ximock.T)

        return cov_
        