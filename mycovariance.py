import configparser
import numpy as np

class CovMat():
    def __init__(self, config_file, input_mocks):
        config = configparser.ConfigParser()
        config.read(config_file)

        self.compute_cov = config["paths"].getboolean("compute_cov")
        self.save_cov = config["paths"].getboolean("save_cov")
        self.cov_file = config["paths"]["cov_file"]
        self.input_mocks = input_mocks
                
        self.nmock, self.Rcov = self.get_cov()
        
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
                ximock[i] = np.loadtxt(mocks[i], usecols=(1, ), unpack=True)
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
        