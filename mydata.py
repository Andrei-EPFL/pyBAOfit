import sys
import configparser
import numpy as np

def get_index(sd, fitmin, fitmax):
    '''Get indices of the data to be fitted.
    Arguments:
      sd: s bins for the data;
      fitmin: minimum s for the fitting;
      fitmax: maximum s for the fitting.
    Return: [index_min, index_max, num_index].'''
    imin = imax = -1
    for i in range(sd.size):
        if imin == -1 and sd[i] >= fitmin:
            imin = i
        if imax == -1 and sd[i] > fitmax:
            imax = i
            break

    if imin == -1 or imax == -1:
        print('ERROR: cannot find the fitting range in data.', file=sys.stderr)
        sys.exit(1)

    nidx = imax - imin
    if nidx < 1:
        print('ERROR: cannot find enough bins for fitting.', file=sys.stderr)
        sys.exit(1)

    return [imin, imax, nidx]

class XiData():
    def __init__(self, config_file, input_data):
        config = configparser.ConfigParser()
        config.read(config_file)
        
        self.fit_smin = config['params'].getfloat('fit_smin')
        self.fit_smax = config['params'].getfloat('fit_smax')
        self.input_data = input_data
        
        self.sd, self.xid, self.ndbin = self.get_xid()

        print('STATUS: Get the indices of data for fitting.')
        self.imin, self.imax, self.nidx = get_index(self.sd, self.fit_smin, self.fit_smax)
        

    def get_xid(self):    
        print('STATUS: Read the input 2PCF to be fitted.')
        try:
            sd, xid = np.loadtxt(self.input_data, usecols=(0, 1), unpack=True)
            ndbin = sd.size
        except:
            print('ERROR: cannot read the input data.', file=sys.stderr)
            sys.exit(1)
        return [sd, xid, ndbin]
