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
    def __init__(self, config_file, args):
        config = configparser.ConfigParser()
        config.read(config_file)
        
        input_data = args.input_data
        if input_data is None:
            input_data = config["paths"]["input_data"]

        fit_smin = args.fit_smin
        if fit_smin is None:
            fit_smin = config['params'].getfloat('fit_smin')

        fit_smax = args.fit_smax
        if fit_smax is None:
            fit_smax = config['params'].getfloat('fit_smax')

        min_s_index = args.min_s_index
        if min_s_index is None:
            min_s_index = config['params'].getint('min_s_index')

        self.input_data = input_data
        self.fit_smin = fit_smin
        self.fit_smax = fit_smax
        self.min_s_index = min_s_index                
        
        self.sd, self.xid, self.ndbin = self.get_xid()

        print('STATUS: Get the indices of data for fitting.')
        self.imin, self.imax, self.nidx = get_index(self.sd, self.fit_smin, self.fit_smax)

    def get_xid(self):    
        print('STATUS: Read the input 2PCF to be fitted.')
        try:
            sdt, xidt = np.loadtxt(self.input_data, usecols=(0, 1), unpack=True)
            sd = sdt[self.min_s_index: ]
            xid = xidt[self.min_s_index: ]

            ndbin = sd.size
        except:
            print('ERROR: cannot read the input data.', file=sys.stderr)
            sys.exit(1)
        return [sd, xid, ndbin]
