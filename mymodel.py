import sys
import configparser
import numpy as np
import numpy.fft as fft
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.special import spherical_jn, loggamma

def windowfn(x, dlnxleft=0.46, dlnxright=0.46):
    xmin = min(x)
    xmax = max(x)
    xleft = np.exp(np.log(xmin) + dlnxleft)
    xright = np.exp(np.log(xmax) - dlnxright)
    w = np.zeros_like(x)
    w[(x > xleft) & (x < xright)] = 1

    il = (x < xleft) & (x > xmin)
    ir = (x > xright) & (x < xmax)

    rl = (x[il] - xmin) / (xleft - xmin)
    rr = (xmax - x[ir]) / (xmax - xright)
    w[il] = rl - np.sin(np.pi * 2 * rl) / (2 * np.pi)
    w[ir] = rr - np.sin(np.pi * 2 * rr) / (2 * np.pi)
    return w

def calc_Mellnu(tt, alpha, q=0):
    n = q - 1 - 1j * tt
    intjlttn = 2**(n-1) * np.sqrt(np.pi) * \
            np.exp(loggamma((1+n)/2.0) - loggamma((2-n)/2.0))
    A = alpha**(1j * tt - q)
    return A * intjlttn

def calc_phi(pk, k0, N, L, q):
    k = k0 * np.exp(np.arange(0, N) * 2 * np.pi / L)
    P = pk(k)
    
    kpk = (k / k0)**(3-q) * P * windowfn(k)
    phi = np.conj(fft.rfft(kpk)) / L
    phi *= windowfn(k)[len(k) - len(phi):]
    return phi
        
class XiModel():
    '''
    This is a class of the 2PCF model.
    It computes the required model and returns it.
    '''
    def __init__(self, config_file, args, cosmoparams):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.cosmoparams = cosmoparams
        
        self.kmin = args.kmin
        if self.kmin is None:
            self.kmin = config['params'].getfloat('kmin')
        
        self.kmax = args.kmax
        if self.kmax is None:
            self.kmax = config['params'].getfloat('kmax')
        
        self.num_lnk_bin = args.num_lnk_bin
        if self.num_lnk_bin is None:
            self.num_lnk_bin = config['params'].getint('num_lnk_bin')
        
        self.model = args.model
        if self.model is None:
            self.model = config["params"]["model"]

        self.method = args.method
        if self.method is None:
            self.method = config['params']['method']
        
        self.k_norm = config['params'].getfloat('k_norm')
        self.k_interp = config['params'].getboolean('k_interp')
        self.list_params = self.xi_model_params()
        print("INFO: the required model is %s" % self.model)
        print("INFO: the required method is %s" % self.method)
        print(("INFO: The parameters of the model are ["+', '.join(['%s']*len(self.list_params))+"]") % tuple(self.list_params))
        
        self.npoly = config['params'].getint('npoly')

        self.smin = config['params'].getfloat('smin')
        self.smax = config['params'].getfloat('smax')
        self.num_s_bin = config['params'].getint('num_s_bin')
        self.damp_a = config['params'].getfloat('damp_a')
        
        self.input_plin = config['paths']['input_plin']
        self.pnw_run = config['paths'].getboolean('pnw_run')
        
        if(self.pnw_run==False):
            self.input_pnw = config['paths']['input_pnw']
        
        input_pvoid = args.input_pvoid
        if input_pvoid is None:
            input_pvoid = config["paths"]["input_pvoid"]
        self.input_pvoid = input_pvoid

        self.k, self.Plin = self.pk_lin()
        self.k2 = self.k**2
        self.Pnw = self.pk_nw()

        if('pvoid' in self.model): self.Pvoid = self.pk_void()

        self.eka2 = np.exp(-self.k2 * self.damp_a**2) * 0.5 / np.pi**2
        self.sm = np.linspace(self.smin, self.smax, self.num_s_bin)
        self.nkbin = self.k.size
        self.nsbin = self.sm.size

        self.j0 = np.zeros([self.nsbin, self.nkbin])
        for i in range(self.nsbin):
            self.j0[i, :] = spherical_jn(0, self.sm[i] * self.k)
        
    def pk_lin(self):
        try:
            k0, Pk0 = np.loadtxt(self.input_plin, comments='#', unpack=True, usecols=(0, 1))
        except:
            print('ERROR: cannot read the linear power spectrum file.', file=sys.stderr)
            sys.exit(1)

        if self.k_interp == False:
            return [k0, Pk0]   
        
        lnk = np.linspace(np.log(self.kmin), np.log(self.kmax), self.num_lnk_bin)
        k = np.exp(lnk)
        fint = interp1d(np.log(k0), np.log(Pk0), kind='cubic')
        Pk = np.exp(fint(lnk))
        return [k, Pk]

    def pk_void(self):
        try:
            kv, Pv = np.loadtxt(self.input_pvoid, comments='#', unpack=True, usecols=(0, 1))
        except:
            print('ERROR: cannot read the linear power spectrum of voids file.', file=sys.stderr)
            sys.exit(1)

        fv = interp1d(np.log(kv), kv*(Pv), kind='cubic')
        Pvoid = (fv(np.log(self.k))/self.k)
        return Pvoid

    def pk_nw(self):
        '''Compute the non-wiggle matter P(k) with the Eisenstein & Hu 1998 formulae.
        Arguments: array of k and the linear P(k).
        Return: P_nw (k).'''
        if self.pnw_run == True:
            print("INFO: The cosmological parameters for the nw P(k) are" + str(self.cosmoparams))
            Omega_m = self.cosmoparams["Omega_m"]
            Omega_b = self.cosmoparams["Omega_b"]
            h = self.cosmoparams["h"]
            Tcmb = self.cosmoparams["Tcmb"]
            ns= self.cosmoparams["ns"]

            Omh2 = Omega_m * h**2
            Obh2 = Omega_b * h**2
            Ofac = Omega_b / Omega_m
            # Eq. 26
            s = 44.5 * np.log(9.83 / Omh2) / np.sqrt(1 + 10 * Obh2**0.75)
            # Eq. 31
            alpha = 1 - 0.328*np.log(431*Omh2)*Ofac + 0.38*np.log(22.3*Omh2)*Ofac**2
            # Eq. 30
            Gamma = Omega_m * h * (alpha + (1 - alpha) / (1 + (0.43 * self.k * s)**4))
            # Eq. 28
            q = self.k * (Tcmb / 2.7)**2 / Gamma
            # Eq. 29
            L0 = np.log(2 * np.e + 1.8 * q)
            C0 = 14.2 + 731.0 / (1 + 62.5 * q)
            T0 = L0 / (L0 + C0 * q**2)
            Pnw = T0**2 * self.k**ns
        else:  # pnw_run == False      
            try:
                knw, Pnw0 = np.loadtxt(self.input_pnw, comments='#', unpack=True, usecols=(0, 1))
            except:
                print('ERROR: cannot read the linear nw power spectrum file.', file=sys.stderr)
                sys.exit(1)

            fnw = interp1d(np.log(knw), np.log(Pnw0), kind='cubic')
            Pnw = np.exp(fnw(np.log(self.k)))
        # Re-normalize the non-wiggle P(k) with the amplitudes at k < k_norm.
        #np.savetxt('./Eisenstein_nw.txt', np.array([self.k,Pnw]).T)
        
        A = np.mean(self.Plin[self.k < self.k_norm] / Pnw[self.k < self.k_norm])
        Pnw = Pnw * A
        
        #np.savetxt('./Normalized_Eisenstein_nw.txt', np.array([self.k,Pnw]).T)
        
        return Pnw
    
    def xicalc(self, pk, r0=1e-4):
        '''Arguments:
            pk: callable
            N: number of grids for FFT
            kmin, kmax: k range
            r0: minimum r value (~1/kmax)
        '''
        qnu = 1.95
        N = self.num_lnk_bin
        N2 = int(N / 2) + 1
        k0 = self.kmin
        G = np.log(self.kmax / self.kmin)
        alpha = k0 * r0
        L = 2 * np.pi * N / G
        tt = np.arange(0, N2) * 2 * np.pi / G
        rr = r0 * np.exp(np.arange(0, N) * (G / N))
        prefac = k0**3 / (np.pi * G) * (rr / r0)**(-qnu)
        
        Mellnu = calc_Mellnu(tt, alpha, qnu)
        phi = calc_phi(pk, k0, N, L, qnu)
        xi = prefac * fft.irfft(phi * Mellnu, N) * N
        return rr, xi

    def xi_model(self, params):
        '''Compute the template correlation function.
        Arguments:
            k, Plin, Pnw: arrays for the linear power spectra;
            Prt: the ratio between the void and linear non-wiggle power spectra;
            nbin: s bins;
            Snl: the BAO damping factor;
        Return: xi_model.'''
        Snl = params[1]

        Pm = ((self.Plin - self.Pnw) * np.exp(-0.5 * self.k2 * Snl**2) + self.Pnw)
        if(self.model == 'galaxy'):
            Pm = Pm * 1.
        elif(self.model == 'pvoid'):       
            Pm = Pm * (self.Pvoid / self.Pnw)
        elif(self.model == 'parab'):
            c = params[2]
            Pm = Pm * (1 + c * self.k2)    
        else:
            print('ERROR: The model %s does not exist! Exiting the code.' %(self.model))
            print('OPTIONS: galaxy, pvoid, parab')
            sys.exit(1)

        if(self.method == 'FFTlog'):
            Pint = interp1d(np.log(self.k), np.log(Pm), kind='cubic')
            Pkfn = lambda k: np.exp(Pint(np.log(k)))
        elif (self.method == 'FFTlinlog'):
            Pint = interp1d(np.log(self.k), self.k*Pm, kind='cubic')
            Pkfn = lambda k: Pint(np.log(k))/k
        elif (self.method == 'fast'):
            Pm *= self.k2 * self.k * self.eka2
            xim = np.zeros(self.nsbin)
            if self.k_interp == True:
                for i in range(self.nsbin):
                    xim[i] = np.sum(Pm * self.j0[i,:]) * (np.log(self.k[1]) - np.log(self.k[0]))
            else:
                for i in range(self.nsbin):
                    xim[i] = simps(Pm * self.j0[i,:], np.log(self.k))
            return xim
        else:
            print('ERROR: The method %s does not exist! Exiting the code.' %(self.method))
            print('OPTIONS: FFTlinlog, FFTlog')
            sys.exit(1)

        s0, xi0 = self.xicalc(Pkfn, self.sm[0])
        Xint = interp1d(s0, xi0*s0**2, kind='cubic')
        xi = Xint(self.sm) / self.sm**2
        return xi
    
    def xi_model_2(self, params):
        '''Compute the template correlation function.
        Arguments:
            k, Plin, Pnw: arrays for the linear power spectra;
            Prt: the ratio between the void and linear non-wiggle power spectra;
            nbin: s bins;
            Snl: the BAO damping factor;
        Return: xi_model.'''
        B = params[0]
        Snl = params[1]

        Pm = ((self.Plin - self.Pnw) * np.exp(-0.5 * self.k2 * Snl**2) + self.Pnw)
        if(self.model == 'galaxy'):
            Pm = Pm * 1.
        elif(self.model == 'pvoid'):       
            Pm = Pm * (self.Pvoid / self.Pnw)
        elif(self.model == 'parab'):
            c = params[2]
            Pm = Pm * (1 + c * self.k2)    
        else:
            print('ERROR: The model %s does not exist! Exiting the code.' %(self.model))
            print('OPTIONS: galaxy, pvoid, parab')
            sys.exit(1)

        if(self.method == 'FFTlog'):
            Pint = interp1d(np.log(self.k), np.log(Pm), kind='cubic')
            Pkfn = lambda k: np.exp(Pint(np.log(k)))
        elif (self.method == 'FFTlinlog'):
            Pint = interp1d(np.log(self.k), self.k*Pm, kind='cubic')
            Pkfn = lambda k: Pint(np.log(k))/k
        elif (self.method == 'fast'):
            Pm *= self.k2 * self.k * self.eka2
            xim = np.zeros(self.nsbin)
            if self.k_interp == True:
                for i in range(self.nsbin):
                    xim[i] = np.sum(Pm * self.j0[i,:]) * (np.log(self.k[1]) - np.log(self.k[0]))
            else:
                for i in range(self.nsbin):
                    xim[i] = simps(Pm * self.j0[i,:], np.log(self.k))
            return xim
        else:
            print('ERROR: The method %s does not exist! Exiting the code.' %(self.method))
            print('OPTIONS: FFTlinlog, FFTlog')
            sys.exit(1)

        s0, xi0 = self.xicalc(Pkfn, self.sm[0])
        Xint = interp1d(s0, xi0*s0**2, kind='cubic')
        xi = Xint(self.sm) / self.sm**2
        return self.sm, xi * B**2
            
    def xi_model_params(self):
        if(self.model == 'pvoid' or self.model == 'galaxy'):
            print('INFO: Using the model: %s.' %(self.model))
            return ['alpha', 'B', 'Snl']
        if(self.model == 'parab'):
            print('INFO: Using the model: %s.' %(self.model))
            return ['alpha', 'B', 'Snl', 'c']
        
        print('ERROR: The model %s does not exist! Exiting the code.' %(self.model))
        sys.exit(1)
