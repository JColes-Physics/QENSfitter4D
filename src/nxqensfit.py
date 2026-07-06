import numpy as np
import math
import lmfit
from lmfit.models import GaussianModel, VoigtModel, QuadraticModel, ExponentialGaussianModel, LorentzianModel, ExponentialModel, ConstantModel
from lmfit import minimize, Parameters, report_fit, fit_report, Model, CompositeModel
from scipy.special import wofz, erfc
from nexusformat.nexus import NXdata, NXentry, NXfield, nxopen
from lmfit.lineshapes import s2, tiny
from nexusformat.nexus import (NeXusError, NXdata, NXentry, NXfield, NXlink,
                               nxgetconfig, nxopen, nxsetconfig)
from pathlib import Path
import logging


class NXQENS:
    def __init__(self, x,y,mask=None,g1=False,v1=False,v2=False,decay=False,g1frac=0.01,v1frac=0.01,v2frac=0.01,sigma=0.1,center=0,decayval=0,fixsigma=False,fixcenter=False,v2gamm=None,method='Powell'):
        """
        Function for fitting QENS linewidths
        ----------
        x,y : float arrays
            Data x,y values
        mask : float array
            Data to not fit (data points corresponding to 1 in the mask do not factor in the fit)
        g1,v1,v2,decay : bool
            Boolian decisions for including gauss, voigt(s), and exponential decay contributions to fit
        g1frac,v1frac,v2frac : float
            Float values indicating proportion of LMFIT guess amplitudes to force before fitting
        sigma : float
            Gaussian sigma parameter
        center : float
            Peak center position
        decayval : float
            Value used by exponential decay in the event there is an exponential decay contribution
        fixsigma,fixcenter : bool
            State if given sigma or center values should be fixed (True) or allowed to vary (False)
        v2gamm : float
            If not None, will pass v2gamm as fixed variable for second Voigt Gamma parameter
        """
        self.x = x
        self.y = y
        self.mask = mask
        self.g1 = g1
        self.v1 = v1
        self.v2 = v2
        self.decay = decay
        self.g1frac = g1frac
        self.v1frac = v1frac
        self.v2frac = v2frac
        self.decay = decay
        self.sigma = sigma
        self.center = center
        self.decayval = decayval
        self.fixsigma = fixsigma
        self.fixcenter = fixcenter
        self.v2gamm = v2gamm
        self.method = method


        self._logger = None


    @property
    def logger(self):
        """Log file handler."""
        if self._logger is None:
            self._logger = logging.getLogger(
                f"{self.sample}_{self.scan}['entry']")
            self._logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                "%(asctime)s %(name)-12s: %(message)s",
                datefmt='%Y-%m-%d %H:%M:%S')
            fileHandler = logging.FileHandler(
                self.task_directory.joinpath('nxlogger.log'))
            fileHandler.setFormatter(formatter)
            self._logger.addHandler(fileHandler)
        return self._logger

    def fit_data(self):
        """
        Perform the fit using selected parameters and method.
        """
        x = self.x
        y = self.y
        mask = self.mask
        g1 = self.g1
        v1 = self.v1
        v2 = self.v2
        decay = self.decay
        g1frac = self.g1frac
        v1frac = self.v1frac
        v2frac = self.v2frac
        decay = self.decay
        sigma = self.sigma
        center = self.center
        decayval = self.decayval
        fixsigma = self.fixsigma
        fixcenter = self.fixcenter
        v2gamm = self.v2gamm
        method = self.method

        weights = np.ones(mask.shape)
        weights *= np.logical_not(mask)


        # Inititialize QENS Functions according to boolean user inputs
        if decay:
            gauss = ExponentialGaussianModel(prefix='g0_')
            gbkg = ExponentialGaussianModel(prefix='gbkg_')
            if v1:
                voigt1 = LorentzianModel(prefix='v1_')
                gauss1 = ExponentialGaussianModel(prefix='g1_')
            if v2:
                voigt2 = LorentzianModel(prefix='v2_')
                gauss2 = ExponentialGaussianModel(prefix='g2_')
        
        if not decay:
            if g1:
                gauss = GaussianModel(prefix='g0_')
            if v1:
                voigt1 = VoigtModel(prefix='v1_')
            if v2:
                voigt2 = VoigtModel(prefix='v2_')

        residual = y
        # Establish Gaussian initalization parameters, and boundary conditions
        if g1:
            gaussParams = gauss.guess(y,x)
            gaussParams['g0_amplitude'].set(
                min=0,
                max=y.max()*1.2,
                vary=True,
                value=gaussParams['g0_amplitude']*g1frac,
                )
            if gaussParams['g0_amplitude']<0:
                gaussParams['g0_amplitude'].set(value=-1*gaussParams['g0_amplitude'])
            gaussParams['g0_center'].set(value=center,
                                         vary=(not fixcenter))
            gaussParams['g0_sigma'].set(value=sigma,
                                        vary=(not fixsigma))
            if decay:
                gaussParams['g0_gamma'].set(value=decayval)
            residual -= gauss.eval(params=gaussParams,x=x)

                            
        # Establish Voigt 1 initialization parameters, boundary conditions, and relationships with function g0 if initialized
        if v1:
            voigt1Params = voigt1.guess(residual,x)
            voigt1Params['v1_amplitude'].set(
                value=voigt1Params['v1_amplitude']*v1frac,
                vary=True,
                min=0,
                max=y.max()*1.2,
            )
            if voigt1Params['v1_amplitude'].value<0:
                voigt1Params['v1_amplitude'].set(value=-1*voigt1Params['v1_amplitude'])
            residual -= voigt1.eval(params=voigt1Params,x=x)
            if not decay:
                voigt1Params['v1_gamma'].set(expr='',vary=True,min=tiny)
                if not g1:
                    voigt1Params['v1_sigma'].set(value=sigma,
                                                vary=(not fixsigma))
                    voigt1Params['v1_center'].set(value=center,
                                            vary=(not fixcenter))
                if g1:
                    voigt1Params['v1_sigma'].set(expr='g0_sigma')
                    voigt1Params['v1_center'].set(expr='g0_center')

            if decay:
                gauss1Params = gauss1.guess(y,x)
                voigt1Params['v1_center'].set(expr='g1_center')
                if g1:
                    for param in gauss1Params:
                        gauss1Params[param].set(expr='g0_'+param[3:])
                if not g1:
                    gauss1Params['g1_gamma'].set(value=decayval)
                    gauss1Params['g1_amplitude'].set(min=0,
                                                     max=y.max(),
                                                     vary=True,
                                                     value=gauss1Params['g1_amplitude']*g1frac,
                                                     )
                    if gauss1Params['g1_amplitude']<0:
                        gauss1Params['g1_amplitude'].set(value=-1*gauss1Params['g1_amplitude'])
                    gauss1Params['g1_center'].set(value=center,
                                                  vary=(not fixcenter))
                    gauss1Params['g1_sigma'].set(value=sigma,
                                                 vary=(not fixsigma))

                    
        
        # Establish Voigt 2 initialization parameters, boundary conditions, and relationships with functions g0 and v1 if inizialized
        if v2:
            voigt2Params = voigt2.guess(residual,x)
            voigt2Params['v2_amplitude'].set(
                value=voigt2Params['v2_amplitude']*v2frac,
                vary=True,
                min=0,
                max=y.max()*1.2,
                )
            if voigt2Params['v2_amplitude'].value<0:
                voigt2Params['v2_amplitude'].set(value=-1*voigt2Params['v2_amplitude'])
            residual -= voigt2.eval(params=voigt2Params,x=x)
            if not decay:
                if v2gamm == None:
                    voigt2Params['v2_gamma'].set(expr='',vary=True,min=tiny)
                elif v2gamm != None:
                    voigt2Params['v2_gamma'].set(expr='',value=v2gamm,vary=False)
                if g1:
                    voigt2Params['v2_sigma'].set(expr='g0_sigma')
                    voigt2Params['v2_center'].set(expr='g0_center')
                if (not g1) & v1:
                    voigt2Params['v2_sigma'].set(expr='v1_sigma')
                    voigt2Params['v2_center'].set(expr='v1_center')
                if (not g1) & (not v1):
                    voigt2Params['v2_sigma'].set(value=sigma,
                                                vary=(not fixsigma))
                    voigt2Params['v2_center'].set(value=center,
                                    vary=(not fixcenter))
            if decay:
                gauss2Params = gauss2.guess(y,x)
                voigt2Params['v2_center'].set(expr='g2_center')
                if v2gamm != None:
                    voigt2Params['v2_sigma'].set(expr='',value=v2gamm,vary=False)
                if g1:
                    for param in gauss2Params:
                        gauss2Params[param].set(expr='g0_'+param[3:])
                if (not g1) & v1:
                    for param in gauss2Params:
                        gauss2Params[param].set(expr='g1_'+param[3:])
                if (not g1) & (not v1):
                    gauss2Params['g2_gamma'].set(value=decayval)
                    gauss2Params['g2_amplitude'].set(min=0,
                                                     max=y.max(),
                                                     vary=True,
                                                     value=gauss2Params['g2_amplitude']*g1frac,
                                                     )
                    if gauss2Params['g2_amplitude']<0:
                        gauss2Params['g2_amplitude'].set(value=-1*gauss2Params['g2_amplitude'])
                    gauss2Params['g2_center'].set(value=center,
                                                  vary=(not fixcenter))
                    gauss2Params['g2_sigma'].set(value=sigma,
                                                 vary=(not fixsigma))


        
        const = QuadraticModel(prefix='bkg_')
        constParams = const.guess(residual,x)
        constParams['bkg_a'].set(value=0,vary=False)
        constParams['bkg_b'].set(value=0,vary=False)
        if constParams['bkg_c'].value < 0 or constParams['bkg_c'].value > (min(y[y!=0])):
            constParams['bkg_c'].set(value=tiny)
        constParams['bkg_c'].set(
            min=0.00,
            #max=min(min(y[y!=0]),1)
            )


        #initialize model using functional constant background
        if not decay:
            model = const
            params = constParams
        if decay:
            gbkgParams = gbkg.guess(y,x)
            if v2:
                for param in gbkgParams:
                        gbkgParams[param].set(expr='g2_'+param[5:])
            if v1:
                for param in gbkgParams:
                        gbkgParams[param].set(expr='g1_'+param[5:])
            if g1:
                for param in gbkgParams:
                        gbkgParams[param].set(expr='g0_'+param[5:])
            elif not (v2 or v1 or g1):
                gbkgParams['g2_gamma'].set(value=decayval)
                gbkgParams['gbkg_amplitude'].set(min=0,
                                                    max=y.max(),
                                                    vary=True,
                                                    value=gbkgParams['gbkg_amplitude']*g1frac,
                                                    )
                if gbkgParams['gbkg_amplitude']<0:
                    gbkgParams['gbkg_amplitude'].set(value=-1*gbkgParams['gbkg_amplitude'])
                gbkgParams['gbkg_center'].set(value=center,
                                                vary=(not fixcenter))
                gbkgParams['gbkg_sigma'].set(value=sigma,
                                                vary=(not fixsigma))

            model = CompositeModel(gbkg,const,self.convolve)
            params = constParams + gbkgParams


        # Create full model
        if decay:
            if g1:
                model += gauss
                params += gaussParams
            if v1:
                model += CompositeModel(gauss1,voigt1,self.convolve)
                params += gauss1Params 
                params += voigt1Params
            if v2:
                model += CompositeModel(gauss2,voigt2,self.convolve)
                params += gauss2Params 
                params += voigt2Params
        if not decay:
            if g1:
                model += gauss
                params += gaussParams
            if v1:
                model += voigt1
                params += voigt1Params
            if v2:
                model += voigt2
                params += voigt2Params

        # Fit results, once with prefered fitting method, and second time to establish fitting errors and QOF
        result = model.fit(y,params,x=x,nan_policy='propagate',method=method,weights=weights) 
        result = model.fit(y,result.params,x=x,nan_policy='propagate',method='leastsq',weights=weights)
        
        return result
    
    def convolve(self, dat, kernel):
        """simple convolution of two arrays"""
        npts = min(len(dat), len(kernel))
        pad = np.ones(npts)
        tmp = np.concatenate((pad*dat[0], dat, pad*dat[-1]))
        out = np.convolve(tmp, kernel, mode='valid')
        noff = int((len(out) - npts) / 2)
        return (out[noff:])[:npts]
    
    def fit_data_unpickled(self):
        result = self.fit_data()

        if self.g1:
            gauss_amplitude_value = result.params['g0_amplitude'].value
            gauss_amplitude_stderr = result.params['g0_amplitude'].stderr
            sigma_value = result.params['g0_sigma'].value
            sigma_stderr = result.params['g0_sigma'].stderr
            center_value = result.params['g0_center'].value
            center_stderr = result.params['g0_center'].stderr
            if self.decay:
                decay_value = result.params['g0_gamma'].value
                decay_value = result.params['g0_gamma'].stderr
        if not self.g1:
            gauss_amplitude_value = 0
            gauss_amplitude_stderr = 0

        if self.v1:
            v1_amplitude_value = result.params['v1_amplitude'].value
            v1_amplitude_stderr = result.params['v1_amplitude'].stderr
            if not self.decay:
                v1_gamma_value = result.params['v1_gamma'].value
                v1_gamma_stderr = result.params['v1_gamma'].stderr
                if not self.g1:
                    sigma_value = result.params['v1_sigma'].value
                    sigma_stderr = result.params['v1_sigma'].stderr
                    center_value = result.params['v1_center'].value
                    center_stderr = result.params['v1_center'].stderr
            elif self.decay:
                v1_gamma_value = result.params['v1_sigma'].value
                v1_gamma_stderr = result.params['v1_sigma'].stderr
                if not self.g1:
                    sigma_value = result.params['g1_sigma'].value
                    sigma_stderr = result.params['g1_sigma'].stderr
                    center_value = result.params['g1_center'].value
                    center_stderr = result.params['g1_center'].stderr
                    decay_value = result.params['g1_gamma'].value
                    decay_stderr = result.params['g1_gamma'].value
        if not self.v1:
            v1_amplitude_value = 0
            v1_amplitude_stderr = 0
            v1_gamma_value = 0
            v1_gamma_stderr = 0


        if self.v2:
            v2_amplitude_value = result.params['v2_amplitude'].value
            v2_amplitude_stderr = result.params['v2_amplitude'].stderr
            if not self.decay:
                v2_gamma_value = result.params['v2_gamma'].value
                v2_gamma_stderr = result.params['v2_gamma'].stderr
                if (not self.g1) and (not self.v1):
                    sigma_value = result.params['v2_sigma'].value
                    sigma_stderr = result.params['v2_sigma'].stderr
                    center_value = result.params['v2_center'].value
                    center_stderr = result.params['v2_center'].stderr
            elif self.decay:
                v2_gamma_value = result.params['v2_sigma'].value
                v2_gamma_stderr = result.params['v2_sigma'].stderr
                if (not self.g1) and (not self.v1):
                    sigma_value = result.params['g2_sigma'].value
                    sigma_stderr = result.params['g2_sigma'].stderr
                    center_value = result.params['g2_center'].value
                    center_stderr = result.params['g2_center'].stderr
                    decay_value = result.params['g2_gamma'].value
                    decay_stderr = result.params['g2_gamma'].value
        if not self.v2:
            v2_amplitude_value = 0
            v2_amplitude_stderr = 0
            v2_gamma_value = 0
            v2_gamma_stderr = 0

        if not self.decay:
            decay_value = 0
            decay_stderr = 0

        bkg_c_value = result.params['bkg_c'].value
        bkg_c_stderr = result.params['bkg_c'].stderr

        output = (
            gauss_amplitude_value,
            gauss_amplitude_stderr,
            v1_amplitude_value,
            v1_amplitude_stderr,
            v1_gamma_value,
            v1_gamma_stderr,
            v2_amplitude_value,
            v2_amplitude_stderr,
            v2_gamma_value,
            v2_gamma_stderr,
            decay_value,
            decay_stderr,
            bkg_c_value,
            bkg_c_stderr,
            sigma_value,
            sigma_stderr,
            center_value,
            center_stderr,
            result.rsquared,
            result.redchi,
        )
        return output
