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
    def __init__(self, root, scan, center=None, sigma=None, qfit=False, expfit=False, emax=None, emin=None, incoherentsig = None, method='powell'):
        """
        Initialize the QENS.
        
        Parameters
        ----------
        sigma : float
            Gaussian sigma parameter
        center : float
            Peak center position
        qfit : bool
            Include quadratic background fitting
        expfit : bool
            Use exponential-Gaussian model instead of standard Gaussian
        emax : double
            Maximum energy offset
        emin : double
            Minimum energy offset
        method : str
            Fitting method for lmfit (default: 'powell')
        incoherentsig : double
            Sigma value for any prefit incoherent sigma
        """
        self.sigma = sigma
        self.center = center
        self.qfit = qfit
        self.expfit = expfit
        self.emax = emax
        self.emin = emin
        self.method = method
        self.root = root
        self.entry = root.entry
        self.scan = self.entry[scan]
        self.directory = Path(
            self.entry['transform'].nxsignal.nxfilename).parent
        self.task_directory = self.directory.parent.parent/'tasks'
        self.sample = self.directory.parent.name
        self.incoherent = incoherentsig


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

 # =====================================================================
    # Utility Methods
    # =====================================================================
    
    @staticmethod
    def ErrorModel(x):
        """
        Complementary error function model.
        
        Parameters
        ----------
        x : array-like
            Input values
            
        Returns
        -------
        array-like
            erfc(x)
        """
        return erfc(x)
    
    @staticmethod
    def incoherentmodel(x):
        """
        Complementary error function model.
        
        Parameters
        ----------
        x : array-like
            Input values
            
        Returns
        -------
        array-like
            erfc(x)
        """
        return erfc(x)
    
    @staticmethod
    def convolve(dat, kernel):
        """
        Simple convolution of two arrays with padding.
        
        Pads the data at edges to preserve array length and applies 
        convolution symmetrically.
        
        Parameters
        ----------
        dat : array-like
            Input data array
        kernel : array-like
            Convolution kernel
            
        Returns
        -------
        array-like
            Convolved array with same length as input
        """
        npts = min(len(dat), len(kernel))
        pad = np.ones(npts)
        tmp = np.concatenate((pad * dat[0], dat, pad * dat[-1]))
        out = np.convolve(tmp, kernel, mode='valid')
        noff = int((len(out) - npts) / 2)
        return out[noff:noff + npts]
    
    @staticmethod
    def invexpgauss(x, amplitude=1, center=0, sigma=1.0, gamma=1.0):
        """
        Inverse exponential-Gaussian function.
        
        Parameters
        ----------
        x : array-like
            x values
        amplitude : float
            Peak amplitude
        center : float
            Peak center
        sigma : float
            Gaussian width
        gamma : float
            Exponential decay parameter
            
        Returns
        -------
        array-like
            Function values
        """
        gss = gamma * sigma * sigma
        arg1 = gamma * (x - center + gss / 2.0)
        arg2 = (center + gss - x) / max(tiny, sigma * np.sqrt(2))
        return amplitude * (gamma / 2) * np.exp(arg1) * erfc(arg2)
    
    # =====================================================================
    # Model Creation Methods
    # =====================================================================
    
    def _create_gaussian_model(self):
        """Create Gaussian model based on expfit setting."""
        if self.expfit:
            return ExponentialGaussianModel(prefix='g1_')
        else:
            return GaussianModel(prefix='g1_')
    
    def _create_voigt_model(self):
        """Create Voigt model (Lorentzian for expfit case)."""
        if self.expfit:
            return LorentzianModel(prefix='v1_')
        else:
            return VoigtModel(prefix='v1_')
    
    def _create_quadratic_model(self):
        """Create quadratic background model."""
        return QuadraticModel(prefix='q1_')
    
    # =====================================================================
    # Parameter Setup Methods
    # =====================================================================
    
    def _setup_gaussian_params(self, gauss, y, x):
        """
        Initialize and configure Gaussian model parameters.
        
        Parameters
        ----------
        gauss : lmfit Model
            Gaussian model instance
        y : array-like
            Data to fit
        x : array-like
            x values
            
        Returns
        -------
        lmfit Parameters
            Configured parameters
        """
        params = gauss.guess(y, x)
        params['g1_amplitude'].set(
            value=params['g1_amplitude'] * 0.4, 
            vary=True, 
            min=0
        )
        params['g1_center'].set(value=self.center, vary=True)
        params['g1_sigma'].set(value=self.sigma, vary=False)
        
        if self.expfit:
            params['g1_gamma'].set(min=3, max=6)
        
        return params
    
    def _setup_voigt_params(self, voigt, y_residual, x, params_gauss):
        """
        Initialize and configure Voigt model parameters.
        
        Parameters
        ----------
        voigt : lmfit Model
            Voigt model instance
        y_residual : array-like
            Residual data after removing Gaussian
        x : array-like
            x values
        params_gauss : lmfit Parameters
            Gaussian parameters for reference
            
        Returns
        -------
        lmfit Parameters
            Configured parameters
        """
        params = voigt.guess(y_residual, x)
        params['v1_center'].set(expr='g1_center')
        params['v1_amplitude'].set(min=0)
        
        if not self.expfit:
            params['v1_sigma'].set(expr='g1_sigma')
            params['v1_gamma'].set(expr='', vary=True, min=0)
        
        return params
    
    def _setup_quadratic_params(self, quad, y_residual, x):
        """
        Initialize and configure quadratic background parameters.
        
        Parameters
        ----------
        quad : lmfit Model
            Quadratic model instance
        y_residual : array-like
            Residual data
        x : array-like
            x values
            
        Returns
        -------
        lmfit Parameters
            Configured parameters
        """
        params = quad.guess(y_residual, x)
        params['q1_c'].set(min=0)
        params['q1_a'].set(value=0, max=0, min=-1, vary=False)
        params['q1_b'].set(value=0, vary=False)
        return params
    
    
    # =====================================================================
    # Data Extraction Methods
    # =====================================================================
    
    def _extract_data(self, scan, H, K, L, cuberad):
        """
        Extract data from scan around specified coordinates.
        
        Parameters
        ----------
        scan : array-like
            4D scan data (E, H, K, L)
        H, K, L : float
            Coordinates
        cuberad : int
            Radius of cube around center
            
        Returns
        -------
        tuple
            (extracted data, validity flag)
        """
        center = [L, K, H]
        Emax = self.emax
        Emin = self.emin
        
        data = scan[
            Emin:Emax,
            (center[0]-cuberad):(center[0]+cuberad),
            (center[1]-cuberad):(center[1]+cuberad),
            (center[2]-cuberad):(center[2]+cuberad)
        ].sum((1, 2, 3))
        
        # Validate data quality
        valid = data.data[:20].min() < data.data[-20:].min()
        
        return data, valid
    
    def _extract_xy_data(self, data):
        """
        Extract x and y arrays from data object.
        
        Parameters
        ----------
        data : object
            Data object with E centers and weighted_data methods
            
        Returns
        -------
        tuple
            (x values, y values)
        """
        x = data['E'].centers().nxvalue
        y = data.weighted_data().nxsignal
        
        if self.expfit:
            y = np.flip(y)
            x = np.flip(x)
        
        return x, y
    
    # =====================================================================
    # Main Fitting Methods
    # =====================================================================
    
    def fit_full(self, scan, H, K, L, cuberad):
        """
        Fit complete model (Gaussian + Voigt ± Quadratic).
        
        Parameters
        ----------
        scan : array-like
            4D scan data
        H, K, L : float
            Coordinates
        cuberad : int
            Radius around center
            
        Returns
        -------
        tuple
            (result, data, x, y, valid)
            where result is the lmfit FitResult object
        """
        # Extract data
        data, valid = self._extract_data(scan, H, K, L, cuberad)
        x, y = self._extract_xy_data(data)
        
        # Create models
        gauss = self._create_gaussian_model()
        voigt = self._create_voigt_model()
        
        # Setup Gaussian parameters
        params_gauss = self._setup_gaussian_params(gauss, y, x)
        
        # Setup Voigt parameters
        y_voigt_residual = y - gauss.eval(params=params_gauss, x=x)
        params_voigt = self._setup_voigt_params(voigt, y_voigt_residual, x, params_gauss)
        
        # Setup quadratic if requested
        if self.qfit:
            quad = self._create_quadratic_model()
            y_quad_residual = y - (
                gauss.eval(params=params_gauss, x=x) +
                voigt.eval(params=params_voigt, x=x)
            )
            params_quad = self._setup_quadratic_params(quad, y_quad_residual, x)
        
        # Build composite model
        if self.qfit:
            if self.expfit:
                gauss2 = ExponentialGaussianModel(prefix='g2_')
                params_gauss2 = self._setup_gaussian2_params(gauss2, y, x, params_gauss)
                model = CompositeModel(gauss2, voigt, self.convolve) + gauss + quad
                params = params_gauss + params_gauss2 + params_voigt + params_quad
            else:
                model = gauss + voigt + quad
                params = params_gauss + params_voigt + params_quad
        else:
            model = gauss + voigt
            params = params_gauss + params_voigt

        mask = np.ones_like(x)
        mask = y>0
        
        # Perform fit
        result = model.fit(y, params, x=x, nan_policy='propagate', method=self.method, weights=mask)
        result = model.fit(y, result.params, x=x, nan_policy='propagate', method = 'leastsq', weights=mask)
        
        return result, data, x, y, valid
    
    def fit_gauss(self, scan, H, K, L, cuberad):
        """
        Fit simple model (Gaussian ± Quadratic, no Voigt).
        
        Used for extracting signal without full spectral decomposition.
        
        Parameters
        ----------
        scan : array-like
            4D scan data
        H, K, L : float
            Coordinates
        cuberad : int
            Radius around center
            
        Returns
        -------
        tuple
            (result, data, x, y, valid)
        """
        # Extract data
        data, valid = self._extract_data(scan, H, K, L, cuberad)
        x, y = self._extract_xy_data(data)
        
        
        # Create models
        gauss = self._create_gaussian_model()
        
        # Setup Gaussian parameters
        params_gauss = self._setup_gaussian_params(gauss, y, x)
        
        # Setup quadratic if requested
        if self.qfit:
            quad = self._create_quadratic_model()
            y_quad_residual = y - gauss.eval(params=params_gauss, x=x)
            params_quad = self._setup_quadratic_params(quad, y_quad_residual, x)
            model = gauss + quad
            params = params_gauss + params_quad
        else:
            model = gauss
            params = params_gauss
        
        mask = np.ones_like(x)
        mask = y>0
        
        # Perform fit
        result = model.fit(y, params, x=x, nan_policy='propagate', method=self.method, weights=mask)
        result = model.fit(y, result.params, x=x, nan_policy='propagate', method = 'leastsq', weights=mask)
        
        return result, data, x, y
    
    # =====================================================================
    # Configuration Methods
    # =====================================================================
    
    def set_parameters(self, sigma=None, center=None, qfit=None, expfit=None, 
                       emax=None, emin=None, method=None):
        """
        Update fitter parameters.
        
        Parameters
        ----------
        sigma : float, optional
        center : float, optional
        qfit : bool, optional
        expfit : bool, optional
        emax : int, optional
        emin : int, optional
        method : str, optional
        """
        if sigma is not None:
            self.sigma = sigma
        if center is not None:
            self.center = center
        if qfit is not None:
            self.qfit = qfit
        if expfit is not None:
            self.expfit = expfit
        if emax is not None:
            self.emax = emax
        if emin is not None:
            self.emin = emin
        if method is not None:
            self.method = method

