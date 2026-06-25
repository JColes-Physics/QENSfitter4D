import numpy as np
from nexusformat.nexus import NeXusError, NXdata, NXfield
from nexusformat.nexus.tree import centers

from nexpy.gui.utils import display_message, format_float, load_models, report_error
from nexpy.gui.widgets import (GridParameters, NXDialog, NXCheckBox, NXColorBox, NXComboBox, NXLabel,
                      NXLineEdit, NXMessageBox, NXPanel, NXPushButton,
                      NXrectangle, NXScrollArea, NXTab, NXDoubleSpinBox)
from nexpy.gui.pyqt import QtCore, QtWidgets, QtGui
from nexpy.gui.plotview import NXPlotView, linestyles

from nexusformat.nexus import (NeXusError, NXdata, NXentry, NXfield, NXlink, NXprocess,
                               nxgetconfig, nxopen, nxsetconfig)

from lmfit.models import GaussianModel, VoigtModel, QuadraticModel, ExponentialGaussianModel, LorentzianModel, ConstantModel
from lmfit import minimize, Parameters, report_fit, fit_report, Model, CompositeModel
from lmfit.lineshapes import s2, tiny
from lmfit import __version__ as lmfit_version

import inspect
import re
from itertools import cycle

def get_models():
    """
    Return a dictionary of LMFIT models.

    This function returns a dictionary of LMFIT models, including those
    defined in the LMFIT package and those defined in the
    ``nexpy.models`` package. Additional models can also be defined in
    the ``~/.nexpy/models`` directory or in another installed package,
    which declares the entry point ``nexpy.models``. The models are
    returned as a dictionary where the keys are the names of the models
    and the values are the classes defining the models.
    """
    from lmfit.models import lmfit_models
    models = lmfit_models
    if 'Expression' in models:
        del models['Expression']
    if 'Gaussian-2D' in models:
        del models['Gaussian-2D']

    nexpy_models = load_models()

    for model in nexpy_models:
        try:
            models.update(
                dict((n.strip('Model'), m)
                for n, m in inspect.getmembers(nexpy_models[model],
                                               inspect.isclass)
                if issubclass(m, Model) and n != 'Model'))
        except ImportError:
            pass

    return models


all_models = get_models()

def get_methods():
    """Return a dictionary of minimization methods in LMFIT."""
    methods = {'leastsq': 'Levenberg-Marquardt',
               'least_squares': 'Least-Squares minimization, '
                                'using Trust Region Reflective method',
               'differential_evolution': 'differential evolution',
               'nelder': 'Nelder-Mead',
               'lbfgsb': ' L-BFGS-B',
               'powell': 'Powell',
               'cg': 'Conjugate-Gradient',
               'newton': 'Newton-CG',
               'cobyla': 'Cobyla',
               'bfgs': 'BFGS',
               'tnc': 'Truncated Newton',
               'trust-ncg': 'Newton-CG trust-region',
               'trust-exact': 'nearly exact trust-region',
               'trust-krylov': 'Newton GLTR trust-region',
               'trust-constr': 'trust-region for constrained optimization',
               'dogleg': 'Dog-leg trust-region',
               'slsqp': 'Sequential Linear Squares Programming'}
    return methods


all_methods = get_methods()

def show_dialog(parent=None):
    """Entry point called when menu item is clicked."""
    try:
        dialog = initparamswindow()
        dialog.show()
    except NeXusError as error:
        report_error("Pameter Initialization", error)


class initparamswindow(NXDialog):
    """Create a dialog window inheriting from NeXpy's NXDialog."""
    
    def __init__(self, parent=None):
        super(initparamswindow, self).__init__(parent)
        
        # Select the NeXus entry to work with
        self.select_data()

        # Create parameters with a grid layout
        self.parameters = GridParameters()

        self.parameters.add('H',0.,'H')
        self.parameters.add('K',0.,'K')
        self.parameters.add('L',0.,'L')

        # self.parameters.add('hkl',(0.0,0.0,0.0),'HKL')
        
        self.parameters.add('dQ', .1, 'Q Box Diameter')
        self.parameters.add('Emin', -1., 'Minimum Energy (meV)')
        self.parameters.add('Emax', 1., 'Maximum Energy (meV)')

        self.parameters2 = GridParameters()
        self.parameters2.add('sigma',0.01, 'Resolution width (σ)',False)
        self.parameters2.add('center', 0.0, 'Expected Center (meV)',False)

        # Create conditional parameters for energy range cutting        
        self.energy_cut_params = GridParameters()        
        self.energy_cut_params.add('energy_cut_min', -0.0, 'Energy Cut Min (meV)')        
        self.energy_cut_params.add('energy_cut_max', 0.0, 'Energy Cut Max (meV)')

        layout_list = [self.entry_layout,
                       self.parameters.grid(),
                       self.parameters2.grid()]


        layout_list.append(self.checkboxes(
            ("fit_resolution_checkbox", "Fit Only Resolution", False),
            ("fit_gaussian_checkbox", "Use Exponential Gaussian", False),
            ("cut_energy_checkbox", "Cut Energy Band", False),
            ("fit_background_checkbox","Fit constant background?", False)))
        layout_list.append(self.checkboxes(
            ("over_plot_checkbox", "Overplot", False),
            ("comp_plot_checkbox","Plot Components", False)))

        self.energy_cut_layout = self.energy_cut_params.grid()
        layout_list.append(self.energy_cut_layout)

        # Create array of command buttons
        layout_list.append(self.action_buttons(('Reset', self.reset_dialog),
                                            ('Run Test Fit',self.open_fit_window),
                                            ('Plot Data', self.plot_data),
                                            ('Plot Fit', self.plot_fit),
                                            ('Save', self.save_data)))
        layout_list.append(self.close_buttons())


        self.set_layout(*layout_list)
        self.setWindowTitle('Setting Inititial Conditions')
        self.reset_dialog

        #self.set_layout(self.Qgrid, self.close_buttons(save=True))
   
    def on_cut_energy_changed(self):
        """Handle cut energy checkbox state change - toggle visibility of energy cut parameters."""
        is_checked = self.cut_energy_checkbox.isChecked()
        self.energy_cut_layout.parent().setVisible(is_checked)

    @property
    def H(self):
        return self.parameters['H'].value
    
    @property
    def K(self):
        return self.parameters['K'].value
    
    @property
    def L(self):
        return self.parameters['L'].value
    
    @property
    def dQ(self):
        return self.parameters['dQ'].value
    
    @property
    def Emin(self):
        return self.parameters['Emin'].value

    @property
    def Emax(self):
        return self.parameters['Emax'].value
    @property
    def sigma(self):
        return abs(self.parameters2['sigma'].value)
    
    @property
    def center(self):
        return self.parameters2['center'].value



    @property
    def fit_resolution_only(self):
        return self.checkbox['fit_resolution_checkbox'].isChecked()
    
    @property
    def fit_with_gaussian(self):
        return self.checkbox['fit_gaussian_checkbox'].isChecked()
    
    @property
    def cut_energy(self):
        return self.checkbox['cut_energy_checkbox'].isChecked()
    
    @property
    def fit_background(self):
        return self.checkbox['fit_background_checkbox'].isChecked()
    


    @property
    def over_plot(self):
        return self.checkbox['over_plot_checkbox'].isChecked()
    
    @property
    def comp_plot(self):
        return self.checkbox['comp_plot_checkbox'].isChecked()
    


    @property
    def energy_cut_min(self):
        return self.energy_cut_params['energy_cut_min'].value
    
    @property
    def energy_cut_max(self):
        return self.energy_cut_params['energy_cut_max'].value


    def open_fit_window(self):
        #dialog2 = FitTab(self.selected_data.title.nxvalue,self.selected_data)
        dialog2 = SpectralFittingWidget(nxdata=self.selected_data,nxentry=self.entry)
        dialog2.show()

    def plot_data(self):
        if self.over_plot:
            self.convert_QE().oplot()
        if not self.over_plot:
            self.convert_QE().plot()

    def plot_fit(self):
        if self.over_plot:
            self.convert_QE().oplot()
        if not self.over_plot:
            self.convert_QE().plot()
        self.fit_QENS().oplot()

    def convert_QE(self):
        return self.selected_data[self.Emin:self.Emax,
                         self.L-(self.dQ/2):self.L+(self.dQ/2),
                         self.K-(self.dQ/2):self.K+(self.dQ/2),
                         self.H-(self.dQ/2):self.H+(self.dQ/2)].sum((1,2,3))

    def test_fit(self):
        """Convert S(phi,eps) to S(Q,eps)"""
        data = self.convert_QE
        x = data.E.centers().nxvalue
        try:
            y = data.weighted_data().nxsignal
        except:
            y = data.nxsignal.nxvalue
        sigma = self.sigma
        cent = self.center
        resolution = self.fit_resolution_only
        expfit = self.fit_with_gaussian
        qfit = self.fit_background
        cutbool = self.cut_energy
        method = 'powell'

        mask = y!=0
        y = y[mask]
        x = x[mask]
        weights = weights[mask]
        
        
        if len(y) == 0:
            print(len(y))
            return 0,data,x,y,False

        if expfit:
            gauss = ExponentialGaussianModel(prefix='g1_')
            y = np.flip(y)
            x = np.flip(x)
            if qfit:
                quad = QuadraticModel(prefix='q1_')
        if not expfit:
            gauss = GaussianModel(prefix='g1_')
            if qfit:
                quad = QuadraticModel(prefix='q1_') 
        if not resolution:
            voigt = VoigtModel(prefix='v1_')



        paramsGauss = gauss.guess(y, x)
        paramsGauss['g1_amplitude'].set(value=paramsGauss['g1_amplitude']*.4,vary=True)
        paramsVoigt = voigt.guess(y - gauss.eval(params=paramsGauss, x=x), x)

        if qfit:
            paramsQuad = quad.guess(y - (gauss.eval(params=paramsGauss, x=x) +
                                        voigt.eval(params=paramsVoigt, x=x)), x)
            paramsQuad['q1_c'].set(min=0)#,max=y.min()) 
            if not expfit:
                paramsQuad['q1_a'].set(value=0, max=0, min=-1, vary=False)
                paramsQuad['q1_b'].set(value=0,vary=False)     
        
        paramsGauss['g1_center'].set(value=cent, vary=True)#, min = -0.1, max = 0.1)
        paramsGauss['g1_sigma'].set(value=sigma, vary=False)
        paramsGauss['g1_amplitude'].set(min=0)
        if not resolution:
            paramsVoigt['v1_center'].set(expr='g1_center')
            paramsVoigt['v1_amplitude'].set(min=0)
            paramsVoigt['v1_sigma'].set(expr='g1_sigma')
            paramsVoigt['v1_gamma'].set(expr='', vary=True, min=0)

        if qfit:
            model = gauss + voigt + quad
            params = paramsGauss + paramsVoigt + paramsQuad
        else:
            model = gauss + voigt
            params = paramsGauss + paramsVoigt

        result = model.fit(y, params, x=x, nan_policy='propagate', method=method,weights=cmask)
        result = model.fit(y, result.params, x=x, nan_policy='propagate', method='leastsq',weights=cmask)

        return result

    def fit_QE(self):
        self.fitresult = self.convert_QE()
    
    def save_data(self):
        self.entry['sqe'] = self.convert_QE()

    def reset_dialog(self):
        try:
            self.init_condition = self.entry.QENSfit_conditions

            self.parameters['dQ'].set((self.init_condition['dQ']))        
            self.parameters['Emin'].set((self.init_condition['Emin']))        
            self.parameters['Emax'].set((self.init_condition['Emax']))
            self.parameters2['sigma'].set((self.init_condition['sigma']))
            self.parameters2['center'].set((self.init_condition['center']))              
            
            self.checkbox['fit_resolution_checkbox'].setChecked(bool(self.init_condition['fit_resolution_only']))        
            self.checkbox['fit_gaussian_checkbox'].setChecked(bool(self.init_condition['fit_exponential_gaussian']))        
            self.checkbox['cut_energy_checkbox'].setChecked(bool(self.init_condition['cut_energy']))                
            
            self.energy_cut_params['energy_cut_min'].set((self.init_condition['energy_cut_min']))        
            self.energy_cut_params['energy_cut_max'].set((self.init_condition['energy_cut_max']))
        except:
            self.parameters['dQ'].set(0.01)        
            self.parameters['Emin'].set(-0.5)        
            self.parameters['Emax'].set(0.5) 
            self.parameters2['sigma'].set(0.01)
            self.parameters2['center'].set(0.0)               
            
            self.checkbox['fit_resolution_checkbox'].setChecked(False)        
            self.checkbox['fit_gaussian_checkbox'].setChecked(False)        
            self.checkbox['cut_energy_checkbox'].setChecked(False)                
            
            self.energy_cut_params['energy_cut_min'].set(-0.0)        
            self.energy_cut_params['energy_cut_max'].set(0.0)
        
    def save_data(self):
        if 'QENSfit_conditions' in self.entry:
            del self.entry['QENSfit_conditions']
        self.entry['QENSfit_conditions'] = NXprocess(dQ=self.dQ,
                                                     Emax=self.Emax,
                                                     Emin=self.Emin,
                                                     sigma=self.sigma,
                                                     center=self.center,
                                                     fit_resolution_only=self.fit_resolution_only,
                                                     fit_exponential_gaussian=self.fit_with_gaussian,
                                                     cut_energy=self.cut_energy,
                                                     energy_cut_min=self.energy_cut_min,
                                                     energy_cut_max=self.energy_cut_max,
                                                     )
        #self.logger.info(f"'QENSfit_conditions' added to entry")


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
    
    def get_parent_subentries(self):
        try:
            return self.entry.NXdata
        except Exception:
            pass
        return[]




"""
Spectral Fitting GUI Widget for NeXpy

A GUI widget for fitting spectral data using lmfit with support for:
- Voigt profiles (1 or 2)
- Exponential decay
- Energy spectra flipping
- Interactive parameter adjustment
- Multiple fitting methods
- Data range cutting
"""

class SpectralFittingWidget(NXDialog):
    """
    A GUI widget for fitting spectral data with interactive parameter adjustment.
    
    Attributes:
        nxdata: NXdata object to be fitted
        fit_results: Dictionary storing fit results
        parameters: lmfit Parameters object
        model: lmfit Model object
    """
    
    def __init__(self, parent=None, nxdata=None, nxentry=None):
        """
        Initialize the spectral fitting widget.
        
        Parameters
        ----------
        parent : QtWidgets.QWidget, optional
            Parent widget
        nxdata : NXdata, optional
            Initial NXdata object for fitting
        """
        super().__init__(parent)
        self.nxdata = nxdata
        self.nxentry = nxentry
        self.fit_results = None
        self.parameters = Parameters()
        self.model = None
        self.data_mask = None
        
        self.setWindowTitle('Spectral Fitting Tool')
        self.init_ui()
        self.reset_all()
    
    # ==================== UI INITIALIZATION ====================
    
    def init_ui(self):
        """Initialize the main user interface."""
        main_layout = QtWidgets.QVBoxLayout()
        
        # Add all sections in order
        main_layout.addWidget(self.create_input_section())
        main_layout.addWidget(self.create_separator())
        
        main_layout.addWidget(self.create_fit_options_section())
        main_layout.addWidget(self.create_separator())
        
        main_layout.addWidget(self.create_parameters_section())
        main_layout.addWidget(self.create_separator())
        
        main_layout.addWidget(self.create_fitting_method_section())
        main_layout.addWidget(self.create_separator())
        
        main_layout.addWidget(self.create_results_section())
        main_layout.addWidget(self.create_separator())
        
        main_layout.addWidget(self.create_action_buttons_section())
        
        # Add stretch at the end to compress layout
        main_layout.addStretch()
        
        self.setLayout(main_layout)
        
        # Set reasonable window size
        self.resize(900, 1200)
    
    # ==================== SECTION 1: INPUT PARAMETERS ====================
    
    def create_input_section(self):
        """Create section for H, K, L, dQ, dE input."""
        group = QtWidgets.QGroupBox('Input Parameters')
        layout = QtWidgets.QGridLayout()
        
        # Input labels and spin boxes
        self.h_input = NXDoubleSpinBox()
        self.k_input = NXDoubleSpinBox()
        self.l_input = NXDoubleSpinBox()
        self.dq_input = NXDoubleSpinBox()
        self.emin_input = NXDoubleSpinBox()
        self.emax_input = NXDoubleSpinBox()

        self.h_input.setSingleStep(0.1)
        self.k_input.setSingleStep(0.1)
        self.l_input.setSingleStep(0.1)
        self.dq_input.setSingleStep(0.1)
        self.emin_input.setSingleStep(0.1)
        self.emax_input.setSingleStep(0.1)

        
        layout.addWidget(QtWidgets.QLabel('H:'), 0, 0)
        layout.addWidget(self.h_input, 0, 1)
        
        layout.addWidget(QtWidgets.QLabel('K:'), 0, 2)
        layout.addWidget(self.k_input, 0, 3)
        
        layout.addWidget(QtWidgets.QLabel('L:'), 0, 4)
        layout.addWidget(self.l_input, 0, 5)
        
        layout.addWidget(QtWidgets.QLabel('dQ (Å⁻¹):'), 0, 6)
        layout.addWidget(self.dq_input, 0, 7)
        
        layout.addWidget(QtWidgets.QLabel('Emin (meV):'), 1, 0)
        layout.addWidget(self.emin_input, 1, 1)

        layout.addWidget(QtWidgets.QLabel('Emax (meV):'), 1, 2)
        layout.addWidget(self.emax_input, 1, 3)
        
        group.setLayout(layout)
        return group
    
    # ==================== SECTION 2: FIT OPTIONS ====================
    
    def create_fit_options_section(self):
        """Create section for selecting fit models and options."""
        group = QtWidgets.QGroupBox('Fit Options')
        layout = QtWidgets.QVBoxLayout()
        
        # Voigt and decay options
        options_layout = QtWidgets.QGridLayout()
        
        self.elastic_cb = QtWidgets.QCheckBox('Elastic Gaussian')
        self.voigt_1_cb = QtWidgets.QCheckBox('1 Voigt')
        self.voigt_2_cb = QtWidgets.QCheckBox('2 Voigt')
        self.decay_cb = QtWidgets.QCheckBox('Exponential Decay')
        self.flip_cb = QtWidgets.QCheckBox('Flip Energy Spectra')
        
        options_layout.addWidget(self.elastic_cb, 0, 0)
        options_layout.addWidget(self.voigt_1_cb, 0, 1)
        options_layout.addWidget(self.voigt_2_cb, 0, 2)
        options_layout.addWidget(self.decay_cb, 0, 3)
        options_layout.addWidget(self.flip_cb, 0, 4)
        
        layout.addLayout(options_layout)
        
        # Two Voigt model selection
        two_voigt_layout = QtWidgets.QHBoxLayout()
        two_voigt_layout.addWidget(QtWidgets.QLabel('2 Voigt Model:'))
        
        self.voigt_model_combo = QtWidgets.QComboBox()
        self.voigt_model_combo.addItems(['None (Fit Independently)', 'Load Model from File'])
        two_voigt_layout.addWidget(self.voigt_model_combo)
        
        self.load_model_btn = QtWidgets.QPushButton('Load Model')
        self.load_model_btn.clicked.connect(self.load_model_from_file)
        two_voigt_layout.addWidget(self.load_model_btn)
        
        layout.addLayout(two_voigt_layout)
        
        # Enable/disable model selection based on 2 Voigt checkbox
        self.voigt_2_cb.stateChanged.connect(
            lambda: self.voigt_model_combo.setEnabled(self.voigt_2_cb.isChecked())
        )
        self.voigt_2_cb.stateChanged.connect(
            lambda: self.load_model_btn.setEnabled(self.voigt_2_cb.isChecked())
        )
        self.voigt_model_combo.setEnabled(False)
        self.load_model_btn.setEnabled(False)
        
        group.setLayout(layout)
        return group
    
    def load_model_from_file():
        pass
    
    # ==================== SECTION 3: PARAMETERS ADJUSTMENT ====================
    
    def create_parameters_section(self):
        """Create section for displaying and adjusting lmfit parameters."""
        group = QtWidgets.QGroupBox('Fit Parameters')

        self.sigma_input = NXDoubleSpinBox()
        self.fix_sigma_cb = QtWidgets.QCheckBox()
        self.center_input = NXDoubleSpinBox()
        self.fix_center_cb = QtWidgets.QCheckBox()
        self.decay_input = NXDoubleSpinBox()

        self.sigma_input.setDecimals(8)
        self.center_input.setDecimals(8)
        self.decay_input.setDecimals(8)

        myFont=QtGui.QFont()
        myFont.setBold(True)

        layout = QtWidgets.QGridLayout()

        layout.addWidget(QtWidgets.QLabel('Variable'), 0, 0)
        layout.addWidget(QtWidgets.QLabel('Value'), 0, 1)
        layout.addWidget(QtWidgets.QLabel('Fix?'), 0, 2)
        
        layout.addWidget(QtWidgets.QLabel('σ (meV)'), 1, 0)
        layout.addWidget(self.sigma_input, 1, 1)
        layout.addWidget(self.fix_sigma_cb, 1, 2)
        
        layout.addWidget(QtWidgets.QLabel('Center (meV)'), 2, 0)
        layout.addWidget(self.center_input, 2, 1)
        layout.addWidget(self.fix_center_cb, 2, 2)
        
        layout.addWidget(QtWidgets.QLabel('decay constant'), 3, 0)
        layout.addWidget(self.decay_input, 3, 1)
        self.decay_input.setEnabled(self.decay_cb.isChecked())
        
        group.setLayout(layout)

        self.decay_cb.stateChanged.connect(
            lambda: self.decay_input.setEnabled(self.decay_cb.isChecked())
        )

        return group
    # def create_parameters_section(self):
    #     """Create section for displaying and adjusting lmfit parameters."""
    #     group = QtWidgets.QGroupBox('Fit Parameters')
    #     layout = QtWidgets.QVBoxLayout()
        
    #     # Create a scroll area for parameters
    #     scroll = QtWidgets.QScrollArea()
    #     scroll.setWidgetResizable(True)
        
    #     self.param_widget = QtWidgets.QWidget()
    #     self.param_layout = QtWidgets.QGridLayout()
    #     self.param_widget.setLayout(self.param_layout)
    #     scroll.setWidget(self.param_widget)
        
    #     layout.addWidget(scroll)
    #     layout.setContentsMargins(0, 0, 0, 0)
        
    #     group.setLayout(layout)
    #     return group
    
    # def update_parameter_display(self):
    #     """Update the parameter adjustment widgets based on current model."""
    #     # Clear existing widgets
    #     while self.param_layout.count():
    #         self.param_layout.takeAt(0).widget().deleteLater()
        
    #     # Add parameter controls for key fitting parameters
    #     # This will be populated after fit initialization
    #     self.param_controls = {}
        
    #     # Example structure (to be populated dynamically):
    #     row = 0
    #     for param_name in ['mu', 'sigma', 'gamma', 'amplitude', 'decay']:
    #         label = QtWidgets.QLabel(param_name)
    #         value_input = NXDoubleSpinBox()
            
    #         self.param_layout.addWidget(label, row, 0)
    #         self.param_layout.addWidget(value_input, row, 1)
            
    #         self.param_controls[param_name] = value_input
    #         row += 1
    
    # ==================== SECTION 4: FITTING METHOD & DATA CUTTING ====================
    
    def create_fitting_method_section(self):
        """Create section for selecting fit method and data range options."""
        group = QtWidgets.QGroupBox('Fitting Options')
        layout = QtWidgets.QVBoxLayout()
        
        # Fitting method
        method_layout = QtWidgets.QHBoxLayout()
        method_layout.addWidget(QtWidgets.QLabel('Fit Method:'))
        
        self.method_combo = QtWidgets.QComboBox()
        self.method_combo.addItems([
            'powell','leastsq', 'least_squares', 'differential_evolution',
            'brute', 'basinhopping', 'ampgo', 'nelder',
        ])
        self.method_combo.setCurrentText('leastsq')
        method_layout.addWidget(self.method_combo)
        
        layout.addLayout(method_layout)
        
        # Data cutting
        cutting_layout = QtWidgets.QHBoxLayout()
        
        self.cut_data_cb = QtWidgets.QCheckBox('Cut Data Range')
        cutting_layout.addWidget(self.cut_data_cb)
        
        # Range input (hidden by default)
        self.cut_min_label = QtWidgets.QLabel('Min:')
        self.cut_min_input = NXDoubleSpinBox()
        self.cut_max_label = QtWidgets.QLabel('Max:')
        self.cut_max_input = NXDoubleSpinBox()
        
        cutting_layout.addWidget(self.cut_min_label)
        cutting_layout.addWidget(self.cut_min_input)
        cutting_layout.addWidget(self.cut_max_label)
        cutting_layout.addWidget(self.cut_max_input)
        
        # Hide range inputs initially
        self.cut_min_label.hide()
        self.cut_min_input.hide()
        self.cut_max_label.hide()
        self.cut_max_input.hide()
        
        # Connect checkbox to show/hide range inputs
        self.cut_data_cb.stateChanged.connect(self.toggle_cut_range)
        
        layout.addLayout(cutting_layout)
        
        group.setLayout(layout)
        return group
    
    def toggle_cut_range(self):
        """Show/hide data range cut inputs based on checkbox state."""
        is_checked = self.cut_data_cb.isChecked()
        self.cut_min_label.setVisible(is_checked)
        self.cut_min_input.setVisible(is_checked)
        self.cut_max_label.setVisible(is_checked)
        self.cut_max_input.setVisible(is_checked)
    
    # ==================== SECTION 5: RESULTS DISPLAY ====================
    
    def create_results_section(self):
        """Create section for displaying fitting results."""
        group = QtWidgets.QGroupBox('Fit Results')
        layout = QtWidgets.QVBoxLayout()
        
        # Results text area
        self.results_text = QtWidgets.QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setMaximumHeight(200)
        #self.results_text.setFont(QtGui.QFont('Courier', 9))
        
        layout.addWidget(self.results_text)
        
        group.setLayout(layout)
        return group
    
    def update_results_display(self, fit_result):
        """
        Update the results text area with fitting results.
        
        Parameters
        ----------
        fit_result : lmfit.model.ModelResult
            Result object from lmfit fit
        """
        if fit_result is None:
            self.results_text.setText('No fit results available.')
            return
        
        result_text = f"""
Fit Summary
===========
Chi-square: {fit_result.chisqr:.6e}
Reduced Chi-square: {fit_result.redchi:.6e}
R-Squared: {fit_result.rsquared:.6e}
AIC: {fit_result.aic:.2f}
BIC: {fit_result.bic:.2f}
Success: {fit_result.success}
Message: {fit_result.message}

Parameters
==========
"""
        # for param_name, param in fit_result.params.items():
        #     result_text += f"{param_name}: {param.value:.6e} ± {param.stderr or 0:.6e}\n"
        result_text += fit_report(fit_result.params)[14:]
        self.results_text.setText(result_text)
        # self.results_text.setText(report_fit(fit_result.params))
    
    # ==================== SECTION 6: ACTION BUTTONS ====================
    
    def create_action_buttons_section(self):
        """Create section with action buttons for fitting and plotting."""
        group = QtWidgets.QGroupBox('Actions')
        layout = QtWidgets.QVBoxLayout()
        
        # Top row of buttons
        button_row1 = QtWidgets.QHBoxLayout()
        
        self.fit_btn = QtWidgets.QPushButton('Fit Data')
        self.fit_btn.clicked.connect(self.fit_data)
        button_row1.addWidget(self.fit_btn)
        
        self.plot_data_btn = QtWidgets.QPushButton('Plot Data')
        self.plot_data_btn.clicked.connect(self.plot_data)
        button_row1.addWidget(self.plot_data_btn)
        
        layout.addLayout(button_row1)
        
        # Second row: checkboxes
        button_row2 = QtWidgets.QHBoxLayout()
        
        self.plot_components_cb = QtWidgets.QCheckBox('Plot Components')
        self.plot_components_cb.setChecked(False)
        button_row2.addWidget(self.plot_components_cb)

        self.plot_cut_data_cb = QtWidgets.QCheckBox('Plot Masked Data')
        self.plot_components_cb.setChecked(False)
        button_row2.addWidget(self.plot_cut_data_cb)
        
        self.overplot_cb = QtWidgets.QCheckBox('Overplot')
        self.overplot_cb.setChecked(False)
        button_row2.addWidget(self.overplot_cb)
        
        layout.addLayout(button_row2)
        
        # Third row of buttons
        button_row3 = QtWidgets.QHBoxLayout()
        
        self.plot_fit_btn = QtWidgets.QPushButton('Plot Fit')
        self.plot_fit_btn.clicked.connect(self.plot_fit)
        button_row3.addWidget(self.plot_fit_btn)
        
        self.reset_btn = QtWidgets.QPushButton('Reset')
        self.reset_btn.clicked.connect(self.reset_all)
        button_row3.addWidget(self.reset_btn)
        
        self.save_btn = QtWidgets.QPushButton('Save Inititial Params')
        self.save_btn.clicked.connect(self.save_init_params)
        button_row3.addWidget(self.save_btn)
        
        layout.addLayout(button_row3)
        
        group.setLayout(layout)

        return group
    
    def plot_data(self):
        data_to_plot = self.select_data()
        if self.plot_cut_data_cb.isChecked():
            data = data_to_plot
            mask = (data[data.axes].nxvalue>self.cut_min_input.value()) & (data[data.axes].nxvalue<self.cut_max_input.value())
            data_to_plot = NXdata(signal=NXfield(np.ma.array(data.nxsignal,mask=mask),name='data'),axes=data[data.axes])
            data_to_plot[data_to_plot.signal].rename(f'[{self.h_input.value()},{self.k_input.value()},{self.l_input.value()}]_cut')
        if self.overplot_cb.isChecked():
            data_to_plot.oplot()
            # plotview.legend()
        if not self.overplot_cb.isChecked():
            data_to_plot.plot()
        pass

    def select_data(self):
        data = self.nxdata[self.emin_input.value():self.emax_input.value(),
                         self.l_input.value()-(self.dq_input.value()/2):self.l_input.value()+(self.dq_input.value()/2),
                         self.k_input.value()-(self.dq_input.value()/2):self.k_input.value()+(self.dq_input.value()/2),
                         self.h_input.value()-(self.dq_input.value()/2):self.h_input.value()+(self.dq_input.value()/2)].sum((1,2,3))
        if data.nxsignal.shape != data[data.axes].shape:
            data = NXdata(data.nxsignal,data[data.axes].centers())
        
        if self.flip_cb.isChecked():
            axes = data[data.axes]
            diff = abs(axes.max())-abs(axes.min())
            data = NXdata(np.flip(data.nxsignal),axes=(data[data.axes]-diff))

        data[data.signal].rename(f'[{self.h_input.value()},{self.k_input.value()},{self.l_input.value()}]')
        return data

    def plot_fit(self):
        data = self.select_data()
        result = self.fit_result
        E = data[data.axes]
        x = np.linspace(E.min(),E.max(),len(E)*4)
        resultgraph = result.best_fit
        components = result.eval_components(x=x)
        nxout = NXdata(signal=resultgraph,axes=self.fit_x,name=f'fit of [{self.h_input.value()},{self.k_input.value()},{self.l_input.value()}]')
        nxout.oplot(linestyle='-',marker=None,color='blue')
        if self.plot_components_cb.isChecked():
            nxout1 = NXentry()
            for model_name, model_value in components.items():
                nxout1[model_name]=NXdata(signal=NXfield(model_value,name=model_name),axes=x,name=model_name)
                if 'g0_' in model_name:
                    nxout1[model_name].oplot(linestyle='--',marker=None,color='green')
                if 'v1_' in model_name:
                    nxout1[model_name].oplot(linestyle='dashdot',marker=None,color='cyan')
                if 'v2_' in model_name:
                    nxout1[model_name].oplot(linestyle='dashdot',marker=None,color='purple')
                if 'bkg' in model_name:
                    nxout1[model_name].oplot(linestyle='dotted',marker=None,color='orange')    
        pass

    def reset_all(self):
        try:
            self.init_condition = self.nxentry['QENSfit_conditions']

            self.dq_input.setValue(self.init_condition['dQ'].nxvalue)       
            self.emax_input.setValue(self.init_condition['Emax'].nxvalue)        
            self.emin_input.setValue(self.init_condition['Emin'].nxvalue)
            self.sigma_input.setValue(self.init_condition['sigma'].nxvalue)
            self.center_input.setValue((self.init_condition['center'].nxvalue))

            self.fix_sigma_cb.setChecked(bool(self.init_condition['sigma_fixed'].nxvalue))
            self.fix_center_cb.setChecked(bool(self.init_condition['center_fixed'].nxvalue))              
            
            self.elastic_cb.setChecked(bool(self.init_condition['gaussian_fit'].nxvalue))        
            self.voigt_1_cb.setChecked(bool(self.init_condition['voigt_1_fit'].nxvalue))        
            self.voigt_2_cb.setChecked(bool(self.init_condition['voigt_2_fit'].nxvalue))
            self.decay_cb.setChecked(bool(self.init_condition['fit_exponential_gauss'].nxvalue)) 

            self.flip_cb.setChecked(bool(self.init_condition['flip_energy_spectra'].nxvalue))               
            
            self.cut_min_input.setValue((self.init_condition['energy_cut_min'].nxvalue))        
            self.cut_max_input.setValue((self.init_condition['energy_cut_max'].nxvalue))
            self.method_combo.setCurrentText(self.init_condition['method'].nxvalue)
        except:
            self.dq_input.setValue(0.1)       
            self.emax_input.setValue(1.)        
            self.emin_input.setValue(-1.)
            self.sigma_input.setValue(0.1)
            self.center_input.setValue(0.)

            self.fix_sigma_cb.setChecked(False)
            self.fix_center_cb.setChecked(False)              
            
            self.elastic_cb.setChecked(False)        
            self.voigt_1_cb.setChecked(False)        
            self.voigt_2_cb.setChecked(False)
            self.decay_cb.setChecked(False)
            self.method_combo.setCurrentText('leastsq')    

            self.flip_cb.setChecked(False)            
            
            self.cut_min_input.setValue(-0.)        
            self.cut_max_input.setValue(0.)
        pass

    def save_init_params(self):
        if 'QENSfit_conditions' in self.nxentry:
            del self.nxentry['QENSfit_conditions']
        self.nxentry['QENSfit_conditions'] = NXprocess(dQ=self.dq_input.value(),
                                                     Emax=self.emax_input.value(),
                                                     Emin=self.emin_input.value(),
                                                     sigma=self.sigma_input.value(),
                                                     sigma_fixed=self.fix_sigma_cb.isChecked(),
                                                     center=self.center_input.value(),
                                                     center_fixed=self.fix_center_cb.isChecked(),
                                                     fit_exponential_gauss=self.decay_cb.isChecked(),
                                                     decay=self.decay_input.value(),
                                                     gaussian_fit=self.elastic_cb.isChecked(),
                                                     voigt_1_fit=self.voigt_1_cb.isChecked(),
                                                     voigt_2_fit=self.voigt_2_cb.isChecked(),
                                                     cut_energy=self.cut_data_cb.isChecked(),
                                                     energy_cut_min=self.cut_min_input.value(),
                                                     energy_cut_max=self.cut_max_input.value(),
                                                     flip_energy_spectra=self.flip_cb.isChecked(),
                                                     method=str(self.method_combo.currentText())
                                                     )
        pass
    
    # ==================== HELPER METHODS ====================
    
    def create_separator(self):
        """Create a horizontal separator line."""
        separator = QtWidgets.QFrame()
        separator.setFrameShape(QtWidgets.QFrame.HLine)
        separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        return separator
    
    # ==================== FITTING METHODS ====================
    def convolve(dat, kernel):
        """simple convolution of two arrays"""
        npts = min(len(dat), len(kernel))
        pad = np.ones(npts)
        tmp = np.concatenate((pad*dat[0], dat, pad*dat[-1]))
        out = np.convolve(tmp, kernel, mode='valid')
        noff = int((len(out) - npts) / 2)
        return (out[noff:])[:npts]
    
    def fit_data(self):
        """
        Perform the fit using selected parameters and method.
        This is a placeholder that will be implemented with actual fitting logic.
        """
        if self.nxdata is None:
            QtWidgets.QMessageBox.warning(self, 'Warning', 'No data loaded.')
            return
        
        # Extract input parameters
        data = self.select_data()

        if self.cut_data_cb.isChecked():
            mask = (data[data.axes].nxvalue>self.cut_min_input.value()) & (data[data.axes].nxvalue<self.cut_max_input.value()) | (data.nxsignal.nxvalue==0)
        if not self.cut_data_cb.isChecked():
            mask = data.nxsignal.nxvalue==0
        y = data.nxsignal.nxvalue
        y = np.nan_to_num(y)
        x = data[data.axes].nxvalue
        self.fit_x = x

        decay = self.decay_cb.isChecked()
        g1 = self.elastic_cb.isChecked()
        v1 = self.voigt_1_cb.isChecked()
        v2 = self.voigt_2_cb.isChecked()

        # Inititialize QENS Functions
        if decay:
            gauss = ExponentialGaussianModel(prefix='g0_')
            if v1:
                voigt1 = LorentzianModel(prefix='v1_')
                gauss1 = ExponentialGaussianModel(prefix='g1_')
            if v2:
                voigt2 = LorentzianModel(prefix='v2_')
                gauss2 = LorentzianModel(prefix='g2_')
        
        if not decay:
            if g1:
                gauss = GaussianModel(prefix='g0_')
            if v1:
                voigt1 = VoigtModel(prefix='v1_')
            if v2:
                voigt2 = VoigtModel(prefix='v2_')

        residual = y
        # Initialize QENS Parameters
        if decay | g1:
            gaussParams = gauss.guess(y,x)
            residual -= gauss.eval(params=gaussParams,x=x)
            gaussParams['g0_amplitude'].set(min=tiny,max=y.max(),vary=True)
            gaussParams['g0_center'].set(value=self.center_input.value(),vary=(not self.fix_center_cb.isChecked()))
            gaussParams['g0_sigma'].set(value=self.sigma_input.value(),vary=(not self.fix_sigma_cb.isChecked()))
            if decay:
                gaussParams['g0_gamma'].set(value=self.decay_input.value())
                if v1:
                    gauss1Params = gauss1.guess(y,x)
                    for param in gauss1Params:
                        gauss1Params[param].set(expr='g0_'+param[3:])
                if v2:
                    gauss2Params = gauss2.guess(y,x)
                    for param in gauss2Params:
                        gauss2Params[param].set(expr='g0'+param[3:])
        
        if v1:
            voigt1Params = voigt1.guess(y,x)
            residual -= voigt1.eval(params=voigt1Params,x=x)
            voigt1Params['v1_amplitude'].set(min=tiny,max=y.max(),vary=True)
            if voigt1Params['v1_amplitude'].value<0:
                voigt1Params['v1_amplitude'].set(value=tiny)
            voigt1Params['v1_center'].set(value=self.center_input.value(),vary=(not self.fix_center_cb.isChecked()))
            if not decay:
                voigt1Params['v1_gamma'].set(expr='',vary=True,min=0)
                voigt1Params['v1_sigma'].set(value=self.sigma_input.value(),vary=(not self.fix_sigma_cb.isChecked()))
        
        if v2:
            voigt2Params = voigt2.guess(y,x)
            residual -= voigt2.eval(params=voigt2Params,x=x)
            voigt2Params['v2_amplitude'].set(min=tiny,max=y.max(),vary=True)
            if voigt2Params['v2_amplitude'].value<0:
                voigt2Params['v2_amplitude'].set(value=tiny)
            voigt2Params['v2_center'].set(value=self.center_input.value(),vary=(not self.fix_center_cb.isChecked()))
            if not decay:
                voigt2Params['v2_gamma'].set(expr='',vary=True,min=tiny)
                voigt2Params['v2_sigma'].set(value=self.sigma_input.value(),vary=(not self.fix_sigma_cb.isChecked()))
        
        const = ConstantModel(prefix='bkg_')
        constParams = const.guess(y,x)
        constParams['bkg_c'].set(min=tiny)

        model = const
        params = constParams

        # Create function relationships and boundary conditions
        if (g1) & (v1 | v2) & (not decay):
            if v1:
                voigt1Params['v1_sigma'].set(expr='g0_sigma')
                voigt1Params['v1_center'].set(expr='g0_center')
            if v2:
                voigt2Params['v2_sigma'].set(expr='g0_sigma')
                voigt2Params['v2_center'].set(expr='g0_center')
        if decay:
            if v1:
                voigt1Params['v1_center'].set(expr='g0_center')
            if v2:
                voigt2Params['v2_center'].set(expr='g0_center')
        
        if (v1 & v2) & ((not decay) | (not g1)):
            voigt2Params['v2_sigma'].set(expr='v1_sigma')
            voigt2Params['v2_center'].set(expr='v1_center')
        

        # Create full model
        if decay:
            # I need an extra gauss function per voigt. Ideally these would be copies but figure it outs
            model += gauss
            params += gaussParams
            if v1:
                model += CompositeModel(gauss1,voigt1,self.convolve)
                params += voigt1Params + gauss1Params
            if v2:
                model += CompositeModel(gauss2,voigt2,self.convolve)
                params += voigt2Params + gauss2Params
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

        weights = np.ones_like(mask)
        weights *= np.logical_not(mask)
        result = model.fit(y,params,x=x,nan_policy='propagate',method=self.method_combo.currentText(),weights=weights) 
        result = model.fit(y,result.params,x=x,nan_policy='propagate',method='leastsq',weights=weights)
        
        self.fit_result = result
        self.update_results_display(self.fit_result)
    