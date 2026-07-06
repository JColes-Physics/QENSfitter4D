import numpy as np
import logging
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


from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, 
                              QDoubleSpinBox, QPushButton, QProgressBar, QGroupBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont

import multiprocessing as mp
from queue import Queue
import threading

from nxqensfit import NXQENS

import inspect
import re
from itertools import cycle



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

        layout_list = [
            self.entry_layout,
            ]


        # Create array of command buttons
        layout_list.append(self.action_buttons(
            ('Test Fits',self.open_fit_window),
            ('Raster-Scan Fits', self.open_raster_window),
            ))
        layout_list.append(self.close_buttons())


        self.set_layout(*layout_list)
        self.setWindowTitle('Setting Inititial Conditions')


    def open_fit_window(self):
        #dialog2 = FitTab(self.selected_data.title.nxvalue,self.selected_data)
        dialog2 = SpectralFittingWidget(nxdata=self.selected_data,nxentry=self.entry)
        dialog2.show()
        pass

    def open_raster_window(self):
        dialog3 = FittingDialog(nxroot=self.root,nxdata=self.selected_data,nxentry=self.entry)
        dialog3.show()
        pass




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
        
        layout.addWidget(QtWidgets.QLabel('dQ (Å⁻¹) (integration cube length):'), 0, 6)
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

        self.g1frac_input = NXDoubleSpinBox()
        self.v1frac_input = NXDoubleSpinBox()
        self.v2frac_input = NXDoubleSpinBox()
        
        # Voigt and decay options
        options_layout = QtWidgets.QGridLayout()
        
        self.elastic_cb = QtWidgets.QCheckBox('Elastic Gaussian')
        self.voigt_1_cb = QtWidgets.QCheckBox('1 Voigt')
        self.voigt_2_cb = QtWidgets.QCheckBox('2 Voigt')
        self.decay_cb = QtWidgets.QCheckBox('Exponential Decay')
        self.flip_cb = QtWidgets.QCheckBox('Flip Energy Spectra')
        
        options_layout.addWidget(self.elastic_cb, 0, 0)
        options_layout.addWidget(self.g1frac_input, 0, 1)
        options_layout.addWidget(self.voigt_1_cb, 0, 2)
        options_layout.addWidget(self.v1frac_input, 0, 3)
        options_layout.addWidget(self.voigt_2_cb, 0, 4)
        options_layout.addWidget(self.v2frac_input, 0, 5)
        options_layout.addWidget(self.decay_cb, 0, 6)
        options_layout.addWidget(self.flip_cb, 0, 7)
        
        layout.addLayout(options_layout)
        
        # Two Voigt model selection
        two_voigt_layout = QtWidgets.QHBoxLayout()
        two_voigt_layout.addWidget(QtWidgets.QLabel('2 Voigt Model:'))
        
        self.voigt_model_combo = QtWidgets.QComboBox()
        fileslist = ['None (Fit Independently)']
        for item in self.nxentry.entries.values(): 
            if item.nxclass=='NXdata':
                fileslist.append(str(item))
        self.voigt_model_combo.addItems(fileslist)
        two_voigt_layout.addWidget(self.voigt_model_combo)
        
        self.load_model_btn = QtWidgets.QPushButton('Load Model')
        self.load_model_btn.clicked.connect(self.load_model_from_file)
        two_voigt_layout.addWidget(self.load_model_btn)
        
        layout.addLayout(two_voigt_layout)
        
        self.elastic_cb.stateChanged.connect(
            lambda: self.g1frac_input.setEnabled(self.elastic_cb.isChecked() | self.decay_cb.isChecked())
        )
        self.decay_cb.stateChanged.connect(
            lambda: self.g1frac_input.setEnabled(self.elastic_cb.isChecked() | self.decay_cb.isChecked())
        )

        self.voigt_1_cb.stateChanged.connect(
            lambda: self.v1frac_input.setEnabled(self.voigt_1_cb.isChecked())
        )
        # Enable/disable model selection based on 2 Voigt checkbox
        self.voigt_2_cb.stateChanged.connect(
            lambda: self.voigt_model_combo.setEnabled(self.voigt_2_cb.isChecked())
        )
        self.voigt_2_cb.stateChanged.connect(
            lambda: self.load_model_btn.setEnabled(self.voigt_2_cb.isChecked())
        )
        self.voigt_2_cb.stateChanged.connect(
            lambda: self.v2frac_input.setEnabled(self.voigt_2_cb.isChecked())
        )
        self.voigt_model_combo.setEnabled(False)
        self.load_model_btn.setEnabled(False)
        
        group.setLayout(layout)
        return group
    
    def load_model_from_file(self):
        if str(self.voigt_model_combo.currentText()) != 'None (Fit Independently)':
            self.v2model = self.nxentry[str(self.voigt_model_combo.currentText())]
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
        self.results_text.setFont(QtGui.QFont())
        
        layout.addWidget(self.results_text)
        
        group.setLayout(layout)
        return group
    
    def update_results_display(self, fit_result,weights = np.array([])):
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
        if weights.shape != 0:
            result_text += np.array_str(weights)
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
            if self.cut_data_cb.isChecked():
                cutmin = self.cut_min_input.value()
                cutmax = self.cut_max_input.value()
            if not self.cut_data_cb.isChecked():
                cutmin = 0
                cutmax = 0
            mask = (data[data.axes].nxvalue>cutmin) & (data[data.axes].nxvalue<cutmax)
            data_to_plot = NXdata(signal=NXfield(np.ma.array(data.nxsignal,mask=mask),name='data'),axes=data[data.axes],errors=data[data.errors])
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
            try:
                data = NXdata(data.nxsignal,data[data.axes].centers(),errors=data[data.errors])
            except:
                data = NXdata(data.nxsignal,data[data.axes].centers())
        
        if self.flip_cb.isChecked():
            axes = data[data.axes]
            diff = abs(axes.max())-abs(axes.min())
            try:
                data = NXdata(np.flip(data.nxsignal),axes=(data[data.axes]-diff),errors=np.flip(data[data.errors]))
            except:
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

            self.g1frac_input.setValue(self.init_condition['gaussian_amplitude_fraction'].nxvalue)
            self.v1frac_input.setValue(self.init_condition['voigt_1_amplitude_fraction'].nxvalue)
            self.v2frac_input.setValue(self.init_condition['voigt_2_amplitude_fraction'].nxvalue)
            
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

            self.g1frac_input.setValue(1.)
            self.v1frac_input.setValue(1.)
            self.v2frac_input.setValue(1.)              
            
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
                                                     method=str(self.method_combo.currentText()),
                                                     gaussian_amplitude_fraction = self.v1frac_input.value(),
                                                     voigt_1_amplitude_fraction = self.v1frac_input.value(),
                                                     voigt_2_amplitude_fraction = self.v2frac_input.value(),
                                                     v2_model = str(self.voigt_model_combo.currentText()),
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
    def convolve(self, dat, kernel):
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
        
        # Extract input data and format to prefered formatting
        data = self.select_data()
        y = data.nxsignal.nxvalue
        y = np.nan_to_num(y)
        x = data[data.axes].nxvalue
        self.fit_x = x

        # Find initial data masks for weighting fits
        if self.cut_data_cb.isChecked():
            mask = (data[data.axes].nxvalue>self.cut_min_input.value()) & (data[data.axes].nxvalue<self.cut_max_input.value()) | (data.nxsignal.nxvalue==0) | np.isnan(data.nxsignal.nxvalue)
        if not self.cut_data_cb.isChecked():
            mask = data.nxsignal.nxvalue==0 | np.isnan(data.nxsignal.nxvalue)

        # Establish boolean defined fitting conditions
        decay = self.decay_cb.isChecked()
        g1 = self.elastic_cb.isChecked()
        v1 = self.voigt_1_cb.isChecked()
        v2 = self.voigt_2_cb.isChecked()
        g1frac = self.g1frac_input.value()
        v1frac = self.v1frac_input.value()
        v2frac = self.v2frac_input.value()
        sigma = self.sigma_input.value()
        center = self.center_input.value()
        decayval = self.decay_input.value()
        fixsigma = self.fix_sigma_cb.isChecked()
        fixcenter = self.fix_center_cb.isChecked()
        if str(self.voigt_model_combo.currentText()) == 'None (Fit Independently)':
            v2gamm = None
        elif str(self.voigt_model_combo.currentText()) != 'None (Fit Independently)':
            v2gamm = self.v2model[self.l_input.value(),self.k_input.value(),self.h_input.value()]
        method = self.method_combo.currentText()


        qensfit = NXQENS(x,y,mask,g1,v1,v2,decay,g1frac,v1frac,v2frac,sigma,center,decayval,fixsigma,fixcenter,v2gamm,method)
        result = qensfit.fit_data()

        self.fit_result = result

        self.update_results_display(self.fit_result)




class FittingWorker(QThread):
    """Worker thread that manages the fitting process with multithreading."""
    
    # Signals
    progress_updated = pyqtSignal(int)  # Progress percentage
    result_ready = pyqtSignal(dict)     # Individual fit result: {coords, result_data}
    fitting_complete = pyqtSignal()
    fitting_failed = pyqtSignal(str)    # Error message
    
    def __init__(self, h_range, k_range, l_range, data, num_cores, 
                 init_conditions):
        """
        Args:
            h_range: tuple (min, max, step) for H
            k_range: tuple (min, max, step) for K
            l_range: tuple (min, max, step) for L
            num_cores: number of worker threads
            fitting_function: callable that performs the fit
            *fitting_args: additional args to pass to fitting_function
        """
        super().__init__()
        self.h_range = h_range
        self.k_range = k_range
        self.l_range = l_range
        self.num_cores = num_cores
        self.data = data
        self.init_conditions = init_conditions
        
        params={param:self.init_conditions[param].nxvalue for param in self.init_conditions}
        
        self.__dict__.update(params)

        self._logger = None

        self.is_running = True
        self._lock = threading.Lock()

    @property
    def logger(self):
        """Log file handler."""
        if self._logger is None:
            self._logger = logging.getLogger(
                f"Multi-core-drifting")
            self._logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                "%(asctime)s %(name)-12s: %(message)s",
                datefmt='%Y-%m-%d %H:%M:%S')
            fileHandler = logging.FileHandler(
                self.task_directory.joinpath('nxlogger.log'))
            fileHandler.setFormatter(formatter)
            self._logger.addHandler(fileHandler)
        return self._logger
    
    def stop(self):
        """Signal the worker to stop."""
        with self._lock:
            self.is_running = False
    
    def run(self):
        """Main fitting loop with multithreaded Q-space sampling."""
        try:
            # Generate all (H, K, L) coordinates
            h_vals = np.arange(self.h_range[0], self.h_range[1] + tiny, 
                                self.h_range[2])
            k_vals = np.arange(self.k_range[0], self.k_range[1] + tiny, 
                                self.k_range[2])
            l_vals = np.arange(self.l_range[0], self.l_range[1] + tiny, 
                                self.l_range[2])
            
            # Create mesh of all coordinate combinations
            hkl_coords = np.array(np.meshgrid(h_vals, k_vals, l_vals, indexing='ij'))
            hkl_coords = hkl_coords.reshape(3, -1).T  # Shape: (N, 3)
            
            total_points = len(hkl_coords)
            result_queue = Queue()
            
            # Start worker thread pool
            worker_threads = []
            for i in range(self.num_cores):
                t = threading.Thread(
                    target=self._worker_thread,
                    args=(result_queue, hkl_coords, i)
                )
                t.daemon = True
                t.start()
                worker_threads.append(t)
            
            # Collect results from queue
            processed = 0
            while processed < total_points:
                with self._lock:
                    if not self.is_running:
                        self.fitting_failed.emit("Fitting cancelled by user")
                        return
                
                try:
                    result = result_queue.get(timeout=1)
                    if result is not None:
                        self.result_ready.emit(result)
                        processed += 1
                        progress = int((processed / total_points) * 100)
                        self.progress_updated.emit(progress)
                except:
                    # Timeout or empty queue, check if threads are alive
                    if not any(t.is_alive() for t in worker_threads):
                        break
            
            # Wait for all threads to finish
            for t in worker_threads:
                t.join(timeout=5)
            
            self.fitting_complete.emit()
            
        except Exception as e:
            self.logger.info(f"Error during fitting: {str(e)}")
    
    def _worker_thread(self, result_queue, hkl_coords, thread_id):
        """Worker thread function - processes HKL points."""
        num_points = len(hkl_coords)
        points_per_thread = (num_points + self.num_cores - 1) // self.num_cores
        start_idx = thread_id * points_per_thread
        end_idx = min(start_idx + points_per_thread, num_points)
        
        try:
            self.logger.info(f"Thread {thread_id} starting: processing points {start_idx} to {end_idx}")
        except:
            pass

        for idx in range(start_idx, end_idx):
            with self._lock:
                if not self.is_running:
                    return
            
            h, k, l = hkl_coords[idx]
            
            try:
                if self.dQ > 0:
                    data = self.data[self.Emin:self.Emax,
                                l-(self.dQ/2):l+(self.dQ/2),
                                k-(self.dQ/2):k+(self.dQ/2),
                                h-(self.dQ/2):h+(self.dQ/2)].sum((1,2,3))
                if self.dQ <= 0:
                    data = self.data[self.Emin:self.Emax,
                                    l,k,h]
                if data.nxsignal.shape != data[data.axes].shape:
                    try:
                        data = NXdata(data.nxsignal,data[data.axes].centers(),errors=data[data.errors])
                    except:
                        data = NXdata(data.nxsignal,data[data.axes].centers())
                
                if self.flip_energy_spectra:
                    axes = data[data.axes]
                    diff = abs(axes.max())-abs(axes.min())
                    try:
                        data = NXdata(np.flip(data.nxsignal),axes=(data[data.axes]-diff),errors=np.flip(data[data.errors]))
                    except:
                        data = NXdata(np.flip(data.nxsignal),axes=(data[data.axes]-diff))

                if (np.isnan(data.nxsignal).sum() > data.nxsignal.shape[0]/4) or ((data.nxsignal==0).sum() > data.nxsignal.shape[0]/4):
                    result_queue.put(None)
                    pass

                y = data.nxsignal.nxvalue
                y = np.nan_to_num(y)
                x = data[data.axes].nxvalue
                self.fit_x = x

                # Find initial data masks for weighting fits
                if self.cut_energy:
                    mask = (data[data.axes].nxvalue>self.cut_min_input.value()) & (data[data.axes].nxvalue<self.cut_max_input.value()) | (data.nxsignal.nxvalue==0) | np.isnan(data.nxsignal.nxvalue)
                if not self.cut_energy:
                    mask = data.nxsignal.nxvalue==0 | np.isnan(data.nxsignal.nxvalue)

                if self.v2_model == 'None (Fit Independently)':
                    v2gamm = None
                elif self.v2_model != 'None (Fit Independently)':
                    v2gamm = self.nxentry[self.v2_model][l,k,h]


                # Call fitting function with the dQ step size
                fit_result = NXQENS(
                    x, y, mask,
                    g1=self.gaussian_fit,
                    v1=self.voigt_1_fit,
                    v2=self.voigt_2_fit,
                    decay=self.fit_exponential_gauss,
                    g1frac=self.gaussian_amplitude_fraction,
                    v1frac=self.voigt_1_amplitude_fraction,
                    v2frac=self.voigt_2_amplitude_fraction,
                    sigma=self.sigma,
                    center=self.center,
                    decayval=self.decay,
                    fixsigma=self.sigma_fixed,
                    fixcenter=self.center_fixed,
                    v2gamm=v2gamm,
                    method=self.method
                ).fit_data()
                
                result_dict = {
                    'coords': (float(h), float(k), float(l)),
                    'result': fit_result
                }
                result_queue.put(result_dict)

                # try:
                #     self.logger.info(f"Thread {thread_id} finished")
                # except:
                #     pass
                
            except Exception as e:
                result_dict = {
                    'coords': (float(h), float(k), float(l)),
                    'result': None
                }
                # Log error but continue processing
                #self.logger.info(f"Error fitting ({h}, {k}, {l}): {str(e)}")
                result_queue.put(result_dict)


class FittingDialog(NXDialog):
    """Dialog for configuring and running multithreaded fitting."""
    
    def __init__(self,nxroot=None, parent=None, nxdata=None, nxentry=None):
        """
        Args:
            parent: Parent widget
            nxdata: NXdata passed by previous window
            nxentry: NXentry passed by previous window
        """
        super().__init__(parent=parent)
        
        try:
            self.nxroot = nxroot
        except:
            self.nxroot = None
        self.nxdata = nxdata
        self.nxentry = nxentry
        self.init_conditions = self.nxentry.QENSfit_conditions
        self.sample = nxentry.sample
        self.scan = nxdata
        
        # Storage for results (thread-safe via Qt signals)
        self.results_matrix = {}  # Dictionary: (h, k, l) -> fit_result
        
        # Worker thread
        self.worker = None
        self.is_fitting = False

        self._logger = None
        
        self.init_ui()

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
    
    def init_ui(self):
        """Initialize the dialog UI."""
        self.setWindowTitle("Advanced Fitting Configuration")
        data = self.nxdata
        # Main layout
        layout = QVBoxLayout()
        
        # --- HKL Range Group ---
        hkl_group = QGroupBox("HKL Range")
        hkl_layout = QHBoxLayout()
        
        # H range
        hkl_layout.addWidget(QLabel("H min:"))
        self.h_min = QDoubleSpinBox()
        self.h_min.setRange(data[data.axes[3]].min(), data[data.axes[3]].max())
        self.h_min.setValue(data[data.axes[3]].min())
        self.h_min.setSingleStep(abs(data[data.axes[3]][1]-data[data.axes[3]][0]))
        hkl_layout.addWidget(self.h_min)
        
        hkl_layout.addWidget(QLabel("H max:"))
        self.h_max = QDoubleSpinBox()
        self.h_max.setRange(data[data.axes[3]].min(), data[data.axes[3]].max())
        self.h_max.setValue(data[data.axes[3]].max())
        self.h_max.setSingleStep(abs(data[data.axes[3]][1]-data[data.axes[3]][0]))
        hkl_layout.addWidget(self.h_max)
        
        # K range
        hkl_layout.addWidget(QLabel("K min:"))
        self.k_min = QDoubleSpinBox()
        self.k_min.setRange(data[data.axes[2]].min(), data[data.axes[2]].max())
        self.k_min.setValue(data[data.axes[2]].min())
        self.k_min.setSingleStep(abs(data[data.axes[2]][1]-data[data.axes[2]][0]))
        hkl_layout.addWidget(self.k_min)
        
        hkl_layout.addWidget(QLabel("K max:"))
        self.k_max = QDoubleSpinBox()
        self.k_max.setRange(data[data.axes[2]].min(), data[data.axes[2]].max())
        self.k_max.setValue(data[data.axes[2]].max())
        self.k_max.setSingleStep(abs(data[data.axes[2]][1]-data[data.axes[2]][0]))
        hkl_layout.addWidget(self.k_max)
        
        # L range
        hkl_layout.addWidget(QLabel("L min:"))
        self.l_min = QDoubleSpinBox()
        self.l_min.setRange(data[data.axes[1]].min(), data[data.axes[1]].max())
        self.l_min.setValue(data[data.axes[1]].min())
        self.l_min.setSingleStep(abs(data[data.axes[1]][1]-data[data.axes[1]][0]))
        hkl_layout.addWidget(self.l_min)
        
        hkl_layout.addWidget(QLabel("L max:"))
        self.l_max = QDoubleSpinBox()
        self.l_max.setRange(data[data.axes[1]].min(), data[data.axes[1]].max())
        self.l_max.setValue(data[data.axes[1]].max())
        self.l_max.setSingleStep(abs(data[data.axes[1]][1]-data[data.axes[1]][0]))        
        hkl_layout.addWidget(self.l_max)


        hkl_group.setLayout(hkl_layout)
        layout.addWidget(hkl_group)
        
        # --- dQ and Q-step Display Group ---
        dq_group = QGroupBox("Q-Space Sampling")
        dq_layout = QHBoxLayout()
        
        dq_layout.addWidget(QLabel("dQ (step size):"))
        self.dq_input = QDoubleSpinBox()
        self.dq_input.setRange(0.001, 100.0)
        self.dq_input.setValue(self.init_conditions.dQ.nxvalue/2)
        self.dq_input.setSingleStep(0.001)
        self.dq_input.setDecimals(4)
        self.dq_input.valueChanged.connect(self.update_dq_display)
        dq_layout.addWidget(self.dq_input)
        
        dq_layout.addWidget(QLabel("dQ_integration / dQ_step:"))
        ratio = (self.init_conditions['dQ']) / self.dq_input.value()
        self.dq_display = QLabel(f"{ratio:.4f}")
        font = QFont()
        font.setPointSize(11)
        font.setBold(True)
        self.dq_display.setFont(font)
        dq_layout.addWidget(self.dq_display)
        
        dq_layout.addStretch()
        dq_group.setLayout(dq_layout)
        layout.addWidget(dq_group)
        
        # --- Cores Configuration ---
        cores_group = QGroupBox("Computation")
        cores_layout = QHBoxLayout()
        
        cores_layout.addWidget(QLabel("Number of cores:"))
        self.cores_input = QSpinBox()
        max_cores = mp.cpu_count()
        self.cores_input.setRange(1, max_cores)
        self.cores_input.setValue(max_cores - 1)
        cores_layout.addWidget(self.cores_input)
        
        cores_layout.addWidget(QLabel(f"(System has {max_cores} cores)"))
        cores_layout.addStretch()
        cores_group.setLayout(cores_layout)
        layout.addWidget(cores_group)
        
        # --- Control Buttons ---
        button_layout = QHBoxLayout()
        
        self.fit_button = QPushButton("Begin Fitting")
        self.fit_button.clicked.connect(self.start_fitting)
        button_layout.addWidget(self.fit_button)
        
        self.cancel_button = QPushButton("Cancel Fitting")
        self.cancel_button.clicked.connect(self.cancel_fitting)
        self.cancel_button.setEnabled(False)
        button_layout.addWidget(self.cancel_button)
        
        layout.addLayout(button_layout)
        
        # --- Progress Bar ---
        progress_layout = QHBoxLayout()
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)
        
        self.close_button = QPushButton("Close")
        self.close_button.clicked.connect(self.accept)
        progress_layout.addWidget(self.close_button)
        
        layout.addLayout(progress_layout)
        
        self.setLayout(layout)
        self.resize(900, 300)
    
    def update_dq_display(self, dq_value):
        """Update the dQ / reference display dynamically."""
        # Use 1.0 as default reference; you can pass a custom reference value
        reference = 1#self.data_dict.get('dq_reference', 1.0)
        ratio = (self.init_conditions['dQ']) / (dq_value)
        self.dq_display.setText(f"{ratio:.4f}")
    
    def start_fitting(self):
        """Start the fitting worker thread."""
        # if not self.fitting_function:
        #     nexpy.report.error("No fitting function provided")
        #     return
        
        # Validate inputs
        if self.h_min.value() > self.h_max.value():
            raise Exception("H min must be <= H max")
            return
        if self.k_min.value() > self.k_max.value():
            raise Exception("K min must be <= K max")
            return
        if self.l_min.value() > self.l_max.value():
            raise Exception("L min must be <= L max")
            return
        
        # Disable fitting controls
        self.fit_button.setEnabled(False)
        self.h_min.setEnabled(False)
        self.h_max.setEnabled(False)
        self.k_min.setEnabled(False)
        self.k_max.setEnabled(False)
        self.l_min.setEnabled(False)
        self.l_max.setEnabled(False)
        self.dq_input.setEnabled(False)
        self.cores_input.setEnabled(False)
        
        self.cancel_button.setEnabled(True)
        self.is_fitting = True
        
        # Create and start worker
        h_range = (self.h_min.value(), self.h_max.value(), self.dq_input.value())
        k_range = (self.k_min.value(), self.k_max.value(), self.dq_input.value())
        l_range = (self.l_min.value(), self.l_max.value(), self.dq_input.value())
        
        self.worker = FittingWorker(
            h_range, k_range, l_range,
            self.nxdata,
            self.cores_input.value(),
            self.init_conditions
        )
        
        # Connect signals
        self.worker.progress_updated.connect(self.on_progress)
        self.worker.result_ready.connect(self.on_result_ready)
        self.worker.fitting_complete.connect(self.on_fitting_complete)
        self.worker.fitting_failed.connect(self.on_fitting_failed)
        
        self.worker.start()
    
    def on_progress(self, progress):
        """Update progress bar."""
        self.progress_bar.setValue(progress)
    
    def on_result_ready(self, result):
        """Store a fitting result (thread-safe via Qt signal)."""
        coords = result['coords']
        fit_result = result['result']
        self.results_matrix[coords] = fit_result

        #self.logger.info(f"Result stored: {coords} -> {fit_result}")    
        #self.logger.info(f"Total results so far: {len(self.results_matrix)}")
        
        # Optionally log or process the result
        # print(f"Fitted {coords}: {fit_result}")

    def cancel_fitting(self):
        self.worker.stop()
        self.is_fitting=False

        # Re-enable controls
        self.fit_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.close_button.setEnabled(True)
        self.h_min.setEnabled(True)
        self.h_max.setEnabled(True)
        self.k_min.setEnabled(True)
        self.k_max.setEnabled(True)
        self.l_min.setEnabled(True)
        self.l_max.setEnabled(True)
        self.dq_input.setEnabled(True)
        self.cores_input.setEnabled(True)
    
    def on_fitting_complete(self):
        """Handle fitting completion."""
        self.is_fitting = False
        self.progress_bar.setValue(100)

        self.logger.info(f"=== FITTING COMPLETE ===")    
        self.logger.info(f"Total points in results_matrix: {len(self.results_matrix)}")
        #nexpy.report.info(f"Fitting complete! Processed {len(self.results_matrix)} points.")
        
        # Re-enable controls
        self.fit_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.close_button.setEnabled(True)
        self.h_min.setEnabled(True)
        self.h_max.setEnabled(True)
        self.k_min.setEnabled(True)
        self.k_max.setEnabled(True)
        self.l_min.setEnabled(True)
        self.l_max.setEnabled(True)
        self.dq_input.setEnabled(True)
        self.cores_input.setEnabled(True)

        self.store_results()
    
    def on_fitting_failed(self, error_msg):
        """Handle fitting error."""
        self.is_fitting = False
        
        # Re-enable controls
        self.fit_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.close_button.setEnabled(True)
        self.h_min.setEnabled(True)
        self.h_max.setEnabled(True)
        self.k_min.setEnabled(True)
        self.k_max.setEnabled(True)
        self.l_min.setEnabled(True)
        self.l_max.setEnabled(True)
        self.dq_input.setEnabled(True)
        self.cores_input.setEnabled(True)

        raise Exception(error_msg)
    
    def store_results(self):
        
        if not self.results_matrix:
            self.logger.info("No points fit, no data saved")
            return

        h_range = (self.h_min.value(), self.h_max.value(), self.dq_input.value())
        k_range = (self.k_min.value(), self.k_max.value(), self.dq_input.value())
        l_range = (self.l_min.value(), self.l_max.value(), self.dq_input.value())

        h_vals = np.arange(h_range[0], h_range[1] + tiny, 
                            h_range[2])
        k_vals = np.arange(k_range[0], k_range[1] + tiny, 
                            k_range[2])
        l_vals = np.arange(l_range[0], l_range[1] + tiny, 
                            l_range[2])
        size = (len(l_vals),len(k_vals),len(h_vals))

        H = NXfield(h_vals,name='Qh')
        K = NXfield(k_vals,name='Qk')
        L = NXfield(l_vals,name='Ql')
        
        #initialize nxdata sets to save
        import random
        randpoint,randresult=random.choice(list(self.results_matrix.items()))
        params_to_save = [param for param in randresult.params if randresult.params[param].vary]

        nxsave = NXentry()
        for param in params_to_save:
            nxsave[param] = NXdata(signal=np.zeros(size),axes=(L,K,H),errors=np.zeros(size))
        nxsave['redchi'] = NXdata(signal=np.zeros(size),axes=(L,K,H))
        nxsave['rsquared'] = NXdata(signal=np.zeros(size),axes=(L,K,H))
        nxsave['fitted_parameters'] = self.init_conditions

        for idx,result in self.results_matrix.items():
            h,k,l = idx
            nxsave['redchi'].signal[l,k,h] = result.redchi
            nxsave['rsquared'].signal[l,k,h] = result.rsquared
            for param in params_to_save:
                nxsave[param].signal[l,k,h] = result.params[param].value
                nxsave[param].signal_errors[l,k,h] = result.params[param].stderr

        if self.root:
            self.root.unlock()
            if 'QENSfit_results' in self.nxentry:
                del self.nxentry['QENSfit_results']
            self.nxentry['QENSfit_results'] = nxsave

            self.logger.info('Results saved to QENSfit_results')
            self.root.lock()
        elif not self.root:
            try:
                if 'QENSfit_results' in self.nxentry:
                    del self.nxentry['QENSfit_results']
                self.nxentry['QENSfit_results'] = nxsave

                self.logger.info('Results saved to QENSfit_results')
            except:
                self.logger.info('Failed to save results')
            

        
