import numpy as np
from nexpy.gui.dialogs import GridParameters, NXDialog
from nexpy.gui.plotview import NXPlotView
from nexpy.gui.pyqt import QtCore
from nexpy.gui.utils import display_message, is_file_locked, report_error
from nexpy.gui.widgets import NXCheckBox, NXLabel
from nexusformat.nexus import NeXusError, NXdata, NXfield
from nexpy.api.nexus import *
from nexpy.gui.importdialog import NXImportDialog
from nxqens.nxquensfit import NXQENS


def show_dialog(parent=None):
    """Entry point called when menu item is clicked."""
    try:
        dialog = MyDialog()
        dialog.show()
    except NeXusError as error:
        report_error("Pameter Initialization", error)


class MyDialog(NXDialog):
    """Create a dialog window inheriting from NeXpy's NXDialog."""
    
    def __init__(self, parent=None):
        super(MyDialog, self).__init__(parent)
        
        # Select the NeXus entry to work with
        self.select_entry()
        self.fit = NXQENS
        
        # Create parameters with a grid layout
        self.parameters = GridParameters()
        self.parameters.add(
            'Ei', self.entry['instrument/monochromator/energy'],
            'Incident Energy')
        self.parameters.add('dQ', self.round(np.sqrt(self.Ei/2)/50), 'Q Step')
        self.parameters.add('dE', self.round(self.Ei/50), 'Energy Step')
        self.set_layout(self.entry_layout,
                        self.parameters.grid(),
                        self.action_buttons(('Plot', self.plot_data),
                                            ('Save', self.save_data)),
                        self.close_buttons())

        self.setWindowTitle('Setting Inititial Conditions')

        self.set_layout(self.Qgrid, self.close_buttons(save=True))
        
        # Add buttons and connect them to handler methods
        self.action_button.clicked.connect(self.perform_action)
        
    def perform_action(self):
        """Handler method for processing data."""
        param1 = self.parameters['param1']  # Read parameter values
        # Your data processing logic here
        pass

    def read_parameters(self):
        self.L1 = - self.entry['sample/distance']
        self.L2 = self.entry['instrument/detector/distance'].average()
        self.m1 = self.entry['monitor1']
        self.t_m1 = self.m1.moment()
        self.d_m1 = self.entry['monitor1/distance']

    def plot_data(self):
        self.convert_QE().plot() 

    def convert_QE(self):
        """Convert S(phi,eps) to S(Q,eps)"""

        self.read_parameters()

        Ei = self.Ei
        dQ = self.dQ
        dE = self.dE

        signal = self.entry['data'].nxsignal
        pol = centers(self.entry['data/polar_angle'], signal.shape[0])
        tof = centers(self.entry['data/time_of_flight'], signal.shape[1])
        en = self.convert_tof(tof)

        idx_max = min(np.where(np.abs(en-0.75*Ei) < 0.1)[0])

        en = en[:idx_max]

        data = signal.nxdata[:, :idx_max]
        if self.entry['data'].nxerrors:
            errors = self.entry['data'].nxerrors.nxdata[:]

        Q = np.zeros((len(pol), len(en)))
        E = np.zeros((len(pol), len(en)))

        for i in range(0, len(pol)):
            p = pol[i]
            Q[i, :] = np.array(np.sqrt((2*Ei - en - 2*np.sqrt(Ei*(Ei-en))
                                       * np.cos(p*np.pi/180.0))/2.0721))
            E[i, :] = np.array(en)

        s = Q.shape
        Qin = Q.reshape(s[0]*s[1])
        Ein = E.reshape(s[0]*s[1])
        datain = data.reshape(s[0]*s[1])
        if self.entry['data'].nxerrors:
            errorsin = errors.reshape(s[0]*s[1])

        qmin = Q.min()
        qmax = Q.max()
        emin = E.min()
        emax = E.max()
        NQ = int((qmax-qmin)/dQ) + 1
        NE = int((emax-emin)/dE) + 1
        Qb = np.linspace(qmin, qmax, NQ)
        Eb = np.linspace(emin, emax, NE)
        # histogram and normalize
        norm, nbin = np.histogramdd((Ein, Qin), bins=(Eb, Qb))
        hist, hbin = np.histogramdd((Ein, Qin), bins=(Eb, Qb), weights=datain)
        if self.entry['data'].nxerrors:
            histe, hbin = np.histogramdd((Ein, Qin), bins=(Eb, Qb),
                                         weights=errorsin * errorsin)
            histe = histe**0.5
            err = histe/norm

        Ib = NXfield(hist/norm, name='S(Q,E)')

        Qb = NXfield(Qb[:-1]+dQ/2., name='Q')
        Eb = NXfield(Eb[:-1]+dE/2., name='E')

        result = NXdata(Ib, (Eb, Qb))
        if self.entry.data.nxerrors:
            result.errors = NXfield(err)
        return result
