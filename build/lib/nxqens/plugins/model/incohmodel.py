import numpy as np
from nexpy.gui.dialogs import GridParameters, NXDialog
from nexpy.gui.plotview import NXPlotView
from nexpy.gui.pyqt import QtCore
from nexpy.gui.utils import display_message, is_file_locked, report_error
from nexpy.gui.widgets import NXCheckBox, NXLabel
from nexusformat.nexus import NeXusError, NXdata, NXfield

def show_dialog():
    """Entry point called when menu item is clicked."""
    try:
        dialog = MyDialog()
        dialog.show()
    except NeXusError as error:
        report_error("My Plugin Action", error)


class MyDialog(NXDialog):
    """Create a dialog window inheriting from NeXpy's NXDialog."""
    
    def __init__(self, parent=None):
        super(MyDialog, self).__init__(parent)
        
        # Select the NeXus entry to work with
        self.select_entry()
        
        # Create parameters with a grid layout
        self.parameters = GridParameters()
        self.parameters.add('param1', 100, 'Parameter 1')
        self.parameters.add('param2', 50, 'Parameter 2')
        
        # Add buttons and connect them to handler methods
        self.action_button.clicked.connect(self.perform_action)
        
    def perform_action(self):
        """Handler method for processing data."""
        param1 = self.parameters['param1']  # Read parameter values
        # Your data processing logic here
        pass
