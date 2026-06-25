from . import paraminit

from . import incohmodel

def plugin_menu():
    """Return the menu name and list of actions for the plugin."""
    menu = 'QENS'  # Top-level menu name in NeXpy
    actions = []
    actions.append(('Initialize Paramters', paraminit.show_dialog))
    actions.append(('Incoherent Model', incohmodel.show_dialog))
    return menu, actions
