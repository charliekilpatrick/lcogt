"""
lcogt package initialization.

"""

# Imports for signal and log handling
import os
import sys
import signal
import warnings

def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s (%s:%s)\n' % (category.__name__, message, os.path.split(filename)[1], lineno)

warnings.formatwarning = short_warning


# Set version
__version__ = '1.0.0'
