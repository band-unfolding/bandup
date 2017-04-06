import warnings
import os

def _new_showwarning(message, category = UserWarning, filename = '', lineno = -1):
    warn_msg = 'WARNING (%s, line %d): %s'%(os.path.basename(filename), lineno, 
                                           message)
    print(warn_msg)
warnings.showwarning = _new_showwarning
