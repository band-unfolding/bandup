import warnings
import os
import sys

class WarningError(UserWarning):
    pass

def _new_showwarning(message, category = UserWarning, filename = '', lineno = -1):
    if(category==WarningError):
        warn_msg = 'ERROR '
    else:
        warn_msg = 'WARNING '
    warn_msg += '(%s, line %d): %s'%(os.path.basename(filename), lineno, 
                                           message)
    if(category==WarningError):
        warn_msg += '\nCannot continue. Stopping now.'
    print(warn_msg)
    if(category==WarningError):
        sys.exit(1)
warnings.showwarning = _new_showwarning

    
