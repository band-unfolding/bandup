import os
import sys

try:
    bandup_dir = os.environ['BANDUPDIR']
except(KeyError):
    msg = 'TEST'
    msg = ('The variable "BANDUPDIR" is not defined in your environment.\n' +
           18*' '+"Please set it yo BandUP's main directory.\n"+
           18*' '+'    > This is the directory where "build.sh" is located.\n'+
           18*' '+'Example:\n'+
           28*' '+'export BANDUPDIR=$HOME/codes/BandUP\n'+
           18*' '+'To have this automatically set next time you open a terminal, \n'+
           18*' '+'you can add the same command to, e.g., the file "~/.bash_profile".'
          )
    raise EnvironmentError(msg)

user_home = os.path.expanduser('~')
plot_path = os.path.join(bandup_dir, "utils", "post_unfolding", "plot")
plot_script = os.path.join(plot_path, "plot_unfolded_EBS_BandUP.py")
sys.path.insert(0, os.path.dirname(plot_script))

working_dir = os.getcwd()
BandUp_exe = os.path.join(bandup_dir, "BandUP_bin", "BandUP.x")

bandup_raw = {'symmetry_avgd':'unfolded_EBS_symmetry-averaged.dat',
              'not_symmetry_avgd':'unfolded_EBS_not-symmetry_averaged.dat'}
