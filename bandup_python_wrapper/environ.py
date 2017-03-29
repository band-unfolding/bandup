import os
import sys

try:
    bandup_dir = os.environ['BANDUPDIR']
    #bandup_folder = os.path.join(user_home, "repos", "git", "my_codes", "BandUP")
except(KeyError):
    print('Environment variable "BANDUPDIR" not set. Please set it to \n'+
          "BandUP's main directory " + '(the one where "build.sh" is located).')
    print 'Stopping now.'
    sys.exit(1)

user_home = os.path.expanduser('~')
plot_script = os.path.join(bandup_dir, "utils", "post_unfolding", "plot", 
                           "plot_unfolded_EBS_BandUP.py")
sys.path.insert(0, os.path.dirname(plot_script))

working_dir = os.getcwd()
BandUp_exe = os.path.join(bandup_dir, "BandUP_bin", "BandUP.x")

bandup_raw = {'symmetry_avgd':'unfolded_EBS_symmetry-averaged.dat',
              'not_symmetry_avgd':'unfolded_EBS_not-symmetry_averaged.dat'}
