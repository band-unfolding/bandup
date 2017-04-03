import os
import shutil
import sys
from bandup_python_wrapper.environ import default_pre_unfolding_inputs_dir

def mkdir(path, ignore_existing=False):
    try:
        os.makedirs(path)
    except(OSError):
        if(os.path.isdir(path) and not ignore_existing): raise
        if(not os.path.isdir(path)): raise
def rmdir(path):
    try:
        shutil.rmtree(path)
    except(OSError):
        if(os.path.isdir(path)):
            raise
def rmfile(path):
    try:
        os.remove(path)
    except(OSError):
        if(os.path.isfile(path)):
            raise

def continuation_lines(fname, marker='&'):
    """Reads data from an opened file object taking into account line continuation marks.

    This routine reads a the data from an input file object and returns its contents 
    taking into account the presence of line continuation marks. 

    The code in this routine is an adaptation of the code available at
    http://stackoverflow.com/questions/16480495/read-a-file-with-line-continuation-
    characters-in-python

    """
    continuation_flines = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.rstrip()
            while line.endswith(marker):
                line = line[:-1] + next(f).rstrip()
            continuation_flines.append(line + '\n')
    return continuation_flines


def get_efermi_fpath(args):
    # The Fermi energy is taken only from the self consistent calculation. 
    # Please ensure convergence w.r.t. k-point mesh for your self-consistent calcs!
    if(args.qe):
        in_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                   os.listdir(args.self_consist_calc_dir) if fname.endswith('.in')][0]
        with open(in_file, 'r') as f:
            in_file_lines = f.readlines()
        for line in in_file_lines:
            if('prefix' in line):
                qe_prefix = line.split('=').split()[0]
        efermi_file = os.path.join(args.self_consist_calc_dir, "%.out"%(qe_prefix))
    elif(args.abinit):
        files_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                      os.listdir(args.self_consist_calc_dir) if 
                      fname.endswith('.files')][0]
        with open(files_file, 'r') as f:
            files_file_lines = f.readlines()
        for line in files_file_lines:
            if('.out' in line):
                efermi_file = os.path.join(args.self_consist_calc_dir,line.strip())
                continue
    elif(args.castep):
        efermi_file = [os.path.join(args.self_consist_calc_dir,fname) for fname in 
                       os.listdir(args.self_consist_calc_dir) if 
                       fname.endswith('.bands')][0]
    else:
        efermi_file = os.path.join(args.self_consist_calc_dir,'OUTCAR')


def get_efermi(args):
    fpath = get_efermi_fpath(args)
    efermi = None

    if(args.castep):
        with open(efermi_file, 'r') as f:
            flines = f.readlines()
            for line in flines:
                if("Fermi energy" in line):
                    efermi_au = float(line.split()[-1])
                    continue
        efermi = efermi_au / physical_constants["electron volt-hartree relationship"][0]

    return efermi


def create_bandup_input(args):
    # PC, SC and PC-KPT files
    input_file_args_var_names = ['pc_file', 'sc_file', 'pckpts_file']
    for arg_name in input_file_args_var_names:
        using_default = False
        if(not '-%s'%(arg_name) in sys.argv[1:]):
            using_default = True
        if(using_default):
            default_fname = args.default_values[arg_name]
            fpath = os.path.join(default_pre_unfolding_inputs_dir, default_fname)
            new_fpath = os.path.join(args.results_dir, default_fname)
            rmfile(new_fpath)
            os.symlink(fpath, new_fpath)

    # Wavefunction symlinks
    wf_fname = None
    wf_file = None
    if(args.castep):
        wf_fname = "%s.orbitals"%(args.seed)
    elif(args.abinit and (not '-wf_file' in sys.argv[1:])):
        with open(args.files_file, 'r') as f:
            flines = f.readlines()
        try:
            wf_fname = flines[3].strip()
        except(IndexError):
            pass
    elif(args.qe):
        raise Exception('Not yet implemmented!')
    elif(not '-wf_file' in sys.argv[1:]):
        wf_fname = args.default_values['wf_file']
    if(wf_fname is not None): 
        wf_file = os.path.join(args.wavefunc_calc_dir, wf_fname)
    if(wf_file is not None):
        wf_file_symlink = os.path.join(args.results_dir, wf_fname)
        rmfile(wf_file_symlink)
        os.symlink(wf_file, wf_file_symlink)
