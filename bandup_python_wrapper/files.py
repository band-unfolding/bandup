import os
import shutil

def mkdir(path):
    try:
        os.makedirs(path)
    except(OSError):
        if(not os.path.isdir(path)):
            raise
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


