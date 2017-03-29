import os
import shutil

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
def rmdir(path):
    try:
        shutil.rmtree(path)
    except(OSError):
        if os.path.isdir(path):
            raise
def rmfile(path):
    try:
        os.remove(path)
    except(OSError):
        if os.path.isdir(path):
            raise

def continuation_lines(fin):
    """ Routine obtained from

    http://stackoverflow.com/questions/16480495/read-a-file-with-line-continuation-
    characters-in-python

    """
    for line in fin:
        #line = line.rstrip('\n').rstrip()
        line = line.strip()
        #while line.endswith('\\'):
        while line.endswith('&'):
            #line = line[:-1] + next(fin).rstrip('\n')
            line = line[:-1] + next(fin).strip()
        yield line
