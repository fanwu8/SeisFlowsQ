
import os
import subprocess
import sys
import numpy as np

from collections import defaultdict
from os.path import abspath, join, exists
from seisflows.tools import msg, unix
from seisflows.tools.tools import iterable


def call_solver(mpiexec, executable, output='solver.log'):
    """ Calls MPI solver executable

      A less complicated version, without error catching, would be
      subprocess.call(mpiexec +' '+ executable, shell=True)
    """
    try:
        f = open(output,'w')
        subprocess.check_call(
            mpiexec +' '+ executable,
            shell=True,
            stdout=f)
    except subprocess.CalledProcessError as err:
        print(msg.SolverError % (mpiexec +' '+ executable))
        sys.exit(-1)
    except OSError:
        print(msg.SolverError % (mpiexec +' '+ executable))
        sys.exit(-1)
    finally:
        f.close()


class Minmax(defaultdict):
    """ Keeps track of min,max values of model or kernel
    """
    def __init__(self):
        super(Minmax, self).__init__(lambda: [+np.inf, -np.inf])

    def update(self, keys, vals):
        for key, val in _zip(keys, vals):
            if min(val) < self.dict[key][0]:
                self.dict[key][0] = min(val)
            if max(val) > self.dict[key][1]:
                self.dict[key][1] = max(val)

    def __call__(self, key):
        return self.dict[key]


class Container(defaultdict):
    """ Dictionary-like object for holding models or kernels
    """
    def __init__(self):
        super(Container, self).__init__(lambda: [])
        self.minmax = Minmax()


class Writer(object):
    """ Utility for appending values to text files
    """
    def __init__(self, path='./output.stat'):
        self.path = abspath(path)
        try:
            os.mkdir(path)
        except:
            raise IOError

        self.__call__('step_count', 0)

    def __call__(self, filename, val):
        fullfile = join(self.path, filename)
        with open(fullfile, 'a') as f:
            f.write('%e\n' % val)





def getpar(key, file='DATA/Par_file', sep='=', cast=str):
    """ Reads parameter from text file
    """
    val = None
    with open(file, 'r') as f:
        # read line by line
        for line in f:
            if line.find(key) == 0:
                # read key
                key, val = _split(line, sep)
                if not key:
                    continue
                # read val
                val, _ = _split(val, '#')
                val.strip()
                break

    if val:
        if cast == float:
            val = val.replace('d', 'e')
        return cast(val)

    else:
        print('Not found in parameter file: %s\n' % key)
        raise Exception


def setpar(key, val, filename='DATA/Par_file', path='.', sep='='):
    """ Writes parameter to text file
    """

    val = str(val)

    # read line by line
    with open(path +'/'+ filename, 'r') as file:
        lines = []
        for line in file:
            if line.find(key) == 0:
                # read key
                key, _ = _split(line, sep)
                # read comment
                _, comment = _split(line, '#')
                n = len(line) - len(key) - len(val) - len(comment) - 2
                # replace line
                if comment:
                    line = _merge(key, sep, val, ' '*n, '#', comment)
                else:
                    line = _merge(key, sep, str(val), '\n')
            lines.append(line)

    # write file
    with open(path +'/'+ filename, 'w') as file:
        file.writelines(lines)


def setpararray(key, val, filename='DATA/Par_file', path='/', sep='='):
    """ Writes parameter to SPECFEM parfile
    """

    #val = str(val)

    # read line by line
    with open(path +'/'+ filename, 'r') as file:
        lines = []
        k = 0
        for line in file:
            if line.find(key) == 0:
                # read key
                key, _ = _split(line, sep)
                # read comment
                _, comment = _split(line, '#')
                n = len(line) - len(key) - len(str(val)) - len(comment) - 2
                # replace line
                if comment:
                    line = _merge(key, sep, str(val[0]), ' '*n, '#', comment)
                else:
                    line = _merge(key, sep, str(val[k]), '\n')
                    k+=1
            lines.append(line)

    # write file
    with open(path +'/'+ filename, 'w') as file:
        file.writelines(lines)



### utility functions

def _split(str, sep):
    n = str.find(sep)
    if n >= 0:
        return str[:n], str[n + len(sep):]
    else:
        return str, ''


def _merge(*parts):
    return ''.join(parts)


def _zip(keys, vals):
    return list(zip(iterable(keys), iterable(vals)))
