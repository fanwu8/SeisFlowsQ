import re
import sys
from glob import glob
from os import getcwd
from os.path import abspath, exists, join

import numpy as np
import matplotlib
import pylab
import scipy.interpolate


def read_fortran(filename):
    """ Reads Fortran style binary data and returns a numpy array.
    """
    with open(filename, 'rb') as f:
        # read size of record
        f.seek(0)
        n = np.fromfile(f, dtype='int32', count=1)[0]

        # read contents of record
        f.seek(4)
        v = np.fromfile(f, dtype='float32')

    return v[:-1]


def mesh2grid(v, x, z):
    """ Interpolates from an unstructured coordinates (mesh) to a structured
        coordinates (grid)
    """
    lx = x.max() - x.min()
    lz = z.max() - z.min()
    nn = v.size
    mesh = _stack(x, z)

    nx = int(np.around(np.sqrt(nn * lx / lz)))
    nz = int(np.around(np.sqrt(nn * lz / lx)))
    dx = lx / nx
    dz = lz / nz

    # construct structured grid
    x = np.linspace(x.min(), x.max(), nx)
    z = np.linspace(z.min(), z.max(), nz)
    X, Z = np.meshgrid(x, z)
    grid = _stack(X.flatten(), Z.flatten())

    # interpolate to structured grid
    V = scipy.interpolate.griddata(mesh, v, grid, 'linear')

    # workaround edge issues
    if np.any(np.isnan(V)):
        W = scipy.interpolate.griddata(mesh, v, grid, 'nearest')
        for i in np.where(np.isnan(V)):
            V[i] = W[i]

    return np.reshape(V, (int(nz), int(nx))), x, z


def exist_files(names):
    """Wrapper for os.path.exists"""
    for name in names:
        if not name:
            return False
        elif not isinstance(name, str):
            raise TypeError
        elif not exists(name):
            return False
    else:
        return True


def get_coord_files(dir):
    for files in [
        [abspath(join(dir, subdir, 'proc000000_x.bin')),
         abspath(join(dir, subdir, 'proc000000_z.bin'))]
        for subdir in [
            '.',
            '../model_init',
            '../../model_init',
            '../../../model_init',
            '../../../../model_init',
            '../output/model_init',
            '../../output/model_init',
            '../../../output/model_init',
            '../../../../output/model_init',
        ]]:
        if exist_files(files):
            return files
    else:
        raise Exception("coordinate files not found")


def get_data_files(dir):
    files = []
    for file in glob(dir + '/' + 'proc*.bin*'):
        if 'jacobian' in file or \
                'ibool' in file or \
                'x.bin' in file or \
                'z.bin' in file:
            continue
        files += [file]
    return sorted(files)


def get_title(file):
    pattern = '.*proc.*_(.*).bin'
    return re.match(pattern, file).group(1)


def _stack(*args):
    return np.column_stack(args)


if __name__ == '__main__':
    """ Plots data on 2-D unstructured mesh
      Original code by Ryan Modrak:
      http://tigress-web.princeton.edu/~rmodrak/visualize/plot2d

      Can be used to plot models or kernels

      SYNTAX
          spplot  <model_dir>
          e.g. spplot projects/example_checker_sh/model_true
    """

    # get coornidate files
    x_coords_file, z_coords_file = get_coord_files(sys.argv[1])
    x = read_fortran(x_coords_file)
    z = read_fortran(z_coords_file)

    files = get_data_files(sys.argv[1])

    for database_file in files:
        # component name
        if len(sys.argv) >= 3 and sys.argv[2] not in database_file:
            continue

        # read database file
        v = read_fortran(database_file)

        # check mesh dimensions
        assert x.shape == z.shape == v.shape, 'Inconsistent mesh dimensions.'

        # interpolate to uniform rectangular grid
        V, X, Z = mesh2grid(v, x, z)

        # display figure
        pylab.pcolor(X, Z, V)
        locs = np.arange(X.min(), X.max() + 1, (X.max() - X.min()) / 5)
        pylab.xticks(locs)
        pylab.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        locs = np.arange(Z.min(), Z.max() + 1, (Z.max() - Z.min()) / 5)
        pylab.yticks(locs)
        pylab.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        pylab.colorbar()
        pylab.xlabel('x')
        pylab.ylabel('z')
        pylab.title(get_title(database_file))
        pylab.gca().invert_yaxis()
        pylab.set_cmap('seismic')

        pylab.show()