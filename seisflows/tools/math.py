
from copy import copy

import os

import numpy as np
import scipy.signal as _signal
import scipy.interpolate as _interp
from scipy.signal import hilbert as analytic
from scipy import ndimage


def gauss2(X, Y, mu, sigma, normalize=True):
    """ Evaluates Gaussian over points of X,Y
    """
    # evaluates Gaussian over X,Y
    D = sigma[0, 0]*sigma[1, 1] - sigma[0, 1]*sigma[1, 0]
    B = np.linalg.inv(sigma)
    X = X - mu[0]
    Y = Y - mu[1]
    Z = B[0, 0]*X**2. + B[0, 1]*X*Y + B[1, 0]*X*Y + B[1, 1]*Y**2.
    Z = np.exp(-0.5*Z)

    if normalize:
        Z *= (2.*np.pi*np.sqrt(D))**(-1.)

    return Z


def backtrack2(f0, g0, x1, f1, b1=0.1, b2=0.5):
    """ Safeguarded parabolic backtrack
    """
    # parabolic backtrack
    x2 = -g0*x1**2/(2*(f1-f0-g0*x1))

    # apply safeguards
    if x2 > b2*x1:
        x2 = b2*x1
    elif x2 < b1*x1:
        x2 = b1*x1
    return x2


def backtrack3(f0, g0, x1, f1, x2, f2):
    """ Safeguarded cubic backtrack
    """
    raise NotImplementedError


def polyfit2(x, f):
    # parabolic fit
    i = np.argmin(f)
    p = np.polyfit(x[i-1:i+2], f[i-1:i+2], 2)

    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print(-1)
        raise Exception()


def lsq2(x, f):
    # parabolic least squares fit
    p = np.polyfit(x, f, 2)
    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print(-1)
        raise Exception()


def angle(x,y):
    xy = dot(x,y)
    xx = dot(x,x)
    yy = dot(y,y)
    return np.arccos(xy/(xx*yy)**0.5)


def dot(x,y):
    return np.dot(
        np.squeeze(x),
        np.squeeze(y))


def hilbert(w):
    return np.imag(analytic(w))


infinity = np.inf



### finite difference

def nabla(V, h=[]):
    """ Returns sum of first-order spatial derivatives of a function defined on
        a 2D rectangular grid; generalizes Laplacian
    """
    W = np.zeros(V.shape)

    if h==[]:
       h = np.ones((V.ndim, 1))

    # interior
    W[1:-1,1:-1] += (V[1:-1,2:] - V[1:-1,:-2])/(2.*h[0])
    W[1:-1,1:-1] += (V[2:,1:-1] - V[:-2,1:-1])/(2.*h[1])

    # top/bottom edges
    W[0,1:-1] = (V[1,1:-1] - V[0,1:-1])/h[1] + (V[0,2:] - V[0,:-2])/(2.*h[0])
    W[-1,1:-1] = (V[-1,1:-1] - V[-2,1:-1])/h[1] + (V[-1,2:] - V[-1,:-2])/(2.*h[0])

    # left/right edges
    W[1:-1,0] = (V[2:,0] - V[:-2,0])/(2.*h[1]) + (V[1:-1,1] - V[1:-1,0])/h[0]
    W[1:-1,-1] = (V[2:,-1] - V[:-2,-1])/(2.*h[1]) + (V[1:-1,-1] - V[1:-1,-2])/h[0]

    # corners
    W[0,0] = (V[1,0] - V[0,0])/h[1] + (V[0,1] - V[0,0])/h[0]
    W[0,-1] = (V[1,-1] - V[0,-1])/h[1] + (V[0,-2] - V[0,-1])/h[0]
    W[-1,0] = (V[-2,0] - V[-1,0])/h[1] + (V[-1,1] - V[-1,0])/h[0]
    W[-1,-1] = (V[-1,-1] - V[-2,-1])/h[1] + (V[-1,-1] - V[-1,-2])/h[0]

    return W



def nabla2(V, h=[]):
    """ Returns sum of second-order spatial derivatives of a function defined on
        a 2D rectangular grid; generalizes Laplacian
    """
    W = np.zeros(V.shape)

    if h==[]:
       h = np.ones((V.ndim, 1))

    # interior
    W[1:-1,1:-1] += (V[1:-1,2:] -2.*V[1:-1,1:-1] + V[1:-1,:-2])/h[0]**2
    W[1:-1,1:-1] += (V[2:,1:-1] -2.*V[1:-1,1:-1] + V[:-2,1:-1])/h[1]**2

    # left/right edges
    W[0,1:-1] = W[1,1:-1]
    W[-1,1:-1] = W[-2,1:-1]

    # top/bottom edges
    W[0,1:-1] = W[1,1:-1]
    W[-1,1:-1] = W[-2,1:-1]

    # corners
    W[0,0] = (W[0,1] + W[1,0])/2
    W[0,-1] = (W[0,-2] + W[1,-1])/2
    W[-1,0] = (W[-1,1] + W[-2,0])/2
    W[-1,-1] = (W[-1,-2] + W[-2,-1])/2

    return W


def grad(V, h=[]):
    """ Evaluates derivatives on a 2D rectangular grid
    """
    ny, nx = V.shape

    X = np.zeros((ny, nx))
    Y = np.zeros((ny, nx))

    if h==[]:
       h = np.ones((V.ndim, 1))

    # interior
    X[:,1:-1] = (V[:,2:] - V[:,:-2])/(2.*h[0])
    Y[1:-1,:] = (V[2:,:] - V[:-2,:])/(2.*h[1])

    # left/right edges
    X[:,0] = (V[:,1] - V[:,0])/h[1]
    X[:,-1] = (V[:,-1] - V[:,-2])/h[1]

    # top/bottom edges
    Y[0,:] = (V[1,:] - V[0,:])/h[0]
    Y[-1,:] = (V[-1,:] - V[-2,:])/h[0]

    return X,Y


def tv(Z, h=[], epsilon=1.e-6):
    nrow = Z.shape[0]
    ncol = Z.shape[1]

    Zx = (Z[:,1:] - Z[:,:-1])/h[0]
    Zy = (Z[1:,:] - Z[:-1,:])/h[1]

    top = np.zeros((nrow, ncol))
    bot = np.zeros((nrow, ncol))

    top[:,1:] += Zx
    top[1:,:] += Zy
    top[ :,-1] += Zx[:,-1]
    top[-1, :] += Zy[-1,:]

    bot[:,1:] += Zx**2
    bot[1:,:] += Zy**2
    bot[ :,-1] += Zx[:,-1]**2
    bot[-1, :] += Zy[-1,:]**2

    return top/(bot + epsilon*bot.max())**0.5

def get_mute_ratio(radius, maxradius, minratio):
    # ratio = np.zeros(len(radius))

    ratio = 1 - (1-minratio)*np.cos(np.pi/(2*maxradius)*radius)

    ratio[np.where(radius>maxradius)] = 1
    return ratio

def mysmooth(v, x, z, nx, nz, smo, xlim=None, zlim=None, topolim=None):
    """ Interpolates from an unstructured coordinates (mesh) to a structured
        coordinates (grid)
    """
    if xlim is None:
        xmax = x.max()
        xmin = x.min()
        lx = xmax - xmin
    else:
        xmax = xlim[1]
        xmin = xlim[0]
        lx = xmax - xmin

    if zlim is None:
        zmax = z.max()
        zmin = z.min()
        lz = zmax - zmin
    else:
        zmax = zlim[1]
        zmin = zlim[0]
        lz = zmax - zmin

    mesh = _stack(x, z)

    dx = lx / (nx-1)
    dz = lz / (nz-1)

    # construct structured grid
    x_grid = np.linspace(xmin, xmax, nx)
    if topolim is not None:
        z1 = np.linspace(zmin, topolim, nz)
        z2 = np.linspace(topolim, zmax, nz * 2)
        z_grid = np.hstack((z1, z2))
    else:
        z_grid = np.linspace(zmin, zmax, nz)

    X, Z = np.meshgrid(x_grid, z_grid)
    grid = _stack(X.flatten(), Z.flatten())

    # interpolate to structured grid
    v_tmp = myinterpolate(mesh, v, grid)
    V = np.reshape(v_tmp, (len(z_grid), len(x_grid)))

    print(np.max(V))
    print(np.min(V))
    print('dx:%f'%dx)
    print('dz:%f'%dz)

    V = ndimage.gaussian_filter(V,sigma=smo)

    mesh1 = _stack(X.flatten(), Z.flatten())
    grid1 = _stack(x, z)
    v1 = V.flatten()
    v_new = myinterpolate(mesh1, v1, grid1)

    return v_new

def myinterpolate(mesh,v, grid):
    V = _interp.griddata(mesh, v, grid, 'linear')
    if np.any(np.isnan(V)):
        W = _interp.griddata(mesh, v, grid, 'nearest')
        for i in np.where(np.isnan(V)):
            V[i] = W[i]

    return V

def _stack(*args):
    return np.column_stack(args)