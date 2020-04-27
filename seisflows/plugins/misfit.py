
# used by the PREPROCESS class and specified by the MISFIT parameter


import sys
import numpy as np
import cmath
from scipy.signal import hilbert as _analytic
from scipy.fftpack import fft, fftfreq
from seisflows.tools.array import loadnpy

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

def Waveform(syn, obs, nt, dt):
    # waveform difference
    wrsd = syn-obs
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope(syn, obs, nt, dt, eps=0.05):
    # envelope difference
    # (Yuan et al 2015, eq 9)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    ersd = esyn-eobs
    return np.sqrt(np.sum(ersd*ersd*dt))


def InstantaneousPhase(syn, obs, nt, dt, eps=0.05):
    # instantaneous phase 
    # from Bozdag et al. 2011

    r = np.real(_analytic(syn))
    i = np.imag(_analytic(syn))
    phi_syn = np.arctan2(i,r)

    r = np.real(_analytic(obs))
    i = np.imag(_analytic(obs))
    phi_obs = np.arctan2(i,r)

    phi_rsd = phi_syn - phi_obs
    return np.sqrt(np.sum(phi_rsd*phi_rsd*dt))


def Traveltime(syn, obs, nt, dt):
    cc = abs(np.convolve(obs, np.flipud(syn)))
    return (np.argmax(cc)-nt+1)*dt


def TraveltimeInexact(syn, obs, nt, dt):
    # much faster but possibly inaccurate
    it = np.argmax(syn)
    jt = np.argmax(obs)
    return (jt-it)*dt


def Amplitude(syn, obs, nt, dt):
    # cross correlation amplitude
    cc = abs(np.convolve(obs, np.flipud(syn)))
    ioff = np.argmax(cc)-nt+1
    # print(ioff)
    if ioff < 0:
        wrsd = syn[-ioff:] - obs[:ioff]
        syn0 = syn[-ioff:]
        obs0 = obs[:ioff]
    elif ioff == 0:
        syn0 = syn[:]
        obs0 = obs[:]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]
        syn0 = syn[:-ioff]
        obs0 = obs[ioff:]

    A_obs = np.sqrt(np.sum(obs0*obs0*dt))
    A_syn = np.sqrt(np.sum(syn0*syn0*dt))

    if A_obs < 1e-20:
        return 0
    if A_syn < 1e-20:
        return 0


    return A_syn/A_obs - 1
    # return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope2(syn, obs, nt, dt, eps=0.):
    # envelope amplitude ratio
    # (Yuan et al 2015, eq B-1)
    cc = abs(np.convolve(obs, np.flipud(syn)))
    ioff = np.argmax(cc)-nt+1
    # print(ioff)
    if ioff < 0:
        wrsd = syn[-ioff:] - obs[:ioff]
        syn0 = syn[-ioff:]
        obs0 = obs[:ioff]
    elif ioff == 0:
        syn0 = syn[:]
        obs0 = obs[:]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]
        syn0 = syn[:-ioff]
        obs0 = obs[ioff:]

    esyn = abs(_analytic(syn0))
    eobs = abs(_analytic(obs0))
    ersd = esyn-eobs
    return np.sqrt(np.sum(ersd*ersd*dt))



def Envelope3(syn, obs, nt, dt, eps=0.):
    # envelope cross-correlation lag
    # (Yuan et al 2015, eqs B-4)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    return Traveltime(esyn, eobs, nt, dt)


def InstantaneousPhase2(syn, obs, nt, dt, eps=0.):
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))

    esyn1 = esyn + eps*max(esyn)
    eobs1 = eobs + eps*max(eobs)

    diff = syn/esyn1 - obs/eobs1

    return np.sqrt(np.sum(diff*diff*dt))



def Displacement(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')

def Velocity(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')

def Acceleration(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')


def Phase2_se(syn, nt, dt,ft_obs, freq_mask):
    # waveform difference in the frequency domain, considering orthogonal frequencies
    nstep = len(syn)
    wadj = 0.0 #np.zeros(nstep)
    ntpss = PAR.NTPSS
    #create a frequential mask
    m = loadnpy(PATH.ORTHO + '/freq_idx')
    ft_syn = fft(syn[-ntpss:])[m]
    obs = ft_syn / ft_obs

    phase = np.vectorize(cmath.phase)
    phase_obs = phase(obs)

    wadj = (freq_mask * phase_obs**2).sum()
    return wadj

def GCE(syn,obs,nt,dt):
    if np.max(np.abs(obs)) < 1e-18:
        # print(_np.max(_np.abs(obs)))
        return 0
    if np.max(np.abs(syn)) < 1e-18:
        # print(_np.max(_np.abs(syn)))
        return 0


    syn_n = syn / np.sqrt(np.sum(syn*syn*dt))
    obs_n = obs / np.sqrt(np.sum(obs*obs*dt))
    return np.sqrt(-np.sum(syn_n*obs_n*dt) + 1 + 5e-16)

#TODO: make mesh denser, increase time to 3s