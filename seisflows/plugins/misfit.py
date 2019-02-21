
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
    ioff = (np.argmax(cc)-nt+1)*dt
    if ioff <= 0:
        wrsd = syn[ioff:] - obs[:-ioff]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope2(syn, obs, nt, dt, eps=0.):
    # envelope amplitude ratio
    # (Yuan et al 2015, eq B-1)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    raise NotImplementedError


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


def Phase_freq2(syn, nt, dt,ft_obs,sff_freq,sff_freq_true, freq_mask):
    # waveform difference in the frequency domain, considering orthogonal frequencies
    nstep = len(syn)
    wadj = 0.0 #np.zeros(nstep)
    period = PAR.PERIOD
    freq_min = PAR.BW_L
    freq_max = PAR.BW_H 
    #create a frequential mask
    freq  = fftfreq(period,dt)
    m = loadnpy(PATH.ORTHO + '/freq_idx')
    ft_syn = fft(syn[-period:])[m]
    obs = ft_syn /  ( (ft_obs) * sff_freq_true / sff_freq )

    phase = np.vectorize(cmath.phase)
    phase_obs = phase(obs)
#    phase_syn = phase(ft_syn)
    wadj = (freq_mask * phase_obs**2).sum() 
    return wadj