
# used by the PREPROCESS class and specified by the MISFIT parameter



import sys
import numpy as _np
import cmath
from scipy.signal import hilbert as _analytic
from scipy.fftpack import fft, ifft, fftfreq

from seisflows.tools.array import loadnpy
from seisflows.plugins import misfit
from seisflows.tools.math import hilbert as _hilbert

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

### adjoint traces generators


def Waveform(syn, obs, nt, dt):
    # waveform difference
    # (Tromp et al 2005, eq 9)
    wadj = syn - obs
    return wadj


def Envelope(syn, obs, nt, dt, eps=0.05):
    # envelope difference
    # (Yuan et al 2015, eq 16)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    etmp = (esyn - eobs)/(esyn + eps*esyn.max())
    wadj = etmp*syn - _np.imag(_analytic(etmp*_np.imag(_analytic(syn))))
    return wadj


def InstantaneousPhase(syn, obs, nt, dt, eps=0.05):
    # instantaneous phase 
    # (Bozdag et al 2011, eq 27)
    r = _np.real(_analytic(syn))
    i = _np.imag(_analytic(syn))
    phi_syn = _np.arctan2(i,r)

    r = _np.real(_analytic(obs))
    i = _np.imag(_analytic(obs))
    phi_obs = _np.arctan2(i,r)

    phi_rsd = phi_syn - phi_obs
    esyn = abs(_analytic(syn))
    emax = max(esyn**2.)

    wadj = phi_rsd*_np.imag(_analytic(syn))/(esyn**2. + eps*emax) + \
           _np.imag(_analytic(phi_rsd * syn/(esyn**2. + eps*emax)))

    return wadj


def Traveltime(syn, obs, nt, dt):
    # cross correlation traveltime
    # (Tromp et al 2005, eq 45)
    wadj = _np.zeros(nt)
    wadj[1:-1] = (syn[2:] - syn[0:-2])/(2.*dt)
    wadj *= 1./(sum(wadj*wadj)*dt)
    wadj *= misfit.Traveltime(syn,obs,nt,dt)
    return wadj


def TraveltimeInexact(syn, obs, nt, dt):
    # must faster but possibly inaccurate
    wadj = _np.zeros(nt)
    wadj[1:-1] = (syn[2:] - syn[0:-2])/(2.*dt)
    wadj *= 1./(sum(wadj*wadj)*dt)
    wadj *= misfit.TraveltimeInexact(syn,obs,nt,dt)
    return wadj


def Amplitude(syn, obs, nt, dt):
    # cross correlation amplitude
    wadj = 1./(sum(syn*syn)*dt) * syn
    wadj *= misfit.Amplitude(syn,obs,nt,dt)
    return wadj




def Envelope2(syn, obs, nt, dt, eps=0.):
    # envelope amplitude ratio
    # (Yuan et al 2015, eqs B-2, B-3)
    raise NotImplementedError


def Envelope3(syn, obs, nt, dt, eps=0.):
    # envelope lag
    # (Yuan et al 2015, eqs B-2, B-5)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))

    erat = _np.zeros(nt)
    erat[1:-1] = (esyn[2:] - esyn[0:-2])/(2.*dt)
    erat[1:-1] /= esyn[1:-1]
    erat *= misfit.Envelope3(syn, obs, nt, dt)

    wadj = -erat*syn + _hilbert(erat*_hilbert(esyn))
    return wadj


def InstantaneousPhase2(syn, obs, nt, dt, eps=0.):
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))

    esyn1 = esyn + eps*max(esyn)
    eobs1 = eobs + eps*max(eobs)
    esyn3 = esyn**3 + eps*max(esyn**3)

    diff1 = syn/(esyn1) - obs/(eobs1)
    diff2 = _hilbert(syn)/esyn1 - _hilbert(obs)/eobs1

    part1 = diff1*_hilbert(syn)**2/esyn3 - diff2*syn*_hilbert(syn)/esyn3
    part2 = diff1*syn*_hilbert(syn)/esyn3 - diff2*syn**2/esyn3

    wadj = part1 + _hilbert(part2)
    return wadj



### migration

def Displacement(syn, obs, nt, dt):
    return obs

def Velocity(syn, obs, nt, dt):
    adj[1:-1] = (obs[2:] - obs[0:-2])/(2.*dt)
    return adj

def Acceleration(syn, obs, nt, dt):
    adj[1:-1] = (-obs[2:] + 2.*obs[1:-1] - obs[0:-2])/(2.*dt)
    return adj



def Phase_freq2(syn,nt,dt,ft_obs,sff_freq,sff_freq_true,freq_mask):
    # waveform difference
    # (Tromp et al 2005, eq 9)
    wadj = _np.zeros(len(syn))
    period = PAR.PERIOD

    m = loadnpy(PATH.ORTHO + '/freq_idx')
    freq_loc = loadnpy(PATH.ORTHO + '/freq')
    
    ft_syn = fft(syn[-period:])[m]
    obs = ft_syn / ( ft_obs * sff_freq_true / sff_freq )

    phase = _np.vectorize(cmath.phase)

    phase_obs = phase(obs)
    phase_syn = phase(ft_syn)
    #amp_syn = abs(ft_syn)/ (abs(freq_loc)/500.0)
    inv_fft = _np.zeros(period,dtype=complex)
    inv_fft[m] = (freq_mask * phase_obs  ) *  _np.vectorize(_np.complex)(-_np.sin(phase_syn),_np.cos(phase_syn))
    wadj[-period:] = - 10e15*ifft(inv_fft).real
    #repeat the periodic signal
    j=1
    while ((j+1)*period < len(syn) ):
      wadj[ -(j+1)*period : -j*period ] = wadj[-period:]
      j+=1
    if (j==1):
      wadj[:-period] = wadj[-(len(syn)-period):]
    else:
      wadj[:-(j-1)*period] = wadj[-(len(syn)-(j-1)*period):]

    return wadj