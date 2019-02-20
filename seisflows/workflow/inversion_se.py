
import sys
import numpy as np
import random

from os.path import join

from seisflows.tools import unix
from seisflows.workflow.inversion import inversion

from scipy.fftpack import fft, fftfreq
from seisflows.tools.array import loadnpy, savenpy
# from seisflows.tools.seismic import setpararray

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']


class inversion_se(inversion):
    """ Waveform inversion with source encoding
    """

    def check(self):
        super().check()

        # get random source
        if 'RANDOM_OVER_IT' not in PAR:
            setattr(PAR, 'RANDOM_OVER_IT', 1)

        # increase frequency over iterations
        if 'FREQ_INCREASE_PER_IT' not in PAR:
            setattr(PAR, 'FREQ_INCREASE_PER_IT', 0)

        # maximum frequency shift over iterations
        if 'MAX_FREQ_SHIFT' not in PAR:
            setattr(PAR, 'MAX_FREQ_SHIFT', None)

        # number of frequency per event
        if 'NFREQ_PER_EVENT' not in PAR:
            setattr(PAR, 'NFREQ_PER_EVENT', 1)

        # default number of super source
        if 'NSRC' not in PAR:
            setattr(PAR, 'NSRC', 1)


    def setup(self):
        super().setup()

        unix.mkdir(join(PATH.FUNC, 'residuals'))
        unix.mkdir(join(PATH.GRAD, 'residuals'))

    def initialize(self):
        """ Prepares for next model update iteration
        """
        self.write_model(path=PATH.GRAD, suffix='new')

        if PAR.RANDOM_OVER_IT or optimize.iter == 1:
            self.get_random_frequencies()

        print('Generating synthetics')
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.GRAD)

        self.write_misfit(path=PATH.GRAD, suffix='new')

    def clean(self):
        super().clean()

        unix.mkdir(join(PATH.FUNC, 'residuals'))
        unix.mkdir(join(PATH.GRAD, 'residuals'))

    def get_random_frequencies(self):
        """ Randomly assign a unique frequency for each source
        """
        period = PAR.PERIOD
        dt = PAR.DT
        nt = PAR.NT
        nrec = PAR.NREC
        nevt = PAR.NEVT
        nfpe = PAR.NFREQ_PER_EVENT
        nsrc = nevt * nfpe

        # get the number of relevant frequencies
        freq_min = float(PAR.BW_L)
        freq_max = float(PAR.BW_H)

        # read data processed py ortho
        freq_idx = loadnpy(PATH.ORTHO + '/freq_idx')
        freq = loadnpy(PATH.ORTHO + '/freq')
        ft_stf = loadnpy(PATH.ORTHO + '/ft_stf')
        ft_obs = loadnpy(PATH.ORTHO + '/ft_obs')
        
        nfreq = len(freq_idx)
        ntrace = ft_obs.shape[3]

        # declaring arrays
        ft_obs_se = np.zeros((nfreq, nrec, ntrace), dtype=complex)
        ft_stf_se = np.zeros((nfreq), dtype=complex)
        ft_stf_se_sinus = np.zeros((nfreq), dtype=complex)
        ft_stf_sinus = np.zeros((nfreq, nsrc), dtype=complex)
        
        # frequency processing
        # TODO freq_mask
        freq_mask_se = np.ones((nfreq, nrec))
        freq_shift = (optimize.iter - 1) * PAR.FREQ_INCREASE_PER_IT
        if PAR.MAX_FREQ_SHIFT != None:
            freq_shift = min(freq_shift, PAR.MAX_FREQ_SHIFT)

        # random frequency
        freq_range = np.linspace(freq_min + freq_shift, freq_max + freq_shift, nsrc + 1)[:-1]
        freq_thresh = (freq_max - freq_min) / nsrc / 20
        rdm_idx = random.sample(range(0, nsrc), nsrc)
        freq_rdm = freq_range[rdm_idx]

        # get individual frequency
        for isrc in range(nsrc):
            f0 = freq_rdm[isrc]
            T = 2 * np.pi * dt * np.linspace(0, nt - 1, nt) * f0
            stf_sinus = 1000 * np.sin(T)
            ft_stf_sinus[:, isrc] = fft(stf_sinus[-period:])[freq_idx]

        # encode frequencies
        n=0
        for ifpe in range(nfpe):
            for ievt in range(nevt):
                ife = ifpe * nevt + ievt
                for ifreq in range(nfreq):
                    # FIXME frequency threshold
                    if abs(abs(freq_rdm[ife]) - abs(freq[ifreq])) < freq_thresh:
                        n+=1
                        ft_obs_se[ifreq, :, :]  = ft_obs[ifreq, ievt, :, :]
                        # TODO freq_mask
                        ft_stf_se[ifreq] = ft_stf[ifreq, ievt]
                        ft_stf_se_sinus[ifreq]  = ft_stf_sinus[ifreq, rdm_idx[ife]]
        # assert that random frequency is a subset of ferquency band
        if (n != nsrc * 2):
            print('Warning: discrete frequency mismatch')

        savenpy(PATH.ORTHO +'/ft_obs_se', ft_obs_se)
        savenpy(PATH.ORTHO +'/ft_stf_se', ft_stf_se)
        savenpy(PATH.ORTHO +'/ft_stf_se_sinus', ft_stf_se_sinus)
        savenpy(PATH.ORTHO +'/freq_mask_se', freq_mask_se)

        for source_name in solver.source_names_all[1:]:
            src = PATH.SPECFEM_DATA + '/' + solver.source_prefix +'_'+ source_name
            dst = PATH.SOLVER + '/000000/DATA/' + solver.source_prefix
            unix.cat(src, dst)
            # setpararray('f0', f_s, filename= PATH.SOLVER + '/000000/DATA/SOURCE')