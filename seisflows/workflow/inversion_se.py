
import sys
import numpy as np
import random

from os.path import join

from seisflows.tools import unix
from seisflows.workflow.inversion import inversion

from scipy.fftpack import fft, fftfreq
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.seismic import setpar, setpararray

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
        sff_obs = loadnpy(PATH.ORTHO + '/sff_obs')
        ft_obs = loadnpy(PATH.ORTHO + '/ft_obs')
        
        nfreq = len(freq_idx)
        # ntrace = ft_obs.shape[3]

        # declaring arrays
        ft_obs_se = np.zeros((nfreq, nrec), dtype=complex)
        sff_obs_se = np.zeros((nfreq), dtype=complex)
        sff_syn_se = np.zeros((nfreq), dtype=complex)
        sff_syn = np.zeros((nfreq), dtype=complex)
        
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

        # assign frequencies
        stf_filenames = [None] * nsrc
        for ifpe in range(nfpe):
            for ievt in range(nevt):
                isrc = ifpe * nevt + ievt
                f0 = freq_rdm[isrc]

                # get sinus source time function
                T = 2 * np.pi * dt * np.linspace(0, nt - 1, nt) * f0
                sinus = 1000 * np.sin(T)
                sff_syn = fft(sinus[-period:])[freq_idx]

                # find and encode matching frequencies
                nmatch = 0
                for ifreq in range(nfreq):
                    if abs(abs(f0) - abs(freq[ifreq])) < freq_thresh:
                        nmatch += 1
                        ft_obs_se[ifreq, :]  = ft_obs[ifreq, ievt, :]
                        # TODO freq_mask
                        sff_obs_se[ifreq] = sff_obs[ifreq, ievt]
                        sff_syn_se[ifreq] = sff_syn[ifreq]

                # assert that source encoding frequency is a subset of ferquency band
                if nmatch != 2:
                    print('Warning: descrete frequency is not a subset of frequency band')

                filename = PATH.SOLVER + '/000000/DATA/STF_' + str(ievt) + '_' + str(ifpe)
                stf_filenames[isrc] = filename

                if optimize.iter == 1:
                    stf_syn = np.zeros([nt, 2])
                    stf_syn[:, 0] = T
                    stf_syn[:, 1] = sinus
                    np.savetxt(filename, stf_syn)


        savenpy(PATH.ORTHO +'/ft_obs_se', ft_obs_se)
        savenpy(PATH.ORTHO +'/sff_obs_se', sff_obs_se)
        savenpy(PATH.ORTHO +'/sff_syn_se', sff_syn_se)
        savenpy(PATH.ORTHO +'/freq_mask_se', freq_mask_se)

        dst = PATH.SOLVER + '/000000/DATA/' + solver.source_prefix
        unix.rm(dst)
        for ifpe in range(nfpe):
            for ievt in range(nevt):
                source_name = solver.source_names_all[ievt]
                src = PATH.SPECFEM_DATA + '/' + solver.source_prefix +'_'+ source_name
                unix.cat(src, dst)

        setpararray('time_function_type', np.ones(nsrc).astype(int) * 8, filename= dst)
        setpararray('f0', freq_rdm, filename= dst)
        setpararray('name_of_source_file', stf_filenames, filename= dst)
        
        # set number of sources fo solver
        if optimize.iter == 1:
            setpar('NSOURCES', nsrc, 'DATA/Par_file', PATH.SOLVER + '/000000')
