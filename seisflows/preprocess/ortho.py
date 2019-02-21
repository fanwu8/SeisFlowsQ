
import sys
import numpy as np
import obspy
from os.path import join

from seisflows.tools import msg, unix
from seisflows.tools.tools import exists, getset
from seisflows.config import ParameterError, custom_import

from seisflows.plugins import adjoint, misfit, readers, writers
from seisflows.tools import signal
from seisflows.tools.array import loadnpy, savenpy

from scipy.fftpack import fft, fftfreq

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class ortho(custom_import('preprocess', 'base')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        super().check()

        # scratch directory
        if 'ORTHO' not in PATH:
            setattr(PATH, 'ORTHO', join(PATH.SCRATCH, 'ortho'))


    def setup(self):
        """ Sets up data preprocessing machinery
        """
        super().setup()

        # prepare scratch directory
        unix.mkdir(PATH.ORTHO)

        # get data file names from solver
        solver = sys.modules['seisflows_solver']

        nevt = PAR.NEVT
        period = PAR.PERIOD
        dt = PAR.DT
        nrec = PAR.NREC
        # ntrace = len(solver.data_filenames)

        # get the number of relevant frequencies
        freq_min = float(PAR.BW_L)
        freq_max = float(PAR.BW_H)
        
        #create a mask on relevant frequencies
        freq_full = fftfreq(period, dt)
        freq_thresh = (freq_max - freq_min) / PAR.NFREQ_PER_EVENT / nevt / 20
        freq_idx = np.squeeze(np.where((freq_min <= abs(freq_full)) & (abs(freq_full) < freq_max - freq_thresh)))
        freq = freq_full[freq_idx]
        nfreq = len(freq_idx)
        print('Number of frequencies considered: ' +str(nfreq)+' / '+str(len(freq_full)))

        # converts time data to Fourier domain
        ft_stf = np.zeros((nfreq, nevt), dtype=complex)
        ft_obs = np.zeros((nfreq, nevt, nrec), dtype=complex) # TODO ntrace

        for isrc in range(nevt):
            source_name = solver.source_names_all[isrc]
            stf_file = solver.stf_files_all[isrc]
            with open(stf_file) as f:
                lines = f.readlines()
                stf_obs = []
                for line in lines:
                    stf_obs.append(float(line.split()[1]))

            ft_stf[:, isrc] = fft(stf_obs, n=period)[freq_idx]
            # for itrace in range(ntrace):
            #     trace = self.reader(PATH.DATA + '/' + source_name, solver.data_filenames[itrace])
            #     for irec in range(nrec):
            #         ft_obs[:, isrc, irec, itrace] = fft(trace[irec].data, n=period)[freq_idx]
            for irec in range(nrec):
                trace = self.reader(PATH.DATA + '/' + source_name, solver.data_filenames[0])
                ft_obs[:, isrc, irec] = fft(trace[irec].data, n=period)[freq_idx]
        
        self.save('freq_idx', freq_idx)
        self.save('freq', freq)
        self.save('ft_stf', ft_stf)
        self.save('ft_obs', ft_obs)
 
    def prepare_eval_grad(self, path='.',wat=True):
        """ Prepares solver for gradient evaluation by writing residuals and
          adjoint traces

          INPUT
            PATH - directory containing observed and synthetic seismic data
        """
        solver = sys.modules['seisflows_solver']
        for filename in solver.data_filenames:
            obs = self.reader(path+'/'+'traces/obs', filename)
            syn = self.reader(path+'/'+'traces/syn', filename)
            nt, dt, _ = self.get_time_scheme(syn)

            if PAR.MISFIT:
                self.write_residuals(path, syn, obs)
            if wat:
              self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs, filename)
              if PAR.ATTENUATION =='yes':
                self.write_adjoint_traces(path+'/'+'traces/adj_att', syn, obs, filename,att='Yes')

    def write_residuals(self, path, syn, obs):
        """ Computes residuals from observations and synthetics

          INPUT
            PATH - location residuals will be written
            SYN - obspy Stream object containing synthetic data
            OBS - obspy Stream object containing observed data
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        residuals = []
        # TODO freq_mask = np.loadtxt('/data1/etienneb/freq_mask.txt')
        ft_obs_se = self.load('ft_obs_se')
        ft_stf_se = self.load('ft_stf_se')
        ft_stf_se_sinus = self.load('ft_stf_se_sinus')
        freq_mask = self.load('freq_mask_se')
        
        for ii in range(nn):
            residuals.append(self.misfit(syn[ii].data, nt, dt,ft_obs_se[:,ii],ft_stf_se,ft_stf_se_sinus,freq_mask[:,ii]))
        #residuals.append(self.misfit(syn, nt, dt,ft_obs_se,freq_mask))
        
        filename = path+'/'+'residuals'
        if exists(filename):
            residuals.extend(list(np.loadtxt(filename)))

        np.savetxt(filename, residuals)

    def write_adjoint_traces(self, path, syn, obs, channel,att=""):
        """ Writes "adjoint traces" required for gradient computation
         (overwrites synthetic data in the process)

          INPUT
            PATH - location "adjoint traces" will be written
            SYN - obspy Stream object containing synthetic data
            OBS - obspy Stream object containing observed data
            CHANNEL - channel or component code used by writer
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)
        
        adj = syn.copy()
        if att =='Yes' :
          self.adjoint = getattr(adjoint, PAR.MISFIT + '_att')
        else :
          self.adjoint = getattr(adjoint, PAR.MISFIT)
        #freq_mask = np.loadtxt('/data1/etienneb/freq_mask.txt')
        
        ft_obs_se = self.load('ft_obs_se')
        ft_stf_se = self.load('ft_stf_se')
        ft_stf_se_sinus = self.load('ft_stf_se_sinus')
        freq_mask = self.load('freq_mask_se')
        
        for ii in range(nn):
            adj[ii].data = self.adjoint(syn[ii].data, nt, dt,ft_obs_se[:,ii],ft_stf_se,ft_stf_se_sinus,freq_mask[:,ii])
        #self.adjoint(syn, adj, nt, dt,ft_obs_se,freq_mask)
        adj = self.apply_filter(adj,dt)
        
        #subset = np.random.choice([i for i in range(nn)],nn-nn/3)
        #for ii in range(nn-nn/3):
        #   adj[subset[ii]].data = np.zeros(len(adj[0].data)) 

        for tr in adj:
          tr.taper(0.005, type='hann')
        self.writer(adj, path, channel)

    def load(self, filename):
        # reads vectors from disk
        return loadnpy(PATH.ORTHO +'/'+filename)

    def save(self, filename, array):
        # writes vectors to disk
        savenpy(PATH.ORTHO +'/'+filename, array)