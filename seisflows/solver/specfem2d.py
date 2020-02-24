
import sys
from os.path import basename, join, dirname, relpath
from glob import glob

import numpy as np

from glob import glob
from seisflows.plugins.solver.specfem2d import smooth_legacy
from seisflows.tools.seismic import getpar, setpar

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.seismic import call_solver
from seisflows.tools.tools import exists
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class specfem2d(custom_import('solver', 'base')):
    """ Python interface for SPECFEM2D

      See base class for method descriptions
    """
    if PAR.MATERIALS == 'LegacyAcoustic':
        parameters = []
        parameters += ['vs']


    def check(self):
        """ Checks parameters and paths
        """
        super(specfem2d, self).check()

        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception

        # check data format
        if 'FORMAT' not in PAR:
            raise Exception()

        if PAR.FORMAT != 'su':
            raise Exception()


    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('NSTEP', cast=int)
        dt = getpar('DT', cast=float)
        #f0 = getpar('f0', file='DATA/SOURCE', cast=float)

        if nt != PAR.NT:
            if self.taskid == 0: print("WARNING: nt != PAR.NT")
            setpar('NSTEP', PAR.NT)

        if dt != PAR.DT:
            if self.taskid == 0: print("WARNING: dt != PAR.DT")
            setpar('deltat', PAR.DT)

        #if f0 != PAR.F0:
        #    if self.taskid == 0: print "WARNING: f0 != PAR.F0"
        #    setpar('f0', PAR.F0, filename='DATA/SOURCE')

        if self.mesh_properties.nproc != PAR.NPROC:
            if self.taskid == 0:
                print('Warning: mesh_properties.nproc != PAR.NPROC')

        if 'MULTIPLES' in PAR:
            if PAR.MULTIPLES:
                setpar('absorbtop', '.false.')
            else:
                setpar('absorbtop', '.true.')


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.cwd)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.false.')

        call_solver(system.mpiexec(), 'bin/xmeshfem2D',output='mesher.log')
        call_solver(system.mpiexec(), 'bin/xspecfem2D')

        if PAR.FORMAT in ['SU', 'su']:
            src = glob('OUTPUT_FILES/*.su')
            dst = 'traces/obs'
            unix.mv(src, dst)

        if PAR.SAVETRACES:
            self.export_traces(PATH.OUTPUT+'/'+'traces/obs')


    def initialize_adjoint_traces(self):
        super(specfem2d, self).initialize_adjoint_traces()

        # work around SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj*/*.su')
            unix.rename('.su', '.su.adj', files)

        # work around SPECFEM2D's requirement that all components exist,
        # even ones not in use
        if PAR.FORMAT in ['SU', 'su']:
            unix.cd(self.cwd +'/'+ 'traces/adj')
            for channel in ['x', 'y', 'z', 'p']:
                src = 'U%s_file_single.su.adj' % PAR.CHANNELS[0]
                dst = 'U%s_file_single.su.adj' % channel           ##4 adjoint components are the same.
                if not exists(dst):
                    unix.cp(src, dst)
            if PAR.ATTENUATION == 'yes':
                unix.cd(self.cwd + '/' + 'traces/adj_att')
                for channel in ['x', 'y', 'z', 'p']:
                    src = 'U%s_file_single.su.adj' % PAR.CHANNELS[0]
                    dst = 'U%s_file_single.su.adj' % channel  ##4 adjoint components are the same.
                    if not exists(dst):
                        unix.cp(src, dst)



    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.cwd)

        assert(exists(model_path))
        self.check_mesh_properties(model_path)

        src = glob(join(model_path, '*'))
        dst = join(self.cwd, 'DATA')
        unix.cp(src, dst)

        if self.taskid == 0:
            self.export_model(PATH.OUTPUT +'/'+ model_name)


    ### low-level solver interface

    def forward(self, path='traces/syn'):
        """ Calls SPECFEM2D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        call_solver(system.mpiexec(), 'bin/xmeshfem2D')

        call_solver(system.mpiexec(), 'bin/xspecfem2D')

        if PAR.FORMAT in ['SU', 'su']:
            filenames = glob('OUTPUT_FILES/*.su')
            unix.mv(filenames, path)


    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        # hack to deal with different SPECFEM2D name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj/*.su')
            unix.rename('.su', '.su.adj', files)

        call_solver(system.mpiexec(), 'bin/xspecfem2D')

    def adjoint_att(self):
        """ Calls SPECFEM2D adjoint solver

        """

        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj_att', 'SEM')

        # hack to deal with different SPECFEM2D name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj_att/*.su')
            unix.rename('.su', '.su.adj', files)

        call_solver(system.mpiexec(), 'bin/xspecfem2D')


    def rename_kernels(self):
        """ Works around conflicting kernel filename conventions
        """
        #files = []
        #files += glob('*proc??????_alpha_acoustic_kernel.bin')
        #unix.rename('alpha_acoustic', 'vp', files)
        files = []
        files += glob('*proc??????_alpha_kernel.bin')
        unix.rename('alpha', 'vp', files)
        files = []
        files += glob('*proc??????_c_acoustic_kernel.bin')
        unix.rename('c_acoustic', 'vp', files)
     #   unix.rename('c_acoustic', 'Qkappa', files)
        files = []
        files += glob('*proc??????_rhop_acoustic_kernel.bin')
        unix.rename('rhop_acoustic', 'rho', files)

        files = []
        files += glob('*proc??????_beta_kernel.bin')
        unix.rename('beta', 'vs', files)

    def export_att_kernel(self, path):
        unix.cd(self.kernel_databases)

        # work around conflicting name conventions
        files = []
        files += glob('*proc??????_kappa_kernel.bin')
        unix.rename('kappa', 'Qkappa', files)

        files = []
        files += glob('*proc??????_mu_kernel.bin')
        unix.rename('mu', 'Qmu', files)

        src = glob('*Q*_kernel.bin')
        dst = join(path, 'kernels', self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)


    def initialize_solver_directories(self):
        """ Creates directory structure expected by SPECFEM3D, copies 
          executables, and prepares input files. Executables must be supplied 
          by user as there is currently no mechanism for automatically
          compiling from source.
        """
        unix.mkdir(self.cwd)
        unix.cd(self.cwd)

        # create directory structure
        unix.mkdir('bin')
        unix.mkdir('DATA')
        unix.mkdir('OUTPUT_FILES')

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')
        if PAR.ATTENUATION == 'yes' :
          unix.mkdir('traces/adj_att')

        unix.mkdir(self.model_databases)
        unix.mkdir(self.kernel_databases)

        # copy exectuables
        src = glob(PATH.SPECFEM_BIN +'/'+ '*')
        dst = 'bin/'
        unix.cp(src, dst)

        # copy input files
        #src = glob(PATH.SPECFEM_DATA +'/'+ '*')
        #dst = 'DATA/'
        #unix.cp(src, dst)

        unix.cp(PATH.SPECFEM_DATA +'/Par_file','DATA/')
        unix.cp(PATH.SPECFEM_DATA +'/interfaces.dat','DATA/')
        unix.cp(PATH.SPECFEM_DATA +'/STATIONS','DATA/')

        src = PATH.SPECFEM_DATA + '/' + self.source_prefix +'_'+ self.source_name
        dst = 'DATA/' + self.source_prefix
        unix.cp(src, dst)

        if self.stf_file:
            src = self.stf_file
            rel = relpath(src, dirname(PATH.SPECFEM_DATA))
            if (rel[0:3] != '../'):
                dst = dirname(rel)
                unix.mkdir(dst)
                unix.cp(src, dst)
            

        #src = 'DATA/STATIONS'+'_'+ self.source_name
        #dst = 'DATA/STATIONS'
        #unix.cp(src, dst)

        self.check_solver_parameter_files()



    ### file transfer utilities

    def import_model(self, path):
        src = glob(path +'/'+ 'model/*')
        dst = join(self.cwd, 'DATA/')
        unix.cp(src, dst)

    def export_model(self, path):
        unix.mkdir(path)
        src = glob(join(self.cwd, 'DATA/*.bin'))
        dst = path
        unix.cp(src, dst)


    def check_stf_files(self):
        """ Get path of source time function files
        """
        stf_files = []
        for source_name in self.source_names_all:
            src = PATH.SPECFEM_DATA + '/' + self.source_prefix +'_'+ source_name
            stf_file = getpar('name_of_source_file', src).lstrip().rstrip()
            
            if stf_file[0:2] == './':
                # stf_file = join(PATH.SPECFEM_DATA,stf_file[2:])
                stf_file = dirname(PATH.SPECFEM_DATA) + '/' + stf_file[2:]

            if not exists(stf_file):
                stf_file = None
            


            stf_files.append(stf_file)

        self._stf_files = stf_files[:PAR.NTASK]
        self._stf_files_all = stf_files


    @property
    def data_filenames(self):
        if PAR.CHANNELS:
            if PAR.FORMAT in ['SU', 'su']:
               filenames = []
               for channel in PAR.CHANNELS:
                   filenames += ['U%s_file_single.su' % channel]
               return filenames

        else:
            unix.cd(self.cwd)
            unix.cd('traces/obs')

            if PAR.FORMAT in ['SU', 'su']:
                return glob('U?_file_single.su')

    @property
    def model_databases(self):
        return join(self.cwd, 'DATA')

    @property
    def kernel_databases(self):
        return join(self.cwd, 'OUTPUT_FILES')

    @property
    def source_prefix(self):
        return 'SOURCE'

    # workaround for older versions of SPECFEM2D,
    # which lacked a smoothing utility
    if not exists(PATH.SPECFEM_BIN+'/'+'xsmooth_sem') or \
        len(glob(PATH.MODEL_INIT+'/'+'proc*_NSPEC_ibool.bin')) == 0 or \
        len(glob(PATH.MODEL_INIT+'/'+'proc*_jacobian.bin')) == 0:
        smooth = staticmethod(smooth_legacy)
