
import sys
import numpy as np

from os.path import join
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']


class base(object):
    """ Regularization, smoothing, sharpening, masking and related operations
      on models or gradients
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'RATIO' not in PAR:
            setattr(PAR, 'RATIO', 1.)

        if 'MASK' not in PATH:
            setattr(PATH, 'MASK', None)

        if PATH.MASK:
            assert exists(PATH.MASK)


    def setup(self):
        """ Placeholder for initialization or setup tasks
        """
        pass


    def write_gradient(self, path):
        """
        Combines contributions from individual sources and material parameters
        to get the gradient, and optionally applies user-supplied scaling

        :input path: directory from which kernels are read and to which
                     gradient is written
        """
        if not exists(path):
            raise Exception

        # because processing operations can be quite expensive, they must be
        # run through the HPC system interface; processing does not involve
        # embarassingly parallel tasks, we use system.run_single instead of 
        # system.run
        system.run_single('postprocess', 'process_kernels',
                 path=path+'/kernels',
                 parameters=solver.parameters)

        gradient = solver.load(
            path+'/'+'kernels/sum', suffix='_kernel')

        # merge into a single vector
        gradient = solver.merge(gradient)

        # convert to absolute perturbations, log dm --> dm
        # see Eq.13 Tromp et al 2005
        gradient *= solver.merge(solver.load(path +'/'+ 'model'))

        if PATH.MASK:
            # to scale the gradient, users can supply "masks" by exactly
            # mimicking the file format in which models stored
            mask = solver.merge(solver.load(PATH.MASK))

            # while both masking and preconditioning involve scaling the
            # gradient, they are fundamentally different operations:
            # masking is ad hoc, preconditioning is a change of variables;
            # see Modrak & Tromp 2016 GJI
            solver.save(solver.split(gradient),
                        path +'/'+ 'gradient_nomask',
                        parameters=solver.parameters,
                        suffix='_kernel')

            solver.save(solver.split(gradient*mask),
                        path +'/'+ 'gradient',
                        parameters=solver.parameters,
                        suffix='_kernel')

        else:
            solver.save(solver.split(gradient),
                        path +'/'+ 'gradient',
                        parameters=solver.parameters,
                        suffix='_kernel')



    def process_kernels(self, path='',iter=0, parameters=[]):
        """ Combines contributions from individual sources and performs any 
         required processing steps

          INPUT
              PATH - directory containing sensitivity kernels
              PARAMETERS - list of material parameters to be operated on
        """
        if not exists(path):
            raise Exception

        if not parameters:
            parameters = solver.parameters

        solver.combine(
               input_path=path,
               output_path=path+'/'+'sum',
               parameters=parameters)

        if PATH.MASK:
            # apply mask
            g = solver.merge(solver.load(path+'/'+'sum',suffix='_kernel'))
            g *= solver.merge(solver.load(PATH.MASK))
            solver.save(solver.split(g), path+'/'+'sum',suffix='_kernel')

        smo = PAR.SMOOTH * PAR.RATIO**(iter)       
        
        if PAR.SMOOTH > 0:
          src = path+'/'+'sum'
          dst = path+'/'+'sum_nosmooth' 
          unix.mv(src, dst)
          solver.smooth(
                   input_path=path+'/'+'sum_nosmooth',
                   output_path=path+'/'+'sum',
                   parameters=parameters,
                   span=smo)


