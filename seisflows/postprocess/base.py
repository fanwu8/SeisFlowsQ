
import sys
import numpy as np

from glob import glob
from os.path import join
from seisflows.tools import unix, math
from seisflows.tools.tools import exists
from seisflows.config import ParameterError
from seisflows.plugins import solver_io

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']


class base(object):
    """ Gradient postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        # check parameters
        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'KERNELTYPE' not in PAR:
            setattr(PAR, 'KERNELTYPE', 'Relative')

        if 'PRECOND' not in PAR:
            setattr(PAR, 'PRECOND', False)

        if 'MUTESRC' not in PAR:
            setattr(PAR, 'MUTESRC', False)

        if 'MYSMOOTH' not in PAR:
            setattr(PAR, 'MYSMOOTH', False)

        if 'FILL' not in PAR:
            setattr(PAR, 'FILL', None)

        # check paths
        if 'MASK' not in PATH:
            setattr(PATH, 'MASK', None)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)

        if PATH.MASK:
            assert exists(PATH.MASK)


    def setup(self):
        """ Can be used to perform any required initialization or setup tasks
        """
        pass


    def write_gradient(self, path, iterf=0):
        """ Writes gradient of objective function

          Combines and processes contributions to the gradient from individual
          sources

          INPUT
              PATH - directory containing output of adjoint simulation
              iterf - number of iterations
        """
        if not exists(path):
            raise Exception

        system.run('postprocess', 'process_kernels',
                 hosts='head',
                 path=path+'/kernels',iter=iterf,
                 parameters=solver.parameters)

        if 'Qmu' in solver.parameters:
            Qmu_tmp = solver_io.fortran_binary._read(path + '/' + 'model/proc000000_Qmu.bin')
            grad_tmp = solver_io.fortran_binary._read(path + '/' + 'kernels/sum/proc000000_Qmu_kernel.bin')
            grad = grad_tmp * Qmu_tmp / PAR.COEF
            solver_io.fortran_binary._write(grad, path + '/' + 'kernels/sum/proc000000_Qmu_kernel.bin')


        gradient = solver.merge(solver.load(
                 path +'/'+ 'kernels/sum',
                 suffix='_kernel'))

        self.save(gradient, path)

        if PAR.KERNELTYPE=='Relative':
            # convert from relative to absolute perturbations

            gradient /= solver.merge(solver.load(path +'/'+ 'model'))
            self.save(gradient, path, backup='relative')


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
            gradient = solver.merge(solver.load(path+'/'+'sum',suffix='_kernel'))
            gradient *= solver.merge(solver.load(PATH.MASK))
            solver.save(solver.split(gradient), path+'/'+'sum',
                        suffix='_kernel',parameters=parameters)




        smo = PAR.SMOOTH * PAR.RATIO**(iter)       

        if PAR.SMOOTH > 0 and PAR.MYSMOOTH == False:
          src = path+'/'+'sum'
          dst = path+'/'+'sum_nosmooth'
          unix.mv(src, dst)
          solver.smooth(
                   input_path=path+'/'+'sum_nosmooth',
                   output_path=path+'/'+'sum',
                   parameters=parameters,
                   span=smo)

        if PAR.SMOOTH > 0 and PAR.MYSMOOTH == True:
            src = path + '/' + 'sum'
            dst = path + '/' + 'sum_nosmooth'
            unix.mv(src, dst)

            gradient = solver.load(path + '/' + 'sum_nosmooth', suffix='_kernel')
            files = []
            x_coords_file = PATH.MODEL_INIT + '/proc000000_x.bin'
            z_coords_file = PATH.MODEL_INIT + '/proc000000_z.bin'
            x = solver_io.fortran_binary._read(x_coords_file)
            z = solver_io.fortran_binary._read(z_coords_file)
            nx = PAR.NX
            nz = PAR.NZ
            smo = PAR.SMOOTH
            fill = PAR.FILL
            for para in parameters:
                gradient[para][0] = math.mysmooth(gradient[para][0],x,z,nx,nz,smo,fill)
            solver.save(gradient, path + '/' + 'sum', suffix='_kernel', parameters=parameters)


        if PATH.MASK:
            # apply mask
            gradient = solver.merge(solver.load(path+'/'+'sum',suffix='_kernel'))
            gradient *= solver.merge(solver.load(PATH.MASK))
            solver.save(solver.split(gradient), path+'/'+'sum',
                        suffix='_kernel',parameters=parameters)


        ### kernel for Q or inv of Q?
        # ker = solver.load(path+'/'+'sum',suffix='_kernel')
        # model = solver.load(path + '/../model/')
        # for key in solver.parameters:
        #     if key in ['Qkappa', 'Qmu']:
        #         for iproc in range(solver.mesh_properties.nproc):
        #             ker[key][iproc] /= -model[key][iproc]
        #             solver.io.write_slice(ker[key][iproc],path + '/' + 'sum',key+'_kernel',iproc)


            # for iproc in range(solver.mesh_properties.nproc):

        # for iproc in range(self.mesh_properties.nproc):
        #     for key in parameters:
        #         self.io.write_slice(
        #             dict[key][iproc], path, prefix+key+suffix, iproc)
        # solver.save(ker, path + '/' + 'sum', suffix='_kernel', fillin=False)

        # ker3 = solver.load(path+'/'+'sum',suffix='_kernel')

        print("finish")


    def save(self, gradient, path='', parameters=[], backup=None):
        """ Convience function for saving dictionary representation of 
          gradient
        """
        if not exists(path):
            raise Exception

        if not parameters:
            parameters = solver.parameters

        if backup:
            src = path +'/'+ 'gradient'
            dst = path +'/'+ 'gradient_'+backup
            unix.mv(src, dst)

        solver.save(solver.split(gradient),
                    path +'/'+ 'gradient',
                    parameters=parameters,
                    suffix='_kernel')


