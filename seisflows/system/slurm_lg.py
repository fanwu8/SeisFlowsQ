
import os
import math
import sys
import time

from os.path import abspath, basename, join
from subprocess import check_output
from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.tools import call, findpath, saveobj, timestamp
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class slurm_lg(custom_import('system', 'base')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      Intermediate files are written to a global scratch path PATH.SCRATCH,
      which must be accessible to all compute nodes.

      Optionally, users can provide a local scratch path PATH.LOCAL if each
      compute node has its own local filesystem.

      For important additional information, please see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """

    def check(self):
        """ Checks parameters and paths
        """
        # name of job
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', basename(abspath('.')))

        # time allocated for workflow in minutes
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        # time allocated for each individual task in minutes
        if 'TASKTIME' not in PAR:
            setattr(PAR, 'TASKTIME', 15.)

        # number of tasks
        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        # number of cores per task
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # limit on number of concurrent tasks
        if 'NTASKMAX' not in PAR:
            setattr(PAR, 'NTASKMAX', PAR.NTASK)

        # number of cores per node
        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # how to invoke executables
        if 'MPIEXEC' not in PAR:
            setattr(PAR, 'MPIEXEC', 'srun')


        # optional additional SLURM arguments
        if 'SLURMARGS' not in PAR:
            setattr(PAR, 'SLURMARGS', '')

        # optional environment variable list VAR1=val1,VAR2=val2,...
        if 'ENVIRONS' not in PAR:
            setattr(PAR, 'ENVIRONS', '')

        # level of detail in output messages
        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # where job was submitted
        if 'WORKDIR' not in PATH:
            setattr(PATH, 'WORKDIR', abspath('.'))

        # where output files are written
        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', PATH.WORKDIR+'/'+'output')

        # where temporary files are written
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', PATH.WORKDIR+'/'+'scratch')

        # where system files are written
        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', PATH.SCRATCH+'/'+'system')

        # optional local scratch path
        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)


    def submit(self, workflow):
        """ Submits workflow
        """
        # create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # create output directories
        unix.mkdir(PATH.OUTPUT)
        unix.mkdir(PATH.WORKDIR+'/'+'output.slurm')

        self.checkpoint()

        # prepare sbatch arguments
        call('sbatch '
                + '%s ' % PAR.SLURMARGS
                + '--job-name=%s ' % PAR.TITLE
                + '--output %s ' % (PATH.WORKDIR+'/'+'output.log')
                + '--ntasks-per-node=%d ' % PAR.NODESIZE
                + '--nodes=%d ' % 1
                + '--time=%d ' % PAR.WALLTIME
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, method, hosts='all', **kwargs):
        """ Executes the following task:
              classname.method(*args, **kwargs)
        """
        self.checkpoint()

        self.save_kwargs(classname, method, kwargs)
        jobs = self.submit_job_array(classname, method, hosts)
        while True:
            # wait a few seconds before checking again
            time.sleep(5)

            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return


    def mpiexec(self):
        """ Specifies MPI executable used to invoke solver
        """
        return 'srun '


    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        try:
            return int(os.getenv('SLURM_ARRAY_TASK_ID'))
        except:
            raise Exception("TASK_ID environment variable not defined.")


    ### job array methods

    def submit_job_array(self, classname, method, hosts='all'):
        """ Submits job array and returns associated job ids
        """
        # submit job array
        cmd = self.job_array_cmd(classname, method, hosts)
        stdout = check_output(cmd, shell=True)

        # construct job id list
        id = stdout.split()[-1].strip().decode()
        if hosts=='all':
            tasks = range(PAR.NTASK)
            jobs = [id+'_'+str(task) for task in tasks]
        else:
            jobs = [id+'_0']
        return jobs


    def job_array_cmd(self, classname, method, hosts):
        return ('sbatch '
                + '%s ' % PAR.SLURMARGS
                + '--job-name=%s ' % PAR.TITLE
                + '--nodes=%d ' % math.ceil(PAR.NPROC/float(PAR.NODESIZE))
                + '--ntasks-per-node=%d ' % PAR.NODESIZE
                + '--ntasks=%d ' % PAR.NPROC
                + '--time=%d ' % PAR.TASKTIME
                + self.job_array_args(hosts)
                + findpath('seisflows.system') +'/'+ 'wrappers/run '
                + PATH.OUTPUT + ' '
                + classname + ' '
                + method + ' ' 
                + PAR.ENVIRONS)


    def job_array_args(self, hosts):
        if hosts == 'all':
            args = ('--array=%d-%d ' % (0,PAR.NTASK-1%PAR.NTASKMAX)
                   +'--output %s ' % (PATH.WORKDIR+'/'+'output.slurm/'+'%A_%a'))

        elif hosts == 'head':
            args = ('--array=%d-%d ' % (0,0)
                   +'--output=%s ' % (PATH.WORKDIR+'/'+'output.slurm/'+'%j'))

        else:
            raise KeyError('Bad keyword argument: system.run: hosts')

        return args


    def job_array_status(self, classname, method, jobs):
        """ Determines completion status of one or more jobs
        """
        states = []
        for job in jobs:
            state = self._query(job)
            if state in ['TIMEOUT']:
                print(msg.TimoutError % (classname, method, job, PAR.TASKTIME))
                sys.exit(-1)
            elif state in ['FAILED', 'NODE_FAIL']:
                print(msg.TaskError_SLURM % (classname, method, job))
                sys.exit(-1)
            elif state in ['COMPLETED']:
                states += [1]
            else:
                states += [0]

        isdone = all(states)

        return isdone, jobs


    def _query(self, job):
        """ Queries job state from SLURM database
        """
        stdout = check_output(
            'sacct -n -o jobid,state -j '+ job.split('_')[0],
            shell=True)

        state = ''
        lines = stdout.decode().strip().split('\n')
        for line in lines:
            if line.split()[0]==job:
                state = line.split()[1]
        return state


    def save_kwargs(self, classname, method, kwargs):
        kwargspath = join(PATH.OUTPUT, 'kwargs')
        kwargsfile = join(kwargspath, classname+'_'+method+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

