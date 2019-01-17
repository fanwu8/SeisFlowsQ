
import sys

from seisflows.config import save

PAR = sys.modules['seisflows_parameters']

class base(object):
    """ Abstract base class
    """

    def check(self):
        """ Checks parameters and paths
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def submit(self):
        """ Submits workflow
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def run(self, classname, method, *args, **kwargs):
        """ Executes the following task:
              classname.method(*args, **kwargs)
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        raise NotImplementedError('Must be implemented by subclass.')



    def checkpoint(self):
        """ Writes information to disk so workflow can be resumed in case of
          interruption        
        """
        save()


    def progress(self, taskid):
        """ Provides status update
        """
        if PAR.NTASK > 1:
            print(' task ' + '%02d of %02d' % (taskid+1, PAR.NTASK))

