
WORKFLOW='inversion' # inversion, migration, modeling
SOLVER='specfem2d'   # specfem2d, specfem3d
SYSTEM='slurm_sm'      # serial, pbs, slurm
OPTIMIZE='LBFGS'     # base
PREPROCESS='base'  # base
POSTPROCESS='base'   # base

MISFIT='Waveform'
MATERIALS='Elastic'
DENSITY='Constant'


# WORKFLOW
BEGIN=1                 # first iteration
END=50                  # last iteration
NREC=500                # number of receivers
NSRC=32                 # number of sources
SAVEGRADIENT=1          # save gradient how often


# PREPROCESSING
FORMAT='su'             # data file format
CHANNELS='p'            # data channels


# OPTIMIZATION
PRECOND=None
STEPMAX=10              # maximum trial steps
STEPTHRESH=0.1          # step length safeguard


# POSTPROCESSING
SMOOTH=5.               # smoothing radius


# SOLVER
NT=7500                 # number of time steps
DT=9.0e-4               # time step
F0=4.0                  # dominant frequency


# SYSTEM
NTASK=NSRC              # number of tasks
NPROC=1                 # processors per task
WALLTIME=500            # walltime

