
## change
WORKFLOW='inversion'# inversion, migration
SOLVER='specfem2d'      # specfem2d, specfem3d
## change
SYSTEM='multicore'         # serial, pbs, slurm
OPTIMIZE='LBFGS'        # NLCG, LBFGS
PREPROCESS='base'       # base
POSTPROCESS='base'      # base

MISFIT='Waveform'
MATERIALS='Elastic'
DENSITY='Constant'
ATTENUATION='yes'
ATTSIMU='yes'
INVPARA=['vs','Qmu']
FREQREF=0.084
COEF=10000


# WORKFLOW
BEGIN=1                 # first iteration
END=20                   # last iteration
NREC=132                # number of receivers
## change
NSRC=25                 # number of sources
SAVEGRADIENT=1          # save gradient how often
SAVERESIDUALS=1
SAVETRACES=1
SAVEKERNELS=1

# PREPROCESSING
FORMAT='su'             # data file format
CHANNELS=['z']          # data channels
NORMALIZE=0             # normalize
BANDPASS=0              # bandpass
MUTE=0                  # mute direct arrival
# FILTER = 'Bandpass'
# FREQMIN=0.05                # low frequency corner
# FREQMAX=0.2                # high frequency corner
MUTECONST=0.            # mute constant
MUTESLOPE=0.            # mute slope


# POSTPROCESSING
SMOOTH=20000.           # smoothing radius
RATIO=0.98              # reduce strength over iterations
SCALE=6.0e6             # scaling factor


# OPTIMIZATION
PRECOND=None            # preconditioner type
STEPMAX=10              # maximum trial steps
STEPTHRESH=0.1          # step length safeguard


# SOLVER
NT=4000                 # number of time steps
DT=0.06                 # time step
F0=0.084                # dominant frequency


# SYSTEM
NTASK=NSRC                # must satisfy 1 <= NTASK <= NSRC
NPROC=1                 # processors per task
NPROCMAX=25
WALLTIME=120
