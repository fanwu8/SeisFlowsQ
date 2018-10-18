
WORKFLOW='inversion'      # inversion, migration
SOLVER='specfem2d'        # specfem2d, specfem3d
SYSTEM='tiger_sm'         # serial, parallel, pbs, slurm
OPTIMIZE='LBFGS'          # base, newton
PREPROCESS='base'         # base
POSTPROCESS='base'        # base

MISFIT='Waveform'
MATERIALS='LegacyAcoustic'
DENSITY='Constant'


# WORKFLOW
BEGIN=1                   # first iteration
END=50                    # last iteration
NREC=32                   # number of receivers
NSRC=16                   # number of sources
SAVEGRADIENT=1            # save gradient how often


# PREPROCESSING
FORMAT='su'               # data file format
CHANNELS='y'              # data channels
NORMALIZE=0               # normalize
BANDPASS=0                # bandpass
MUTE=0                    # mute direct arrival
FREQLO=0.                 # low frequency corner
FREQHI=0.                 # high frequency corner
MUTECONST=0.              # mute constant
MUTESLOPE=0.              # mute slope


# POSTPROCESSING
SMOOTH=20.                # smoothing radius


# OPTIMIZATION
PRECOND=0                 # preconditioner flag
STEPMAX=10                # maximum trial steps
STEPTHRESH=0.1            # step length safeguard
ADHOCFACTOR=7.92e+00      # scaling factor



# SOLVER
NT=8000                   # number of time steps
DT=0.25                   # time step
F0=0.02                   # dominant frequency


# SYSTEM
NTASK=NSRC                # must satisfy 1 <= NTASK <= NSRC
NPROC=1                   # processers per source
WALLTIME=3000             # walltime

