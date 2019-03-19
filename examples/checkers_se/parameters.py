
WORKFLOW='inversion_se' # inversion, migration
SOLVER='specfem2d'      # specfem2d, specfem3d
SYSTEM='serial'         # serial, pbs, slurm
OPTIMIZE='LBFGS'        # NLCG, LBFGS
PREPROCESS='ortho'      # base
POSTPROCESS='base'      # base

MISFIT='Phase2_se'
MATERIALS='LegacyAcoustic'
DENSITY='Constant'
ATTENUATION='no'


# WORKFLOW
BEGIN=1                 # first iteration
END=10                  # last iteration
NREC=132                # number of receivers
NSRC=1                  # number of sources
SAVEGRADIENT=1          # save gradient how often
SAVERESIDUALS=1

# SOURCE ENCODING
NEVT=25                 # number of encoded sources
FREQ_MIN=0.02           # minimium frequency of interest
FREQ_MAX=0.52               # maximium frequency of interest
RANDOM_OVER_IT=1        # reassign frequencies at each iteration
NFREQ_PER_EVENT=9       # number of frequencies per source
# NTPSS=7500
# NTPSS is the number of timesteps after steady state,
# which is determined by frequency range and length of timestep:
# NTPSS_DF = 1/((FREQ_MAX-FREQ_MIN)/NEVT/NFREQ_PER_EVENT)/DT
# if NTPSS is set here, an error will be thrown if NTPSS != NTPSS_DF


# PREPROCESSING
FORMAT='su'             # data file format
CHANNELS=['y']          # data channels
NORMALIZE=0             # normalize
BANDPASS=0              # bandpass
MUTE=0                  # mute direct arrival
FREQLO=0.               # low frequency corner
FREQHI=0.               # high frequency corner
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
NT=12300                # number of time steps
DT=0.06                 # time step
F0=0.084                # dominant frequency


# SYSTEM
NTASK=1                 # must satisfy 1 <= NTASK <= NSRC
NPROC=1                 # processors per task
