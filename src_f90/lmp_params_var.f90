! Input file to generate LAMMPS file for SIC systems
! Use in conjunction with lammps_inp.f90 and ran_numbers.f90
! Change the required parameters
! Monomers can contain 1 or 2 blobs
MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file - inputs from python script

  INTEGER, PARAMETER :: M_poly      = py_mwchain     ! MW of polymerized
  ! chain (VEC + anion)
  INTEGER, PARAMETER :: N_poly      = py_nchains     ! # of
  ! polymerized chains
  REAL,    PARAMETER :: an_poly_rat = py_fracanions  ! Ratio b/w anion
  ! and VEC monomers: domain = [0,1) in polymerized chains 
  REAL,    PARAMETER :: frac_unpoly = py_unpoly_frac ! Unpolymerized
  ! VEC monomer fraction: domain = [0,1)

  REAL,    PARAMETER :: charge_poly = py_chargblob   ! Charge on polymer
  ! blob
  INTEGER, PARAMETER :: CG_per_mon  = py_cgpermon ! No. of CG blobs per VEC 
  REAL,    PARAMETER :: density     = py_density  ! System density


! Calculations based on input parameters for system initialization

  ! Ions
  INTEGER, PARAMETER :: ideal_an_per_ch = INT(an_poly_rat&
       &*REAL(M_poly)) ! number of anions per chain
  INTEGER, PARAMETER :: N_anions  = N_poly*ideal_an_per_ch ! Tot. #
  ! of anions in the system
  INTEGER, PARAMETER :: N_cations = N_anions ! Tot. # of cations
  INTEGER, PARAMETER :: N_ions    = N_anions + N_cations ! Tot. ions

  ! VECs (or uncharged moieties)
  INTEGER, PARAMETER :: VEC_per_ch   = INT(REAL(M_poly*REAL(1 -&
       & an_poly_rat)) ! Number of VEC per chain 
  INTEGER, PARAMETER :: T_poly_VEC   = VEC_per_ch*N_poly ! Tot. # of
  ! polymerized VECs in the system
  INTEGER, PARAMETER :: T_VEC_mons   = INT(REAL(T_poly_VEC)/REAL(1&
       &-frac_unpoly)) ! Tot. # of VEC mons
  INTEGER, PARAMETER :: T_unpoly_VEC = T_VEC_mons-T_poly_VEC ! Tot. #
  ! of unpolymerized mons

  ! Number of blobs
  INTEGER, PARAMETER :: blob_per_ch = CG_per_mon*VEC_per_ch +&
       & ideal_an_per_ch ! # of blobs per chain
  INTEGER, PARAMETER :: tot_poly_blobs = N_poly*blob_per_ch ! tot. #
  ! of blobs from polymerized system - VEC + anion
  INTEGER, PARAMETER :: unpoly_VEC_blobs = CG_per_mon*T_unpoly_VEC ! #
  ! of unpolymerized VEC blobs
  INTEGER, PARAMETER :: totblobs  = N_cations + tot_poly_blobs +&
       & unpoly_VEC_blobs ! Total number of blobs (LAMMPS atoms)
  INTEGER, PARAMETER :: nmols     = N_poly + T_unpoly_VEC ! Tot. # of
  ! LAMMPS molecules
  LOGICAL, PARAMETER :: unwrapped = .true. ! unwrapped coordinates
  LOGICAL, PARAMETER :: ifort     = .true. ! intel compiler

! Box/Particle details

  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox
  INTEGER :: npolyatoms

! Flags for creating the data file

  INTEGER, PARAMETER :: atomic = 0
  INTEGER, PARAMETER :: triblock = 0
  INTEGER, PARAMETER :: stretched = 0
  INTEGER, PARAMETER :: numatomtypes = CG_per_mon + 2
  INTEGER, PARAMETER :: numbondtypes = CG_per_mon + 2
  INTEGER, PARAMETER :: numangltypes = 0! 2*CG_per_mon
  INTEGER, PARAMETER :: numdihdtypes = 0
  INTEGER, PARAMETER :: bondtype = 1
  INTEGER, PARAMETER :: angltype = 1
  INTEGER, PARAMETER :: dihdtype = 0
  INTEGER, PARAMETER :: outfile  = 17
  INTEGER, PARAMETER :: nbonds = N_poly*(blob_per_ch-1)
  INTEGER, PARAMETER :: nangls = 0
  INTEGER, PARAMETER :: ndihds = 0

! Global Arrays involved in creating data file
  
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rxyz, uxyz
  REAL,ALLOCATABLE,DIMENSION(:) :: charge
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aidvals,ixyz
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: topo_bond
  
! Character Arrays for creating the data file name

  CHARACTER (LEN = 12) :: f_char
  CHARACTER (LEN = 60 ):: datafile

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS