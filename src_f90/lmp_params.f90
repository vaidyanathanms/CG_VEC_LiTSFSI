! Input file to generate LAMMPS file for SIC systems
! Use in conjunction with lammps_inp.f90 and ran_numbers.f90
! Change the required parameters
! Monomers can contain 1 or 2 blobs
MODULE PARAMS

!  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file

  INTEGER, PARAMETER :: N_poly      = 100  ! No of polymer chains
  INTEGER, PARAMETER :: M_poly      = 30  ! MW of polymer chain
  INTEGER, PARAMETER :: CG_per_mon  = 2   ! No. of CG blobs per VEC monomer
  REAL,    PARAMETER :: an_poly_rat = 0.2 ! Ratio b/w VEC and STSFI monomers
  REAL,    PARAMETER :: charge_poly = 0.25 ! Charge on polymer blob
  REAL,    PARAMETER :: density = 0.8
  INTEGER, PARAMETER :: ideal_an_per_ch = INT(an_poly_rat*M_poly)
  INTEGER, PARAMETER :: N_anions    = N_poly*ideal_an_per_ch
  INTEGER, PARAMETER :: N_cations   = N_anions
  INTEGER, PARAMETER :: N_ions = N_anions + N_cations
  INTEGER, PARAMETER :: totpart = N_ions + N_poly*M_poly
  INTEGER, PARAMETER :: blob_per_ch = CG_per_mon*M_poly+ideal_an_per_ch
  INTEGER, PARAMETER :: totblobs = N_cations + N_poly*blob_per_ch
  
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

!!$  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS
