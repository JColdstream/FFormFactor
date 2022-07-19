! real kind precision parameters
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)

! pi
real, parameter :: pi = 3.14159265

logical :: first_read = .TRUE.
integer(8) :: trajstatus, frame

! wavevectors
integer(8) :: nq
real(dp) :: logqmin, logqmax

real(dp) :: pqnrm, rg
real(dp) :: scattering_sum = 0.0

real(dp), dimension(:), allocatable :: rx, ry, rz, p, q, b
real(dp), dimension(:,:), allocatable :: atom, k

real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax, lx, ly, lz

! index and types for atoms
integer(8), dimension(:), allocatable :: nindex, atomtype, typesanalyse, molanalyse

! trajectory numbers, timestep, number of total frames, number of frames to skip, number of frames to calculate, 
! number of atoms to analyse (nmol*molsize), number of atom types
integer(8) :: timestep, nframe, nskip, nframecalc, nanalyse, ntypes, lastframe

! centre of geometry variables
real(dp), dimension(3) :: cog

! system numbers
integer(8) :: natom, nmol, molsize, ntotal

! calculate centre of mass
integer(8) :: calc_com

! input, trajectory, and mass file names
character(LEN=256) :: inputfile, trajfile, scatteringlengths
character(LEN=256) :: strtypesanalyse

! storing strings in the LAMMPS dump file headers
character(LEN=256), dimension(4) :: headertext
