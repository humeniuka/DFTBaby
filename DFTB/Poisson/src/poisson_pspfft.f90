program poisson
  use PSPFFT

  implicit none
  include 'mpif.h'

  integer :: &
       Error, &
       nProcs, &
       nProcsRoot
  integer, dimension(3) :: &
       nCellsPerProcess, &
       nTotalCells
  real(KR), dimension(3) :: &
       CellWidth
  type(ArrayReal_3D_Base), dimension(1) :: &
       SourceSolution
  type(PSPFFT_Form), pointer :: &
       PS
  integer :: &
       nx,ny,nz
  integer, parameter :: &
       fh_data=10, fh_header=11
  real(KR) :: x0,y0,z0, dx,dy,dz
  call MPI_INIT(Error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, Error)
  WRITE(*,*) "nProcs = ", nProcs
  nProcsRoot = nProcs**(1.0_KR/3) + 0.5_KR
  ! read header
  OPEN(unit=fh_header, form='formatted', file='source.dat.header')
  READ(fh_header,*) nx,ny,nz, x0,y0,z0, dx,dy,dz
  CLOSE(fh_header)
  ! allocate space for source distribution
  CellWidth = (/ dx,dy,dz /)
  nCellsPerProcess = (/ nx,ny,nz /)
  nTotalCells = (/ nx,ny,nz /)
  allocate( SourceSolution(1)%Data(nx,ny,nz) )
  ! read source data from file
  OPEN(unit=fh_data, form='unformatted', access='direct', &
       recl=sizeof(SourceSolution(1)%Data), &
       file='source.dat')
  READ(fh_data,rec=1) SourceSolution(1)%Data
  CLOSE(unit=fh_data)
  ! solve Poisson equation
  call Create( PS , CellWidth , nTotalCells , MPI_COMM_WORLD)
  call Solve ( PS , SourceSolution )
  ! potential is returned in SourceSolution
  ! write it to a file called 'potential.dat'
  ! write header
  OPEN(unit=fh_header, form='formatted', &
       file='potential.dat.header')
  WRITE(fh_header,*) nx,ny,nz, x0,y0,z0, dx,dy,dz
  CLOSE(fh_header)
  ! write data
  OPEN(unit=fh_data, form='unformatted', access='direct', &
       recl=sizeof(SourceSolution(1)%Data), &
       file='potential.dat')
  WRITE(fh_data,rec=1) SourceSolution(1)%Data
  CLOSE(unit=fh_data)

  call Destroy ( PS )
  deallocate( SourceSolution(1)%Data )
  call MPI_FINALIZE(Error)
end program poisson

