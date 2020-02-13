!--------------------------------------------------
! read and write data defined on rectangular grids
! to binary files
!--------------------------------------------------

SUBROUTINE write_grid(x,y,z, f, nx,ny,nz, filename)
  !
  ! Purpose
  ! -------
  ! writes grid data to a binary file
  !
  ! Parameters
  ! ----------
  ! x,y,z   : vectors with x,y, and z-position of grid points
  ! f       : grid data, f[i,j,k] = f(x[i],y[j],z[k])
  ! nx,ny,nz: number of grid points on each axis, nx,ny,nz >= 2
  ! filename: path to binary file
  !
  ! Returns
  ! -------
  ! nothing
  !
  ! See Also
  ! --------
  ! read_grid
  !
  
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE
  ! INPUT
  DOUBLE PRECISION, INTENT(IN) :: x(nx),y(ny),z(nz), f(nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  CHARACTER (len=*), INTENT(IN) :: filename
  ! LOCAL
  INTEGER, PARAMETER :: fh_data = 10            ! file handle
  INTEGER, PARAMETER :: fh_header = 11     ! file handle for header
  DOUBLE PRECISION :: x0,y0,z0  ! origin of grid
  DOUBLE PRECISION :: dx,dy,dz  ! grid spacing

  x0 = x(1)
  y0 = y(1)
  z0 = z(1)
  dx=x(2)-x(1)
  dy=y(2)-y(1)
  dz=z(2)-z(1)
  ! write header
  OPEN(unit=fh_header, form='formatted', &
       file=filename//'.header')
  WRITE(fh_header,*) nx,ny,nz, x0,y0,z0, dx,dy,dz
  CLOSE(fh_header)
  ! write data
  OPEN(unit=fh_data, form='unformatted', access='direct', &
       recl=sizeof(f), &
       file=filename)
  WRITE(fh_data,rec=1) f
  CLOSE(unit=fh_data)
END SUBROUTINE write_grid

SUBROUTINE read_grid(filename, x,y,z, f, nx,ny,nz)
  !
  ! Purpose
  ! -------
  ! writes grid data to a binary file
  !
  ! Parameters
  ! ----------
  ! filename: path to binary file
  ! x,y,z   : vectors with x,y, and z-position of grid points.
  !           Theses should be the same vectors used in ``write_grid``
  ! nx,ny,nz: number of grid points on each axis, nx,ny,nz >= 2
  !
  ! Returns
  ! -------
  ! f       : grid data, f[i,j,k] = f(x[i],y[j],z[k])
  !
  ! See Also
  ! --------
  ! write_grid
  !
  
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE
  ! INPUT
  CHARACTER (len=*), INTENT(IN) :: filename
  DOUBLE PRECISION, INTENT(IN) :: x(nx),y(ny),z(nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: f(nx,ny,nz)
  ! LOCAL
  INTEGER, PARAMETER :: fh_data = 10     ! file handle for data
  INTEGER, PARAMETER :: fh_header = 11 ! file handle for header
  DOUBLE PRECISION :: mx,my,mz  ! dummy variables with same meaning as nx,ny,nz
  DOUBLE PRECISION :: x0,y0,z0  ! origin of grid
  DOUBLE PRECISION :: dx,dy,dz  ! grid spacing
  ! read header
  OPEN(unit=fh_header, form='formatted', file=filename//'.header')
  READ(fh_header,*) mx,my,mz, x0,y0,z0, dx,dy,dz
  CLOSE(fh_header)
  ! read data
  OPEN(unit=fh_data, form='unformatted', access='direct', &
       recl=sizeof(f), &
       file=filename)
  READ(fh_data,rec=1) f
  CLOSE(unit=fh_data)
END SUBROUTINE read_grid
