double precision function bcoulomb(a,zg,t,zion)
  use constants_mod

  implicit none

  integer,intent(in) :: zg,zion
  double precision,intent(in) :: a,t
  double precision :: e2a,kt

  e2a = e_c*e_c/(1e-4*a)
  kt = k_b*t

  if (zg .eq. 0) then
     bcoulomb = 1. + sqrt(pi*zion**2*e2a/(2*kt))
  else if (zg*zion .gt. 0) then
     bcoulomb = exp(-zg*zion*e2a/kt)
  else if (zg*zion .lt. 0) then
     bcoulomb = (1 - zg*zion*e2a/kt)
  else if (zion .eq. 0) then
     bcoulomb = 1.
  end if

end function bcoulomb

double precision function coulombfac(a,nz,zg,fz,t,zion)
  use constants_mod

  implicit none

  integer,intent(in) :: nz,zg(nz),zion
  double precision,intent(in) :: a,fz(nz),t
  double precision :: bcoulomb
  integer :: i

  coulombfac = 0.

  do i=1,nz
     coulombfac = coulombfac + fz(i)*bcoulomb(a,zg(i),t,zion)
  end do

end function coulombfac
