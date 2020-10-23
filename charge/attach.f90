double precision function stick(a,z,zmin)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: a,zmin
  double precision :: le,nc,al

  le = 10.
  al = 1e4*a/le
  nc = 468 * (a/1e-3)**3

  if (z .eq. 0) then
     stick = 0.5 * (1. - exp(-al)) / (1. + exp(20 - nc))
  else if (z .gt. 0) then
     stick = 0.5 * (1. - exp(-al))
  else if (z .lt. 0) then
     if (z .gt. zmin) then
        stick = 0.5 * (1. - exp(-al)) / (1. + exp(20 - nc))
     else
        stick = 0.
     end if
  end if

end function stick

double precision function jbar(a,z,q,t)
  use constants_mod

  implicit none

  integer,intent(in) :: z,q
  double precision,intent(in) :: a,t
  integer :: vee
  double precision :: qq,tau,theta

  qq = q*e_c
  tau = 1e-4*a * k_b*t / qq**2
  vee = z/q

  if (vee .eq. 0) then
     jbar = 1. + sqrt(pi/(2*tau))
  else if (vee .lt. 0) then
     jbar = (1. - vee/tau) * (1. + sqrt(2./(tau - 2*vee)))
  else if (vee .gt. 0) then
     theta = vee/(1. + vee**(-0.5))
     jbar = (1. + 1./sqrt(4*tau + 3*vee))**2 * exp(-theta/vee)
  end if

end function jbar
