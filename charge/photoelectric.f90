double precision function ionpot(a,z,w)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: a,w
  double precision :: aang,e2a

  aang = 1e4*a
  e2a = e_c*e_c/(1e-4*a)

  ionpot = w + (z+0.5)*e2a
  ionpot = ionpot + (z+2)*e2a * 0.3 / aang

end function ionpot

double precision function elecaff(a,z,w,e)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: a,w,e
  double precision :: e2a

  e2a = e_c*e_c/(1e-4*a)

  elecaff = w - e + (z-0.5)*e2a

end function elecaff

double precision function photoyield(hv,z,a,w,e,krad)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: hv,a,w,e,krad
  double precision :: e2a,aang
  double precision :: emin,hvpet,theta,le,la
  double precision :: alpha,beta
  double precision :: elow,ehigh
  double precision :: ionpot
  double precision :: y0,y1,y2

  aang = 1e4*a
  e2a = e_c*e_c/(1e-4*a)
  le = 10.
  la = hc/(hv*fourpi*krad)

  if (z .lt. 0) then
     emin = -(z+1) * e2a / (1 + (27/aang)**0.75)
  end if

  if (z .ge. -1) then
     hvpet = ionpot(a,z,w)
  else
     hvpet = ionpot(a,z,w) + emin
  end if
  
  if (z .ge. 0) then
     theta = hv - hvpet + (z+1)*e2a
  else
     theta = hv - hvpet
  end if

  y0 = 0.5*theta/(w + 5*theta)

  alpha = aang/la + aang/le
  beta = aang/la

  y1 = (beta/alpha)**2 * (alpha**2 - 2*alpha + 2 - 2*exp(-alpha))/ &
       & (beta**2 - 2*beta + 2 - 2*exp(-beta))

  if (z .lt. 0) then
     elow = emin
     ehigh = emin + hv - hvpet
  else
     elow = -(z+1)*e2a
     ehigh = hv - hvpet
  end if

  if (z .ge. 0) then
     y2 = ehigh**2 * (ehigh - 3*elow) / (ehigh - elow)**3
  else
     y2 = 1.
  end if

  photoyield = y2 * min(y0*y1,1.)

end function photoyield

double precision function lampet(a,z,w)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: a,w
  double precision :: aang,e2a,emin,ionpot,hvpet

  aang = 1e4*a
  e2a = e_c*e_c/(1e-4*a)
  
  if (z .lt. 0) then
     emin = -(z+1) * e2a / (1 + (27/aang)**0.75)
  end if

  if (z .ge. -1) then
     hvpet = ionpot(a,z,w)
  else
     hvpet = ionpot(a,z,w) + emin
  end if

  lampet = hc/hvpet*1e4

end function lampet

double precision function photdetach(hv,z,a,w,e)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: hv,a,w,e
  double precision :: elecaff
  double precision :: emin,hvpdt,dele,x
  double precision :: aang,e2a

  aang = 1e4*a
  e2a = e_c*e_c/(1e-4*a)
  dele = 3*ev

  emin = -(z+1) * e2a / (1 + (27/aang)**0.75)
  hvpdt = elecaff(a,z+1,w,e) + emin

  x = (hv - hvpdt)/dele

  photdetach = 1.2e-17 * abs(z) * x/(1+x**2/3)**2

  if (z .ge. 0) photdetach = 0.

end function photdetach

double precision function lampdt(a,z,w,e)
  use constants_mod

  implicit none

  integer,intent(in) :: z
  double precision,intent(in) :: a,w,e
  double precision :: aang,e2a,dele,emin,elecaff,hvpdt

  aang = 1e4*a
  e2a = e_c*e_c/(1e-4*a)
  dele = 3*ev

  emin = -(z+1) * e2a / (1 + (27/aang)**0.75)
  hvpdt = elecaff(a,z+1,w,e) + emin
  
  lampdt = hc/hvpdt*1e4

end function lampdt

double precision function JMat(lambda)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda
  double precision :: Jbb

  JMat = 0.0d0

  if (lambda .lt. 0.0912) then
     JMat = 0.0d0
  else if (lambda .lt. 0.110) then
     JMat = 3.0690d4*lambda**3.4172
  else if (lambda .lt. 0.134) then
     JMat = 1.627d1
  else if (lambda .lt. 0.250) then
     JMat = 0.566d0*lambda**(-1.6678)
  else if (lambda .ge. 0.250) then
     JMat = 1d-14*Jbb(lambda,7.5d3) + 1d-13*Jbb(lambda,4.0d3) &
          & + 4d-13*Jbb(lambda,3.0d3)
  end if

  Jmat = pi*Jmat

end function JMat

double precision function Jbb(lambda,T)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,T
  double precision :: lambdacm,invlambda,boltzfac

  Jbb = 0.

  lambdacm = lambda*1d-4 ! convert to cm
  
  invlambda = 1d0/lambdacm
  boltzfac = 1d0/(k_b*T)

  Jbb = 2d0*hc*c_s*invlambda**5 / (exp(hc*invlambda*boltzfac) - 1d0)

end function Jbb
