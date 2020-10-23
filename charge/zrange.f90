integer function getzmin(a)
  use constants_mod

  implicit none

  double precision,intent(in) :: a
  double precision :: aang,uv

  aang = 1e4*a
  uv = 2.5 + 0.07*aang + 8./aang

  getzmin = uv/14.4 * aang
  getzmin = max(int(-getzmin) + 1,-5)

end function getzmin

integer function getzmax(a,wf)
  use constants_mod

  implicit none

  double precision,intent(in) :: a,wf
  double precision :: hvmax,aang

  aang = 1e4*a
  hvmax = 13.6*ev

  getzmax = (hvmax - wf)/(14.4*ev) * aang + 0.5 - 0.3/aang
  getzmax = getzmax/(1 + 0.3/aang)
  getzmax = min(int(getzmax),5)

end function getzmax
