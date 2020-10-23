subroutine calcmie(nWav,lambda,nrad,krad,a,Q)
  use constants_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,lambda(nWav),nrad(nWav),krad(nWav)
  double precision,intent(out) :: Q(nWav)
  double precision :: x,qext,qsca,gg
  double complex :: index
  integer :: i

  do i=1,nWav
     x = twopi*a/lambda(i)
     index = cmplx(nrad(i),krad(i))
     call bhmie(x,index,qext,qsca,gg)
     Q(i) = qext-qsca
     if (Q(i) .le. 0.0d0) Q(i) = 0.0d0
  end do

end subroutine calcmie

subroutine readnk(nkfile,nWav,lambda,nrad,krad)

  implicit none

  character(len=100),intent(in) :: nkfile
  integer,intent(in) :: nWav
  double precision,intent(out) :: lambda(nWav),nrad(nWav),krad(nWav)
  character(len=3) :: type
  integer :: i

  type = nkfile(len(trim(nkfile))-2:len(trim(nkfile)))

  open(unit=1,file=trim(nkfile))

  if (type .eq. 'lnk') then
     do i=1,nWav
        read(1,*) lambda(i),nrad(i),krad(i)
     end do
  else if (type .eq. 'vnk') then
     do i=nWav,1,-1
        read(1,*) lambda(i),nrad(i),krad(i)
        lambda(i) = 1.0d4/lambda(i)
     end do
  else
     do i=1,nWav
        lambda(i) = 0.0d0
        nrad(i) = 0.0d0
        krad(i) = 0.0d0
     end do
  end if

  close(unit=1)

end subroutine readnk
