subroutine readfactor(coulombname,nq,aq,fq)
  use constants_mod

  implicit none

  character(len=100),intent(in) :: coulombname
  integer,intent(in) :: nq
  double precision,intent(out) :: aq(nq),fq(nq)
  integer :: i

  open(unit=40,file='input/'//trim(coulombname),status='old')

  do i=1,nq
     read(40,*) aq(i),fq(i)
  end do

  close(unit=40)

end subroutine readfactor

subroutine calcfactor(nsizes,agrain,qgrain,nq,aq,fq)
  use constants_mod

  implicit none

  integer,intent(in) :: nsizes,nq
  double precision,intent(in) :: agrain(nsizes),aq(nq),fq(nq)
  double precision,intent(out) :: qgrain(nsizes)
  integer :: i
  double precision :: q1,q2

  do i=1,nsizes
     call interpolate(agrain(i),qgrain(i),aq,fq,nq)
  end do

end subroutine calcfactor
