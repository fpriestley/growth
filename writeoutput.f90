subroutine openfiles(name)

  implicit none

  character(len=100) :: name
  

  open(unit=31,file='output/'//trim(name)//'.depletion.out',status='replace')
  open(unit=32,file='output/'//trim(name)//'.sizes.out',status='replace')

end subroutine openfiles

subroutine closefiles()

  implicit none

  close(unit=31)
  close(unit=32)

end subroutine closefiles

subroutine writedepletion(time,fgas,fdust,mdust,gtimescale)

  implicit none

  double precision,intent(in) :: time,fdust,fgas,mdust,gtimescale

  write(31,"(ES10.3,2(2X,ES27.20),2(2X,ES12.3))") time/1e6,fgas,fdust,mdust,gtimescale/1e6

end subroutine writedepletion

subroutine writesizes(nsizes,ngrain)

  implicit none

  integer,intent(in) :: nsizes
  double precision,intent(in) :: ngrain(nsizes)

  write(32,"(1000ES12.3)") ngrain

end subroutine writesizes
