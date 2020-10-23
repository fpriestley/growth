FC = /usr/bin/gfortran

SOURCE = main.f90 constants_mod.f90 sizedist.f90 integrate.f90 depletion.f90 grate.f90 writeoutput.f90 qfac.f90 countlines.f90 interpolate.f90 shatter.f90

main: $(SOURCE) constants_mod.o
	$(FC) -o growth $(SOURCE)

constants_mod.o: constants_mod.f90
	$(FC) -c constants_mod.f90

.PHONY: clean

clean:
	/bin/rm *.o *.mod growth
