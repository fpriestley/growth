FC = /usr/bin/gfortran

SOURCE = main.f90 zrange.f90 photoelectric.f90 calcmie.f90 bhmie.f90 integrate.f90 interpolate.f90 countlines.f90 attach.f90 coulomb.f90

main: $(SOURCE) constants_mod.o
	$(FC) -o charge $(SOURCE)

constants_mod.o: constants_mod.f90
	$(FC) -c constants_mod.f90

.PHONY: clean

clean:
	/bin/rm *.o *.mod growth
