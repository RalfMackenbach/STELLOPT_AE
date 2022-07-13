#GCC10 Macports
FC=gfortran
FFLAGS=  -g -fbacktrace -fbounds-check -fallow-argument-mismatch # -fcheck=all -Wall -Wextra -Wconversion -pedantic -O0
LDFLAGS=
AR=ar -ruvs

#Marconi
#FC=ifort
#FFLAGS = -g -traceback -O0 -xCORE-AVX512 -qopt-zmm-usage=high
#LDFLAGS=
#AR=ar -ruvs

OBJ=AE_proxy.o \
make_periodic.o \
bounce_wells.o \
inner_bounce_trapz.o \
left_bounce_trapz.o \
right_bounce_trapz.o \
trapz.o \
w_bounce.o \
zero_cross_idx.o \
zero_cross_count.o \
bounce_average.o \
AE_total.o \
AE_integrand.o \
trapped_avail_energy_mod.o 

.PHONY: clean run

# Default to makeing libtrap
libtrapAE.a: $(OBJ)
	$(AR) $@ $(OBJ)

AE_proxy.o : trapped_avail_energy_mod.o

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

xaeproxy: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

run: xaeproxy
	@./xaeproxy

clean:
	@rm -rf *.o xaeproxy libtrapAE.a
