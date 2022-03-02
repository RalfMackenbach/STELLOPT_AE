
FC=gfortran
FFLAGS=  -g -fbacktrace -fbounds-check -fallow-argument-mismatch # -fcheck=all -Wall -Wextra -Wconversion -pedantic -O0
LDFLAGS=

OBJ=AE_proxy.o \
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
AE_integrand.o


%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

xaeproxy: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

run: xaeproxy
	@./xaeproxy

clean:
	@rm -rf *.o xaeproxy
