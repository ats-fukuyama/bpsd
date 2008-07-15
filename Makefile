### $Id$ ###
include ../task/make.header

.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .a

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

SRCSFIXED = bpsd_flags.f bpsd_types.f bpsd_types_internal.f bpsd_subs.f \
            bpsd_shot.f bpsd_device.f bpsd_species.f \
            bpsd_equ1D.f bpsd_metric1D.f bpsd_plasmaf.f bpsd_base.f
SRCSFREE = bpsd_kinds.f90 bpsd_constants.f90

OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = bpsdlib.a ../task/lib/tasklib.a

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR)

all : bpsdlib.a

libs: ../task/lib/tasklib.a

../task/lib/tasklib.a:
	(cd ../task/lib; make tasklib.a)

bpsdlib.a: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a

veryclean: clean

new:
	-mkdir ../bpsdnew
	cp -f Makefile ../bpsdnew
	cp -f *.f ../bpsdnew
	cp -f *.f90 ../bpsdnew
