# Makfile

include ../task/make.header

.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

SRCSFIXED = 
SRCSFREE = bpsd_libfio.f90 \
           bpsd_kinds.f90 bpsd_constants.f90 bpsd_flags.f90 \
           bpsd_types.f90 bpsd_types_internal.f90 bpsd_subs.f90 \
           bpsd_shot.f90 bpsd_device.f90 bpsd_species.f90 \
           bpsd_equ1D.f90 bpsd_metric1D.f90 bpsd_plasmaf.f90 \
           bpsd_trmatrix.f90 bpsd_trsource.f90 \
           bpsd_base.f90

OBJS = $(SRCSFREE:.f90=.o) $(SRCSFIXED:.f=.o) 
OBJO = $(SRCO:.f=.o) 

LIBS = libbpsd.a
MODINCLUDE=-I $(MOD)

.f.o :
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
.f90.o :
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

all : libbpsd.a

libbpsd.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(MOD)/*.mod

veryclean: clean

BPSD_COMMON = bpsd_subs.f90 bpsd_types_internal.f90 bpsd_types.f90 \
	      bpsd_flags.f90 bpsd_constants.f90 bpsd_kinds.f90 
BPSD_TYPES  = bpsd_shot.f90 bpsd_device.f90 bpsd_species.f90 bpsd_equ1D.f90 \
	      bpsd_metric1D.f90 bpsd_plasmaf.f90

bpsd_libfio.o:		bpsd_libfio.f90
bpsd_kinds.o:		bpsd_kinds.f90
bpsd_constants.o: 	bpsd_constants.f90 bpsd_kinds.f90
bpsd_flags.o:		bpsd_flags.f90
bpsd_types.o:		bpsd_types.f90 bpsd_kinds.f90
bpsd_types_internal.o:	bpsd_types_internal.f90 bpsd_kinds.f90
bpsd_subs.o:		bpsd_subs.f90 bpsd_types_internal.f90 bpsd_kinds.f90
bpsd_shot.o:		bpsd_shot.f90 $(BPSD_COMMON)
bpsd_device.o:		bpsd_device.f90 $(BPSD_COMMON)
bpsd_species.o:		bpsd_species.f90 $(BPSD_COMMON)
bpsd_equ1D.o:		bpsd_equ1D.f90 $(BPSD_COMMON)
bpsd_metric1D.o:	bpsd_metric1D.f90 $(BPSD_COMMON)
bpsd_plasmaf.o:		bpsd_plasmaf.f90 $(BPSD_COMMON)
bpsd_trmatrix.o:	bpsd_trmatrix.f90 $(BPSD_COMMON)
bpsd_trsource.o:	bpsd_trsource.f90 $(BPSD_COMMON)
bpsd_base.o:		bpsd_base.f90 bpsd_subs.f90 $(BPSD_TYPES) $(BPSD_COMMON)
