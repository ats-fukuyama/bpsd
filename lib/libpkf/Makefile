### Makefile for wq 2021-11-04 ###

include ../../make.header

#LIB_MTX=$(LIB_MTX_BND)
#LIBX_MTX=$(LIBX_MTX_BND)
#LIB_MTX=$(LIB_MTX_PCG)
#LIBX_MTX=$(LIBX_MTX_PCG)
#LIB_MTX=$(LIB_MTX_MUMPS)
#LIBX_MTX=$(LIBX_MTX_MUMPS)
#LIB_MTX=$(LIBXMTX_KSP)
#LIBX_MTX=$(LIBX_MTX_KSP)

#FFLAGS=$(OFLAGS)
FFLAGS=$(DFLAGS)

OBJDIR   = ./obj

${OBJDIR}/%.o:%.f90
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)
${OBJDIR}/%.o:%.f
	$(FCFIXED) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

MODINCLUDE= -I./$(MOD) -I../$(MOD) -I../../mtxp/$(MOD) \
            -I../../../bpsd/$(MOD)

SRCS = libpgkf.f90 \
       libpdkf.f90 \
       libpdkf2.f90 \

OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS = libpkf.a ../../lib/libtask.a ../../lib/libgrf.a ../../../bpsd/libbpsd.a

all: libs testpkf

libpkf.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

../../bpsd/libbpsd.a:
	(cd ../../../bpsd; make libbpsd.a)
../task/lib/libtask.a:
	(cd ../../lib; make libtask.a)
../task/lib/libgrf.a:
	(cd ../../lib; make libgrf.a)
$(LIB_MTX):
	(cd ../../mtxp; make)

libs:
	(cd ../../../bpsd; make libbpsd.a)
	(cd ../../lib; make libtask.a)
	(cd ../../lib; make libgrf.a)
	(cd ../../mtxp; make)

testpkf: $(LIB_MTX) $(LIBS) $(OBJDIR)/testpkf.o
	$(FCFREE) $(OBJDIR)/testpkf.o $(LIBS) -o $@ $(FFLAGS) $(FLIBS) $(MODINCLUDE)

clean:
	-rm -f core a.out *.o *.a ./#* ./*~
	-rm -f ./mod/*.mod
	-rm -f ./obj/*.o
	-rm -rf ./*.dSYM

veryclean: clean
	-rm -f wq

$(OBJDIR)/testpkf.o: testpkf.f90 libpdkf.f90 libpgkf.f90 libpdkf2.f90
$(OBJDIR)/libpdkf.o: libpdkf.f90
$(OBJDIR)/libpdkf2.o: libpdkf2.f90
$(OBJDIR)/libpgkf.o: libpgkf.f90

