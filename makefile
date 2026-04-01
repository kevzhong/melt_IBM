
#FC = h5pfc -r8 -O3 -fpp -xHost -funroll-loops -module $(OBJDIR) # FOR INTEL
FC = h5pfc -cpp -O3 -C -g -fdefault-real-8 -fdefault-double-8 -funroll-loops -Wno-line-truncation -fallow-argument-mismatch -J $(OBJDIR) # GNU

#FC+= -fcheck=all -fbacktrace # DEBUG FLAGS
################################## UNCONMMENT BELOW FOR DISCOVERER #############################################
## GNU
#LINKS = -lz -lhdf5_fortran -lhdf5
#LINKS+= -L/opt/software/lapack/3/3.12.0-gcc/lib -llapack -ldl -lblas
#LINKS+= -L/opt/software/fftw/3/3.3.10-gcc-openmpi/lib -lfftw3

## INTEL
#FC+=${FFTW3_FFLAGS}
#LINKS = -lz -lhdf5_fortran -lhdf5 -qmkl=sequential#-mkl=sequential

################################## UNCONMMENT BELOW FOR SNELLIUS #############################################
#LINKS = -lfftw3 -lz -lhdf5_fortran -lhdf5 -qmkl=sequential
#LINKS = -lfftw3 -lz -lhdf5_fortran -lhdf5

################################## UNCONMMENT BELOW FOR PERSONAL (MAC) MACHINE (GNU) #############################################
FFTW_PREFIX := $(shell brew --prefix fftw)
LINKS = -L$(FFTW_PREFIX)/lib -lfftw3 -lz -llapack -lblas -ldl -lhdf5_fortran -lhdf5 # for GNU personal machine, ABOVE FOR INTEL

PROGRAM = a.out

# Flag for computing spectra, do "make SPEC=1"
ifdef SPEC
	FC += -DSPEC
endif

FFILES  = auxroutines.f90 cfl.f90 cordin.f90 divg.f90 gcurv.f90  hdf.f90       \
          hdnl1.f90 hdnl2.f90 hdnl3.f90 hdnlte.f90 hit.f90 inirea.f90 init.f90 inqpr.f90  \
          interp.f90 invtr1.f90 invtr2.f90 invtr3.f90 invtrte.f90 matrix_transpose.f90     \
          mpi_routines.f90 mpiauxroutines.f90 papero.f90 phcalc.f90 phini.f90  \
          prcalc.f90 quit.f90 solxi.f90 solxj.f90 solxk.f90     \
          tridiag_periodic.f90 tsch.f90 updvp.f90 inicut.f90 movcut.f90 hdf2.f90    \
          diss.f90 vorticity.f90 spectra.f90 phasefield.f90 \
		  solxi_ib.f90 solxj_ib.f90 solxk_ib.f90

MFILES = param.f90 #KZ: no collisions
OBJDIR = obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)


all: objdir outdir $(PROGRAM)

$(PROGRAM) : $(MOBJS) $(OBJS)
	$(FC) -o $@ $^ $(LINKS) 


$(OBJDIR)/param.o: param.f90
	$(FC) -c -o $@ $<

$(OBJDIR)/geom.o: geom.f90
	$(FC) -c -o $@ $<

$(OBJDIR)/coll.o: coll.f90
	$(FC) -c -o $@ $<

$(OBJDIR)/rayAux.o: rayAux.f90
	$(FC) -c -o $@ $<



$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) -c -o $@ $<       


#-- Output folders 
OUTDIR := flowmov continuation stringdata spectra

clean :
	rm -rf $(OBJDIR)
	rm -rf $(PROGRAM)
	rm -rf *.mod

veryclean :
	rm -rf $(OBJDIR)
	rm -rf $(PROGRAM)
	rm -rf $(OUTDIR)
	rm -rf *.mod

.PHONY: objdir outdir
objdir: ${OBJDIR}
${OBJDIR}:
	mkdir -p ${OBJDIR}

outdir: ${OUTDIR}
${OUTDIR}:
	mkdir -p ${OUTDIR}