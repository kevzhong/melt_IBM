#FC = h5pfc -O0 -r8 -fpp -module $(OBJDIR) -g -traceback -check bounds -debug all -warn all -fpe0 -ftrapuv

FC = h5pfc -r8 -O3 -fpp -module $(OBJDIR)# -traceback  -check bounds -fpe0
#uncomment the following 2 lines for discoverer
#FC+=${FFTW3_FFLAGS}
#LINKS = -lz -lhdf5_fortran -lhdf5 -qmkl=sequential#-mkl=sequential
#Uncomment the following lines for Snellius and personal machine
LINKS = -lfftw3 -lz -lhdf5_fortran -lhdf5 -qmkl=sequential
#export HDF5_FC=mpiifort #KZ: specify which Fortran compiler to use (for personal machine)


PROGRAM = a.out

FFILES  = auxroutines.f90 cfl.f90 cordin.f90 divg.f90 gcurv.f90  hdf.f90       \
          hdnl1.f90 hdnl2.f90 hdnl3.f90 hit.f90 inirea.f90 init.f90 inqpr.f90  \
          interp.f90 invtr1.f90 invtr2.f90 invtr3.f90 invtrte.f90 matrix_transpose.f90     \
          mpi_routines.f90 mpiauxroutines.f90 papero.f90 phcalc.f90 phini.f90  \
          prcalc.f90 quit.f90 solxi.f90 solxj.f90 solxk.f90 stat.f90           \
          tripvmyline.f90 tsch.f90 updvp.f90 inicut.f90 movcut.f90 hdf2.f90    \
          diss.f90 findCMindices.f90 vorticity.f90 injection.f90               \
		  tagging.f90

FFILES += allotri.f90 RigidAuxRoutines.f90 create_geo.f90 findindices.f90      \
         	forc1.f90 forc2.f90 forc3.f90 forctemp.f90 mlsForce.f90 mlsWeight.f90 partaux.f90 \
          initmls.f90 particle.f90 triaux.f90 tri_geo.f90 velforce.f90 cutvol.f90

MFILES = param.f90 geom.f90 coll.f90
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





$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) -c -o $@ $<       


#-- Output folders 
OUTDIR := flowmov continuation 



clean :
	rm -rf $(OBJDIR)
	rm -rf $(PROGRAM)

veryclean :
	rm -rf $(OBJDIR)
	rm -rf $(PROGRAM)
	rm -rf $(OUTDIR)

.PHONY: objdir outdir
objdir: ${OBJDIR}
${OBJDIR}:
	mkdir -p ${OBJDIR}

outdir: ${OUTDIR}
${OUTDIR}:
	mkdir -p ${OUTDIR}
