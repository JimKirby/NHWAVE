#-----------BEGIN MAKEFILE-------------------------------------------------
         SHELL         = /bin/sh
         DEF_FLAGS     = -P -C -traditional 
         EXEC          = mk_nesting
#==========================================================================
#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
         FLAG_1 = -DDOUBLE_PRECISION
#         FLAG_2 = -DPARALLEL
#         FLAG_3 = -DLANDSLIDE
#          FLAG_4 = -DSALINITY
#         FLAG_5 = -DTEMPERATURE
         FLAG_6 = -DCOUPLING
          FLAG_7 = -DSEDIMENT
#         FLAG_8 = -DVEGETATION
#         FLAG_9 = -DINTEL
#--------------------------------------------------------------------------
#  mpi defs 
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp 
         CPPFLAGS = $(DEF_FLAGS)
#         FC       = ifort
         FC        = mpif90
         DEBFLGS  = 
         OPT      = #-g
         CLIB     = 
#==========================================================================
         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib
#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
         CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) $(FLAG_3) \
                   $(FLAG_4) $(FLAG_5) $(FLAG_6)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------
         LIBS  = -L/home/1670/hypre-2.8.0b/src/hypre/lib/ -lHYPRE
         INCS  = -L/home/1670/hypre-2.8.0b/src/hypre/include/
#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPARGS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $*.f90
#	\rm $*.f90
#--------------------------------------------------------------------------
#  NHWAVE Source Code.
#--------------------------------------------------------------------------

MODS  = mod_global.F	mod_util.F 

#MAIN  = mk_input.F
MAIN   = mk_nesting.F

SRCS = $(MODS)  $(MAIN)

OBJS = $(SRCS:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS)
#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod

clobber:	clean
		/bin/rm -f *.f90 *.o mk_input







