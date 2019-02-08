#-----------BEGIN MAKEFILE-------------------------------------------------
         SHELL         = /bin/sh
         DEF_FLAGS     = -P -C -traditional 
         EXEC          = nhwave
#==========================================================================
#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------
         FLAG_1 = -DDOUBLE_PRECISION
         FLAG_2 = -DPARALLEL
#         FLAG_3 = -DLANDSLIDE
#         FLAG_4 = -DSALINITY
#         FLAG_5 = -DTEMPERATURE
#         FLAG_6 = -DBUBBLE
#         FLAG_7 = -DSEDIMENT
#         FLAG_8 = -DVEGETATION
         FLAG_9 = -DINTEL
#	      FLAG_10 = -DBALANCE2D
#	      FLAG_11 = -DOBSTACLE
#	      FLAG_12 = -DTWOLAYERSLIDE
#	      FLAG_13 = -DCORALREEF
#	      FLAG_14 = -DPOROUSMEDIA
         FLAG_15 = -DFROUDE_CAP
#         FLAG_16 = -DCOUPLING
#         FLAG_17 = -DFLUIDSLIDE
#         FLAG_18 = -DLANDSLIDE_COMPREHENSIVE
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
                   $(FLAG_4) $(FLAG_5) $(FLAG_6) $(FLAG_7) $(FLAG_8) $(FLAG_9) \
	           $(FLAG_10) $(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
			   $(FLAG_15) $(FLAG_16)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------
         LIBS  = $(LDFLAGS) -lHYPRE
         INCS  = $(CPPFLAGS)
#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPARGS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
#	\rm $*.f90
#--------------------------------------------------------------------------
#  NHWAVE Source Code.
#--------------------------------------------------------------------------

MODS  = mod_global.F	mod_util.F 

MAIN  = nhwave.F initialize.F two_layer_slide.F source_terms.F waves.F

SRCS = $(MODS)  $(MAIN)

OBJS = $(SRCS:.F=.o) nspcg.o

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)
#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod

clobber:	clean
		/bin/rm -f *.f90 *.o nhwave







