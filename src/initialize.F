!------------------------------------------------------------------------------------------
!
!
!    initialize.F
!
!-----------------------------------------------------------------------------------------
!
!
!    This file is part of NHWAVE.
!
!    This file contains the subroutines:
!
!        (1) wall_time_secs
!        (2) read_input
!        index
!        allocate_variables
!        generate_grid
!        read_bathymetry
!        initial
!
!------------------------------------------------------------------------------------------
!
!   BSD 2-Clause License
!
!   Copyright (c) 2019, NHWAVE Development Group
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
!
!   (1) wall_time_secs
!
!   Calculate current wall time
!
!   Called by: main
!
!   Gangfeng Ma, 09/12/2011
!-----------------------------------------------------------------------------------------
!
    subroutine wall_time_secs(tcurrent)
!
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
!
    end subroutine wall_time_secs
!
!---------------------------------------------------------------------------------------
!
     subroutine read_input
!
!----------------------------------------------------------------------------------------
!
!    (2) read_input
!
!    This subroutine is used to read input.txt
!
!    Called by: main
!
!    Gangfeng Ma, 20/12/2010
!
!----------------------------------------------------------------------------------------
!
     use global
     use input_util

# if defined (INTEL)
     use IFPORT, only: rand
# endif

     implicit none
     character(len=80) :: FILE_NAME
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,  &
                 Fk,Fkdif,Theta_Calc,Wnumy,tmp,tmp1,DFreq, &
                 Freq_Peak,gam,sa,sb,SumInt,A_Jon
     integer :: line,ierr,Iter,i,j,nw,n,nn
     !Added by M.Derakhti for FOCUSED wavemaker
     real(SP) :: Xltmp, Ctmp, xxk,segmasqr
     real(SP), dimension(12000) :: period_i 

     ! log and error file
     open(3,file='log.txt')

     ! read from input.txt
     FILE_NAME='input.txt'
!
!    title
!
     CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)

     IF(ierr==1)THEN

# if defined (PARALLEL)
     if(myid.eq.0) write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# else
     write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# endif

     TITLE='---TEST RUN---'
     ENDIF

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'---- LOG FILE ---'
     if(myid.eq.0) WRITE(3,*)TITLE
     if(myid.eq.0) WRITE(3,*)'--------------input start --------------'
# else
     WRITE(3,*)'---- LOG FILE ---'
     WRITE(3,*)TITLE
     WRITE(3,*)'--------------input start --------------'
# endif
!
!    dimension     
!                                        
     CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
     CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
     CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# else
     WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# endif

# if defined (PARALLEL)
     ! processor number
     CALL GET_INTEGER_VAL(PX,FILE_NAME,'PX',line)
     CALL GET_INTEGER_VAL(PY,FILE_NAME,'PY',line)
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PX= ',PX
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PY= ',PY
     if(PX*PY.ne.NumP) then
       if(myid.eq.0) WRITE(3,'(A6,I5)') 'NumP= ',NumP
       stop
     endif
# endif

     ! grid sizes
     CALL GET_Float_VAL(dx,FILE_NAME,'DX',line)
     CALL GET_Float_VAL(dy,FILE_NAME,'DY',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DX= ',dx
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DY= ',dy
# else
     WRITE(3,'(A4,F8.4)')'DX= ',dx
     WRITE(3,'(A4,F8.4)')'DY= ',dy
# endif

     ! vertical grid option
     call GET_INTEGER_VAL(Ivgrd,FILE_NAME,'IVGRD',line)
     CALL GET_Float_VAL(Grd_R,FILE_NAME,'GRD_R',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
     if(myid.eq.0) WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R
# else
     WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
     WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R
# endif

     ! time step
     CALL GET_Float_VAL(dt_ini,FILE_NAME,'DT_INI',line)
     CALL GET_Float_VAL(dt_min,FILE_NAME,'DT_MIN',line)
     CALL GET_Float_VAL(dt_max,FILE_NAME,'DT_MAX',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max
# else
     WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max
# endif

! result folder                                     
     CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# else
     WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# endif

     ! simulation steps and time
     call GET_INTEGER_VAL(SIM_STEPS,FILE_NAME,'SIM_STEPS',line)
     CALL GET_Float_VAL(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',line)
     CALL GET_Float_VAL(Plot_Start,FILE_NAME,'PLOT_START',line)
     CALL GET_Float_VAL(Plot_Intv,FILE_NAME,'PLOT_INTV',line)
     CALL GET_Float_VAL(Screen_Intv,FILE_NAME,'SCREEN_INTV',line)
	 CALL GET_LOGICAL_VAL(HOTSTART,FILE_NAME,'HOTSTART',line) ! added by Cheng for hot start
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     if(myid.eq.0) WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
	 if(myid.eq.0) WRITE(3,'(A11,L4)')'HOTSTART= ',HOTSTART ! added by Cheng for hot start
# else
     WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
	 WRITE(3,'(A11,L4)')'HOTSTART= ',HOTSTART ! added by Cheng for hot start
# endif

! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     CALL GET_STRING_VAL(Slide_File,FILE_NAME,'SLIDE_FILE',line,ierr)
     CALL GET_Float_VAL(SLIDE_DENSITY ,FILE_NAME,'SLIDE_DENSITY',line)
     CALL GET_Float_VAL(SLIDE_VISCOSITY ,FILE_NAME,'SLIDE_VISCOSITY',line)
     CALL GET_Float_VAL(SLIDE_MINTHICK ,FILE_NAME,'SLIDE_MINTHICK',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A20)')'Slide File= ', TRIM(Slide_File)
       if(myid.eq.0) WRITE(3,'(A15,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       if(myid.eq.0) WRITE(3,'(A17,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       if(myid.eq.0) WRITE(3,'(A16,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
# else
       WRITE(3,'(A12,A20)')'Slide File= ', TRIM(Slide_File)
       WRITE(3,'(A15,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       WRITE(3,'(A17,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       WRITE(3,'(A16,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
# endif
# endif

! added by Cheng for deformable slide
# if defined (DEFORMABLESLIDE)
     CALL GET_STRING_VAL(Slide_File,FILE_NAME,'SLIDE_FILE',line,ierr)
	 CALL GET_STRING_VAL(RHEO_OPT,FILE_NAME,'RHEO_OPT',line,ierr)
	 CALL GET_LOGICAL_VAL(NON_HYDRO_SLD,FILE_NAME,'NON_HYDRO_SLD',line)
	 CALL GET_LOGICAL_VAL(DISP_CORR_SLD,FILE_NAME,'DISP_CORR_SLD',line)
	 CALL GET_LOGICAL_VAL(REDU_GRAV_SLD,FILE_NAME,'REDU_GRAV_SLD',line)
	 CALL GET_LOGICAL_VAL(NON_HYDRO_UP,FILE_NAME,'NON_HYDRO_UP',line)
     CALL GET_Float_VAL(SLIDE_DENSITY ,FILE_NAME,'SLIDE_DENSITY',line)
     CALL GET_Float_VAL(SLIDE_VISCOSITY ,FILE_NAME,'SLIDE_VISCOSITY',line)
     CALL GET_Float_VAL(SLIDE_MINTHICK ,FILE_NAME,'SLIDE_MINTHICK',line)
     CALL GET_Float_VAL(SLIDE_GAMMA ,FILE_NAME,'SLIDE_GAMMA',line)
     CALL GET_Float_VAL(SLIDE_CONC ,FILE_NAME,'SLIDE_CONC',line)
     CALL GET_Float_VAL(GRAIN_DENSITY ,FILE_NAME,'GRAIN_DENSITY',line)
     CALL GET_Float_VAL(PhiInt_A ,FILE_NAME,'PhiInt_A',line)
     CALL GET_Float_VAL(PhiBed_A ,FILE_NAME,'PhiBed_A',line)
     CALL GET_Float_VAL(PhiInt_F ,FILE_NAME,'PhiInt_F',line)
     CALL GET_Float_VAL(PhiBed_F ,FILE_NAME,'PhiBed_F',line)
	 CALL GET_Float_VAL(SLIDE_LAMBDA ,FILE_NAME,'SLIDE_LAMBDA',line)
     CALL GET_Float_VAL(SLIDE_INIU ,FILE_NAME,'SLIDE_INIU',line)
	 CALL GET_Float_VAL(SLIDE_INIV ,FILE_NAME,'SLIDE_INIV',line)
	 CALL GET_Float_VAL(SLIDE_INIW ,FILE_NAME,'SLIDE_INIW',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A15)')'Slide File= ', TRIM(Slide_File)
	   if(myid.eq.0) WRITE(3,'(A10,A10)')'RHEO_OPT= ',TRIM(RHEO_OPT)
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'NON_HYDRO_SLD= ',NON_HYDRO_SLD
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'DISP_CORR_SLD= ',DISP_CORR_SLD
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'REDU_GRAV_SLD= ',REDU_GRAV_SLD
	   if(myid.eq.0) WRITE(3,'(A14,L4)')'NON_HYDRO_UP= ',NON_HYDRO_UP
       if(myid.eq.0) WRITE(3,'(A15,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       if(myid.eq.0) WRITE(3,'(A17,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       if(myid.eq.0) WRITE(3,'(A16,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
       if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SLIDE_GAMMA= ', SLIDE_GAMMA
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_CONC= ', SLIDE_CONC
       if(myid.eq.0) WRITE(3,'(A15,F8.2)')'GRAIN_DENSITY= ', GRAIN_DENSITY
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiInt_A= ', PhiInt_A
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiBed_A= ', PhiBed_A
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiInt_F= ', PhiInt_F
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiBed_F= ', PhiBed_F
	   if(myid.eq.0) WRITE(3,'(A14,F8.2)')'SLIDE_LAMBDA= ', SLIDE_LAMBDA
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIU= ', SLIDE_INIU
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIV= ', SLIDE_INIV
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIW= ', SLIDE_INIW
# else
       WRITE(3,'(A12,A15)')'Slide File= ', TRIM(Slide_File)
	   WRITE(3,'(A10,A10)')'RHEO_OPT= ',TRIM(RHEO_OPT)
	   WRITE(3,'(A15,L4)')'NON_HYDRO_SLD= ',NON_HYDRO_SLD
	   WRITE(3,'(A15,L4)')'DISP_CORR_SLD= ',DISP_CORR_SLD
	   WRITE(3,'(A15,L4)')'REDU_GRAV_SLD= ',REDU_GRAV_SLD
	   WRITE(3,'(A14,L4)')'NON_HYDRO_UP= ',NON_HYDRO_UP
       WRITE(3,'(A15,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       WRITE(3,'(A17,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       WRITE(3,'(A16,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
       WRITE(3,'(A13,F8.2)')'SLIDE_GAMMA= ', SLIDE_GAMMA
       WRITE(3,'(A12,F8.2)')'SLIDE_CONC= ', SLIDE_CONC
       WRITE(3,'(A15,F8.2)')'GRAIN_DENSITY= ', GRAIN_DENSITY
       WRITE(3,'(A8,F10.2)')'PhiInt_A= ', PhiInt_A
       WRITE(3,'(A8,F10.2)')'PhiBed_A= ', PhiBed_A
       WRITE(3,'(A8,F10.2)')'PhiInt_F= ', PhiInt_F
       WRITE(3,'(A8,F10.2)')'PhiBed_F= ', PhiBed_F
	   WRITE(3,'(A14,F8.2)')'SLIDE_LAMBDA= ', SLIDE_LAMBDA
       WRITE(3,'(A12,F8.2)')'SLIDE_INIU= ', SLIDE_INIU
       WRITE(3,'(A12,F8.2)')'SLIDE_INIV= ', SLIDE_INIV
       WRITE(3,'(A12,F8.2)')'SLIDE_INIW= ', SLIDE_INIW
# endif
     PhiInt_A = PhiInt_A*pi/180.
     PhiBed_A = PhiBed_A*pi/180.
     PhiInt_F = PhiInt_F*pi/180.
     PhiBed_F = PhiBed_F*pi/180.
# endif

     ! courant number
     CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# else
     WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# endif

     ! viscous number
     CALL GET_Float_VAL(VISCOUS_NUMBER,FILE_NAME,'VISCOUS_NUMBER',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# else
     WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# endif

! added by Cheng for limiting the maximum Froude number
# if defined (FROUDE_CAP)
     CALL GET_Float_VAL(FROUDECAP,FILE_NAME,'FROUDE_CAP',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A13,F8.3)')'FROUDE_CAP= ',FROUDECAP
# else
     WRITE(3,'(A13,F8.3)')'FROUDE_CAP= ',FROUDECAP
# endif

# endif
! end froude cap

     ! minimum depth
     CALL GET_Float_VAL(MinDep,FILE_NAME,'MinDep',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# else
     WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# endif

     ! laminar viscosity
     CALL GET_LOGICAL_VAL(VISCOUS_FLOW,FILE_NAME,'VISCOUS_FLOW',line)
     CALL GET_INTEGER_VAL(IVturb,FILE_NAME,'IVTURB',line)
     CALL GET_INTEGER_VAL(IHturb,FILE_NAME,'IHTURB',line)
     CALL GET_INTEGER_VAL(ProdType,FILE_NAME,'PRODTYPE',line)
     CALL GET_Float_VAL(Visc,FILE_NAME,'VISCOSITY',line)
     CALL GET_Float_VAL(Schmidt,FILE_NAME,'Schmidt',line)
     CALL GET_Float_VAL(Cvs,FILE_NAME,'Cvs',line)
     CALL GET_Float_VAL(Chs,FILE_NAME,'Chs',line)
     CALL GET_LOGICAL_VAL(RNG,FILE_NAME,'RNG',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IVTURB= ',IVturb
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IHTURB= ',IHturb
     if(myid.eq.0) WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     if(myid.eq.0) WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Chs= ',Chs
     if(myid.eq.0) WRITE(3,'(A5,L4)')'RNG= ',RNG
# else
     WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     WRITE(3,'(A7,I2)')'IVTURB= ',IVturb
     WRITE(3,'(A7,I2)')'IHTURB= ',IHturb
     WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     WRITE(3,'(A6,F8.3)')'Chs= ',Chs
     WRITE(3,'(A11,L4)')'RNG= ',RNG
# endif

! added by Cheng for hot start
     IF(HOTSTART)THEN
       CALL GET_STRING_VAL(Eta_HotStart_File,FILE_NAME,'Eta_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(U_HotStart_File,FILE_NAME,'U_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(V_HotStart_File,FILE_NAME,'V_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(W_HotStart_File,FILE_NAME,'W_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(P_HotStart_File,FILE_NAME,'P_HotStart_File',line,ierr)
# if defined (SALINITY)
       CALL GET_STRING_VAL(Sali_HotStart_File,FILE_NAME,'Sali_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(Temp_HotStart_File,FILE_NAME,'Temp_HotStart_File',line,ierr)
# endif
# if defined (LANDSLIDE) || defined (FLUIDSLIDE) || defined (TWOLAYERSLIDE) || defined (LANDSLIDE_COMPREHENSIVE) || defined (DEFORMABLESLIDE)
       CALL GET_STRING_VAL(Depth_HotStart_File,FILE_NAME,'Depth_HotStart_File',line,ierr)
# endif
# if defined (FLUIDSLIDE) || defined (DEFORMABLESLIDE)
       CALL GET_STRING_VAL(Us_HotStart_File,FILE_NAME,'Us_HotStart_File',line,ierr)
	   CALL GET_STRING_VAL(Vs_HotStart_File,FILE_NAME,'Vs_HotStart_File',line,ierr)
# endif
# if defined (DEFORMABLESLIDE)
       CALL GET_STRING_VAL(Ws_HotStart_File,FILE_NAME,'Ws_HotStart_File',line,ierr)
# endif
      IF(VISCOUS_FLOW)THEN
       ! CALL GET_STRING_VAL(Rho_HotStart_File,FILE_NAME,'Rho_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(TKE_HotStart_File,FILE_NAME,'TKE_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(EPS_HotStart_File,FILE_NAME,'EPS_HotStart_File',line,ierr)
      ENDIF
     ENDIF 
! end hotstart

     ! bathymetry     
     CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)
     CALL GET_LOGICAL_VAL(ANA_BATHY,FILE_NAME,'ANA_BATHY',line)
     CALL GET_Float_VAL(DepConst,FILE_NAME,'DepConst',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     if(myid.eq.0) WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY
     if(myid.eq.0) WRITE(3,'(A10,F8.5)')'DepConst= ',DepConst
# else
     WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY
     WRITE(3,'(A10,F8.5)')'DepConst= ',DepConst
# endif

     ! initial conditions
     CALL GET_LOGICAL_VAL(INITIAL_EUVW,FILE_NAME,'INITIAL_EUVW',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'INITIAL_EUVW= ',INITIAL_EUVW
# else
     WRITE(3,'(A14,L4)')'INITIAL_EUVW= ',INITIAL_EUVW
# endif

     ! bottom roughness
     CALL GET_INTEGER_VAL(Ibot,FILE_NAME,'Ibot',line)
     CALL GET_Float_VAL(Cd0,FILE_NAME,'Cd0',line)
     CALL GET_Float_VAL(Zob,FILE_NAME,'Zob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# else
     WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# endif

# if defined (CORALREEF)
     CALL GET_Float_VAL(Creef,FILE_NAME,'Creef',line)
# endif

     ! wind speed/stress
     CALL GET_INTEGER_VAL(Iws,FILE_NAME,'Iws',line)
     if(Iws==1) then
       CALL GET_Float_VAL(WindU,FILE_NAME,'WindU',line)
       CALL GET_Float_VAL(WindV,FILE_NAME,'WindV',line)
     endif

     ! Coriolis
     CALL GET_Float_VAL(slat,FILE_NAME,'slat',line)

     slat = slat*pi/180.0
     fcor = 2.0*7.29e-5*sin(slat)

     ! barotropic or baroclinic
     CALL GET_LOGICAL_VAL(BAROTROPIC,FILE_NAME,'BAROTROPIC',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC
# else
     WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC
# endif

     ! numerical scheme
     CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
     CALL GET_STRING_VAL(TIME_ORDER,FILE_NAME,'TIME_ORDER',line,ierr)
     CALL GET_STRING_VAL(CONVECTION,FILE_NAME,'CONVECTION',line,ierr)
     CALL GET_LOGICAL_VAL(ADV_HLLC,FILE_NAME,'HLLC',line)
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       HIGH_ORDER='SECOND'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# else
     WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# endif
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       TIME_ORDER='THIRD'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# else
     WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# endif

     ! ramp up the simulation
     CALL GET_Float_VAL(TRamp,FILE_NAME,'TRAMP',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp
# else
     WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp
# endif

     ! if non-hydrostatic simulation
     CALL GET_LOGICAL_VAL(NON_HYDRO,FILE_NAME,'NON_HYDRO',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# else
     WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# endif

     ! poisson solver
     CALL GET_INTEGER_VAL(isolver,FILE_NAME,'ISOLVER',line)
     CALL GET_INTEGER_VAL(itmax,FILE_NAME,'ITMAX',line)
     CALL GET_Float_VAL(tol,FILE_NAME,'TOL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     if(myid.eq.0) WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     if(myid.eq.0) WRITE(3,'(A5,E12.3)')'TOL= ',tol
# else
     WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     WRITE(3,'(A5,E12.3)')'TOL= ',tol
# endif

     ! periodic bc
     CALL GET_LOGICAL_VAL(PERIODIC_X,FILE_NAME,'PERIODIC_X',line)
     CALL GET_LOGICAL_VAL(PERIODIC_Y,FILE_NAME,'PERIODIC_Y',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y
# else
     WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y
# endif

     ! boundary type
     CALL GET_INTEGER_VAL(Bc_X0,FILE_NAME,'BC_X0',line)
     CALL GET_INTEGER_VAL(Bc_Xn,FILE_NAME,'BC_Xn',line)
     CALL GET_INTEGER_VAL(Bc_Y0,FILE_NAME,'BC_Y0',line)
     CALL GET_INTEGER_VAL(Bc_Yn,FILE_NAME,'BC_Yn',line)
     CALL GET_INTEGER_VAL(Bc_Z0,FILE_NAME,'BC_Z0',line)
     CALL GET_INTEGER_VAL(Bc_Zn,FILE_NAME,'BC_Zn',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# else
     WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# endif

! wavemaker 
     CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# else
     WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# endif
     IF(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='RIG'  &
         .or.WaveMaker(1:3)=='INT'.or.WaveMaker(1:3)=='FLU'.or.WaveMaker(1:3)=='WAV')THEN
       CALL GET_Float_VAL(Amp_Wave,FILE_NAME,'AMP',line)
       CALL GET_Float_VAL(Per_Wave,FILE_NAME,'PER',line)
       CALL GET_Float_VAL(Dep_Wave,FILE_NAME,'DEP',line)
       CALL GET_Float_VAL(Theta_Wave,FILE_NAME,'THETA',line)
       CALL GET_Float_VAL(Cur_Wave,FILE_NAME,'CUR',line)
       CALL GET_Float_VAL(sd_return,FILE_NAME,'sd_return',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       if(myid.eq.0) WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
       if(myid.eq.0) WRITE(3,'(A11,F6.3)')'sd_return= ',sd_return
# else
       WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
       WRITE(3,'(A11,F6.3)')'sd_return= ',sd_return
# endif

       IF(WaveMaker(1:3)=='INT') then
         CALL GET_Float_VAL(Xsource_West,FILE_NAME,'Xsource_West',line)
         CALL GET_Float_VAL(Xsource_East,FILE_NAME,'Xsource_East',line)
         CALL GET_Float_VAL(Ysource_Suth,FILE_NAME,'Ysource_Suth',line)
         CALL GET_Float_VAL(Ysource_Nrth,FILE_NAME,'Ysource_Nrth',line)
# if defined (PARALLEL)
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# else
         WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West                                                                                    
         WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East                                                                                    
         WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth                                                                                    
         WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# endif
       ENDIF

! test periodicity
       IF(PERIODIC_Y.and.Theta_Wave.ne.Zero) then
         ! find wave number
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length

         Iter = 0
 75      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85      continue

         if(Theta_Wave>Zero) then
           ! find right angle for periodic bc
           tmp = Large       
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         elseif(Theta_Wave<Zero) then
           ! find right angle for periodic bc 
           tmp = Large
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = -asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         endif

# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave         
           write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
         endif
# else
         write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave                                                      
         write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
# endif
         Theta_Wave = Theta_Calc
       ENDIF
     ENDIF

 
!added by M.Derakhti for FOCUSED wavemaker
     IF(WaveMaker(1:3)=='FOC') then
         CALL GET_INTEGER_VAL(nwave,FILE_NAME,'nwave',line)
         CALL GET_INTEGER_VAL(Component_Amp_Type,FILE_NAME,'Component_Amp_Type',line)
         CALL GET_Float_VAL(k_center,FILE_NAME,'k_center',line)
         CALL GET_Float_VAL(f_center,FILE_NAME,'f_center',line)
         CALL GET_Float_VAL(x_breaking,FILE_NAME,'x_breaking',line)
         CALL GET_Float_VAL(t_breaking,FILE_NAME,'t_breaking',line)
         CALL GET_Float_VAL(Slope_group,FILE_NAME,'Slope_group',line)
         CALL GET_Float_VAL(normalized_delta_f,FILE_NAME,'normalized_delta_f',line)
         CALL GET_Float_VAL(depth_comp,FILE_NAME,'depth_comp',line)

# if defined (PARALLEL)
         if(myid.eq.0) WRITE(3,'(A7,I2)')'nwave= ',nwave
         if(myid.eq.0) WRITE(3,'(A20,I2)')'Component_Amp_Type= ',Component_Amp_Type
         if(myid.eq.0) WRITE(3,'(A10,F6.3)')'k_center= ',k_center
         if(myid.eq.0) WRITE(3,'(A10,F6.3)')'f_center= ',f_center
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'x_breaking= ',x_breaking
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'t_breaking= ',t_breaking
         if(myid.eq.0) WRITE(3,'(A13,F6.3)')'Slope_group= ',Slope_group
         if(myid.eq.0) WRITE(3,'(A20,F6.3)')'normalized_delta_f= ',normalized_delta_f
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'depth_comp= ',depth_comp
# else
         WRITE(3,'(A7,I2)')'nwave= ',nwave
         WRITE(3,'(A20,I2)')'Component_Amp_Type= ',Component_Amp_Type
         WRITE(3,'(A10,F6.3)')'k_center= ',k_center
         WRITE(3,'(A10,F6.3)')'f_center= ',f_center
         WRITE(3,'(A12,F6.3)')'x_breaking= ',x_breaking
         WRITE(3,'(A12,F6.3)')'t_breaking= ',t_breaking
         WRITE(3,'(A13,F6.3)')'Slope_group= ',Slope_group
         WRITE(3,'(A20,F6.3)')'normalized_delta_f= ',normalized_delta_f
         WRITE(3,'(A12,F6.3)')'depth_comp= ',depth_comp
# endif

         !calculate f_i based on packet center frequency and band width
         period_i = 0.0
         do nw = 1 , nwave
            f_i (nw) = f_center - 0.5 * normalized_delta_f * f_center + &
              (nw - 1.0) * normalized_delta_f * f_center / (nwave - 1.0)
            period_i (nw)  = 1.0_SP / f_i (nw)
         end do
!calculate the wave number based on dispersion relation using
!Newton Ralphson method,xxk means k
         k_i = 0.0
         do nw = 1 , nwave
            segmasqr = ( 2.0d0 * pi / period_i (nw) ) ** 2
            ctmp = dsqrt ( Grav * depth_comp)
            xltmp = ctmp * period_i ( nw )
            xxk = 2.0d0 * pi / xltmp
            n = 0
 355        Fk = Grav * xxk * tanh ( xxk * depth_comp ) - segmasqr
            if ( abs ( Fk ) .le. 1.0e-6 .or. n .gt. 1000 ) goto 365
            Fkdif = Grav * xxk * depth_comp * ( 1.0d0 - &
            tanh ( xxk * depth_comp) ** 2 ) &
            + Grav * tanh( xxk * depth_comp )
            xxk = xxk - Fk / Fkdif
            n = n+1
            goto 355
 365        continue
            k_i (nw) = xxk
         end do
         !calculate componets amplitude
         Amp_i = 0.0
         If (Component_Amp_Type == 1) then
            do nw = 1 , nwave
               Amp_i (nw) = Slope_group / nwave / k_center
            end do
         else
            do nw = 1 , nwave
               Amp_i (nw) = Slope_group / nwave / k_i(nw)
            end do
         endif
     ENDIF

     ! random wave, read in 2d spectrum
     IF(WaveMaker(5:7)=='SPC') then
       open(14,file='spc2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
         stop
# endif
       endif
       if(NumDir>MaxNumDir) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
         stop
# endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
       close(14)

       ! random phase for each component
       do j = 1,NumFreq
       do i = 1,NumDir
# if defined (INTEL)
         Random_Phs(i,j) = rand()*2.0*pi
# else
         Random_Phs(i,j) = rand(0)*2.0*pi
# endif
       enddo
       enddo
     ENDIF
	 
     ! Added by Cheng Zhang,irregular wave internal wavemaker
     IF(WaveMaker(5:7)=='IRR') then
       open(14,file='irr2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
         stop
# endif
       endif
       if(NumDir>MaxNumDir) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
         stop
# endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
	   do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Random_Phs(i,j)
       enddo
       enddo
       close(14)
     ENDIF

! JONSWAP spectrum
     if((WaveMaker(5:7)=='JON').or.(WaveMaker(5:7)=='TMA')) then
       CALL GET_Float_VAL(Hm0,FILE_NAME,'Hm0',line)
       CALL GET_Float_VAL(Tp,FILE_NAME,'Tp',line)
       CALL GET_Float_VAL(Freq_Min,FILE_NAME,'Freq_Min',line)
       CALL GET_Float_VAL(Freq_Max,FILE_NAME,'Freq_Max',line)
       CALL GET_INTEGER_VAL(NumFreq,FILE_NAME,'NumFreq',line) 
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A5,f6.2)')'Hm0= ', Hm0
       if(myid.eq.0) WRITE(3,'(A4,f6.2)')'Tp= ', Tp
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Min= ', Freq_Min
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Max= ', Freq_Max
       if(myid.eq.0) WRITE(3,'(A9,I5)')'NumFreq= ', NumFreq
# endif

! jonswap spectrum
       gam = 3.3; sa = 0.07; sb = 0.09
       Freq_Peak = 1.0/Tp

       DFreq = (Freq_Max-Freq_Min)/NumFreq
       do i = 1,NumFreq
         Freq(i) = Freq_Min+0.5*DFreq+(i-1)*DFreq

         Per_Wave = 1.0/Freq(i)
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length
       
         Iter = 0
 76      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
                 Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 76
 86      continue

         if(Freq(i)<Freq_Peak) then
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sa**2)
         else
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sb**2)
         endif

         if(WaveMaker(5:7)=='TMA') then
            Jon_Spc(i) = Jon_Spc(i)*tanh(Wave_Number*Dep_Wave)**2/  &
                 (1.0+2*Wave_Number*Dep_Wave/sinh(2.*Wave_Number*Dep_Wave))
         endif
       enddo
         
       ! make sure m0=Hm0**2/16=int S(f)df
       SumInt = Zero
       do i = 1,NumFreq
         SumInt = SumInt+Jon_Spc(i)*DFreq
       enddo
       A_Jon = Hm0**2/16.0/SumInt

       do i = 1,NumFreq
         Jon_Spc(i) = Jon_Spc(i)*A_Jon
         RanPhs(i) = rand()*2.0*pi
       enddo
     endif

! sponge layer
     CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# else
     WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# endif
     IF(SPONGE_ON)THEN
       CALL GET_Float_VAL(Sponge_West_Width,FILE_NAME,'Sponge_West_Width',line)
       CALL GET_Float_VAL(Sponge_East_Width,FILE_NAME,'Sponge_East_Width',line)
       CALL GET_Float_VAL(Sponge_South_Width,FILE_NAME,'Sponge_South_Width',line)
       CALL GET_Float_VAL(Sponge_North_Width,FILE_NAME,'Sponge_North_Width',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
# else
       WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
# endif
     ENDIF

     ! wave average control
     CALL GET_LOGICAL_VAL(WAVE_AVERAGE_ON,FILE_NAME,'WAVE_AVERAGE_ON',line)
     CALL GET_Float_VAL(Wave_Ave_Start,FILE_NAME,'WAVE_AVERAGE_START',line)
     CALL GET_Float_VAL(Wave_Ave_End,FILE_NAME,'WAVE_AVERAGE_END',line)
     CALL GET_INTEGER_VAL(WaveheightID,FILE_NAME,'WaveheightID',line)

# if defined (SALINITY)
     CALL GET_LOGICAL_VAL(INITIAL_SALI,FILE_NAME,'INITIAL_SALI',line)
# endif

# if defined (SEDIMENT)
     ! sediment parameters
     CALL GET_STRING_VAL(Sed_Type,FILE_NAME,'Sed_Type',line,ierr)
     CALL GET_LOGICAL_VAL(BED_LOAD,FILE_NAME,'BED_LOAD',line)
     CALL GET_LOGICAL_VAL(COUPLE_FS,FILE_NAME,'COUPLE_FS',line) 
     CALL GET_INTEGER_VAL(ntyws,FILE_NAME,'ntyws',line)
     CALL GET_Float_VAL(Sedi_Ws,FILE_NAME,'Sedi_Ws',line)
     CALL GET_Float_VAL(SD50,FILE_NAME,'D50',line)
     CALL GET_Float_VAL(Shields_c,FILE_NAME,'Shields_c',line)
     CALL GET_Float_VAL(Af,FILE_NAME,'Af',line)
     CALL GET_Float_VAL(Tau_ce,FILE_NAME,'Tau_ce',line)
     CALL GET_Float_VAL(Tau_cd,FILE_NAME,'Tau_cd',line)
     CALL GET_Float_VAL(Erate,FILE_NAME,'Erate',line)
     CALL GET_Float_VAL(Mud_Visc,FILE_NAME,'Mud_Visc',line)
     CALL GET_Float_VAL(Tim_Sedi,FILE_NAME,'Tim_Sedi',line)
     CALL GET_Float_VAL(MorDt,FILE_NAME,'MorDt',line)
     CALL GET_LOGICAL_VAL(BED_CHANGE,FILE_NAME,'BED_CHANGE',line)
# endif

# if defined (VEGETATION)
     ! vegetation parameters
     CALL GET_STRING_VAL(Veg_Type,FILE_NAME,'Veg_Type',line,ierr)
     CALL GET_Float_VAL(Veg_X0,FILE_NAME,'Veg_X0',line)
     CALL GET_Float_VAL(Veg_Xn,FILE_NAME,'Veg_Xn',line)
     CALL GET_Float_VAL(Veg_Y0,FILE_NAME,'Veg_Y0',line)
     CALL GET_Float_VAL(Veg_Yn,FILE_NAME,'Veg_Yn',line)
     CALL GET_Float_VAL(VegH,FILE_NAME,'VegH',line)
     CALL GET_Float_VAL(VegDens,FILE_NAME,'VegDens',line)
     CALL GET_Float_VAL(StemD,FILE_NAME,'StemD',line)
     CALL GET_Float_VAL(VegDrag,FILE_NAME,'VegDrag',line)
     CALL GET_Float_VAL(Cfk,FILE_NAME,'Cfk',line)
     CALL GET_Float_VAL(Cfe,FILE_NAME,'Cfe',line)
     CALL GET_Float_VAL(VegVM,FILE_NAME,'VegVM',line)
     if(Veg_Type(1:3)=='FLE') then
       CALL GET_Float_VAL(EI,FILE_NAME,'EI',line)
     endif
# endif

# if defined (POROUSMEDIA)
     CALL GET_Float_VAL(alpha_por,FILE_NAME,'alpha_por',line)
     CALL GET_Float_VAL(beta_por,FILE_NAME,'beta_por',line)
     CALL GET_Float_VAL(D50_por,FILE_NAME,'D50_por',line)
     CALL GET_Float_VAL(Por_n,FILE_NAME,'Por_n',line)
     CALL GET_Float_VAL(Por_X0,FILE_NAME,'Por_X0',line)
     CALL GET_Float_VAL(Por_Xn,FILE_NAME,'Por_Xn',line)
     CALL GET_Float_VAL(Por_Y0,FILE_NAME,'Por_Y0',line)
     CALL GET_Float_VAL(Por_Yn,FILE_NAME,'Por_Yn',line)
     CALL GET_Float_VAL(Por_Z0,FILE_NAME,'Por_Z0',line)
     CALL GET_Float_VAL(Por_Zn,FILE_NAME,'Por_Zn',line)
# endif

# if defined (LANDSLIDE)
     ! landslide parameters                           
     CALL GET_STRING_VAL(SlideType,FILE_NAME,'SlideType',line,ierr)
     CALL GET_Float_VAL(SlideT,FILE_NAME,'SlideT',line)
     CALL GET_Float_VAL(SlideL,FILE_NAME,'SlideL',line)
     CALL GET_Float_VAL(SlideW,FILE_NAME,'SlideW',line)
     CALL GET_Float_VAL(SlideAngle,FILE_NAME,'SlideAngle',line)
     CALL GET_Float_VAL(SlopeAngle,FILE_NAME,'SlopeAngle',line)
     CALL GET_Float_VAL(SlideX0,FILE_NAME,'SlideX0',line)
     CALL GET_Float_VAL(SlideY0,FILE_NAME,'SlideY0',line)
     if(SlideType(1:5)=='RIGID') then ! modified by Cheng to identify 2d and 3d landslide
       CALL GET_Float_VAL(SlideUt,FILE_NAME,'SlideUt',line)
       CALL GET_Float_VAL(SlideA0,FILE_NAME,'SlideA0',line)
     endif
     if(trim(SlideType)=='DEFORMABLE') then
       CALL GET_Float_VAL(SlideDens,FILE_NAME,'SlideDens',line)
     endif
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
       CALL GET_LOGICAL_VAL(SLIDE_SHAPE_BOX,FILE_NAME,'SLIDE_SHAPE_BOX',line)

     IF(SLIDE_SHAPE_BOX)THEN  ! multi-slides for all box shape slides

      CALL GET_STRING_VAL(MultiSlidesFile,FILE_NAME,'MultiSlidesFile',line,ierr)

      OPEN(4,FILE=TRIM(MultiSlidesFile))

       READ(4,*)VoidString
       READ(4,*)NumSlides

       ALLOCATE(T_slide(NumSlides))
       ALLOCATE(b_slide(NumSlides))
       ALLOCATE(subaer_center2water0(NumSlides),subaer_center2water(NumSlides))
       ALLOCATE(w_slide(NumSlides))
       ALLOCATE(e_slide(NumSlides))
       ALLOCATE(alpha_slide(NumSlides))
       ALLOCATE(Coulomb_phi(NumSlides))
       ALLOCATE(x0_slide(NumSlides))
       ALLOCATE(y0_slide(NumSlides))
       ALLOCATE(slope_slide(NumSlides))
       ALLOCATE(term_v(NumSlides))
       ALLOCATE(uslide(NumSlides),ac_slide(NumSlides))
         uslide = 0.0
       ALLOCATE(Smax_slide(NumSlides))
       ALLOCATE(SlideStartTime(NumSlides))


       READ(4,*)VoidString
       READ(4,*)(SlideStartTime(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(T_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(b_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(w_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(e_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(subaer_center2water0(I),I=1,NumSlides)
         subaer_center2water=subaer_center2water0

       READ(4,*)VoidString
       READ(4,*)(alpha_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(Coulomb_phi(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(x0_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(y0_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(slope_slide(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(term_v(I),I=1,NumSlides)

       READ(4,*)VoidString
       READ(4,*)(Smax_slide(I),I=1,NumSlides)

      CLOSE(4)

       CALL GET_LOGICAL_VAL(SMFCut,FILE_NAME,'SMFCut',line) !CO
       CALL GET_LOGICAL_VAL(SlideStop,FILE_NAME,'SlideStop',line) !CO

     ELSE
       NumSlides = 1
       ALLOCATE(T_slide(NumSlides))
       ALLOCATE(b_slide(NumSlides))
       ALLOCATE(subaer_center2water0(NumSlides),subaer_center2water(NumSlides))
       ALLOCATE(w_slide(NumSlides))
       ALLOCATE(e_slide(NumSlides))
       ALLOCATE(alpha_slide(NumSlides))
       ALLOCATE(Coulomb_phi(NumSlides))
       ALLOCATE(x0_slide(NumSlides))
       ALLOCATE(y0_slide(NumSlides))
       ALLOCATE(slope_slide(NumSlides))
       ALLOCATE(term_v(NumSlides))
       ALLOCATE(uslide(NumSlides))
       ALLOCATE(Smax_slide(NumSlides))

       CALL GET_Float_VAL(T_slide(1),FILE_NAME,'H_slide',line)
       CALL GET_Float_VAL(b_slide(1),FILE_NAME,'L_slide',line)

       CALL GET_Float_VAL(subaer_center2water0(1),FILE_NAME,'subaer_center2water',line)
       subaer_center2water(1)=subaer_center2water0(1)

       CALL GET_Float_VAL(w_slide(1),FILE_NAME,'W_slide',line)
       CALL GET_Float_VAL(e_slide(1),FILE_NAME,'e_slide',line)
       CALL GET_Float_VAL(alpha_slide(1),FILE_NAME,'Angle_slide',line)
       CALL GET_Float_VAL(Coulomb_phi(1),FILE_NAME,'Coulomb_phi',line)
       CALL GET_Float_VAL(x0_slide(1),FILE_NAME,'X0_slide',line)  !is saved over later CO
       CALL GET_Float_VAL(y0_slide(1),FILE_NAME,'Y0_slide',line)  !is saved over later CO
       CALL GET_Float_VAL(XX0_Slide,FILE_NAME,'X0_slide',line) !CO
       CALL GET_Float_VAL(YY0_Slide,FILE_NAME,'Y0_slide',line) !CO
       CALL GET_Float_VAL(slope_slide(1),FILE_NAME,'Slope_slide',line)
       CALL GET_Float_VAL(term_v(1),FILE_NAME,'TermV_slide',line)
       CALL GET_Float_VAL(sf_slide,FILE_NAME,'sf_slide',line) !CO
       CALL GET_Float_VAL(tf_slide,FILE_NAME,'tf_slide',line) !CO
       CALL GET_LOGICAL_VAL(SMFCut,FILE_NAME,'SMFCut',line) !CO
       CALL GET_LOGICAL_VAL(Slide,FILE_NAME,'Slide',line) !CO
       CALL GET_LOGICAL_VAL(SlideStop,FILE_NAME,'SlideStop',line) !CO
       CALL GET_Float_VAL(Smax_slide(1),FILE_NAME,'Smax_slide',line) !CO
       CALL GET_Float_VAL(ac_slide(1),FILE_NAME,'ac_slide',line) !CO it will be calculate 


       IF(ANA_BATHY)THEN
       CALL GET_Float_VAL(acceleration_lab,FILE_NAME,'ACC_lab',line)
       ENDIF
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide height= ', T_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide length= ', b_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide width = ', w_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide e para= ', e_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide angle = ', alpha_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide init x= ', x0_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide init y= ', y0_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide slope = ', slope_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide term V= ', term_v
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide sf= ', sf_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide tf= ', tf_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'ac_slide= ', ac_slide
       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Smax_slide= ', Smax_slide
# else
       WRITE(3,'(A15,F12.3)')'Slide height= ', T_slide
       WRITE(3,'(A15,F12.3)')'Slide length= ', b_slide
       WRITE(3,'(A15,F12.3)')'Slide width = ', w_slide
       WRITE(3,'(A15,F12.3)')'Slide e para= ', e_slide
       WRITE(3,'(A15,F12.3)')'Slide angle = ', alpha_slide
       WRITE(3,'(A15,F12.3)')'Slide init x= ', x0_slide
       WRITE(3,'(A15,F12.3)')'Slide init y= ', y0_slide
       WRITE(3,'(A15,F12.3)')'Slide slope = ', slope_slide
       WRITE(3,'(A15,F12.3)')'Slide term V= ', term_v
       WRITE(3,'(A15,F12.3)')'Slide sf= ', sf_slide
       WRITE(3,'(A15,F12.3)')'Slide tf= ', tf_slide
       WRITE(3,'(A15,F12.3)')'ac_slide= ', ac_slide
       WRITE(3,'(A15,F12.3)')'Smax_slide= ', Smax_slide

# endif

     ENDIF ! end boxslide

# endif
! end landslide comprehensive

# if defined (TWOLAYERSLIDE)
     ! landslide parameters           
     CALL GET_Float_VAL(SlideT,FILE_NAME,'SlideT',line)
     CALL GET_Float_VAL(SlideL,FILE_NAME,'SlideL',line)
     CALL GET_Float_VAL(SlideW,FILE_NAME,'SlideW',line)
     CALL GET_Float_VAL(SlideAngle,FILE_NAME,'SlideAngle',line)
     CALL GET_Float_VAL(SlopeAngle,FILE_NAME,'SlopeAngle',line)
     CALL GET_Float_VAL(SlideX0,FILE_NAME,'SlideX0',line)
     CALL GET_Float_VAL(SlideY0,FILE_NAME,'SlideY0',line)
     CALL GET_Float_VAL(SlideDens,FILE_NAME,'SlideDens',line)
     CALL GET_Float_VAL(SlideVisc,FILE_NAME,'SlideVisc',line)
     CALL GET_Float_VAL(SlideLambda,FILE_NAME,'SlideLambda',line)
     CALL GET_Float_VAL(SlideIniU,FILE_NAME,'SlideIniU',line)
     CALL GET_Float_VAL(Cf_ul,FILE_NAME,'Cf_ul',line)
     CALL GET_Float_VAL(PhiInt,FILE_NAME,'PhiInt',line)
     CALL GET_Float_VAL(PhiBed,FILE_NAME,'PhiBed',line)

     PhiInt = PhiInt*pi/180.
     PhiBed = PhiBed*pi/180.
# endif
     
     ! whether to consider rheology
     CALL GET_LOGICAL_VAL(RHEOLOGY_ON,FILE_NAME,'RHEOLOGY_ON',line)
     CALL GET_Float_VAL(Yield_Stress,FILE_NAME,'Yield_Stress',line)
     CALL GET_Float_VAL(Plastic_Visc,FILE_NAME,'Plastic_Visc',line)

# if defined (POROUSMEDIA)
     CALL GET_Float_VAL(Carman_Koseny_Const,FILE_NAME,  &
                                 'Carman_Koseny_Const',line)
# endif

     ! if there is external forcing
     CALL GET_LOGICAL_VAL(EXTERNAL_FORCING,FILE_NAME,'EXTERNAL_FORCING',line)
     if(EXTERNAL_FORCING) then
       CALL GET_Float_VAL(Pgrad0,FILE_NAME,'Pgrad0',line)
     endif

     ! probe output
     CALL GET_INTEGER_VAL(NSTAT,FILE_NAME,'NSTAT',line)
     CALL GET_Float_VAL(Plot_Intv_Stat,FILE_NAME,'PLOT_INTV_STAT',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# else
       WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# endif

     if(NSTAT>0) then
       open(15,file='stat.txt',status='old')
       do i = 1,NSTAT
         read(15,*) xstat(i),ystat(i),zstat(i)
       enddo
       close(15)
     endif

     ! output parameters  
     CALL GET_LOGICAL_VAL(OUT_H,FILE_NAME,'OUT_H',line)
     CALL GET_LOGICAL_VAL(OUT_E,FILE_NAME,'OUT_E',line)
     CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'OUT_U',line)
     CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'OUT_V',line)
     CALL GET_LOGICAL_VAL(OUT_W,FILE_NAME,'OUT_W',line)
     CALL GET_LOGICAL_VAL(OUT_P,FILE_NAME,'OUT_P',line)
     CALL GET_LOGICAL_VAL(OUT_K,FILE_NAME,'OUT_K',line)
     CALL GET_LOGICAL_VAL(OUT_D,FILE_NAME,'OUT_D',line)
     CALL GET_LOGICAL_VAL(OUT_S,FILE_NAME,'OUT_S',line)
     CALL GET_LOGICAL_VAL(OUT_C,FILE_NAME,'OUT_C',line)
     CALL GET_LOGICAL_VAL(OUT_B,FILE_NAME,'OUT_B',line)
     CALL GET_LOGICAL_VAL(OUT_A,FILE_NAME,'OUT_A',line)
     CALL GET_LOGICAL_VAL(OUT_F,FILE_NAME,'OUT_F',line)
     CALL GET_LOGICAL_VAL(OUT_T,FILE_NAME,'OUT_T',line)
     CALL GET_LOGICAL_VAL(OUT_G,FILE_NAME,'OUT_G',line)
     CALL GET_LOGICAL_VAL(OUT_I,FILE_NAME,'OUT_I',line)
	 CALL GET_LOGICAL_VAL(OUT_Z,FILE_NAME,'OUT_Z',line) !added by Cheng for varying depth
	 CALL GET_LOGICAL_VAL(OUT_M,FILE_NAME,'OUT_M',line) !added by Cheng for recording Hmax
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_A= ',OUT_A
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_F= ',OUT_F
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_T= ',OUT_T
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_G= ',OUT_G
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_I= ',OUT_I
	 if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_Z= ',OUT_Z !added by Cheng for varying depth
	 if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_M= ',OUT_M !added by Cheng for recording Hmax
# else
     WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     WRITE(3,'(A7,L4)') 'OUT_A= ',OUT_A
     WRITE(3,'(A7,L4)') 'OUT_F= ',OUT_F
     WRITE(3,'(A7,L4)') 'OUT_T= ',OUT_T
     WRITE(3,'(A7,L4)') 'OUT_G= ',OUT_G
     WRITE(3,'(A7,L4)') 'OUT_I= ',OUT_I
	 WRITE(3,'(A7,L4)') 'OUT_Z= ',OUT_Z !added by Cheng for varying depth
	 WRITE(3,'(A7,L4)') 'OUT_M= ',OUT_M !added by Cheng for recording Hmax
# endif

! added by Cheng for nesting
# if defined (COUPLING)
      CALL GET_STRING_VAL(COUPLING_FILE,FILE_NAME,'COUPLING_FILE',line,ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'COUPLING_FILE:', COUPLING_FILE
# else
      WRITE(3,'(A15,A50)')'COUPLING_FILE:', COUPLING_FILE
# endif
# endif

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'--------------input end --------------'
# else
     WRITE(3,*)'--------------input end --------------'
# endif
!
     end subroutine read_input
!
!---------------------------------------------------------------------------------------
!
     subroutine index
!
!----------------------------------------------------------------------------------------
!
!    (3) index
!
!    This subroutine is used to create work index (what's that?)
!
!    Called by: main  
!                                                              
!    Last update: 20/12/2010, Gangfeng Ma       
!                              
!---------------------------------------------------------------------------------------
!
     use global
     implicit none

# if defined (PARALLEL)
     dims(1)=PX
     dims(2)=PY
     periods(1)=.false.
     periods(2)=.false.
     if(PERIODIC_X) periods(1)=.true.
     if(PERIODIC_Y) periods(2)=.true.
     coords(1)=0
     coords(2)=0

     call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims, &
         periods,reorder,comm2d,ier)
     call MPI_CART_COORDS(comm2d,myid,2,coords,ier)

     npx=coords(1)
     npy=coords(2)
 
     call MPI_CART_SHIFT(comm2d,0,1,n_west,n_east,ier)
     call MPI_CART_SHIFT(comm2d,1,1,n_suth,n_nrth,ier)
# else
     PX=1
     PY=1
# endif

     ! local index
     Mloc = Mglob/PX+2*Nghost
     Nloc = Nglob/PY+2*Nghost
     Kloc = Kglob+2*Nghost
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc1 = Kloc+1

     Ibeg = Nghost+1
     Iend = Mloc-Nghost
     Iend1 = Mloc1-Nghost
     Jbeg = Nghost+1
     Jend = Nloc-Nghost
     Jend1 = Nloc1-Nghost
     Kbeg = Nghost+1
     Kend = Kloc-Nghost
     Kend1 = Kloc1-Nghost

# if defined (DEFORMABLESLIDE)
# if defined (PARALLEL)
	 if (npx<(PX-1).and.npy<(PY-1)) then
	   IendC=Iend;JendC=Jend
	   MlocC=Mloc;NlocC=Nloc
	 elseif (npx==(PX-1).and.npy<(PY-1)) then
	   IendC=Iend1;JendC=Jend
       MlocC=Mloc1;NlocC=Nloc
	 elseif (npx<(PX-1).and.npy==(PY-1)) then
	   IendC=Iend;JendC=Jend1
       MlocC=Mloc;NlocC=Nloc1
	 elseif (npx==(PX-1).and.npy==(PY-1)) then
	   IendC=Iend1;JendC=Jend1
       MlocC=Mloc1;NlocC=Nloc1
	 endif
# else
	 IendC=Iend1;JendC=Jend1
     MlocC=Mloc1;NlocC=Nloc1
# endif
# endif
!
     end subroutine index
!
!-----------------------------------------------------------------------------------------
!
     subroutine allocate_variables
!
!----------------------------------------------------------------------------------------
!
!    (4) allocate_variables
! 
!    This subroutine is used to allocate variables
!
!    Called by: main  
!                                                                                
!    Last update: 23/12/2010, Gangfeng Ma    
!                                               
!----------------------------------------------------------------------------------------
!
     use global
     implicit none

     ! one-dimensional vars
     ALLOCATE(x(Mloc1),xc(Mloc),y(Nloc1),yc(Nloc),sig(Kloc1),dsig(Kloc),sigc(Kloc),  &
              Ein_X0(Nloc),Din_X0(Nloc),Ein_Xn(Nloc),Din_Xn(Nloc))

     ! two-dimensional vars
     ALLOCATE(HCG(Mglob,Nglob),Ho(Mloc,Nloc),H(Mloc1,Nloc1),Hc(Mloc,Nloc),Hc0(Mloc,Nloc), &
              Hfx(Mloc1,Nloc),Hfy(Mloc,Nloc1),Hfx0(Mloc1,Nloc),Hfy0(Mloc,Nloc1),D(Mloc,Nloc),Eta0(Mloc,Nloc),Eta00(Mloc,Nloc),  &
              D0(Mloc,Nloc),DeltH(Mloc,Nloc),DelxH(Mloc,Nloc),DelyH(Mloc,Nloc),Eta(Mloc,Nloc),Mask(Mloc,Nloc),  &
              Mask_Struct(Mloc,Nloc),Mask9(Mloc,Nloc),SourceC(Mloc,Nloc),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc),DeltHo(Mloc,Nloc),Uin_X0(Nloc,Kloc),Vin_X0(Nloc,Kloc),  &
              Win_X0(Nloc,Kloc),Uin_Xn(Nloc,Kloc),Uin_Xni(Nloc,Kloc),Uin_Xni0(Nloc,Kloc),Delt2H(Mloc,Nloc),  &
              Vin_Xn(Nloc,Kloc),Win_Xn(Nloc,Kloc),Bc_Prs(Mloc,Nloc),Brks(Mloc,Nloc))
     ALLOCATE(DxL(Mloc1,Nloc),DxR(Mloc1,Nloc),DyL(Mloc,Nloc1),DyR(Mloc,Nloc1), &
              EtaxL(Mloc1,Nloc),EtaxR(Mloc1,Nloc),EtayL(Mloc,Nloc1),EtayR(Mloc,Nloc1), &
              DelxEta(Mloc,Nloc),DelyEta(Mloc,Nloc),DelxD(Mloc,Nloc),DelyD(Mloc,Nloc),Sponge(Mloc,Nloc), &
              Setup(Mloc,Nloc),WaveHeight(Mloc,Nloc),Umean(Mloc,Nloc),Vmean(Mloc,Nloc),Num_Zero_Up(Mloc,Nloc), &
              Emax(Mloc,Nloc),Emin(Mloc,Nloc),WdU(Mloc,Nloc),WdV(Mloc,Nloc),Wsx(Mloc,Nloc),Wsy(Mloc,Nloc), &
			  HeightMax(Mloc,Nloc))

     ! three-dimensional vars

# if defined (LANDSLIDE_COMPREHENSIVE)
! fyshi add base bathymetry for landslide
     ALLOCATE(Dep0(Mloc1,Nloc1),DepC0(Mloc,Nloc))
# endif

     ALLOCATE(U(Mloc,Nloc,Kloc),V(Mloc,Nloc,Kloc),W(Mloc,Nloc,Kloc),Omega(Mloc,Nloc,Kloc1), &
              P(Mloc,Nloc,Kloc1),DU(Mloc,Nloc,Kloc),DV(Mloc,Nloc,Kloc),DW(Mloc,Nloc,Kloc),  &
              U0(Mloc,Nloc,Kloc),V0(Mloc,Nloc,Kloc),W0(Mloc,Nloc,Kloc),  &
              U00(Mloc,Nloc,Kloc),V00(Mloc,Nloc,Kloc),W00(Mloc,Nloc,Kloc),  &
              DU0(Mloc,Nloc,Kloc),DV0(Mloc,Nloc,Kloc),DW0(Mloc,Nloc,Kloc),Uf(Mloc,Nloc,Kloc1), &
              Vf(Mloc,Nloc,Kloc1),Wf(Mloc,Nloc,Kloc1),Cmu(Mloc,Nloc,Kloc),CmuR(Mloc,Nloc,Kloc), &
              Diffxx(Mloc,Nloc,Kloc),Diffxy(Mloc,Nloc,Kloc),Diffxz(Mloc,Nloc,Kloc), &
              Diffyx(Mloc,Nloc,Kloc),Diffyy(Mloc,Nloc,Kloc),Diffyz(Mloc,Nloc,Kloc),Diffzx(Mloc,Nloc,Kloc), &
              Diffzy(Mloc,Nloc,Kloc),Diffzz(Mloc,Nloc,Kloc),DelxSc(Mloc,Nloc,Kloc),DelySc(Mloc,Nloc,Kloc), &
              CmuHt(Mloc,Nloc,Kloc),CmuVt(Mloc,Nloc,Kloc),Rho(Mloc,Nloc,Kloc),Rmean(Mloc,Nloc,Kloc),Tke(Mloc,Nloc,Kloc), &
              Eps(Mloc,Nloc,Kloc),Skl(Mloc,Nloc,Kloc),DTke(Mloc,Nloc,Kloc),DEps(Mloc,Nloc,Kloc),DTke0(Mloc,Nloc,Kloc), &
              DEps0(Mloc,Nloc,Kloc),Prod_s(Mloc,Nloc,Kloc),Prod_b(Mloc,Nloc,Kloc),Lag_Umean(Mloc,Nloc,Kloc), &
              Lag_Vmean(Mloc,Nloc,Kloc),Lag_Wmean(Mloc,Nloc,Kloc),Euler_Umean(Mloc,Nloc,Kloc),Euler_Vmean(Mloc,Nloc,Kloc), &
              Euler_Wmean(Mloc,Nloc,Kloc),DRhoX(Mloc,Nloc,Kloc),DRhoY(Mloc,Nloc,Kloc),ExtForceX(Mloc,Nloc,Kloc), &
              ExtForceY(Mloc,Nloc,Kloc),UpWp(Mloc,Nloc,Kloc),IsMove(Mloc,Nloc,Kloc),Richf(Mloc,Nloc,Kloc), &
              DelxSl(Mloc,Nloc,Kloc1),DelySl(Mloc,Nloc,Kloc1))

     ! fluxes for construction at cell faces    
     ALLOCATE(UxL(Mloc1,Nloc,Kloc),UxR(Mloc1,Nloc,Kloc),VxL(Mloc1,Nloc,Kloc),VxR(Mloc1,Nloc,Kloc), &
              WxL(Mloc1,Nloc,Kloc),WxR(Mloc1,Nloc,Kloc),DUxL(Mloc1,Nloc,Kloc),DUxR(Mloc1,Nloc,Kloc), &
              DVxL(Mloc1,Nloc,Kloc),DVxR(Mloc1,Nloc,Kloc),DWxL(Mloc1,Nloc,Kloc),DWxR(Mloc1,Nloc,Kloc), &
              UyL(Mloc,Nloc1,Kloc),UyR(Mloc,Nloc1,Kloc),VyL(Mloc,Nloc1,Kloc),VyR(Mloc,Nloc1,Kloc), &
              WyL(Mloc,Nloc1,Kloc),WyR(Mloc,Nloc1,Kloc),DUyL(Mloc,Nloc1,Kloc),DUyR(Mloc,Nloc1,Kloc), &
              DVyL(Mloc,Nloc1,Kloc),DVyR(Mloc,Nloc1,Kloc),DWyL(Mloc,Nloc1,Kloc),DWyR(Mloc,Nloc1,Kloc), &
              UzL(Mloc,Nloc,Kloc1),UzR(Mloc,Nloc,Kloc1),VzL(Mloc,Nloc,Kloc1),VzR(Mloc,Nloc,Kloc1), &
              WzL(Mloc,Nloc,Kloc1),WzR(Mloc,Nloc,Kloc1),OzL(Mloc,Nloc,Kloc1),OzR(Mloc,Nloc,Kloc1), &
              SxL(Mloc1,Nloc,Kloc),SxR(Mloc1,Nloc,Kloc),SxS(Mloc1,Nloc,Kloc), &
              SyL(Mloc,Nloc1,Kloc),SyR(Mloc,Nloc1,Kloc),SyS(Mloc,Nloc1,Kloc), &
              ExL(Mloc1,Nloc,Kloc),ExR(Mloc1,Nloc,Kloc),FxL(Mloc1,Nloc,Kloc),FxR(Mloc1,Nloc,Kloc), &
              GxL(Mloc1,Nloc,Kloc),GxR(Mloc1,Nloc,Kloc),HxL(Mloc1,Nloc,Kloc),HxR(Mloc1,Nloc,Kloc), &
              EyL(Mloc,Nloc1,Kloc),EyR(Mloc,Nloc1,Kloc),FyL(Mloc,Nloc1,Kloc),FyR(Mloc,Nloc1,Kloc), &
              GyL(Mloc,Nloc1,Kloc),GyR(Mloc,Nloc1,Kloc),HyL(Mloc,Nloc1,Kloc),HyR(Mloc,Nloc1,Kloc), &
              Ex(Mloc1,Nloc,Kloc),Ey(Mloc,Nloc1,Kloc),Fx(Mloc1,Nloc,Kloc),Fy(Mloc,Nloc1,Kloc), &
              Gx(Mloc1,Nloc,Kloc),Gy(Mloc,Nloc1,Kloc),Hx(Mloc1,Nloc,Kloc),Hy(Mloc,Nloc1,Kloc), &
              Fz(Mloc,Nloc,Kloc1),Gz(Mloc,Nloc,Kloc1),Hz(Mloc,Nloc,Kloc1),DelxU(Mloc,Nloc,Kloc), &
              DelyU(Mloc,Nloc,Kloc),DelzU(Mloc,Nloc,Kloc),DelxV(Mloc,Nloc,Kloc),DelyV(Mloc,Nloc,Kloc), &
              DelzV(Mloc,Nloc,Kloc),DelxW(Mloc,Nloc,Kloc),DelyW(Mloc,Nloc,Kloc),DelzW(Mloc,Nloc,Kloc), &
              DelxDU(Mloc,Nloc,Kloc),DelyDU(Mloc,Nloc,Kloc),DelxDV(Mloc,Nloc,Kloc),DelyDV(Mloc,Nloc,Kloc), &
              DelxDW(Mloc,Nloc,Kloc),DelyDW(Mloc,Nloc,Kloc),DelzO(Mloc,Nloc,Kloc)) 

# if defined (BUBBLE)
     ALLOCATE(Rbg(Mg),DRbg(Mg),Wbg(Mg),Entrain(Mg))
     ALLOCATE(Nbg(Mloc,Nloc,Kloc,Mg),DNbg(Mloc,Nloc,Kloc,Mg),DNbg0(Mloc,Nloc,Kloc,Mg),Vbg(Mloc,Nloc,Kloc))
# endif

# if defined (SEDIMENT)
     ALLOCATE(Conc(Mloc,Nloc,Kloc),DConc(Mloc,Nloc,Kloc),DConc0(Mloc,Nloc,Kloc),  &
              Euler_Smean(Mloc,Nloc,Kloc),Smean(Mloc,Nloc),Qbedx(Mloc,Nloc),Qbedy(Mloc,Nloc),  &
              SSour(Mloc,Nloc),Bed(Mloc,Nloc),Bed0(Mloc,Nloc),Taub(Mloc,Nloc),Sed_X0(Nloc,Kloc),  &
              Sed_Xn(Nloc,Kloc),QbxAvg(Mloc,Nloc),QbyAvg(Mloc,Nloc),TotDep(Mloc,Nloc))
# endif

# if defined (SALINITY)
     ALLOCATE(Sali(Mloc,Nloc,Kloc),DSali(Mloc,Nloc,Kloc),DSali0(Mloc,Nloc,Kloc))
     ALLOCATE(Sin_X0(Nloc,Kloc),Sin_Xn(Nloc,Kloc))
# endif

# if defined (TEMPERATURE)
     ALLOCATE(Temp(Mloc,Nloc,Kloc),DTemp(Mloc,Nloc,Kloc),DTemp0(Mloc,Nloc,Kloc))                                    
     ALLOCATE(Tin_X0(Nloc,Kloc),Tin_Xn(Nloc,Kloc))
# endif

# if defined (VEGETATION)
     ALLOCATE(Prod_v(Mloc,Nloc,Kloc),FVegH(Mloc,Nloc),foliage(Mloc,Nloc),Prod_w(Mloc,Nloc,Kloc), &
              Tke_w(Mloc,Nloc,Kloc),Eps_w(Mloc,Nloc,Kloc),DWke(Mloc,Nloc,Kloc),DWke0(Mloc,Nloc,Kloc), &
              Wke_p(Mloc,Nloc,Kloc))
# endif

# if defined (BALANCE2D)
     ALLOCATE(DUDT2D(Mloc,Nloc),DVDT2D(Mloc,Nloc),DEDX2D(Mloc,Nloc),DEDY2D(Mloc,Nloc),  &
              DPDX2D(Mloc,Nloc),DPDY2D(Mloc,Nloc),DIFFX2D(Mloc,Nloc),DIFFY2D(Mloc,Nloc),  &
              TAUBX2D(Mloc,Nloc),TAUBY2D(Mloc,NLoc),FVEGX2D(Mloc,Nloc),FVEGY2D(Mloc,Nloc))
# endif

# if defined (POROUSMEDIA)
     ALLOCATE(Porosity(Mloc,Nloc,Kloc),HyCk(Mloc,Nloc,Kloc))
# endif

# if defined (OBSTACLE)
     ALLOCATE(set_flag(Mloc,Nloc,Kloc),ObsForceX(Mloc,Nloc,Kloc),ObsForceY(Mloc,Nloc,Kloc), &
              ObsForceZ(Mloc,Nloc,Kloc),set_dist_x(Mloc,Nloc,Kloc),set_dist_y(Mloc,Nloc,Kloc), &
              set_dist_z(Mloc,Nloc,Kloc),set_flag_old(Mloc,Nloc,Kloc))
# if defined (DEFORMABLESLIDE)
     ALLOCATE(ObsForceX_gs(Mloc,Nloc),ObsForceY_gs(Mloc,Nloc),ObsForceZ_gs(Mloc,Nloc), &
              set_dist_x_gs(Mloc,Nloc),set_dist_y_gs(Mloc,Nloc),set_dist_z_gs(Mloc,Nloc), &
              set_flag_gs(Mloc,Nloc),set_flag_gs_old(Mloc,Nloc))
# endif
# endif

# if defined (TWOLAYERSLIDE)
     ALLOCATE(Ha(Mloc,Nloc),Hs(Mloc,Nloc),Ha0(Mloc,Nloc),Ua(Mloc,Nloc),Ua0(Mloc,Nloc),Va(Mloc,Nloc),Va0(Mloc,Nloc), &
              HUa0(Mloc,Nloc),HVa0(Mloc,Nloc),HUa(Mloc,Nloc),HVa(Mloc,Nloc),Maska(Mloc,Nloc),Slambda(Mloc,Nloc),  &
              SorceaX(Mloc,Nloc),SorceaY(Mloc,Nloc),DelxHa(Mloc,Nloc),DelxUa(Mloc,Nloc),DelxVa(Mloc,Nloc), &
              DelxHUa(Mloc,Nloc),DelxHVa(Mloc,Nloc),DelyHa(Mloc,Nloc),DelyUa(Mloc,Nloc),DelyVa(Mloc,Nloc), &
              DelyHUa(Mloc,Nloc),DelyHVa(Mloc,Nloc),Grav_Lx(Mloc,Nloc),Grav_Ly(Mloc,Nloc),Grav_Lz(Mloc,Nloc), &
              ForceX_UL(Mloc,NLoc),ForceY_UL(Mloc,Nloc),Grav_LL(Mloc,Nloc),DelxGravLz(Mloc,Nloc),DelyGravLz(Mloc,Nloc),  &
              wbar(Mloc,Nloc))

     ALLOCATE(HaxL(Mloc1,Nloc),HaxR(Mloc1,Nloc),UaxL(Mloc1,Nloc),UaxR(Mloc1,Nloc),VaxL(Mloc1,Nloc),VaxR(Mloc1,Nloc), &
              HayL(Mloc,Nloc1),HayR(Mloc,Nloc1),UayL(Mloc,Nloc1),UayR(Mloc,Nloc1),VayL(Mloc,Nloc1),VayR(Mloc,Nloc1), &
              HUaxL(Mloc1,Nloc),HUaxR(Mloc1,Nloc),HVaxL(Mloc1,Nloc),HVaxR(Mloc1,Nloc),HUayL(Mloc,Nloc1), &
              HUayR(Mloc,Nloc1),HVayL(Mloc,Nloc1),HVayR(Mloc,Nloc1),Eax(Mloc1,Nloc),Eay(Mloc,Nloc1),EaxL(Mloc1,Nloc),  &
              EaxR(Mloc1,Nloc),EayL(Mloc,Nloc1),EayR(Mloc,Nloc1),Fax(Mloc1,Nloc),FaxL(Mloc1,Nloc),FaxR(Mloc1,Nloc), &
              Fay(Mloc,Nloc1),FayL(Mloc,Nloc1),FayR(Mloc,Nloc1),Gax(Mloc1,Nloc),GaxL(Mloc1,Nloc),GaxR(Mloc1,Nloc), &
              Gay(Mloc,Nloc1),GayL(Mloc,Nloc1),GayR(Mloc,Nloc1),Kap(Mloc,Nloc),KapxL(Mloc1,Nloc),KapxR(Mloc1,Nloc), &
              KapyL(Mloc,Nloc1),KapyR(Mloc,Nloc1),DelxKap(Mloc,Nloc),DelyKap(Mloc,Nloc),SaxL(Mloc1,Nloc), &
              SaxR(Mloc1,Nloc),SaxS(Mloc1,Nloc),SayL(Mloc,Nloc1),SayR(Mloc,Nloc1),SayS(Mloc,Nloc1),GravxL(Mloc1,Nloc), &  
              GravxR(Mloc1,Nloc),GravyL(Mloc,Nloc1),GravyR(Mloc,Nloc1))
# endif

! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     ALLOCATE(Uvs(Mloc,Nloc),Vvs(Mloc,Nloc),DUvs(Mloc,Nloc),DVvs(Mloc,Nloc),Dvs(Mloc,Nloc),DelxH0(Mloc,Nloc),Hvs(Mloc,Nloc), &
              Uvs0(Mloc,Nloc),Vvs0(Mloc,Nloc),DUvs0(Mloc,Nloc),DVvs0(Mloc,Nloc),Dvs0(Mloc,Nloc),DelyH0(Mloc,Nloc), &
			  DelxUvs(Mloc,Nloc),DelxVvs(Mloc,Nloc),DelxDUvs(Mloc,Nloc),DelxDVvs(Mloc,Nloc), &
			  DelyUvs(Mloc,Nloc),DelyVvs(Mloc,Nloc),DelyDUvs(Mloc,Nloc),DelyDVvs(Mloc,Nloc), &
			  DelxHvs(Mloc,Nloc),SrcpvsX(Mloc,Nloc),SrchvsX(Mloc,Nloc) ,SrcdvsX(Mloc,Nloc),  &
			  DelyHvs(Mloc,Nloc),SrcpvsY(Mloc,Nloc),SrchvsY(Mloc,Nloc) ,SrcdvsY(Mloc,Nloc),Maskvs(Mloc,Nloc))
     ALLOCATE(DvsxL(Mloc1,Nloc),DvsxR(Mloc1,Nloc),UvsxL(Mloc1,Nloc),UvsxR(Mloc1,Nloc),VvsxL(Mloc1,Nloc),VvsxR(Mloc1,Nloc), &
              H0fx(Mloc1,Nloc),DUvsxL(Mloc1,Nloc),DUvsxR(Mloc1,Nloc),DVvsxL(Mloc1,Nloc),DVvsxR(Mloc1,Nloc), &
			  HvsxL(Mloc1,Nloc),HvsxR(Mloc1,Nloc),SvsxL(Mloc1,Nloc),SvsxR(Mloc1,Nloc),Evsx(Mloc1,Nloc),Fvsx(Mloc1,Nloc),Gvsx(Mloc1,Nloc), &
			  EvsxL(Mloc1,Nloc),EvsxR(Mloc1,Nloc),FvsxL(Mloc1,Nloc),FvsxR(Mloc1,Nloc),GvsxL(Mloc1,Nloc),GvsxR(Mloc1,Nloc))
     ALLOCATE(DvsyL(Mloc,Nloc1),DvsyR(Mloc,Nloc1),UvsyL(Mloc,Nloc1),UvsyR(Mloc,Nloc1),VvsyL(Mloc,Nloc1),VvsyR(Mloc,Nloc1), &
              H0fy(Mloc,Nloc1),DUvsyL(Mloc,Nloc1),DUvsyR(Mloc,Nloc1),DVvsyL(Mloc,Nloc1),DVvsyR(Mloc,Nloc1), &
			  HvsyL(Mloc,Nloc1),HvsyR(Mloc,Nloc1),SvsyL(Mloc,Nloc1),SvsyR(Mloc,Nloc1),Evsy(Mloc,Nloc1),Fvsy(Mloc,Nloc1),Gvsy(Mloc,Nloc1), &
			  EvsyL(Mloc,Nloc1),EvsyR(Mloc,Nloc1),FvsyL(Mloc,Nloc1),FvsyR(Mloc,Nloc1),GvsyL(Mloc,Nloc1),GvsyR(Mloc,Nloc1))
# endif

! added by Cheng for deformable slide
# if defined (DEFORMABLESLIDE)
     ALLOCATE(Ugs(Mloc,Nloc),Vgs(Mloc,Nloc),Wgs(Mloc,Nloc),DUgs(Mloc,Nloc),DVgs(Mloc,Nloc),DWgs(Mloc,Nloc), &
	          Dgs(Mloc,Nloc),Rhogs(Mloc,Nloc),Hgs(Mloc,Nloc),Wtgs(Mloc,Nloc),Qtgs(Mloc,Nloc), &
	          Ugs0(Mloc,Nloc),Vgs0(Mloc,Nloc),Wgs0(Mloc,Nloc),DUgs0(Mloc,Nloc),DVgs0(Mloc,Nloc),DWgs0(Mloc,Nloc), &
			  Dgs0(Mloc,Nloc),Rhogs0(Mloc,Nloc),Wgs00(Mloc,Nloc),Wbgs(Mloc,Nloc),Qbgs(Mloc,Nloc), &
              DelxUgs(Mloc,Nloc),DelxVgs(Mloc,Nloc),DelxWgs(Mloc,Nloc), &
			  DelyUgs(Mloc,Nloc),DelyVgs(Mloc,Nloc),DelyWgs(Mloc,Nloc), &
			  DelxDUgs(Mloc,Nloc),DelxDVgs(Mloc,Nloc),DelxDWgs(Mloc,Nloc),DelxDgs(Mloc,Nloc), &
			  DelyDUgs(Mloc,Nloc),DelyDVgs(Mloc,Nloc),DelyDWgs(Mloc,Nloc),DelyDgs(Mloc,Nloc), &
			  DelxHgs(Mloc,Nloc),DelxH0(Mloc,Nloc),Delx2H0(Mloc,Nloc),Kx(Mloc,Nloc),Cxgs(Mloc,Nloc), &
			  DelyHgs(Mloc,Nloc),DelyH0(Mloc,Nloc),Dely2H0(Mloc,Nloc),Ky(Mloc,Nloc),Cygs(Mloc,Nloc), &
			  SrcmgsX(Mloc,Nloc),SrcpgsX(Mloc,Nloc),SrctgsX(Mloc,Nloc), &
			  SrcmgsY(Mloc,Nloc),SrcpgsY(Mloc,Nloc),SrctgsY(Mloc,Nloc), &
			  SrctgsZ(Mloc,Nloc),UgsP(Mloc,Nloc),VgsP(Mloc,Nloc),WgsP(Mloc,Nloc),DwDt(Mloc,Nloc), &
			  Cxx(Mloc,Nloc),Cxy(Mloc,Nloc),Cxz(Mloc,Nloc),Cyx(Mloc,Nloc),Cyy(Mloc,Nloc), &
			  Cyz(Mloc,Nloc),Czx(Mloc,Nloc),Czy(Mloc,Nloc),Czz(Mloc,Nloc),Czz0(Mloc,Nloc), &
			  Taxx(Mloc,Nloc),Tayx(Mloc,Nloc),Taxy(Mloc,Nloc),Tayy(Mloc,Nloc),Taxz(Mloc,Nloc),Tayz(Mloc,Nloc), &
			  Tezz(Mloc,Nloc),Pss(Mloc,Nloc),Qsgs(Mloc,Nloc), &
			  Tbxx(Mloc,Nloc),Tbxy(Mloc,Nloc),Tbxz(Mloc,Nloc),Tbyx(Mloc,Nloc),Tbyy(Mloc,Nloc), &
			  Tbyz(Mloc,Nloc),Tbzx(Mloc,Nloc),Tbzy(Mloc,Nloc),Tbzz(Mloc,Nloc),Maskgs(Mloc,Nloc), &
			  PhiInt(Mloc,Nloc),PhiBed(Mloc,Nloc))
     ALLOCATE(DgsxL(Mloc1,Nloc),DgsxR(Mloc1,Nloc),HSgsxL(Mloc1,Nloc),HSgsxR(Mloc1,Nloc),H0fx(Mloc1,Nloc), &
	          UgsxL(Mloc1,Nloc),UgsxR(Mloc1,Nloc),VgsxL(Mloc1,Nloc),VgsxR(Mloc1,Nloc),WgsxL(Mloc1,Nloc),WgsxR(Mloc1,Nloc), &
              DUgsxL(Mloc1,Nloc),DUgsxR(Mloc1,Nloc),DVgsxL(Mloc1,Nloc),DVgsxR(Mloc1,Nloc),DWgsxL(Mloc1,Nloc),DWgsxR(Mloc1,Nloc), &
			  Cxgsx(Mloc1,Nloc),SgsxL(Mloc1,Nloc),SgsxR(Mloc1,Nloc), &
			  Egsx(Mloc1,Nloc),Fgsx(Mloc1,Nloc),Ggsx(Mloc1,Nloc),Hgsx(Mloc1,Nloc), &
			  EgsxL(Mloc1,Nloc),EgsxR(Mloc1,Nloc),FgsxL(Mloc1,Nloc),FgsxR(Mloc1,Nloc), &
			  GgsxL(Mloc1,Nloc),GgsxR(Mloc1,Nloc),HgsxL(Mloc1,Nloc),HgsxR(Mloc1,Nloc), &
			  SrcmgsxL(Mloc1,Nloc),SrcmgsxR(Mloc1,Nloc),SrcpgsxL(Mloc1,Nloc),SrcpgsxR(Mloc1,Nloc))
     ALLOCATE(DgsyL(Mloc,Nloc1),DgsyR(Mloc,Nloc1),HSgsyL(Mloc,Nloc1),HSgsyR(Mloc,Nloc1),H0fy(Mloc,Nloc1), &
	          UgsyL(Mloc,Nloc1),UgsyR(Mloc,Nloc1),VgsyL(Mloc,Nloc1),VgsyR(Mloc,Nloc1),WgsyL(Mloc,Nloc1),WgsyR(Mloc,Nloc1), &
              DUgsyL(Mloc,Nloc1),DUgsyR(Mloc,Nloc1),DVgsyL(Mloc,Nloc1),DVgsyR(Mloc,Nloc1),DWgsyL(Mloc,Nloc1),DWgsyR(Mloc,Nloc1), &
			  Cygsy(Mloc,Nloc1),SgsyL(Mloc,Nloc1),SgsyR(Mloc,Nloc1), &
			  Egsy(Mloc,Nloc1),Fgsy(Mloc,Nloc1),Ggsy(Mloc,Nloc1),Hgsy(Mloc,Nloc1), &
			  EgsyL(Mloc,Nloc1),EgsyR(Mloc,Nloc1),FgsyL(Mloc,Nloc1),FgsyR(Mloc,Nloc1), &
			  GgsyL(Mloc,Nloc1),GgsyR(Mloc,Nloc1),HgsyL(Mloc,Nloc1),HgsyR(Mloc,Nloc1), &
			  SrcmgsyL(Mloc,Nloc1),SrcmgsyR(Mloc,Nloc1),SrcpgsyL(Mloc,Nloc1),SrcpgsyR(Mloc,Nloc1))
     ALLOCATE(QbgsC(MlocC,NlocC),WgsC(MlocC,NlocC))
! poisson solver (for NSPCG use)
     neqnsgs = (IendC-Ibeg+1)*(JendC-Jbeg+1)
     ALLOCATE(Coefgs(neqnsgs,5),JCoefgs(5),Rhsgs(neqnsgs))
!
# endif

! poisson solver (for NSPCG use)
     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))
!
     end subroutine allocate_variables
!
!--------------------------------------------------------------------------------------------------------------------------     
!
     subroutine generate_grid
!
!--------------------------------------------------------------------------------------------------------------------------
!
!    (5) generate_grid
!
!    This subroutine is used to generate grids
!
!    Called by: main
!
!    Last update: 20/12/2010, Gangfeng Ma
!--------------------------------------------------------------------------------------------------------------------------
!
     use global
     implicit none
     integer :: i,j,k
!
!    horizontal grid
!
# if defined (PARALLEL)
     x(Ibeg) = npx*(Mloc-2*Nghost)*dx
# else
     x(Ibeg) = Zero
# endif

     do i = Ibeg+1,Mloc1
     x(i) = x(i-1)+dx
     xc(i-1) = x(i-1)+0.5*dx
     enddo

     do i = Ibeg-1,Ibeg-Nghost,-1
     x(i) = x(i+1)-dx
     xc(i) = x(i+1)-0.5*dx
     enddo

# if defined (PARALLEL)
     y(Jbeg) = npy*(Nloc-2*Nghost)*dy
# else
     y(Jbeg) = Zero
# endif

     do j = Jbeg+1,Nloc1
     y(j) = y(j-1)+dy
     yc(j-1) = y(j-1)+0.5*dy
     enddo

     do j = Jbeg-1,Jbeg-Nghost,-1
     y(j) = y(j+1)-dy
     yc(j) = y(j+1)-0.5*dy
     enddo
!
!    vertical grid
!
     if (Ivgrd==1) then
     do k = 1,Kloc
     dsig(k) = 1.0/float(Kglob)
     enddo
     elseif(Ivgrd==2) then
     dsig(Kbeg) = (Grd_R-1.0)/(Grd_R**float(Kglob)-1.0)
     do k = Kbeg+1,Kend
     dsig(k) = dsig(k-1)*Grd_R
     enddo

     do k = 1,Nghost
     dsig(Kbeg-k) = dsig(Kbeg+k-1)
     enddo

     do k = 1,Nghost
     dsig(Kend+k) = dsig(Kend-k+1)
     enddo
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
     sig(k) = sig(k-1)+dsig(k-1)
     sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
     sig(k) = sig(k+1)-dsig(k)
     sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo
!
     end subroutine generate_grid
!
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!  
     subroutine read_bathymetry
!
!------------------------------------------------------------------------------------------
!
!    (6) read_bathymetry
! 
!    This subroutine is used to read bathymetry 
!                                               
!    Called by: main  
!  
!    Gangfeng Ma, 21/12/2010   
!
!    Needs revision to make input depth grid filename arbitrary
!   
!-------------------------------------------------------------------------------------------
!
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG
! added by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     real(SP), dimension(Mglob,Nglob) :: Dvstmp
# endif
# if defined (DEFORMABLESLIDE)
     real(SP), dimension(Mglob,Nglob) :: Dgstmp
# endif
# if defined (LANDSLIDE)
     real(SP) :: alpha0,L0,T,bl,wl,e,kb,kw,x0,xr,xt,yt,zt,alpha1, &
	             slope0,xl1,xl2,yl1,yl2
     ! added by Cheng for 3D triangle landslide
!     real(SP) :: slope0,xgb,xgn,xg,yg,ztm
# endif
     Hc = Zero ! added by Cheng for initialization
!
!    read bathymetry at grid points
!
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 30.0
         enddo
         enddo
       else
         open(5,file='depth.txt',status='old')
!
!    Need to make file name arbitrary  (kirby, 6/27/16)
!
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

! find permanent dry points (how would this be possible?)
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# endif

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center
       if(ANA_BATHY) then
         do j = 1,Nglob
         do i = 1,Mglob
           HCG(i,j) = 0.5
         enddo
         enddo
       else ! not analytical bathymetry, read from depth file
         open(5,file='depth.txt',status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo
       endif
	   
# if defined (PARALLEL)
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# else
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# endif
     endif

     ! collect data into ghost cells 
     call phi_2D_coll(Hc)

     ! save the initial water depth
     Hc0 = Hc

# if defined (LANDSLIDE)
     ! for landslide problem, put slide over the bathymetry
     if(SlideType(1:8)=='RIGID_2D') then
	   ! 2D landslide (added by Cheng for 2d landslide)
       e = 0.75
       kb = 2.0*acosh(1.0/sqrt(e))/SlideL
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         xt = xc(i)-SlideX0
         zt = SlideT/(1-e)*((1.0/cosh(kb*xt))**2-e)
         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
       enddo
       enddo
       ! cosine shape (added by Cheng to smooth the angle out)
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         xt = xc(i)-SlideX0
!		 zt = zero
!		 if (xt<=SlideL/2.0 .and. xt>=-SlideL/2.0) then
!           zt = SlideT/2.0*(1.0-cos(2.0*pi*(xt-SlideL/2.0)/SlideL))
!		 endif
!         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
!       enddo
!       enddo
     elseif(SlideType(1:8)=='RIGID_3D') then
	   ! 3D landslide
       e = 0.717
       alpha0 = SlideAngle*pi/180.
       kb = 2.0*acosh(1.0/e)/SlideL
       kw = 2.0*acosh(1.0/e)/SlideW
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         xt = (xc(i)-SlideX0)*cos(alpha0)+(yc(j)-SlideY0)*sin(alpha0)
         yt = -(xc(i)-SlideX0)*sin(alpha0)+(yc(j)-SlideY0)*cos(alpha0)
         zt = SlideT/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yt)-e)
         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
       enddo
       enddo
!	   ! 3D triangle landslide
!       slope0 = SlopeAngle*pi/180. 
!!      xgb = SlideL*cos(slope0)
!!	   xgn = SlideL/cos(slope0)
!       xgb = SlideL-0.2
!	   xgn = SlideL+0.2
!	   xg = xgn-xgb
!	   ztm = xgb*tan(slope0)
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         if(xc(i)>=SlideX0.and.xc(i)<=SlideX0+xgb) then
!           zt = (xc(i)-SlideX0)*tan(slope0)
!!		   alpha1 = slope0/2.0*(1.0+cos(pi*(xc(i)-SlideX0-xgb)/xgb)) !cosine curve
!!		   yg = zt/2.0*tan(alpha1)
!           yg = 0.1
!		   if (yc(j)<=SlideW/2.0-yg) then
!		     zt = zt
!		   elseif (yc(j)>SlideW/2.0-yg.and.yc(j)<=SlideW/2.0+yg) then
!!		     zt= zt-(yc(j)-(SlideW/2.0-yg))/tan(alpha1)  ! line
!			 zt= zt/2.0*(1.0+cos(pi*(yc(j)-(SlideW/2.0-yg))/(2.0*yg)))  !cosine curve
!		   else
!             zt = 0.0
!		   endif
!         elseif(xc(i)>SlideX0+xgb .and.xc(i)<=SlideX0+xgn) then
!!		   zt = (xg-(xc(i)-SlideX0-xgb))/tan(slope0)
!		   zt = ztm/2.0*(1.0+cos(pi*(xc(i)-SlideX0-xgb)/xg))
!!		   alpha1 = slope0/2.0*(1.0+cos(pi*(xc(i)-SlideX0-xgb)/xg))
!!		   yg = zt/2.0*tan(alpha1)
!           yg = 0.1
!		   if (yc(j)<=SlideW/2.0-yg) then
!		     zt = zt
!		   elseif (yc(j)>SlideW/2.0-yg.and.yc(j)<=SlideW/2.0+yg) then
!!		     zt= zt-(yc(j)-(SlideW/2.0-yg))/tan(alpha1)
!			 zt= zt/2.0*(1.0+cos(pi*(yc(j)-(SlideW/2.0-yg))/(2.0*yg)))
!		   else
!             zt = 0.0
!		   endif
!         else
!           zt = 0.0
!         endif
!         Hc(i,j) = Hc0(i,j)-max(0.0,zt)
!       enddo
!       enddo
     elseif(SlideType(1:9)=='RIGID_SLP') then
       e = 0.717
       slope0 = SlopeAngle*pi/180.
       kb = 2.0*acosh(1.0/e)/SlideL
       kw = 2.0*acosh(1.0/e)/SlideW
       xl1 = SlideX0-0.5*SlideL*cos(slope0)
       xl2 = SlideX0+0.5*SlideL*cos(slope0)
       yl1 = SlideY0-0.5*SlideW
       yl2 = SlideY0+0.5*SlideW
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=xl1.and.xc(i)<=xl2.and.yc(j)>=yl1.and.yc(j)<=yl2) then
           xt = (xc(i)-SlideX0)/cos(slope0)
           yt = yc(j)-SlideY0
           zt = SlideT/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yt)-e)
           Hc(i,j) = Hc0(i,j)-max(0.0,zt)
         endif
       enddo
       enddo
     endif
# endif

     ! added by Cheng for fluid slide and deformable slide
# if defined (FLUIDSLIDE) || defined (DEFORMABLESLIDE)
# if defined (FLUIDSLIDE)
     Dvs = zero
     Hvs  = zero
     open(5001,file=TRIM(Slide_File),status='old')
     do j = 1,Nglob
       read(5001,*) (Dvstmp(i,j),i=1,Mglob) !slide shape
     enddo
     close(5001)

     do i = Ibeg,Iend
     do j = Jbeg,Jend
# if defined (PARALLEL)
       iglob = npx*(Mloc-2*Nghost)+i-Nghost
       jglob = npy*(Nloc-2*Nghost)+j-Nghost
# else
       iglob = i-Nghost
       jglob = j-Nghost
# endif

       Dvs(i,j) = Dvstmp(iglob,jglob) !slide shape
       if(Dvs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Dvs(i,j) = SLIDE_MINTHICK
	   endif
       Hvs(i,j) = Hc0(i,j)-Dvs(i,j) !chi
     enddo
     enddo
	 
     ! wetting-drying mask for fluid slide
     Maskvs = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if(Dvs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Maskvs(i,j) = 0
       endif
     enddo
     enddo

! interpolate into grid center
     do j=Jbeg,Jend
     do i=Ibeg,Iend
	   if(Dvs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Hc(i,j) = Hc0(i,j) 
	   else
	     Hc(i,j) = Hc0(i,j)-Dvs(i,j) 
	   endif
     enddo
     enddo

     ! collect data into ghost cells
# if defined (PARALLEL)
     call phi_int_exch(Maskvs)
# endif
     call phi_2D_coll(Hc)
     call phi_2D_coll(Dvs) 
     call phi_2D_coll(Hvs)  
# endif

# if defined (DEFORMABLESLIDE)
     Dgs = zero
     Hgs  = zero
     open(5001,file=TRIM(Slide_File),status='old')
     do j = 1,Nglob
       read(5001,*) (Dgstmp(i,j),i=1,Mglob) !slide shape
     enddo
     close(5001)

     do i = Ibeg,Iend
     do j = Jbeg,Jend
# if defined (PARALLEL)
       iglob = npx*(Mloc-2*Nghost)+i-Nghost
       jglob = npy*(Nloc-2*Nghost)+j-Nghost
# else
       iglob = i-Nghost
       jglob = j-Nghost
# endif

       Dgs(i,j) = Dgstmp(iglob,jglob) !slide shape
       if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Dgs(i,j) = SLIDE_MINTHICK
	   endif
       Hgs(i,j) = Hc0(i,j)-Dgs(i,j) !chi
     enddo
     enddo
	 
     ! wetting-drying mask for fluid slide
     Maskgs = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Maskgs(i,j) = 0
       endif
     enddo
     enddo

! interpolate into grid center
     do j=Jbeg,Jend
     do i=Ibeg,Iend
	   if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Hc(i,j) = Hc0(i,j) 
	   else
	     Hc(i,j) = Hc0(i,j)-Dgs(i,j) 
	   endif
     enddo
     enddo
	 
! second derivatives of depth at cell center
     do j = 1,Nloc
     do i = 2,Mloc-1
       Delx2H0(i,j)=((Hc0(i+1,j)-Hc0(i,j))/dx-(Hc0(i,j)-Hc0(i-1,j))/dx)/dx;
     enddo
     enddo
     do j = 1,Nloc
       Delx2H0(1,j) = Delx2H0(2,j)
       Delx2H0(Mloc,j) = Delx2H0(Mloc-1,j)
     enddo
	 
     do j = 2,Nloc-1
     do i = 1,Mloc
        Dely2H0(i,j) = ((Hc0(i,j+1)-Hc0(i,j))/dy-(Hc0(i,j)-Hc0(i,j-1))/dy)/dy;
     enddo
     enddo
     do i = 1,Mloc
       Dely2H0(i,1) = Dely2H0(i,2)
       Dely2H0(i,Nloc) = Dely2H0(i,Nloc-1)
     enddo

     ! collect data into ghost cells
# if defined (PARALLEL)
     call phi_int_exch(Maskgs)
	 call phi_2D_exch(Delx2H0)
	 call phi_2D_exch(Dely2H0)
# endif
     call phi_2D_coll(Hc)
     call phi_2D_coll(Dgs) 
     call phi_2D_coll(Hgs)  
# endif
	 
     ! reconstruct base depth below the fluid slide at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       H0fx(i,j) = 0.5*(Hc0(i-1,j)+Hc0(i,j))
     enddo
     H0fx(1,j) = Hc0(1,j)
     H0fx(Mloc1,j) = Hc0(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       H0fy(i,j) = 0.5*(Hc0(i,j-1)+Hc0(i,j))
     enddo
     H0fy(i,1) = Hc0(i,1)
     H0fy(i,Nloc1) = Hc0(i,Nloc)
     enddo

     ! derivatives of base depth below the fluid slide at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH0(i,j) = (H0fx(i+1,j)-H0fx(i,j))/dx
       DelyH0(i,j) = (H0fy(i,j+1)-H0fy(i,j))/dy
     enddo
     enddo
# if defined (PARALLEL)
	call phi_2D_exch(DelxH0)
	call phi_2D_exch(DelyH0)
# endif

# endif

     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-1000.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hc(1,j)
     Hfx(Mloc1,j) = Hc(Mloc,j)
     enddo
	 !hfx0 is flux of bed, which is different from hfx at wet-dry front (Cheng)
	 Hfx0 = Hfx	

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hc(i,1)
     Hfy(i,Nloc1) = Hc(i,Nloc)
     enddo
     Hfy0 = Hfy

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     end subroutine read_bathymetry


     subroutine read_bathymetry_comprehensive
!------------------------------------------------------ 
!    This subroutine is used to read bathymetry                                                
!    Called by  
!       main    
!    update: 21/12/2010, Gangfeng Ma    
!    update: 04/13/2012, Fengyan Shi
!    fyshi make a standard slide application    
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob,Kslide
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG
!     real(SP), dimension(Mglob,Nglob) :: HCG
# if defined (LANDSLIDE_COMPREHENSIVE)
     real(SP) :: Hd,alpha0,L0,T,bl,wl,e,kb,kw,x0,x1,x2,xr,Hb(Mloc,Nloc),&
              xt,xt1,zt,zt1,ht,yt,Hslide(Mloc1,Nloc1)
     real(SP) :: edge_width
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
   IF(ANA_BATHY)THEN
!      if use analytical slide, 1)need iteration 2) ac is time-dependent
     Hd=1.5  
     T = T_slide(1)
     bl = b_slide(1)
     wl = w_slide(1)
     e = e_slide(1)
     alpha0 = slope_slide(1)*3.1415926/180.0
     kb = 2.0*acosh(1.0/e)/bl
     kw = 2.0*acosh(1.0/e)/wl
     x0=x0_slide(1)
     L0 = x0/cos(alpha0)-T*tan(alpha0)
     xr = 0.0
     x1 = (L0-bl/2.)*cos(alpha0)+xr
     x2 = (L0+bl/2.)*cos(alpha0)+xr

   ELSE

    DO Kslide=1,NumSlides

     alpha_slide(Kslide) = alpha_slide(Kslide)*3.1415926/180.0
     kb_slide = 2.0*acosh(1.0/e_slide(Kslide))/b_slide(Kslide)
     kw_slide = 2.0*acosh(1.0/e_slide(Kslide))/w_slide(Kslide)

     ac_slide(Kslide) = 9.81*sin(slope_slide(Kslide)/180.0*3.1415926)

   ENDDO ! end Kslide

   ENDIF

     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
# if defined (PARALLEL)
       if(myid.eq.0) write(*,*) 'DEPTH_TYPE has to be cell_center'
# else
       write(*,*) 'DEPTH_TYPE has to be cell_center'
# endif  
     endif

# endif

     ! read bathymetry at grid points
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 0.35
         enddo
         enddo
       else
         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

       ! find pernament dry points
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# endif

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center

       if(ANA_BATHY) then

# if defined (LANDSLIDE_COMPREHENSIVE)
         ! base bathymetry                                                                     
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xc(i)<=xr) then
             Hb(i,j) = -(xr-xc(i))*tan(alpha0)
           elseif(xc(i)<=Hd/tan(alpha0)+xr) then
             Hb(i,j) = (xc(i)-xr)*tan(alpha0)
           else
             Hb(i,j) = Hd
           endif
           ! temporarily no runup                                                              
        !   Hb(i,j) = max(0.01,Hb(i,j))
         enddo
         enddo


         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xc(i)<=x1.or.xc(i)>=x2.or.yc(j)>=wl/2.) then
             Hc(i,j) = Hb(i,j)
           else
             iter = 1
             xt = (xc(i)-xr)/cos(alpha0)-L0
             zt = T/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yc(j))-e)
 60          xt1 = ((xc(i)-xr)/cos(alpha0)-zt*tan(alpha0))-L0
             zt1 = T/(1-e)*(1.0/cosh(kb*xt1)/cosh(kw*yc(j))-e)
             if(abs(zt1-zt)/abs(zt)>1.e-8) then
               iter = iter+1
               if(iter>20) write(*,*) 'too many iterations!'
               zt = zt1
               goto 60
             endif
             Hc(i,j) = Hb(i,j)-max(0.0,zt1)/cos(alpha0)
           endif
         enddo
         enddo
# endif
 
       else ! not analytical bathymetry, read from depth file

         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo


# if defined (PARALLEL)
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# else
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# endif

# if defined (LANDSLIDE_COMPREHENSIVE)
         ! add the landslide

      DepC0 = Hc   ! save base bathymetry 



    IF(SLIDE_SHAPE_BOX) THEN  ! box shape
       Hb=Hc

       ! add landslide
! to keep consistent with slide bathy at grid, assuming slide bathy at grid point

     DO Kslide =1,NumSlides
       edge_width=T_slide(Kslide)
       Hslide=0.0_SP

       do j = Jbeg,Jend1
       do i = Ibeg,Iend1
         xt = (x(i)-x0_slide(Kslide))*cos(alpha_slide(Kslide))+(y(j)-y0_slide(Kslide))*sin(alpha_slide(Kslide))
         yt = -(x(i)-x0_slide(Kslide))*sin(alpha_slide(Kslide))+(y(j)-y0_slide(Kslide))*cos(alpha_slide(Kslide))


         IF(abs(xt)<=0.5*b_slide(Kslide)-0.5*edge_width.AND.abs(yt)<=0.5*w_slide(Kslide) &
                         +SlideErrorRange)THEN
            ht=T_slide(Kslide)

         ELSEIF(abs(xt)>0.5*b_slide(Kslide)-0.5*edge_width  &
            .AND.abs(xt)<0.5*b_slide(Kslide)+0.5*edge_width.AND.abs(yt)<=0.5*w_slide(Kslide) &
                         +SlideErrorRange)THEN
            ht=T_slide(Kslide)-T_slide(Kslide)*(abs(xt)-0.5*b_slide(Kslide)+0.5*edge_width)/edge_width
         ELSE
            ht=0.0_SP
         ENDIF
         Hslide(I,J) = max(0.0,ht)
       enddo
       enddo

! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         ht=0.25*(Hslide(I,J)+Hslide(I+1,J)+Hslide(I,J+1) &
                            +Hslide(I+1,J+1))
       if (SMFCut) then
         Hb(i,j) = Hb(i,j)+ht
       else
         Hc(i,j) = Hc(i,j)-ht
       endif
       enddo
       enddo

     ENDDO ! end Kslide

       do j=Jbeg,Jend
       do i=Ibeg,Iend
       DepC0(i,j)=Hb(i,j)
       enddo
       enddo


    ELSE ! for ellipse shape

! to keep consistent with slide at grid, assuming slide bathy at grid, not ceter
         do j = Jbeg,Jend1
         do i = Ibeg,Iend1
           xt = (x(i)-x0_slide(1))*cos(alpha_slide(1))+(y(j)-y0_slide(1))*sin(alpha_slide(1))
           yt = -(x(i)-x0_slide(1))*sin(alpha_slide(1))+(y(j)-y0_slide(1))*cos(alpha_slide(1))
           ht = T_slide(1)/(1-e_slide(1))*(1./cosh(kb_slide*xt)/cosh(kw_slide*yt)-e_slide(1))
           Hslide(I,J) = max(0.0,ht)
         enddo
         enddo

! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         ht=0.25*(Hslide(I,J)+Hslide(I+1,J)+Hslide(I,J+1) &
                            +Hslide(I+1,J+1))
        if (SMFCut) then
           DepC0(i,j) = DepC0(i,j)+ht
        else
           Hc(i,j) = Hc(i,j)-ht
        endif
       enddo
       enddo

     ENDIF ! end box shape

# endif

     endif  ! end analytical 

     endif  ! end cell center

     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-1000.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! collect data into ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces

! gfma fixed the bug 12/20/2011
    do j = 1,Nloc
    do i = 2,Mloc
      Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
    enddo
    Hfx(1,j) = Hc(1,j)
    Hfx(Mloc1,j) = Hc(Mloc,j)
    enddo

    do i = 1,Mloc
    do j = 2,Nloc
      Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
    enddo
    Hfy(i,1) = Hc(i,1)
    Hfy(i,Nloc1) = Hc(i,Nloc)
    enddo


     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo


     end subroutine read_bathymetry_comprehensive
	 
	 
	 
!
!------------------------------------------------------------------------------------------
!
     subroutine initial
!
!--------------------------------------------------------------------------------  
!
!    (7) initial
!
!    This subroutine is used to initialize model run 
!
!    Called by: main 
!
!    Last update: 21/12/2010, Gangfeng Ma  
!
!--------------------------------------------------------------------------------
!
     use global
     implicit none
     integer  :: i,j,k,n,m,nmax,iglob,jglob
     real(SP) :: xsol(80),zsol(80),zmax,xmax,xterp,zterp,tmp,zc(Kloc), &
                 utmp1,wtmp1,utmp2,wtmp2,xk(321,16),zk(321,16),  &
                 uk(321,16),wk(321,16)
     real(SP) :: Ufric,Zlev1,Zlev,mud_dens,Zslide(Mloc,Nloc),conc_slide, &
                 alpha0,L0,T,bl,hs0,ls0,ls1,ls2,lsx,lsx1,hslide,  &                   
                 eslide,xt,yt,zt,kb,kw,Slope,Xslide,SlideX1,SlideX2
     real(SP), dimension(Mglob,Nglob) :: EtaG
     real(SP), dimension(Mglob,Nglob,Kglob) :: UG,VG,WG,SaliG

     ! simulation time
     TIME = Zero
     RUN_STEP = 0
     dt = dt_ini     
     Screen_Count = Zero
     Plot_Count = Zero
     Plot_Count_Stat = Zero

# if defined (SEDIMENT)
     Update_Bed_T = Zero
     nbed = 0
# endif

     Icount = 0
     
     ! working arrays
     D = Zero
     U = Zero
     V = Zero
     W = Zero
     P = Zero
     Omega = Zero
     DU = Zero
     DV = Zero
     DW = Zero
     D0 = Zero
     Eta0 = Zero
     DU0 = Zero
     DV0 = Zero
     DW0 = Zero
     Uf = Zero
     Vf = Zero
     Wf = Zero
     Rho = Rho0
     
     ! source terms
     SourceC = Zero
     SourceX = Zero
     SourceY = Zero

     ! fluxes
     DxL = Zero
     DxR = Zero
     DyL = Zero
     DyR = Zero
     UxL = Zero
     UxR = Zero
     UyL = Zero
     UyR = Zero
     UzL = Zero
     UzR = Zero
     VxL = Zero
     VxR = Zero
     VyL = Zero
     VyR = Zero
     VzL = Zero
     VzR = Zero
     WxL = Zero
     WxR = Zero
     WyL = Zero
     WyR = Zero
     WzL = Zero
     WzR = Zero
     DUxL = Zero
     DUxR = Zero
     DUyL = Zero
     DUyR = Zero
     DVxL = Zero
     DVxR = Zero
     DVyL = Zero
     DVyR = Zero
     DWxL = Zero
     DWxR = Zero
     DWyL = Zero
     DWyR = Zero
     OzL = Zero
     OzR = Zero
     SxL = Zero
     SxR = Zero
     SxS = Zero
     SyL = Zero
     SyR = Zero
     SyS = Zero
     Ex = Zero
	 ExL = Zero
	 ExR = Zero
     Ey = Zero
     EyL = Zero
	 EyR = Zero
     Fx = Zero
	 FxL = Zero
	 FxR = Zero
     Fy = Zero
     FyL = Zero
	 FyR = Zero
     Fz = Zero
     Gx = Zero
	 GxL = Zero
	 GxR = Zero
     Gy = Zero
	 GyL = Zero
	 GyR = Zero
     Gz = Zero
     Hx = Zero
	 HxL = Zero
	 HxR = Zero
     Hy = Zero
	 HyL = Zero
	 HyR = Zero
     Hz = Zero
     EtaxL = Zero
     EtaxR = Zero
     EtayL = Zero
     EtayR = Zero
     DelxEta = Zero
     DelyEta = Zero
     DeltH = Zero
     DeltHo = Zero
     Delt2H = Zero
     DelxD = Zero
     DelyD = Zero
     DelxU = Zero
     DelyU = Zero
     DelzU = Zero
     DelxV = Zero
     DelyV = Zero
     DelzV = Zero
     DelxW = Zero
     DelyW = Zero
     DelzW = Zero
     DelxDU = Zero
     DelyDU = Zero
     DelxDV = Zero
     DelyDV = Zero
     DelxDW = Zero
     DelyDW = Zero
     DelzO = Zero
     Sponge = One
     Cmu = Visc
     CmuHt = Zero
     CmuVt = Zero
     CmuR  = Zero
     Richf = Zero

     Diffxx = Zero
     Diffxy = Zero
     Diffxz = Zero
     Diffyx = Zero
     Diffyy = Zero
     Diffyz = Zero
     Diffzx = Zero
     Diffzy = Zero
     Diffzz = Zero

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     Ein_X0 = Zero
     Din_X0 = Zero
     Uin_Xn = Zero
	 Uin_Xni = Zero
     Uin_Xni0 = Zero
     Vin_Xn = Zero
     Win_Xn = Zero
     Ein_Xn = Zero
     Din_Xn = Zero   

     HeightMax = Zero ! added by cheng for recording Hmax
	 
     Setup = Zero
     WaveHeight = Zero
     Umean = Zero
     Vmean = Zero
     Num_Zero_Up = 0
     Emax = -1000.
     Emin = 1000.

     WdU = Zero
     WdV = Zero

     Lag_Umean = Zero
     Lag_Vmean = Zero
     Lag_Wmean = Zero
     Euler_Umean = Zero
     Euler_Vmean = Zero
     Euler_Wmean = Zero

# if defined (SEDIMENT)
     Smean = Zero
     Euler_Smean = Zero
     Qbedx = Zero
     Qbedy = Zero
     QbxAvg = Zero
     QbyAvg = Zero
# endif

     ExtForceX = Zero
     ExtForceY = Zero

     ! baroclinic terms
     DRhoX = Zero
     DRhoY = Zero

# if defined (BUBBLE)
     Nbg = Zero
     DNbg = Zero
     DNbg0 = Zero
     Vbg = Zero
# endif

# if defined (SEDIMENT)
     Conc = Zero
     DConc = Zero
     DConc0 = Zero
     Bed = Zero
     TotDep = Zero
     Sed_X0 = Zero
     Sed_Xn = Zero
# endif

# if defined (SALINITY)
     Sali = Zero
     DSali = Zero
     DSali0 = Zero
     Sin_X0 = Zero
     Sin_Xn = Zero
# endif

# if defined (TEMPERATURE)
     Temp = 10.
# endif

# if defined (VEGETATION)
     FVegH = VegH
     foliage = one
# endif

# if defined (TWOLAYERSLIDE)
     Ha0 = Zero
     Ua0 = Zero
     Va0 = Zero
     Ha = Zero
     Ua = Zero
     Va = Zero
     HUa = Zero
     HVa = Zero
     HUa0 = Zero
     HVa0 = Zero
     HaxL = Zero
     HaxR = Zero
     UaxL = Zero
     UaxR = Zero
     VaxL = Zero
     VaxR = Zero
     HayL = Zero
     HayR = Zero
     UayL = Zero
     UayR = Zero
     VayL = Zero
     VayR = Zero
     HUaxL = Zero
     HUaxR = Zero
     HVaxL = Zero
     HVaxR = Zero
     HUayL = Zero
     HUayR = Zero
     HVayL = Zero
     HVayR = Zero
     Eax = Zero
     Eay = Zero
     EaxL = Zero
     EaxR = Zero
     EayL = Zero
     EayR = Zero
     Fax = Zero
     Fay = Zero
     FaxL = Zero
     FaxR = Zero
     FayL = Zero
     FayR = Zero
     Gax = Zero
     Gay = Zero
     GaxL = Zero
     GaxR = Zero
     GayL = Zero
     GayR = Zero
     SaxL = Zero
     SaxR = Zero
     SaxS = Zero
     SayL = Zero
     SayR = Zero
     SayS = Zero
     wbar = Zero

     Kap = one
     DelxKap = Zero
     DelyKap = Zero

     SorceaX = Zero
     SorceaY = Zero

     ForceX_UL = Zero
     ForceY_UL = Zero

     Grav_LL = Grav

     Maska = 1  
   
     call specify_slide
# endif

! add by Cheng for fluid slide
# if defined (FLUIDSLIDE)
     Uvs = Zero; Vvs = Zero; DUvs = Zero; DVvs = Zero
     Uvs0 = Zero;Vvs0 = Zero;DUvs0 = Zero;DVvs0 = Zero
     DelxUvs = Zero;DelxVvs = Zero;DelxDUvs = Zero;DelxDVvs = Zero
	 DelyUvs = Zero;DelyVvs = Zero;DelyDUvs = Zero;DelyDVvs = Zero
	 DelxHvs = Zero;SrcpvsX = Zero;SrchvsX = Zero; SrcdvsX = Zero
	 DelyHvs = Zero;SrcpvsY = Zero;SrchvsY = Zero; SrcdvsY = Zero
     DvsxL = Zero;DvsxR = Zero;UvsxL = Zero;UvsxR = Zero;VvsxL = Zero;VvsxR = Zero
	 DUvsxL = Zero;DUvsxR = Zero;DVvsxL = Zero;DVvsxR = Zero
	 HvsxL = Zero;HvsxR = Zero;SvsxL = Zero;SvsxR = Zero;Evsx = Zero;Fvsx = Zero;Gvsx = Zero
	 EvsxL = Zero;EvsxR = Zero;FvsxL = Zero;FvsxR = Zero;GvsxL = Zero;GvsxR = Zero
     DvsyL = Zero;DvsyR = Zero;UvsyL = Zero;UvsyR = Zero;VvsyL = Zero;VvsyR = Zero
	 DUvsyL = Zero;DUvsyR = Zero;DVvsyL = Zero;DVvsyR = Zero
	 HvsyL = Zero;HvsyR = Zero;SvsyL = Zero;SvsyR = Zero;Evsy = Zero;Fvsy = Zero;Gvsy = Zero
	 EvsyL = Zero;EvsyR = Zero;FvsyL = Zero;FvsyR = Zero;GvsyL = Zero;GvsyR = Zero
# endif

! add by Cheng for deformable slide
# if defined (DEFORMABLESLIDE)
     do j = 1,Nloc
     do i = 1,Mloc
       if(Maskgs(i,j)==1) then
         Ugs(i,j) = SLIDE_INIU
		 Vgs(i,j) = SLIDE_INIV
		 Wgs(i,j) = SLIDE_INIW
	   else
	     Ugs(i,j) = Zero
		 Vgs(i,j) = Zero
		 Wgs(i,j) = Zero
       endif
       DUgs(i,j) = Ugs(i,j)*Dgs(i,j)
       DVgs(i,j) = Vgs(i,j)*Dgs(i,j)
	   DWgs(i,j) = Wgs(i,j)*Dgs(i,j)
     enddo
     enddo
     Ugs0=Ugs;Vgs0=Vgs;Wgs0=Wgs;DUgs0=DUgs;DVgs0=DVgs;DWgs0=DWgs;Wgs00=Wgs0
	 Wtgs=Zero;Wbgs=Zero;Qtgs=Zero;Qbgs=Zero
     DelxUgs=Zero;DelxVgs=Zero;DelxWgs=Zero
	 DelxDUgs=Zero;DelxDVgs=Zero;DelxDWgs=Zero
	 DelxDgs=Zero;DelxHgs=Zero
	 DelyUgs=Zero;DelyVgs=Zero;DelyWgs=Zero
	 DelyDUgs=Zero;DelyDVgs=Zero;DelyDWgs=Zero
	 DelyDgs=Zero;DelyHgs=Zero
     Kx=Zero;SrcmgsX=Zero;SrcpgsX=Zero;SrctgsX=Zero
     Ky=Zero;SrcmgsY=Zero;SrcpgsY=Zero;SrctgsY=Zero;SrctgsZ=Zero
	 UgsP=Zero;VgsP=Zero;WgsP=Zero;DwDt=Zero
     Cxx=Zero;Cxy=Zero;Cxz=Zero;Cyx=Zero;Cyy=Zero;Cyz=Zero;Czx=Zero;Czy=Zero;Czz=Zero;Czz0=Zero
	 Taxx=Zero;Tayx=Zero;Taxy=Zero;Tayy=Zero;Taxz=Zero;Tayz=Zero;
	 Tbxx=Zero;Tbxy=Zero;Tbxz=Zero;Tbyx=Zero;Tbyy=Zero;Tbyz=Zero;Tbzx=Zero;Tbzy=Zero;Tbzz=Zero
	 Tezz=Zero;Pss=Zero;Qsgs=Zero
     DgsxL=Zero;DgsxR=Zero;HSgsxL=Zero;HSgsxR=Zero
	 UgsxL=Zero;UgsxR=Zero;VgsxL=Zero;VgsxR=Zero;WgsxL=Zero;WgsxR=Zero
	 DUgsxL=Zero;DUgsxR=Zero;DVgsxL=Zero;DVgsxR=Zero;DWgsxL=Zero;DWgsxR=Zero
	 Cxgsx=Zero;SgsxL=Zero;SgsxR=Zero;Egsx=Zero;Fgsx=Zero;Ggsx=Zero;Hgsx=Zero
	 EgsxL=Zero;EgsxR=Zero;FgsxL=Zero;FgsxR=Zero
	 GgsxL=Zero;GgsxR=Zero;HgsxL=Zero;HgsxR=Zero
	 SrcmgsxL=Zero;SrcmgsxR=Zero;SrcpgsxL=Zero;SrcpgsxR=Zero
     DgsyL=Zero;DgsyR=Zero;HSgsyL=Zero;HSgsyR=Zero
	 UgsyL=Zero;UgsyR=Zero;VgsyL=Zero;VgsyR=Zero;WgsyL=Zero;WgsyR=Zero
	 DUgsyL=Zero;DUgsyR=Zero;DVgsyL=Zero;DVgsyR=Zero;DWgsyL=Zero;DWgsyR=Zero
	 Cygsy=Zero;SgsyL=Zero;SgsyR=Zero;Egsy=Zero;Fgsy=Zero;Ggsy=Zero;Hgsy=Zero
	 EgsyL=Zero;EgsyR=Zero;FgsyL=Zero;FgsyR=Zero
	 GgsyL=Zero;GgsyR=Zero;HgsyL=Zero;HgsyR=Zero
	 SrcmgsyL=Zero;SrcmgsyR=Zero;SrcpgsyL=Zero;SrcpgsyR=Zero
	 QbgsC=Zero;WgsC=Zero
	 if(trim(RHEO_OPT)=='VISCOUS') then
	    Rhogs = SLIDE_DENSITY
	 elseif(trim(RHEO_OPT)=='GRANULAR') then
	    Rhogs=Zero
		PhiInt=Zero;PhiBed=Zero
	    do j = Jbeg,Jend
        do i = Ibeg,Iend
!		  if(Mask(i,j)==1) then
	        Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC + Rho0*(1.0-SLIDE_CONC)
!		  else
!		    Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC
!		  endif
		  if(Mask(i,j)==1) then
	        PhiInt(i,j) = PhiInt_F
			PhiBed(i,j) = PhiBed_F
		  else
	        PhiInt(i,j) = PhiInt_A
			PhiBed(i,j) = PhiBed_A
		  endif
		enddo
		enddo
		call phi_2D_coll(Rhogs)
		call phi_2D_coll(PhiInt)
		call phi_2D_coll(PhiBed)
	 endif
	 Rhogs0=Rhogs
	 call flux_coeff_gs
# endif

# if defined (OBSTACLE)
     ObsForceX=Zero;ObsForceY=Zero;ObsForceZ=Zero
	 set_dist_x=Zero;set_dist_y=Zero;set_dist_z=Zero
	 set_flag=0;set_flag_old=0;
# if defined (DEFORMABLESLIDE)
     ObsForceX_gs=Zero;ObsForceY_gs=Zero;ObsForceZ_gs=Zero
	 set_dist_x_gs=Zero;set_dist_y_gs=Zero;set_dist_z_gs=Zero
	 set_flag_gs=0;set_flag_gs_old=0;
# endif
# endif


     Tke = Zero
     Eps = Zero
     DTke = Zero
     DEps = Zero
     DTke0 = Zero
     DEps0 = Zero

     !added by m.derakhti            
     Wsx = Zero
     Wsy = Zero
	 
# if defined (VEGETATION)
     Tke_w = Zero
     Eps_w = Zero
     DWke = Zero
     DWke0 = Zero
# endif

     ! wave breaking mask
     Brks = 0

     ! pressure boundary
     Bc_Prs = Zero

     ! initial surface elevation (user-specified)
     Eta = Zero

# if defined (POROUSMEDIA)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(xc(i)<=35.0) then
         Eta(i,j) = 0.0
       elseif(xc(i)<=38.0) then
         Eta(i,j) = (xc(i)-35.0)*0.1
       else
         Eta(i,j) = 0.3
       endif
     enddo
     enddo
# endif

     if(IVturb==20.or.IVturb==30) then
       call RandomU
     endif

     if(INITIAL_EUVW) then

        ! initial condition for eta
        open(21,file='eta0.txt',status='old')
        do j = 1,Nglob
          read(21,*) (EtaG(i,j),i=1,Mglob)
        enddo
      
        ! initial condition for U, V, W
        open(22,file='uvw0.txt',status='old')
        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (UG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (VG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (WG(i,j,k),i=1,Mglob)
        enddo
        enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Eta(i,j) = EtaG(iglob,jglob)
         do k = Kbeg,Kend
           U(i,j,k) = UG(iglob,jglob,k-Nghost)
           V(i,j,k) = VG(iglob,jglob,k-Nghost)
           W(i,j,k) = WG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         Eta(i,j) = EtaG(iglob,jglob)
         do k = Kbeg,Kend
           U(i,j,k) = UG(iglob,jglob,k-Nghost)
           V(i,j,k) = VG(iglob,jglob,k-Nghost)
           W(i,j,k) = WG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo
# endif        

!       ! solitary wave from Tanaka solution
!       open(21,file='soliton.dat')
!       do n = 1,80
!         read(21,*) i,xsol(n),zsol(n),tmp,tmp
!       enddo
!       close(21)
!
!       ! find the peak location                                                   
!       zmax = -1.0e+10
!       do n = 1,80
!         if(zsol(n)>zmax) then
!           zmax = zsol(n)
!           xmax = xsol(n)
!         endif
!       enddo
!
!       ! move the peak to x = 3.0m                                                
!       do n = 1,80
!         xsol(n) = xsol(n)+8.0-xmax
!       enddo
!
!       ! interpolate into computational grid
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         do n = 2,80
!           if(xc(i)>=xsol(n-1).and.xc(i)<xsol(n)) then
!             xterp = (xc(i)-xsol(n-1))/(xsol(n)-xsol(n-1))
!             Eta(i,j) = (1.0-xterp)*zsol(n-1)+xterp*zsol(n)
!           endif
!         enddo
!       enddo
!       enddo
!
!       open(22,file='plotuv.dat')
!       do n = 1,16
!       do m = 1,321
!         read(22,*) k,xk(m,n),zk(m,n),tmp,uk(m,n),wk(m,n)
!         xk(m,n) = xk(m,n)+8.0
!       enddo
!       enddo
!       close(22)
!       
!       do j = Jbeg,Jend
!       do k = Kbeg,Kend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         zc(k) = (1.0+Eta(i,j))*sigc(k)-1.0
!
!         do n = 2,16
!         do m = 2,321
!           if(xc(i)>=xk(m-1,n).and.xc(i)<xk(m,n).and.  &
!                    zc(k)>=zk(m,n-1).and.zc(k)<zk(m,n)) then
!             xterp = (xc(i)-xk(m-1,n))/(xk(m,n)-xk(m-1,n))
!             zterp = (zc(k)-zk(m,n-1))/(zk(m,n)-zk(m,n-1)) 
!             utmp1 = (1.0-xterp)*uk(m-1,n-1)+xterp*uk(m,n-1)
!             wtmp1 = (1.0-xterp)*wk(m-1,n-1)+xterp*wk(m,n-1)
!             utmp2 = (1.0-xterp)*uk(m-1,n)+xterp*uk(m,n)
!             wtmp2 = (1.0-xterp)*wk(m-1,n)+xterp*wk(m,n)
!             
!             U(i,j,k) = (1.0-zterp)*utmp1+zterp*utmp2
!             W(i,j,k) = (1.0-zterp)*wtmp1+zterp*wtmp2
!           endif
!         enddo
!         enddo
!       enddo
!       enddo
!       enddo
! 100   continue
     endif    

     ! wetting-drying mask
     ! Mask: 1 - wet; 0 - dry
     ! Mask_Struct: 0 - permanent dry point
     ! Mask9: mask for itself and 8 elements around
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo
     Mask = Mask*Mask_Struct
	 
     ! collect data into ghost cells (move from before above loop to after it by Cheng)
     call phi_2D_coll(Eta)
     Eta0 = Eta

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)
# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     ! total water depth and flux
     D = max(Hc+Eta, MinDep)

     call vel_bc
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
# endif

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     if(VISCOUS_FLOW) then
       ! initial seeding values for turbulence
!       Tke_min = 0.5*(1.4e-3)**2
!       Eps_min = 0.09*Tke_min**2/(0.1*Visc)
!       Tke_min = 1.e-12
!       Eps_min = 0.09*Tke_min**2/(1.e-4*Visc)
       Tke_min = 1.e-9
       Eps_min = 1.e-9
       if (RNG) then
          Cmut_min = 8.5e-2*Tke_min**2/Eps_min
       else
          Cmut_min = 9.0e-2*Tke_min**2/Eps_min
       endif
       do k = 1,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Tke(i,j,k) = Tke_min
         Eps(i,j,k) = Eps_min
         CmuHt(i,j,k) = Cmut_min
         CmuVt(i,j,k) = Cmut_min
         DTke(i,j,k) = D(i,j)*Tke(i,j,k)*Mask(i,j)
         DEps(i,j,k) = D(i,j)*Eps(i,j,k)*Mask(i,j)
       enddo
       enddo
       enddo
     endif

# if defined (SALINITY)

     ! read initial salinity field from an input file
     if(INITIAL_SALI) then
       open(23,file='sali0.txt',status='old')
       do k = 1,Kglob
       do j = 1,Nglob
         read(23,*) (SaliG(i,j,k),i=1,Mglob)
       enddo
       enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         do k = Kbeg,Kend
           Sali(i,j,k) = SaliG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         do k = Kbeg,Kend
           Sali(i,j,k) = SaliG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo
# endif
     else
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         Sali(i,j,k) = (1.0-sigc(k))*0.01
       enddo
       enddo
       enddo
     endif

     ! boundary conditions
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       Sin_X0(j,k) = Sali(Ibeg,j,k)
       Sin_Xn(j,k) = Sali(Iend,j,k)
     enddo
     enddo

     call sali_bc
# if defined (PARALLEL)
     call phi_3D_exch(Sali)
# endif

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 1,Kloc
       DSali(i,j,k) = D(i,j)*Sali(i,j,k)
     enddo
     enddo
     enddo     

     DSali0 = DSali

     call eval_dens
# endif

# if defined (LANDSLIDE)
     if(trim(SlideType)=='DEFORMABLE') then
# if defined (SEDIMENT)
       SlideAngle = SlideAngle*pi/180.
       conc_slide = (SlideDens-Rho0)/(SRho-Rho0)
       ls1 = SlideX0-0.5*SlideL*cos(SlideAngle)
       ls2 = SlideX0+0.5*SlideL*cos(SlideAngle)

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         lsx = (xc(i)-SlideX0)/cos(SlideAngle)
         if(xc(i)>=ls1.and.xc(i)<=ls2) then
           hslide = SlideT*(1.0-(2.0*lsx/SlideL)**2)
           ! The following slide height is an approximation.
           Zslide(i,j) = hslide/cos(SlideAngle)
         else
           Zslide(i,j) = 0.0
         endif

         do k = Kbeg,Kend
            Zlev = sigc(k)*Hc(i,j)
            if(Zlev<=Zslide(i,j)) then
              Conc(i,j,k) = conc_slide
            else
              Conc(i,j,k) = 0.0
            endif
            DConc(i,j,k) = D(i,j)*Conc(i,j,k)
         enddo
       enddo
       enddo

!       ! lab case
!       do i = Ibeg,Iend
!       do j = Jbeg,Jend
!         if(xc(i)>=1.0.and.xc(i)<=1.65) then
!           Zslide(i,j) = xc(i)-1.0
!         else
!           Zslide(i,j) = 0.0
!         endif
!
!         do k = Kbeg,Kend
!            Zlev = sigc(k)*Hc(i,j)
!            if(xc(i)>=1.0.and.xc(i)<=1.65) then
!              if(Zlev<=Zslide(i,j)) then
!                Conc(i,j,k) = conc_slide
!              else
!                Conc(i,j,k) = 0.0
!              endif
!            endif
!            DConc(i,j,k) = D(i,j)*Conc(i,j,k)
!         enddo
!       enddo
!     enddo

       call sedi_bc

# if defined (PARALLEL)
       call phi_3D_exch(DConc)
# endif

       DConc0 = DConc

       do i = 1,Mloc
       do j = 1,Nloc
       do k = 1,Kloc
         Conc(i,j,k) = DConc(i,j,k)/D(i,j)
       enddo
       enddo
       enddo

       ! update mixture density
       call eval_dens

       ! update mixture viscosity
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 1,Kloc
         if(Conc(i,j,k)<0.1) then
           Cmu(i,j,k) = Visc
         else
           Cmu(i,j,k) = Mud_Visc
         endif
       enddo
       enddo
       enddo

       ! rheology
       IsMove = 1
       if(RHEOLOGY_ON) then
         do i = 1,Mloc
         do j = 1,Nloc
         do k = 1,Kloc
           if(Conc(i,j,k)>=0.1) then
             CmuR(i,j,k) = 1.e+10
             IsMove(i,j,k) = 0
           else
             CmuR(i,j,k) = 0.0
             IsMove(i,j,k) = 1
           endif
         enddo
         enddo
         enddo
       endif
# endif
     endif
# endif

    ! added by Cheng for hot start
     if(HOTSTART)then
       call hot_start
     endif

!
!    SSP Runge-Kutta method parameters
!
     if(TIME_ORDER(1:3)=='THI') then
       It_Order = 3
       ALPHA(1) = 0.0
       ALPHA(2) = 3.0/4.0
       ALPHA(3) = 1.0/3.0
       BETA(1) = 1.0
       BETA(2) = 1.0/4.0
       BETA(3) = 2.0/3.0
     elseif(TIME_ORDER(1:3)=='SEC') then
       It_Order = 2
       ALPHA(1) = 0.0
       ALPHA(2) = 1.0/2.0
       BETA(1) = 1.0
       BETA(2) = 1.0/2.0
     else
       It_Order = 1
       ALPHA(1) = 0.0
       BETA(1) = 1.0
     endif
!
!    sponge layer
!
     if(SPONGE_ON) then
       call calculate_sponge
     endif
!    
     end subroutine initial
!
!------------------------------------------------------------------------------------------------------
!
